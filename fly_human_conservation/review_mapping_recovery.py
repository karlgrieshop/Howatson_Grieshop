"""Review excluded mappings; propose evidence-based changes without applying them."""

import argparse
import csv
import gzip
import json
import os
import platform
import re
import subprocess
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path

import prepare_references as prep


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "fly_human_conservation"
ORIGINAL = BASE / "mapping_4935338"
HGNC_URL = (
    "https://storage.googleapis.com/public-download-files/"
    "hgnc/tsv/tsv/hgnc_complete_set.txt"
)
HGNC_HELP = "https://www.genenames.org/download/archive/"
FIXED_RELEASE = "human GRCh38.p14 GCA_000001405.29 Ensembl genebuild 2026-04"
TABLE_COLUMNS = [
    "FlyBaseID", "category", "original_reason", "selected_HGNC_ID",
    "selected_human_symbol", "selected_DIOPT_score", "original_candidate_ENSGs",
    "original_one_to_one", "source_disease_records", "decision",
    "proposed_HGNC_ID", "proposed_ENSG_ID", "proposed_protein_ID",
    "evidence_source", "evidence_release", "evidence", "residual_uncertainty",
]


def fetch(url, path, json_response=False):
    headers = {"Accept": "application/json"} if json_response else {}
    request = urllib.request.Request(url, headers=headers)
    with urllib.request.urlopen(request, timeout=120) as response:
        if response.status != 200:
            raise ValueError(f"Incomplete HTTP response: {url}")
        with path.open("xb") as handle:
            while chunk := response.read(1024 * 1024):
                handle.write(chunk)
        expected = response.headers.get("Content-Length")
        modified = response.headers.get("Last-Modified", "")
    if not path.stat().st_size or (expected and path.stat().st_size != int(expected)):
        raise ValueError(f"Empty/incomplete source retrieval: {url}")
    return {
        "path": str(path), "url": url, "retrieved_utc": prep.now(),
        "last_modified": modified, "bytes": path.stat().st_size,
        "sha256": prep.sha256(path),
    }


def read_inputs():
    audit = prep.read_rows(
        ORIGINAL / "tables/gene_mapping_audit.csv", ["FlyBaseID", "status", "HGNC_ID"]
    )
    selected_rows = prep.read_rows(
        ORIGINAL / "tables/selected_pairs.csv", ["FlyBaseID", "human_symbol"]
    )
    selected = {row["FlyBaseID"]: row for row in selected_rows}
    if len(selected) != len(selected_rows) or len({row["FlyBaseID"] for row in audit}) != len(audit):
        raise ValueError("Duplicate original gene/pair audit")
    excluded = [row for row in audit if row["status"] != "ready_for_alignment"]
    if len(audit) != 9092 or len(excluded) != 657:
        raise ValueError("Original eligible/excluded universes changed")
    source_rows = prep.read_rows(
        ROOT / "Data/outputs/sexspecific_disease_gene_data.csv",
        ["FlyBaseID", "Human.Symbol", "source_type", "source_file"],
    )
    sources = defaultdict(list)
    for row in source_rows:
        sources[row["FlyBaseID"]].append(row)
    return audit, excluded, selected, sources


def initial_rows(excluded, selected, sources):
    rows = []
    for original in excluded:
        gene = original["FlyBaseID"]
        pair = selected.get(gene, {})
        records = sources.get(gene, []) if original["status"] == "no_supplied_orthology_pair" else []
        reason = original["status"]
        if reason == "HGNC_ambiguous":
            uncertainty = (
                "Need authoritative HGNC-to-ENSG assignment corroborated by fixed-release "
                "gene annotation/direct xrefs; protein availability alone is insufficient."
            )
        elif reason == "HGNC_unmapped" and original["HGNC_ID"] == "HGNC:0":
            uncertainty = (
                "HGNC:0 is a placeholder, not a valid unique partner identifier. "
                "Trace the original human symbol/LOC identifier; do not change partner "
                "or redefine the original full-graph one-to-one flag."
            )
        elif reason == "no_supplied_orthology_pair":
            uncertainty = (
                "No supplied pair or DIOPT score; original raw SSDO source files are "
                "absent here. A sourced human symbol is a lead, not a ranked pair. "
                "Any added pair needs evidence and separate approval; one-to-one unassigned."
            )
        else:
            uncertainty = (
                "Inspect fixed-release biotype/CDS/peptide records; do not synthesize "
                "a protein, change release, or substitute another orthologue."
            )
        rows.append({
            "FlyBaseID": gene, "category": original["category"],
            "original_reason": reason, "selected_HGNC_ID": original["HGNC_ID"],
            "selected_human_symbol": pair.get("human_symbol", ""),
            "selected_DIOPT_score": original["DIOPT_score"],
            "original_candidate_ENSGs": original["human_gene_candidates"],
            "original_one_to_one": original["one_to_one"] if pair else "unassigned_no_pair",
            "source_disease_records": json.dumps(records, sort_keys=True),
            "decision": "pending_evidence_review",
            "proposed_HGNC_ID": "", "proposed_ENSG_ID": "", "proposed_protein_ID": "",
            "evidence_source": (
                "mapping_4935338/tables/gene_mapping_audit.csv;selected_pairs.csv"
                + (";Data/outputs/sexspecific_disease_gene_data.csv" if records else "")
            ),
            "evidence_release": FIXED_RELEASE,
            "evidence": "Original selected partner and exclusion preserved.",
            "residual_uncertainty": uncertainty,
        })
    return rows


def hgnc_documents(path):
    records = prep.read_rows(
        path, ["hgnc_id", "symbol", "status", "ensembl_gene_id", "locus_group", "locus_type"]
    )
    by_id, by_symbol = {}, defaultdict(list)
    for row in records:
        if row["hgnc_id"] in by_id:
            raise ValueError(f"Duplicate authoritative HGNC ID: {row['hgnc_id']}")
        by_id[row["hgnc_id"]] = row
        by_symbol[row["symbol"]].append(row)
    return by_id, by_symbol


def official_gene(document):
    value = document.get("ensembl_gene_id", "")
    return value if re.fullmatch(r"ENSG[0-9]+", value) else ""


def propose(row, document, annotation, representative, direct_ids):
    """Use authoritative identity + annotation, never protein availability alone."""
    row = row.copy()
    gene = official_gene(document)
    identity = {
        key: document.get(key, "") for key in (
            "hgnc_id", "symbol", "status", "locus_group", "locus_type",
            "ensembl_gene_id", "location", "date_modified", "prev_symbol", "alias_symbol",
        )
    }
    row["evidence"] = json.dumps({
        "HGNC": identity, "fixed_reference_gene": annotation,
        "direct_reference_HGNC_ids": sorted(direct_ids),
        "representative": representative,
    }, sort_keys=True)
    if document.get("status") != "Approved" or not gene:
        row["decision"] = "retain_excluded_authoritative_mapping_incomplete"
        return row
    if row["selected_HGNC_ID"] != "HGNC:0" and document["hgnc_id"] != row["selected_HGNC_ID"]:
        raise ValueError("Review attempted to change a selected HGNC partner")
    if row["selected_HGNC_ID"] == "HGNC:0" and document["symbol"] != row["selected_human_symbol"]:
        raise ValueError("Placeholder recovery requires an exact authoritative partner symbol")
    row["proposed_HGNC_ID"] = document["hgnc_id"]
    row["proposed_ENSG_ID"] = gene
    row["decision"] = "manual_review_required"
    if document.get("locus_group") != "protein-coding gene":
        row["decision"] = "retain_excluded_HGNC_nonprotein_locus"
        row["residual_uncertainty"] = (
            "Authoritative HGNC classifies the selected locus as non-protein-coding. "
            "No protein substitution or historical release change authorized."
        )
        return row
    if not annotation:
        row["residual_uncertainty"] = "HGNC-designated ENSG is absent from fixed genes.gtf.bgz."
        return row
    if (
        row["original_reason"] == "HGNC_ambiguous"
        and gene not in row["original_candidate_ENSGs"].split(";")
    ):
        row["residual_uncertainty"] = "Current HGNC ENSG is outside original candidate set; manual history review required."
        return row
    names = {document["symbol"]}
    names.update(item for item in document.get("prev_symbol", "").split("|") if item)
    if annotation["gene_name"] not in names:
        row["residual_uncertainty"] = "HGNC/previous approved symbol does not match fixed gene_name."
        return row
    chromosome = annotation["chromosome"]
    if chromosome not in {str(n) for n in range(1, 23)} | {"X", "Y", "MT"}:
        row["residual_uncertainty"] = "Selected gene lies outside main reference chromosomes; locus review required."
        return row
    location = document.get("location", "")
    if not location or not re.match(r"^" + re.escape(chromosome) + r"(?:[pq]|$)", location):
        row["residual_uncertainty"] = "Current HGNC cytogenetic location does not corroborate fixed chromosome."
        return row
    if (
        annotation["gene_biotype"] != "protein_coding"
        or document["hgnc_id"] not in direct_ids
        or direct_ids != {document["hgnc_id"]}
    ):
        row["residual_uncertainty"] = "Fixed-reference biotype/direct HGNC assignment is absent or conflicting."
        return row
    if not representative:
        row["decision"] = "retain_excluded_no_valid_fixed_protein"
        row["residual_uncertainty"] = "Identity is supported but no valid representative exists under unchanged protein rules."
        return row
    row["proposed_protein_ID"] = representative["protein_id"]
    row["decision"] = "propose_recovery_authoritative_identity"
    row["residual_uncertainty"] = (
        "Proposal only: current HGNC designates this ENSG; fixed annotation agrees in "
        "approved/previous symbol, chromosome, coding biotype and direct HGNC xref. "
        "Protein validity is a subsequent eligibility check, not the identity criterion. "
        "Current HGNC evidence is newer than reference genebuild; requires review before application."
    )
    return row


def scan_gene_metadata(path, wanted):
    genes = {}
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise ValueError("Invalid reference GTF")
            if fields[2] != "gene":
                continue
            attr = prep.attributes(fields[8])
            gene = prep.stable(attr.get("gene_id", ""))
            if gene not in wanted:
                continue
            record = {
                "gene_id": gene, "gene_name": attr.get("gene_name", ""),
                "gene_biotype": attr.get("gene_biotype", ""),
                "chromosome": fields[0], "start": fields[3], "end": fields[4],
                "strand": fields[6], "source": fields[1],
            }
            if gene in genes:
                raise ValueError(f"Duplicate candidate gene locus in fixed GTF: {gene}")
            genes[gene] = record
    return genes


def scan_xrefs(path, genes, symbols):
    direct = defaultdict(set)
    selected = []
    loc_identifiers = {
        symbol.removeprefix("LOC") for symbol in symbols
        if re.fullmatch(r"LOC[0-9]+", symbol)
    }
    with gzip.open(path, "rt") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"gene_stable_id", "xref_id", "xref_label", "db_name", "info_type"}
        if not required <= set(reader.fieldnames or []):
            raise ValueError("Fixed xref schema changed")
        for row in reader:
            gene = prep.stable(row["gene_stable_id"])
            is_loc = (
                row["db_name"] == "NCBI gene (formerly Entrezgene)"
                and row["xref_id"] in loc_identifiers
            )
            if gene not in genes and row["xref_label"] not in symbols and not is_loc:
                continue
            if row["info_type"] == "DIRECT" and row["db_name"] == "HGNC Symbol":
                direct[gene].add(row["xref_id"])
            if row["db_name"] in {
                "HGNC Symbol", "NCBI gene (formerly Entrezgene)", "MIM gene",
            }:
                selected.append(row)
    return direct, selected


def verify_references(reference_source=None):
    reference_source = reference_source or BASE
    sources = prep.read_rows(
        reference_source / "provenance/reference_manifest.tsv",
        ["path", "bytes", "sha256", "verified", "provider_checksum"],
    )
    for row in sources:
        path = reference_source / row["path"]
        if path.stat().st_size != int(row["bytes"]) or prep.sha256(path) != row["sha256"]:
            raise ValueError(f"Reference differs from successful preparation: {path}")
        relative = row["path"].removeprefix("references/")
        if relative in prep.EXPECTED:
            if row["verified"] != "True" or row["provider_checksum"] != prep.EXPECTED[relative]:
                raise ValueError(f"Original reference checksum validation missing: {relative}")
            prep.verify_provider_checksum(reference_source / "references", relative)
    return sources


def category_counts(audit, reviews):
    by_gene = {row["FlyBaseID"]: row for row in reviews}
    rows = []
    for category in prep.CATEGORIES:
        group = [row for row in audit if row["category"] == category]
        initial = [row for row in group if row["status"] == "ready_for_alignment"]
        candidates = [
            row for row in group if row["FlyBaseID"] in by_gene
            and by_gene[row["FlyBaseID"]]["decision"] == "propose_recovery_authoritative_identity"
        ]
        one_to_one = [row for row in group if row["one_to_one"] == "True"]
        rows.append({
            "category": category, "canonical_eligible": len(group),
            "original_measured": len(initial), "original_excluded": len(group) - len(initial),
            "proposed_additional_mappings_not_applied": len(candidates),
            "projected_mapped_if_approved": len(initial) + len(candidates),
            "remaining_unresolved_or_unsupported": len(group) - len(initial) - len(candidates),
            "original_full_graph_one_to_one_eligible": len(one_to_one),
            "original_one_to_one_measured": sum(row["one_to_one"] == "True" for row in initial),
            "original_one_to_one_excluded": sum(row["status"] != "ready_for_alignment" for row in one_to_one),
            "proposed_one_to_one_additional": sum(row["one_to_one"] == "True" for row in candidates),
        })
    return rows


def full_review(output, rows, reference_source=None):
    reference_source = reference_source or BASE
    references = verify_references(reference_source)
    source_manifest = [fetch(HGNC_URL, output / "downloads/hgnc_complete_set.tsv")]
    by_id, by_symbol = hgnc_documents(output / "downloads/hgnc_complete_set.tsv")
    documents = {}
    for row in rows:
        hgnc = row["selected_HGNC_ID"]
        if hgnc == "HGNC:0":
            matches = by_symbol.get(row["selected_human_symbol"], [])
            if len(matches) == 1:
                documents[row["FlyBaseID"]] = matches[0]
        elif hgnc in by_id:
            documents[row["FlyBaseID"]] = by_id[hgnc]
    wanted = {
        gene for row in rows for gene in row["original_candidate_ENSGs"].split(";") if gene
    } | {official_gene(doc) for doc in documents.values() if official_gene(doc)}
    gtf = reference_source / "references/human/genes.gtf.bgz"
    genes = scan_gene_metadata(gtf, wanted)
    symbols = {row["selected_human_symbol"] for row in rows if row["selected_HGNC_ID"] == "HGNC:0"}
    direct, xrefs = scan_xrefs(reference_source / "references/human/xref.tsv.gz", wanted, symbols)
    prep.write_rows(output / "tables/candidate_reference_xrefs.csv", xrefs, list(xrefs[0]) if xrefs else [
        "gene_stable_id", "transcript_stable_id", "protein_stable_id", "xref_id",
        "xref_label", "description", "db_name", "info_type", "source", "ensembl_identity", "xref_identity",
    ])
    proteins, canonical, coding, header = prep.annotation(gtf, wanted, "human")
    reps, protein_audit = prep.representative_proteins(
        reference_source / "references/human/pep.fa.bgz", proteins, canonical, coding, "human",
    )
    prep.write_rows(output / "tables/candidate_proteins.csv", protein_audit, [
        "species", "gene_id", "protein_id", "canonical", "status", "length_aa", "selected",
    ])
    candidate_rows = []
    for gene in sorted(wanted):
        candidate_rows.append({
            "ENSG_ID": gene, "annotation": json.dumps(genes.get(gene, {}), sort_keys=True),
            "direct_HGNC_ids": ";".join(sorted(direct[gene])),
            "representative": json.dumps(reps.get(gene, {}), sort_keys=True),
        })
    prep.write_rows(output / "tables/candidate_gene_evidence.csv", candidate_rows, [
        "ENSG_ID", "annotation", "direct_HGNC_ids", "representative",
    ])
    used = {doc["hgnc_id"]: doc for doc in documents.values()}
    if used:
        all_columns = sorted(set().union(*(doc.keys() for doc in used.values())))
        prep.write_rows(
            output / "tables/authoritative_HGNC_records.csv",
            [{key: doc.get(key, "") for key in all_columns} for _, doc in sorted(used.items())],
            all_columns,
        )
    resolved = []
    for row in rows:
        row = row.copy()
        row["evidence_source"] += ";" + HGNC_URL + ";genes.gtf.bgz;xref.tsv.gz"
        row["evidence_release"] += ";HGNC current snapshot retrieved " + source_manifest[0]["retrieved_utc"]
        doc = documents.get(row["FlyBaseID"])
        if doc:
            gene = official_gene(doc)
            row = propose(row, doc, genes.get(gene, {}), reps.get(gene, {}), direct[gene])
        elif row["original_reason"] == "HGNC_unmapped":
            symbol = row["selected_human_symbol"]
            matches = [
                r for r in xrefs if r["xref_label"] == symbol or (
                    re.fullmatch(r"LOC[0-9]+", symbol)
                    and r["db_name"] == "NCBI gene (formerly Entrezgene)"
                    and r["xref_id"] == symbol.removeprefix("LOC")
                )
            ]
            historical = [
                doc for doc in by_id.values()
                if symbol in set(doc.get("prev_symbol", "").split("|")) | set(doc.get("alias_symbol", "").split("|"))
            ]
            row["evidence"] = json.dumps({
                "original_symbol": symbol,
                "fixed_xrefs_matching_symbol": matches,
                "current_HGNC_alias_or_previous_symbol_leads": historical,
            }, sort_keys=True)
            row["decision"] = "manual_identifier_history_review_required"
            row["residual_uncertainty"] += (
                " No unique current approved HGNC symbol match. Do not equate an alias "
                "or LOC label to another paralogue without authoritative history."
            )
        resolved.append(row)
    prep.write_rows(output / "provenance/reference_checks.csv", references, list(references[0]))
    return resolved, source_manifest


def main(preflight=False, run_id=None, reference_source=None):
    run_id = prep.local_run_id(run_id)
    reference_source = Path(reference_source).resolve() if reference_source else BASE
    if not reference_source.is_relative_to(ROOT):
        raise ValueError("Reference source must be within the repository")
    name = f"recovery_review_local_{run_id}"
    output = BASE / name
    if output.exists():
        raise FileExistsError(f"Refusing to overwrite {output}")
    audit, excluded, selected, sources = read_inputs()
    input_paths = [
        ORIGINAL / "tables" / name for name in (
            "gene_mapping_audit.csv", "selected_pairs.csv", "all_pairs.csv",
            "hgnc_mapping_audit.csv", "excluded_human_xrefs.csv", "category_attrition.csv",
        )
    ] + [
        ROOT / "Data/outputs/sexspecific_disease_gene_data.csv",
        Path(__file__), BASE / "prepare_references.py", BASE / "review_mapping_recovery.sh",
    ]
    inputs = [
        {"path": str(path.relative_to(ROOT)), "bytes": path.stat().st_size,
         "sha256": prep.sha256(path)} for path in input_paths
    ]
    output.mkdir()
    for folder in ("tables", "provenance", "downloads"):
        (output / folder).mkdir()
    prep.write_rows(output / "provenance/input_manifest.csv", inputs, list(inputs[0]))
    with (output / "provenance/execution.json").open("x") as handle:
        json.dump({
            "started_utc": prep.now(), "preflight": preflight, "run_id": run_id,
            "python": os.sys.version, "platform": platform.platform(),
            "command": [os.sys.executable, *os.sys.argv],
            "wrapper_command": os.environ.get("WORKFLOW_WRAPPER_COMMAND"),
            "reference_source": str(reference_source.relative_to(ROOT)),
            "fixed_reference_release": FIXED_RELEASE,
            "HGNC_source_documentation": HGNC_HELP,
            "git_commit": subprocess.run(
                ["git", "-C", str(ROOT), "rev-parse", "HEAD"],
                check=True, capture_output=True, text=True,
            ).stdout.strip(),
            "biological_changes_applied": False,
        }, handle, indent=2)
        handle.write("\n")
    rows = initial_rows(excluded, selected, sources)
    if preflight:
        remote = []
        for hgnc in ("HGNC:6698", "HGNC:13605"):
            path = output / "provenance" / (hgnc.replace(":", "_") + ".json")
            url = "https://rest.genenames.org/fetch/hgnc_id/" + hgnc
            remote.append(fetch(url, path, json_response=True))
            response = json.loads(path.read_text())["response"]
            if response["numFound"] != 1 or response["docs"][0]["hgnc_id"] != hgnc:
                raise ValueError(f"Unexpected HGNC response for {hgnc}")
            doc = response["docs"][0]
            for row in rows:
                if row["selected_HGNC_ID"] != hgnc:
                    continue
                row["evidence_source"] += ";" + url
                row["evidence_release"] += ";HGNC live retrieved " + remote[-1]["retrieved_utc"]
                row["evidence"] = json.dumps(doc, sort_keys=True)
                row["proposed_HGNC_ID"] = hgnc
                row["proposed_ENSG_ID"] = doc.get("ensembl_gene_id", "")
                if doc.get("locus_group") == "pseudogene":
                    row["decision"] = "retain_excluded_HGNC_nonprotein_locus"
                    row["residual_uncertainty"] = "HGNC classifies selected FBXL9P as pseudogene; fixed annotation inspection pending. No historical protein or alternate orthologue authorized."
                else:
                    row["decision"] = "authoritative_ENSG_lead_pending_fixed_annotation"
    else:
        rows, remote = full_review(output, rows, reference_source)
    prep.write_rows(output / "tables/per_gene_resolution.csv", rows, TABLE_COLUMNS)
    counts = category_counts(audit, rows)
    prep.write_rows(output / "tables/category_attrition_review.csv", counts, list(counts[0]))
    for entry in remote:
        entry["path"] = str(Path(entry["path"]).relative_to(output))
    prep.write_rows(output / "provenance/source_manifest.csv", remote, [
        "path", "url", "retrieved_utc", "last_modified", "bytes", "sha256",
    ])
    for row in inputs:
        if prep.sha256(ROOT / row["path"]) != row["sha256"]:
            raise ValueError(f"Review input changed during execution: {row['path']}")
    output_manifest = [
        {"path": str(path.relative_to(output)), "bytes": path.stat().st_size,
         "sha256": prep.sha256(path)}
        for path in sorted(output.rglob("*")) if path.is_file()
    ]
    prep.write_rows(
        output / "provenance/output_manifest.csv",
        output_manifest, ["path", "bytes", "sha256"],
    )
    result = {
        "completed_utc": prep.now(), "mode": "preflight" if preflight else "full_reference_review",
        "reviewed_excluded_fly_genes": len(rows),
        "decisions": dict(Counter(row["decision"] for row in rows)),
        "category_counts": counts,
        "applied_mapping_changes": 0, "alignments_run": 0,
        "output_manifest_sha256": prep.sha256(output / "provenance/output_manifest.csv"),
        "note": "Proposed identities require review/approval; existing outputs are unchanged.",
    }
    with (output / "provenance/REVIEW_COMPLETE.json").open("x") as handle:
        json.dump(result, handle, indent=2)
        handle.write("\n")
    prep.log(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--preflight", action="store_true", help="Compact audits and two bounded HGNC samples only.")
    parser.add_argument("--run-id", help="Unique local run identifier (default: UTC timestamp).")
    parser.add_argument("--reference-source", type=Path,
                        help="Directory containing references/ and provenance/reference_manifest.tsv.")
    args = parser.parse_args()
    main(preflight=args.preflight, run_id=args.run_id, reference_source=args.reference_source)
