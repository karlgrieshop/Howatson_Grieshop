"""Prepare fixed-release references and auditable gene/protein mappings only."""

import argparse
import csv
import gzip
import hashlib
import json
import os
import platform
import re
import shutil
import subprocess
import sys
import urllib.request
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SIDE = ROOT / "fly_human_conservation"
HUMAN_URL = (
    "https://ftp.ebi.ac.uk/pub/ensemblorganisms/"
    "GCA/000/001/405/29/ensembl/2026_04/geneset/"
)
FLY_PEP_URL = (
    "https://ftp.ensembl.org/pub/release-115/"
    "fasta/drosophila_melanogaster/pep/"
)
FLY_GTF_URL = (
    "https://ftp.ensembl.org/pub/release-115/"
    "gtf/drosophila_melanogaster/"
)
FLY_PEP = "Drosophila_melanogaster.BDGP6.54.pep.all.fa.gz"
FLY_GTF = "Drosophila_melanogaster.BDGP6.54.115.gtf.gz"
EXPECTED = {
    "human/pep.fa.bgz": "7b709f5f5e6901086a3b4b5c2c8d534e",
    "human/genes.gtf.bgz": "9a46e3af809a6ffcd4002dd4a02c7295",
    "human/xref.tsv.gz": "be11fab9c0bdbb842304c22cbf6333c4",
    "fly/" + FLY_PEP: "53185 6220",
    "fly/" + FLY_GTF: "10582 6962",
}
INPUTS = [
    "Data/inputs/all_orthologs.csv",
    "Data/inputs/all_disease_orthologs.csv",
    "Data/outputs/sexspecific_disease_ids.csv",
    "Data/outputs/all_disease_ortholog_ids.csv",
    "Data/outputs/all_ortholog_ids.csv",
    "Data/outputs/dmel_gene_reference.csv",
]
CATEGORIES = (
    "SSDOs",
    "Disease orthologues excluding SSDOs",
    "Orthologues excluding disease genes",
)
PAIR_COLUMNS = (
    "Dmel_gene_ID", "Human_gene_HGNC_ID", "DIOPT_score",
    "Dmel_gene_symbol", "Human_gene_symbol",
)
PROTEIN_ALPHABET = set("ACDEFGHIKLMNPQRSTVWYBXZJUO")
RECOVERY_REVIEW_SHA256 = (
    "29f9d66162539f99a1ceceff1e367c03cfbd032c7eb54a10e6e1fca4f31a48fb"
)
XREF_EXCLUSION_COLUMNS = [
    "HGNC_ID", "gene_stable_id", "transcript_stable_id", "protein_stable_id",
    "xref_label", "db_name", "info_type", "source", "exclusion_reason",
]


def now():
    return datetime.now(timezone.utc).isoformat()


def log(message):
    print(f"{now()} {message}", flush=True)


def stable(identifier):
    return identifier.split(".", 1)[0]


def write_rows(path, rows, columns):
    with path.open("x", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=columns, extrasaction="raise",
            delimiter="\t" if path.suffix == ".tsv" else ",",
        )
        writer.writeheader()
        writer.writerows(rows)


def read_rows(path, required):
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t" if path.suffix == ".tsv" else ",")
        if not set(required).issubset(reader.fieldnames or []):
            raise ValueError(f"{path}: missing required columns {required}")
        return list(reader)


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def gene_ids(root, filename):
    rows = read_rows(root / "Data/outputs" / filename, ["FlyBaseID"])
    ids = [row["FlyBaseID"] for row in rows]
    if len(ids) != len(set(ids)) or any(
        not re.fullmatch(r"FBgn[0-9]+", value) for value in ids
    ):
        raise ValueError(f"{filename}: duplicate or invalid FlyBase IDs")
    return set(ids)


def prepare_pairs(root):
    ssdo = gene_ids(root, "sexspecific_disease_ids.csv")
    disease = gene_ids(root, "all_disease_ortholog_ids.csv")
    orthologues = gene_ids(root, "all_ortholog_ids.csv")
    reference = read_rows(
        root / "Data/outputs/dmel_gene_reference.csv",
        ["FlyBaseID", "gene_biotype"],
    )
    coding = {
        row["FlyBaseID"] for row in reference
        if row["gene_biotype"] == "protein_coding"
    }
    if not ssdo <= disease <= orthologues <= coding:
        raise ValueError("Canonical nesting/protein-coding invariant failed")
    categories = {}
    for label, ids in zip(CATEGORIES, (ssdo, disease - ssdo, orthologues - disease)):
        categories.update({gene: label for gene in ids})

    tables = {}
    for name in ("all_orthologs", "all_disease_orthologs"):
        rows = read_rows(root / f"Data/inputs/{name}.csv", PAIR_COLUMNS)
        pairs = {}
        for row in rows:
            fly, human, score = (row[key] for key in PAIR_COLUMNS[:3])
            if not (
                re.fullmatch(r"FBgn[0-9]+", fly)
                and re.fullmatch(r"HGNC:[0-9]+", human)
                and re.fullmatch(r"[0-9]+", score)
            ):
                raise ValueError(f"{name}: invalid pair ID or DIOPT score: {row}")
            key = (fly, human)
            if key in pairs:
                raise ValueError(f"{name}: duplicate pair {key}")
            pairs[key] = row
        tables[name] = pairs
    all_pairs = tables["all_orthologs"]
    disease_pairs = tables["all_disease_orthologs"]
    if not disease_pairs.keys() <= all_pairs.keys():
        raise ValueError("Disease pair table is no longer a subset of all pairs")
    if any(
        int(row["DIOPT_score"]) != int(all_pairs[key]["DIOPT_score"])
        for key, row in disease_pairs.items()
    ):
        raise ValueError("DIOPT scores conflict between input tables")

    by_fly, by_human = defaultdict(list), defaultdict(set)
    for (fly, human), row in all_pairs.items():
        by_fly[fly].append(row)
        by_human[human].add(fly)
    selected = {}
    for fly in sorted(categories):
        if fly not in by_fly:
            continue
        candidates = sorted(
            by_fly[fly],
            key=lambda row: (-int(row["DIOPT_score"]), row["Human_gene_HGNC_ID"]),
        )
        winner = candidates[0]
        human = winner["Human_gene_HGNC_ID"]
        selected[fly] = {
            "FlyBaseID": fly,
            "category": categories[fly],
            "HGNC_ID": human,
            "DIOPT_score": int(winner["DIOPT_score"]),
            "fly_symbol": winner["Dmel_gene_symbol"],
            "human_symbol": winner["Human_gene_symbol"],
            "fly_partner_count": len(candidates),
            "human_partner_count": len(by_human[human]),
            "top_score_tie_count": sum(
                int(row["DIOPT_score"]) == int(winner["DIOPT_score"])
                for row in candidates
            ),
            "one_to_one": len(candidates) == 1 and len(by_human[human]) == 1,
        }
    return categories, tables, selected


def download(url, path):
    log(f"Retrieving {url}")
    request = urllib.request.Request(
        url, headers={"User-Agent": "fly-human-conservation-reference-preparation"}
    )
    with urllib.request.urlopen(request, timeout=120) as response:
        if response.status != 200:
            raise ValueError(f"Expected complete HTTP 200 response: {url}")
        expected_length = response.headers.get("Content-Length")
        modified = response.headers.get("Last-Modified", "")
        with path.open("xb") as handle:
            shutil.copyfileobj(response, handle, length=1024 * 1024)
    if not path.stat().st_size:
        raise ValueError(f"Empty download: {url}")
    if expected_length and path.stat().st_size != int(expected_length):
        raise ValueError(f"Incomplete download: {url}")
    return {
        "url": url, "retrieved_utc": now(), "last_modified": modified,
        "bytes": path.stat().st_size, "sha256": sha256(path),
    }


def reference_downloads():
    return [
        ("human/md5sum.txt", HUMAN_URL + "md5sum.txt"),
        ("human/pep.fa.bgz", HUMAN_URL + "pep.fa.bgz"),
        ("human/genes.gtf.bgz", HUMAN_URL + "genes.gtf.bgz"),
        ("human/xref.tsv.gz", HUMAN_URL + "xref.tsv.gz"),
        ("fly/pep.CHECKSUMS", FLY_PEP_URL + "CHECKSUMS"),
        ("fly/gtf.CHECKSUMS", FLY_GTF_URL + "CHECKSUMS"),
        ("fly/" + FLY_PEP, FLY_PEP_URL + FLY_PEP),
        ("fly/" + FLY_GTF, FLY_GTF_URL + FLY_GTF),
    ]


def verify_provider_checksum(refs, relative):
    path = refs / relative
    if relative.startswith("human/"):
        declared = {
            line.split()[1].lstrip("*"): line.split()[0]
            for line in (refs / "human/md5sum.txt").read_text().splitlines()
            if line.strip()
        }
        published = declared.get(path.name)
        digest = hashlib.md5()
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
        actual = digest.hexdigest()
        method = "MD5"
    else:
        filename = "pep.CHECKSUMS" if path.name == FLY_PEP else "gtf.CHECKSUMS"
        declared = {
            parts[2]: " ".join(parts[:2])
            for line in (refs / "fly" / filename).read_text().splitlines()
            if len(parts := line.split()) == 3
        }
        published = declared.get(path.name)
        checksum = 0
        size = 0
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                size += len(chunk)
                for byte in chunk:
                    checksum = ((checksum >> 1) | ((checksum & 1) << 15))
                    checksum = (checksum + byte) & 0xffff
        actual = f"{checksum} {(size + 1023) // 1024}"
        method = "BSD sum -r (checksum, 1KiB blocks)"
    if published != EXPECTED[relative] or actual != published:
        raise ValueError(
            f"Checksum mismatch: {relative}; pinned={EXPECTED[relative]}, "
            f"provider={published}, downloaded={actual}"
        )
    return {"checksum_method": method, "provider_checksum": published, "verified": True}


def reference_files(side):
    refs = side / "references"
    (refs / "human").mkdir()
    (refs / "fly").mkdir()
    manifest = []
    for relative, url in reference_downloads():
        path = refs / relative
        entry = download(url, path)
        entry["path"] = str(path.relative_to(side))
        entry["checksum_method"] = "SHA256"
        entry["provider_checksum"] = ""
        entry["verified"] = ""
        if relative in EXPECTED:
            entry.update(verify_provider_checksum(refs, relative))
        manifest.append(entry)
        log(f"Retrieved and fingerprinted {relative}")
    write_rows(
        side / "provenance/reference_manifest.tsv",
        manifest,
        ["path", "url", "retrieved_utc", "last_modified", "bytes", "sha256",
         "checksum_method", "provider_checksum", "verified"],
    )
    return refs


def reuse_references(source_side, output_side, inputs):
    old_inputs = read_rows(
        source_side / "provenance/input_manifest.csv", ["path", "bytes", "sha256"]
    )
    old_inputs = [
        row for row in old_inputs
        if row["path"] != "fly_human_conservation/recovery_review_4955170/tables/per_gene_resolution.csv"
    ]
    if (
        len(old_inputs) != len(inputs)
        or {row["path"] for row in old_inputs} != {row["path"] for row in inputs}
    ):
        raise ValueError("Reference reuse requires the original input manifest")
    current = {row["path"]: row for row in inputs}
    for row in old_inputs:
        if (
            row["sha256"] != current[row["path"]]["sha256"]
            or int(row["bytes"]) != current[row["path"]]["bytes"]
        ):
            raise ValueError(f"Input changed since the original run: {row['path']}")
    rows = read_rows(
        source_side / "provenance/reference_manifest.tsv",
        ["path", "url", "retrieved_utc", "last_modified", "bytes", "sha256",
         "checksum_method", "provider_checksum", "verified"],
    )
    expected_urls = {
        "references/" + relative: url for relative, url in reference_downloads()
    }
    if len(rows) != len(expected_urls) or {row["path"] for row in rows} != set(expected_urls):
        raise ValueError("Reference reuse requires a complete original reference manifest")
    for row in rows:
        path = source_side / row["path"]
        if row["url"] != expected_urls[row["path"]]:
            raise ValueError(f"Reference source URL mismatch: {path}")
        if path.stat().st_size != int(row["bytes"]) or sha256(path) != row["sha256"]:
            raise ValueError(f"Reference fingerprint changed since original run: {path}")
        relative = str(path.relative_to(source_side / "references"))
        if relative in EXPECTED:
            if row["verified"] != "True" or row["provider_checksum"] != EXPECTED[relative]:
                raise ValueError(f"Original provider verification missing: {path}")
            row.update(verify_provider_checksum(source_side / "references", relative))
        row["reuse_verified_utc"] = now()
        log(f"Reverified existing {row['path']}; not downloading")
    write_rows(
        output_side / "provenance/reference_manifest.tsv",
        rows, list(rows[0]),
    )
    return source_side / "references"


def verify_approved_reference_snapshot(refs):
    rows = read_rows(SIDE / "provenance/reference_manifest.tsv", ["path", "sha256", "bytes"])
    archived = {row["path"].removeprefix("references/"): row for row in rows}
    for relative in EXPECTED:
        path = refs / relative
        row = archived[relative]
        if path.stat().st_size != int(row["bytes"]) or sha256(path) != row["sha256"]:
            raise ValueError(f"Reference differs from approved recovery snapshot: {path}")


def local_run_id(run_id):
    run_id = run_id or datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    if not re.fullmatch(r"[A-Za-z0-9_-]+", run_id):
        raise ValueError(f"Invalid local run ID: {run_id}")
    return run_id


def hgnc_mapping(path, requested):
    genes = defaultdict(set)
    counts = Counter()
    excluded = []
    with gzip.open(path, "rt") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"gene_stable_id", "xref_id", "db_name", "info_type"}
        if not required <= set(reader.fieldnames or []):
            raise ValueError("Human xref schema lacks required ID/source fields")
        for row in reader:
            hgnc = row["xref_id"]
            if hgnc not in requested or row["db_name"] != "HGNC Symbol":
                continue
            if row["info_type"] != "DIRECT":
                continue
            gene = stable(row["gene_stable_id"])
            if re.fullmatch(r"LRG_[0-9]+", gene):
                excluded.append({
                    "HGNC_ID": hgnc,
                    **{key: row.get(key, "") for key in XREF_EXCLUSION_COLUMNS[1:-1]},
                    "exclusion_reason": "LRG_reference_namespace_not_ENSG",
                })
                continue
            if not re.fullmatch(r"ENSG[0-9]+", gene):
                raise ValueError(f"Unexpected gene namespace in HGNC xref: {row}")
            genes[hgnc].add(gene)
            counts[hgnc] += 1
    return genes, counts, excluded


def attributes(text):
    return {
        key: value for key, value in
        re.findall(r'([A-Za-z0-9_]+) "([^"]*)"', text)
    }


def annotation(path, wanted, species):
    proteins = defaultdict(lambda: defaultdict(set))
    canonical_transcripts = set()
    coding = set()
    header = []
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                header.append(line.rstrip())
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"Invalid GTF record in {path}")
            if fields[2] not in {"gene", "transcript", "CDS"}:
                continue
            attr = attributes(fields[8])
            gene = stable(attr.get("gene_id", ""))
            if gene not in wanted:
                continue
            if attr.get("gene_biotype") == "protein_coding":
                coding.add(gene)
            transcript = stable(attr.get("transcript_id", ""))
            if transcript and 'tag "Ensembl_canonical"' in fields[8]:
                canonical_transcripts.add((gene, transcript))
            if fields[2] == "CDS" and attr.get("protein_id"):
                if not transcript:
                    raise ValueError(f"CDS protein lacks transcript in {path}")
                proteins[gene][stable(attr["protein_id"])].add(transcript)
    expected = (
        "#!genome-build GRCh38.p14" if species == "human"
        else "#!genome-build BDGP6.54"
    )
    if expected not in header:
        raise ValueError(f"{species} reference build header mismatch: {header}")
    if species == "human" and "#!genebuild-last-updated 2026-04" not in header:
        raise ValueError("Human genebuild date mismatch")
    canonical = {
        gene: {
            protein for protein, transcripts in members.items()
            if any((gene, transcript) in canonical_transcripts for transcript in transcripts)
        }
        for gene, members in proteins.items()
    }
    return proteins, canonical, coding, header


def fasta_records(path):
    header, chunks = None, []
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header, chunks = line[1:].strip(), []
            else:
                if header is None:
                    raise ValueError(f"Sequence precedes FASTA header in {path}")
                chunks.append(line.strip())
    if header is not None:
        yield header, "".join(chunks)


def representative_proteins(path, proteins, canonical, coding, species):
    by_id = defaultdict(set)
    for gene, members in proteins.items():
        for protein in members:
            by_id[protein].add(gene)
    records = defaultdict(list)
    for header, sequence in fasta_records(path):
        protein = stable(header.split()[0])
        if protein not in by_id:
            continue
        metadata = header.split(" description:", 1)[0]
        tags = dict(
            part.split(":", 1) for part in metadata.split()[1:] if ":" in part
        )
        gene = stable(tags.get("gene", ""))
        transcript = stable(tags.get("transcript", ""))
        terminal_stop = sequence.endswith("*")
        sequence = sequence[:-1] if terminal_stop else sequence
        reasons = []
        if not sequence:
            reasons.append("empty_sequence")
        if set(sequence) - PROTEIN_ALPHABET:
            reasons.append("invalid_residue_or_internal_stop")
        if len(by_id[protein]) != 1 or gene not in by_id[protein]:
            reasons.append("gene_mapping_conflict")
        elif transcript not in proteins[gene][protein]:
            reasons.append("transcript_mapping_conflict")
        if gene not in coding:
            reasons.append("not_protein_coding_in_annotation")
        records[protein].append({
            "protein_id": protein, "versioned_protein_id": header.split()[0],
            "fasta_gene_id": gene, "transcript_id": transcript,
            "length_aa": len(sequence), "terminal_stop_removed": terminal_stop,
            "nonstandard_residue_count": sum(aa not in "ACDEFGHIKLMNPQRSTVWY" for aa in sequence),
            "status": ";".join(reasons) or "valid",
            "sequence_sha256": hashlib.sha256(sequence.encode("ascii")).hexdigest(),
        })
    representatives, audit = {}, []
    for gene in sorted(proteins):
        candidates = []
        for protein in sorted(proteins[gene]):
            entries = records.get(protein, [])
            if len(entries) != 1:
                audit.append({
                    "species": species, "gene_id": gene, "protein_id": protein,
                    "canonical": protein in canonical.get(gene, set()),
                    "status": "missing_fasta" if not entries else "ambiguous_fasta_records",
                    "length_aa": "", "selected": False,
                })
                continue
            record = entries[0]
            audit.append({
                "species": species, "gene_id": gene, "protein_id": protein,
                "canonical": protein in canonical.get(gene, set()),
                "status": record["status"], "length_aa": record["length_aa"],
                "selected": False,
            })
            if record["status"] == "valid":
                candidates.append(record)
        valid_canonical = [
            row for row in candidates
            if row["protein_id"] in canonical.get(gene, set())
        ]
        if len(canonical.get(gene, set())) == 1 and len(valid_canonical) == 1:
            winner = valid_canonical[0].copy()
            rule = "unique_valid_Ensembl_canonical"
        elif candidates:
            winner = sorted(candidates, key=lambda row: (-row["length_aa"], row["protein_id"]))[0].copy()
            rule = "longest_valid_protein_stable_ID_tiebreak"
        else:
            continue
        winner.update({
            "species": species, "gene_id": gene, "selection_rule": rule,
            "canonical_protein_count": len(canonical.get(gene, set())),
            "valid_protein_count": len(candidates),
        })
        representatives[gene] = winner
    for row in audit:
        winner = representatives.get(row["gene_id"])
        row["selected"] = bool(winner and winner["protein_id"] == row["protein_id"])
    return representatives, audit


def mapping_rows(
    categories, selected, hgnc, fly_proteins, human_proteins,
    approved_recoveries=None,
):
    rows = []
    for fly, category in sorted(categories.items()):
        pair = selected.get(fly)
        human = pair["HGNC_ID"] if pair else ""
        original_mapped = sorted(hgnc.get(human, set()))
        recovery = (approved_recoveries or {}).get(fly)
        mapped = (
            [recovery["proposed_ENSG_ID"]] if recovery else original_mapped
        )
        human_gene = mapped[0] if len(mapped) == 1 else ""
        fly_rep = fly_proteins.get(fly)
        human_rep = human_proteins.get(human_gene)
        reasons = []
        if not pair:
            reasons.append("no_supplied_orthology_pair")
        elif not mapped:
            reasons.append("HGNC_unmapped")
        elif len(mapped) > 1:
            reasons.append("HGNC_ambiguous")
        if not fly_rep:
            reasons.append("no_valid_fly_protein")
        if human_gene and not human_rep:
            reasons.append("no_valid_human_protein")
        rows.append({
            "FlyBaseID": fly, "category": category, "HGNC_ID": human,
            "DIOPT_score": pair["DIOPT_score"] if pair else "",
            "human_gene_candidates": ";".join(mapped),
            "human_gene_candidate_count": len(mapped),
            "human_gene_id": human_gene,
            "fly_protein_id": fly_rep["protein_id"] if fly_rep else "",
            "human_protein_id": human_rep["protein_id"] if human_rep else "",
            "one_to_one": pair["one_to_one"] if pair else False,
            "status": ";".join(reasons) or "ready_for_alignment",
        })
        if approved_recoveries is not None:
            rows[-1].update({
                "original_human_gene_candidates": ";".join(original_mapped),
                "original_human_gene_candidate_count": len(original_mapped),
                "recovery_decision": (
                    recovery["decision"] if recovery else "not_proposed"
                ),
                "recovery_applied": bool(recovery),
            })
    return rows


def approved_recovery_rows(
    path, categories, selected, hgnc, original_rows, human_proteins,
):
    if sha256(path) != RECOVERY_REVIEW_SHA256:
        raise ValueError(f"Recovery review differs from the approved snapshot: {path}")
    rows = read_rows(path, [
        "FlyBaseID", "category", "original_reason", "selected_HGNC_ID",
        "selected_DIOPT_score", "original_candidate_ENSGs", "original_one_to_one",
        "decision", "proposed_HGNC_ID", "proposed_ENSG_ID",
        "proposed_protein_ID", "residual_uncertainty",
    ])
    if len(rows) != 657 or len({row["FlyBaseID"] for row in rows}) != len(rows):
        raise ValueError("Approved recovery review must contain 657 unique fly genes")
    expected_decisions = {
        "propose_recovery_authoritative_identity": 628,
        "manual_identifier_history_review_required": 18,
        "pending_evidence_review": 9,
        "retain_excluded_HGNC_nonprotein_locus": 2,
    }
    if Counter(row["decision"] for row in rows) != expected_decisions:
        raise ValueError("Recovery review decisions differ from the approved 628/29 set")
    original_by_fly = {row["FlyBaseID"]: row for row in original_rows}
    originally_excluded = {
        fly for fly, row in original_by_fly.items()
        if row["status"] != "ready_for_alignment"
    }
    reviewed = {row["FlyBaseID"] for row in rows}
    if reviewed != originally_excluded:
        raise ValueError("Recovery review does not match the original excluded-gene universe")
    proposals = {}
    for row in rows:
        fly = row["FlyBaseID"]
        pair = selected.get(fly)
        original = original_by_fly[fly]
        expected_one_to_one = (
            str(original["one_to_one"]) if pair is not None
            else "unassigned_no_pair"
        )
        if (
            row["category"] != categories[fly]
            or row["selected_HGNC_ID"] != original["HGNC_ID"]
            or row["selected_DIOPT_score"] != str(original["DIOPT_score"])
            or row["original_one_to_one"] != expected_one_to_one
        ):
            raise ValueError(f"Recovery review conflicts with the original selected pair: {fly}")
        if row["decision"] != "propose_recovery_authoritative_identity":
            if row["decision"] == "retain_excluded_HGNC_nonprotein_locus":
                candidates = sorted(hgnc.get(row["selected_HGNC_ID"], set()))
                if (
                    row["proposed_HGNC_ID"] != row["selected_HGNC_ID"]
                    or row["proposed_ENSG_ID"] not in candidates
                ):
                    raise ValueError(f"Non-protein exclusion evidence changed for {fly}")
            elif any(row[key] for key in (
                "proposed_HGNC_ID", "proposed_ENSG_ID", "proposed_protein_ID",
            )):
                raise ValueError(f"Unapproved proposal present for {fly}")
            continue
        if pair is None:
            raise ValueError(f"Recovery proposal has no original selected pair: {fly}")
        candidates = sorted(hgnc.get(pair["HGNC_ID"], set()))
        if (
            row["original_reason"] != "HGNC_ambiguous"
            or len(candidates) < 2
            or row["proposed_HGNC_ID"] != pair["HGNC_ID"]
            or row["proposed_ENSG_ID"] not in candidates
            or row["original_candidate_ENSGs"] != ";".join(candidates)
            or not re.fullmatch(r"ENSP[0-9]+", row["proposed_protein_ID"])
        ):
            raise ValueError(f"Recovery proposal does not resolve its original candidates: {fly}")
        representative = human_proteins.get(row["proposed_ENSG_ID"])
        if not representative or representative["protein_id"] != row["proposed_protein_ID"]:
            raise ValueError(f"Recovery representative differs from fixed-reference selection: {fly}")
        proposals[fly] = row
    if len(proposals) != 628:
        raise ValueError(f"Expected 628 approved recovery proposals, found {len(proposals)}")
    return proposals, rows


def main(reuse_verified_references=False, recovery_review=None, run_id=None, reference_source=None):
    run_id = local_run_id(run_id)
    side = SIDE / f"mapping_local_{run_id}"
    if side.exists():
        raise FileExistsError(f"Refusing to overwrite existing {side}")
    source = Path(reference_source).resolve() if reference_source else SIDE
    if not source.is_relative_to(ROOT):
        raise ValueError("Reference source must be within the repository")
    for relative in INPUTS:
        if not (ROOT / relative).is_file():
            raise FileNotFoundError(ROOT / relative)
    output_dirs = (
        ("tables", "provenance") if reuse_verified_references
        else ("references", "tables", "provenance")
    )
    for name in output_dirs:
        if (side / name).exists():
            raise FileExistsError(f"Refusing to overwrite existing {side / name}")
    categories, tables, selected = prepare_pairs(ROOT)
    side.mkdir()
    for name in output_dirs:
        (side / name).mkdir()
    table_dir = side / "tables"
    canonical_inputs = [
        {"path": relative, "bytes": (ROOT / relative).stat().st_size,
         "sha256": sha256(ROOT / relative)}
        for relative in INPUTS
    ]
    inputs = list(canonical_inputs)
    if recovery_review:
        recovery_path = Path(recovery_review).resolve()
        if not recovery_path.is_relative_to(ROOT):
            raise ValueError("Recovery review must be within the repository")
        inputs.append({
            "path": str(recovery_path.relative_to(ROOT)),
            "bytes": recovery_path.stat().st_size,
            "sha256": sha256(recovery_path),
        })
    write_rows(side / "provenance/input_manifest.csv", inputs, ["path", "bytes", "sha256"])
    execution = {
        "started_utc": now(), "python": sys.version, "platform": platform.platform(),
        "run_id": run_id, "command": [sys.executable, *sys.argv],
        "wrapper_command": os.environ.get("WORKFLOW_WRAPPER_COMMAND"),
        "git_commit": subprocess.run(
            ["git", "-C", str(ROOT), "rev-parse", "HEAD"],
            check=True, capture_output=True, text=True,
        ).stdout.strip(),
        "entry_points_sha256": {
            name: sha256(SIDE / name)
            for name in (
                "prepare_references.py", "prepare_references.sh",
                *(("recover_mapping.sh",) if recovery_review else ()),
            )
        },
        "stage": "reference_preparation_and_mapping_only",
        "reuse_verified_references": reuse_verified_references,
        "reference_manifest_source_directory": str(
            (source if reuse_verified_references else side).relative_to(ROOT)
        ),
    }
    if recovery_review:
        execution["approved_recovery_review"] = inputs[-1]
        execution["recovery_review_sha256"] = RECOVERY_REVIEW_SHA256
    with (side / "provenance/execution.json").open("x") as handle:
        json.dump(execution, handle, indent=2)
        handle.write("\n")
    log(f"Canonical genes={len(categories)}, selected pairs={len(selected)}")
    refs = (
        reuse_references(source, side, canonical_inputs)
        if reuse_verified_references else reference_files(side)
    )
    if recovery_review:
        verify_approved_reference_snapshot(refs)
    requested = {human for fly, human in tables["all_orthologs"]}
    log("Reading full human HGNC xref")
    hgnc, xref_counts, excluded_xrefs = hgnc_mapping(refs / "human/xref.tsv.gz", requested)
    excluded_counts = Counter(row["HGNC_ID"] for row in excluded_xrefs)
    write_rows(
        table_dir / "excluded_human_xrefs.csv", excluded_xrefs, XREF_EXCLUSION_COLUMNS,
    )
    log(f"Audited {len(excluded_xrefs)} LRG xref rows separately from ENSG mappings")
    write_rows(
        table_dir / "hgnc_mapping_audit.csv",
        [
            {"HGNC_ID": human, "direct_ENSG_xref_rows": xref_counts[human],
             "excluded_LRG_xref_rows": excluded_counts[human],
             "gene_count": len(hgnc.get(human, set())),
             "gene_ids": ";".join(sorted(hgnc.get(human, set())))}
            for human in sorted(requested)
        ],
        ["HGNC_ID", "direct_ENSG_xref_rows", "excluded_LRG_xref_rows", "gene_count", "gene_ids"],
    )
    human_genes = set().union(*hgnc.values()) if hgnc else set()
    log("Reading complete fly and human annotation and protein references")
    representatives, protein_audit = {}, []
    for species, gtf, fasta, wanted in [
        ("fly", refs / "fly" / FLY_GTF, refs / "fly" / FLY_PEP, set(categories)),
        ("human", refs / "human/genes.gtf.bgz", refs / "human/pep.fa.bgz", human_genes),
    ]:
        proteins, canonical, coding, header = annotation(gtf, wanted, species)
        with (side / f"provenance/{species}_annotation_header.txt").open("x") as handle:
            handle.write("\n".join(header) + "\n")
        reps, audit = representative_proteins(fasta, proteins, canonical, coding, species)
        representatives[species] = reps
        protein_audit.extend(audit)
        log(f"{species}: requested genes={len(wanted)}, representative proteins={len(reps)}")
    write_rows(
        table_dir / "protein_candidates.csv", protein_audit,
        ["species", "gene_id", "protein_id", "canonical", "status", "length_aa", "selected"],
    )
    rep_columns = [
        "species", "gene_id", "protein_id", "versioned_protein_id", "fasta_gene_id",
        "transcript_id", "length_aa", "terminal_stop_removed",
        "nonstandard_residue_count", "status", "sequence_sha256",
        "selection_rule", "canonical_protein_count", "valid_protein_count",
    ]
    write_rows(
        table_dir / "representative_proteins.csv",
        [row for species in ("fly", "human") for gene, row in sorted(representatives[species].items())],
        rep_columns,
    )
    original_mapped = mapping_rows(
        categories, selected, hgnc, representatives["fly"], representatives["human"]
    )
    approved_recoveries, recovery_rows = ({}, [])
    if recovery_review:
        approved_recoveries, recovery_rows = approved_recovery_rows(
            recovery_path, categories, selected, hgnc, original_mapped,
            representatives["human"],
        )
    mapped = mapping_rows(
        categories, selected, hgnc, representatives["fly"], representatives["human"],
        approved_recoveries if recovery_review else None,
    )
    map_columns = list(mapped[0])
    write_rows(table_dir / "gene_mapping_audit.csv", mapped, map_columns)
    mapped_by_fly = {row["FlyBaseID"]: row for row in mapped}
    selected_rows = []
    for fly, pair in sorted(selected.items()):
        selected_rows.append({**pair, **mapped_by_fly[fly]})
    selected_columns = list(selected_rows[0])
    write_rows(table_dir / "selected_pairs.csv", selected_rows, selected_columns)
    write_rows(
        table_dir / "one_to_one_pairs.csv",
        [row for row in selected_rows if row["one_to_one"]], selected_columns,
    )
    all_pair_rows = []
    disease_pairs = tables["all_disease_orthologs"]
    for (fly, human), original in sorted(tables["all_orthologs"].items()):
        winner = selected.get(fly)
        is_selected = bool(winner and human == winner["HGNC_ID"])
        all_pair_rows.append({
            **original,
            "in_all_orthologs_input": True,
            "in_all_disease_orthologs_input": (fly, human) in disease_pairs,
            "category": categories.get(fly, "outside_canonical_universe"),
            "selected": is_selected,
            "selection_status": (
                "outside_canonical_universe" if fly not in categories
                else "selected" if is_selected
                else "not_highest_score_then_sorted_HGNC"
            ),
            "selected_pair_mapping_status": mapped_by_fly[fly]["status"] if fly in mapped_by_fly else "",
        })
    write_rows(table_dir / "all_pairs.csv", all_pair_rows, list(all_pair_rows[0]))
    distributions = []
    for source, pairs in tables.items():
        counts = Counter(int(row["DIOPT_score"]) for row in pairs.values())
        distributions.extend(
            {"source": source, "category": "all_input_rows", "score": score, "count": n}
            for score, n in sorted(counts.items())
        )
    for category in CATEGORIES:
        counts = Counter(
            row["DIOPT_score"] for row in selected_rows if row["category"] == category
        )
        distributions.extend(
            {"source": "selected_pairs", "category": category, "score": score, "count": n}
            for score, n in sorted(counts.items())
        )
    write_rows(
        table_dir / "diopt_distributions.csv", distributions,
        ["source", "category", "score", "count"],
    )
    attrition = []
    for category in CATEGORIES:
        members = [row for row in mapped if row["category"] == category]
        chosen = [row for row in selected_rows if row["category"] == category]
        row = {
            "category": category, "eligible_genes": len(members),
            "selected_pairs": len(chosen),
            "no_supplied_pair": sum(not row["HGNC_ID"] for row in members),
            "HGNC_unmapped": sum(row["human_gene_candidate_count"] == 0 for row in chosen),
            "HGNC_ambiguous": sum(row["human_gene_candidate_count"] > 1 for row in chosen),
            "missing_valid_fly_protein": sum(not row["fly_protein_id"] for row in members),
            "unique_HGNC_missing_valid_human_protein": sum(
                row["human_gene_candidate_count"] == 1 and not row["human_protein_id"]
                for row in chosen
            ),
            "ready_for_alignment": sum(row["status"] == "ready_for_alignment" for row in members),
            "excluded_before_alignment": sum(row["status"] != "ready_for_alignment" for row in members),
            "top_score_ties": sum(row["top_score_tie_count"] > 1 for row in chosen),
            "one_to_one_selected": sum(row["one_to_one"] for row in chosen),
            "one_to_one_ready": sum(
                row["one_to_one"] and row["status"] == "ready_for_alignment" for row in chosen
            ),
        }
        attrition.append(row)
    if recovery_review:
        expected_recovery_counts = {
            "SSDOs": (558, 104),
            "Disease orthologues excluding SSDOs": (4498, 796),
            "Orthologues excluding disease genes": (4007, 1551),
        }
        observed_recovery_counts = {
            row["category"]: (row["ready_for_alignment"], row["one_to_one_ready"])
            for row in attrition
        }
        if observed_recovery_counts != expected_recovery_counts:
            raise ValueError(
                "Recovered mapping counts differ from the approved projected universes: "
                f"{observed_recovery_counts}"
            )
    write_rows(table_dir / "category_attrition.csv", attrition, list(attrition[0]))
    if recovery_review:
        recovery_columns = [
            "FlyBaseID", "category", "original_reason", "selected_HGNC_ID",
            "selected_DIOPT_score", "original_candidate_ENSGs", "original_one_to_one",
            "decision", "proposed_HGNC_ID", "proposed_ENSG_ID",
            "proposed_protein_ID", "residual_uncertainty",
        ]
        recovery_audit = [
            {
                **{column: row[column] for column in recovery_columns},
                "applied": row["decision"] == "propose_recovery_authoritative_identity",
            }
            for row in recovery_rows
        ]
        write_rows(
            table_dir / "recovery_mapping_audit.csv", recovery_audit,
            [*recovery_columns, "applied"],
        )
    for item in inputs:
        if sha256(ROOT / item["path"]) != item["sha256"]:
            raise RuntimeError(f"Input changed during execution: {item['path']}")
    log("Preparation complete; no alignment, conservation metrics, or statistical tests were run")
    for row in attrition:
        log(json.dumps(row, sort_keys=True))
    completion = {
        "completed_utc": now(), "category_attrition": attrition,
    }
    if recovery_review:
        completion["recovery_proposals_applied"] = len(approved_recoveries)
    output_files = [
            {"path": str(path.relative_to(side)), "bytes": path.stat().st_size,
             "sha256": sha256(path)}
            for path in sorted(side.rglob("*")) if path.is_file()
            and not path.is_relative_to(side / "references")
        ]
    write_rows(
            side / "provenance/output_manifest.csv", output_files,
            ["path", "bytes", "sha256"],
        )
    completion["output_manifest_sha256"] = sha256(
            side / "provenance/output_manifest.csv"
        )
    with (side / "provenance/PREPARATION_COMPLETE.json").open("x") as handle:
        json.dump(completion, handle, indent=2)
        handle.write("\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--reuse-verified-references", action="store_true",
        help="Reverify existing downloads from --reference-source instead of downloading.",
    )
    parser.add_argument(
        "--apply-reviewed-recoveries", type=Path,
        help="Apply only the pinned 628 proposals from the completed per-gene review.",
    )
    parser.add_argument("--run-id", help="Unique local run identifier (default: UTC timestamp).")
    parser.add_argument("--reference-source", type=Path,
                        help="Directory containing references/ and provenance/reference_manifest.tsv.")
    arguments = parser.parse_args()
    main(
        reuse_verified_references=arguments.reuse_verified_references,
        recovery_review=arguments.apply_reviewed_recoveries,
        run_id=arguments.run_id, reference_source=arguments.reference_source,
    )
