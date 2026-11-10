"""Global alignments and conservation metrics for fixed, previously selected pairs."""

import csv
import hashlib
import json
import os
import re
import subprocess
import sys
from collections import Counter
from pathlib import Path

import prepare_references as preparation


ROOT = Path(__file__).resolve().parents[2]
SIDE = ROOT / "side/fly_human_conservation"
MAPPING = SIDE / "mapping_4935338"
STANDARD = frozenset("ACDEFGHIKLMNPQRSTVWY")
METHOD = {
    "mode": "global",
    "substitution_matrix": "BLOSUM62",
    "gap_open_score": -10.0,
    "gap_extend_score": -0.5,
    "score_epsilon": 1e-6,
    "gap_cost": "10 + (gap_length - 1) * 0.5, including terminal gaps",
    "tie_rule": "first optimal traceback from pinned Biopython, fly as target",
    "masking": "nonstandard residues to X for scoring only; originals retained",
    "identity": "100 * exact standard-residue matches / all alignment columns",
    "fly_coverage": "100 * residue-residue columns / full fly protein length",
    "human_coverage": "100 * residue-residue columns / full human protein length",
    "coverage_masking": "masked residues count as present, not as gaps",
}
METRIC_COLUMNS = [
    "FlyBaseID", "category", "HGNC_ID", "human_gene_id", "DIOPT_score",
    "one_to_one", "fly_protein_id", "human_protein_id",
    "fly_selection_rule", "human_selection_rule",
    "fly_length_aa", "human_length_aa", "fly_masked_residues", "human_masked_residues",
    "fly_nonstandard_counts", "human_nonstandard_counts",
    "alignment_score", "alignment_columns", "standard_identical_columns",
    "residue_paired_columns", "standard_residue_paired_columns",
    "fly_gap_columns", "human_gap_columns",
    "identity_percent", "fly_coverage_percent", "human_coverage_percent",
    "fly_sequence_sha256", "human_sequence_sha256", "alignment_path",
]


def create_aligner():
    import Bio
    import numpy
    from Bio import Align
    from Bio.Align import substitution_matrices

    if Bio.__version__ != "1.86" or numpy.__version__ != "1.26.4":
        raise RuntimeError(
            f"Expected Biopython 1.86 / NumPy 1.26.4; got {Bio.__version__} / {numpy.__version__}"
        )
    aligner = Align.PairwiseAligner()
    aligner.mode = METHOD["mode"]
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = METHOD["gap_open_score"]
    aligner.extend_gap_score = METHOD["gap_extend_score"]
    aligner.epsilon = METHOD["score_epsilon"]
    return aligner


def mask_sequence(sequence):
    if not sequence or set(sequence) - preparation.PROTEIN_ALPHABET:
        raise ValueError("Expected a nonempty, validated protein without internal stops or gaps")
    return "".join(aa if aa in STANDARD else "X" for aa in sequence)


def aligned_originals(alignment, fly, human):
    coordinates = alignment.coordinates
    if tuple(coordinates[:, 0]) != (0, 0) or tuple(coordinates[:, -1]) != (len(fly), len(human)):
        raise ValueError("Alignment is not full-length global")
    fly_chunks, human_chunks = [], []
    for column in range(coordinates.shape[1] - 1):
        f0, f1 = (int(value) for value in coordinates[0, column:column + 2])
        h0, h1 = (int(value) for value in coordinates[1, column:column + 2])
        nf, nh = f1 - f0, h1 - h0
        if nf < 0 or nh < 0 or (nf == 0 and nh == 0) or (nf and nh and nf != nh):
            raise ValueError("Unexpected pairwise traceback segment")
        length = max(nf, nh)
        fly_chunks.append(fly[f0:f1] if nf else "-" * length)
        human_chunks.append(human[h0:h1] if nh else "-" * length)
    return "".join(fly_chunks), "".join(human_chunks)


def alignment_metrics(fly_aligned, human_aligned):
    if not fly_aligned or len(fly_aligned) != len(human_aligned):
        raise ValueError("Alignment rows must have the same positive length")
    if (set(fly_aligned) | set(human_aligned)) - (preparation.PROTEIN_ALPHABET | {"-"}):
        raise ValueError("Invalid aligned residue")
    if any(fly == human == "-" for fly, human in zip(fly_aligned, human_aligned)):
        raise ValueError("Double-gap column")
    fly_length = len(fly_aligned.replace("-", ""))
    human_length = len(human_aligned.replace("-", ""))
    if not fly_length or not human_length:
        raise ValueError("Both proteins must contain residues")
    columns = len(fly_aligned)
    paired = sum(fly != "-" and human != "-" for fly, human in zip(fly_aligned, human_aligned))
    standard_paired = sum(
        fly in STANDARD and human in STANDARD
        for fly, human in zip(fly_aligned, human_aligned)
    )
    identical = sum(
        fly == human and fly in STANDARD
        for fly, human in zip(fly_aligned, human_aligned)
    )
    result = {
        "alignment_columns": columns,
        "standard_identical_columns": identical,
        "residue_paired_columns": paired,
        "standard_residue_paired_columns": standard_paired,
        "fly_gap_columns": fly_aligned.count("-"),
        "human_gap_columns": human_aligned.count("-"),
        "identity_percent": 100.0 * identical / columns,
        "fly_coverage_percent": 100.0 * paired / fly_length,
        "human_coverage_percent": 100.0 * paired / human_length,
    }
    if columns != fly_length + human_length - paired:
        raise ValueError("Alignment denominator invariant failed")
    return result


def load_sequences(path, representatives, species):
    wanted = {
        row["protein_id"]: row for row in representatives if row["species"] == species
    }
    if len(wanted) != len(representatives):
        raise ValueError(f"{species}: duplicate or wrong-species representatives")
    found = {}
    for header, raw in preparation.fasta_records(path):
        versioned = header.split()[0]
        protein = preparation.stable(versioned)
        if protein not in wanted:
            continue
        if protein in found:
            raise ValueError(f"Duplicate selected FASTA protein: {protein}")
        row = wanted[protein]
        metadata = dict(
            part.split(":", 1)
            for part in header.split(" description:", 1)[0].split()[1:]
            if ":" in part
        )
        sequence = raw[:-1] if raw.endswith("*") else raw
        if (
            versioned != row["versioned_protein_id"]
            or preparation.stable(metadata.get("gene", "")) != row["gene_id"]
            or preparation.stable(metadata.get("transcript", "")) != row["transcript_id"]
            or len(sequence) != int(row["length_aa"])
            or hashlib.sha256(sequence.encode("ascii")).hexdigest() != row["sequence_sha256"]
            or str(raw.endswith("*")) != row["terminal_stop_removed"]
            or sum(aa not in STANDARD for aa in sequence) != int(row["nonstandard_residue_count"])
        ):
            raise ValueError(f"Selected protein differs from successful preparation: {protein}")
        mask_sequence(sequence)
        found[protein] = sequence
    if found.keys() != wanted.keys():
        raise ValueError(f"Selected proteins missing from FASTA: {sorted(wanted.keys() - found.keys())}")
    return found


def require_unique(rows, field):
    values = [row[field] for row in rows]
    if len(set(values)) != len(values):
        raise ValueError(f"Duplicate {field} in input table")


def load_pairs(mapping):
    audit = preparation.read_rows(
        mapping / "tables/gene_mapping_audit.csv",
        ["FlyBaseID", "category", "status", "fly_protein_id", "human_protein_id"],
    )
    selected = preparation.read_rows(
        mapping / "tables/selected_pairs.csv",
        ["FlyBaseID", "category", "HGNC_ID", "human_gene_id", "one_to_one", "status"],
    )
    reps = preparation.read_rows(
        mapping / "tables/representative_proteins.csv",
        ["species", "gene_id", "protein_id", "sequence_sha256", "selection_rule", "length_aa"],
    )
    require_unique(audit, "FlyBaseID")
    require_unique(selected, "FlyBaseID")
    repmap = {(row["species"], row["gene_id"]): row for row in reps}
    if len(repmap) != len(reps):
        raise ValueError("Duplicate species/gene representatives")
    auditmap = {row["FlyBaseID"]: row for row in audit}
    pairs = []
    for row in selected:
        if row["one_to_one"] not in {"True", "False"} or row["category"] not in preparation.CATEGORIES:
            raise ValueError(f"Invalid sensitivity/category annotation: {row}")
        original = auditmap.get(row["FlyBaseID"])
        if original is None or any(row[key] != value for key, value in original.items()):
            raise ValueError(f"Selected pair disagrees with gene audit: {row['FlyBaseID']}")
        if row["status"] != "ready_for_alignment":
            continue
        if row["human_gene_candidate_count"] != "1":
            raise ValueError(f"Nonunique human mapping marked ready: {row['FlyBaseID']}")
        for species, gene, protein in (
            ("fly", row["FlyBaseID"], row["fly_protein_id"]),
            ("human", row["human_gene_id"], row["human_protein_id"]),
        ):
            if (species, gene) not in repmap or repmap[species, gene]["protein_id"] != protein:
                raise ValueError(f"Representative mismatch: {species} {gene}")
        pairs.append(row)
    if {r["FlyBaseID"] for r in pairs} != {
        r["FlyBaseID"] for r in audit if r["status"] == "ready_for_alignment"
    }:
        raise ValueError("Selected/alignment-ready gene universes disagree")
    complete = json.loads((mapping / "provenance/PREPARATION_COMPLETE.json").read_text())
    for group in complete["category_attrition"]:
        members = [r for r in pairs if r["category"] == group["category"]]
        if len(members) != group["ready_for_alignment"] or sum(
            r["one_to_one"] == "True" for r in members
        ) != group["one_to_one_ready"]:
            raise ValueError(f"Completion-marker universe mismatch: {group['category']}")
    if not pairs:
        raise ValueError("No alignment-ready pairs")
    return sorted(pairs, key=lambda row: row["FlyBaseID"]), repmap, audit


def validated_inputs():
    paths = [
        MAPPING / "tables" / name
        for name in ("selected_pairs.csv", "representative_proteins.csv", "gene_mapping_audit.csv")
    ] + [
        MAPPING / "provenance" / name
        for name in ("PREPARATION_COMPLETE.json", "input_manifest.csv", "reference_manifest.tsv", "execution.json")
    ] + [
        SIDE / name for name in (
            "align_proteins.py", "align_proteins.sh", "alignment_requirements.txt",
            "test_align_proteins.py", "prepare_references.py",
        )
    ]
    for row in preparation.read_rows(
        MAPPING / "provenance/input_manifest.csv", ["path", "bytes", "sha256"]
    ):
        path = ROOT / row["path"]
        if path.stat().st_size != int(row["bytes"]) or preparation.sha256(path) != row["sha256"]:
            raise ValueError(f"Canonical/orthology input changed since preparation: {path}")
        paths.append(path)
    sources = preparation.read_rows(
        MAPPING / "provenance/reference_manifest.tsv",
        ["path", "bytes", "sha256", "verified", "provider_checksum"],
    )
    sequences = {}
    for species, relative in (
        ("fly", "references/fly/" + preparation.FLY_PEP),
        ("human", "references/human/pep.fa.bgz"),
    ):
        entries = [row for row in sources if row["path"] == relative]
        if len(entries) != 1:
            raise ValueError(f"Missing/duplicate reference manifest entry: {relative}")
        row = entries[0]
        path = SIDE / relative
        if (
            row["verified"] != "True"
            or row["provider_checksum"] != preparation.EXPECTED[relative.removeprefix("references/")]
            or path.stat().st_size != int(row["bytes"])
            or preparation.sha256(path) != row["sha256"]
        ):
            raise ValueError(f"Peptide reference differs from verified preparation: {relative}")
        paths.append(path)
        sequences[species] = path
    manifest = [
        {"path": str(path.relative_to(ROOT)), "bytes": path.stat().st_size,
         "sha256": preparation.sha256(path)}
        for path in paths
    ]
    return manifest, sequences


def save_alignment(path, pair, fly, human, fly_aligned, human_aligned, alignment):
    with path.open("x") as handle:
        json.dump({
            "FlyBaseID": pair["FlyBaseID"], "HGNC_ID": pair["HGNC_ID"],
            "fly_protein_id": pair["fly_protein_id"], "human_protein_id": pair["human_protein_id"],
            "fly_original": fly, "human_original": human,
            "fly_scoring": mask_sequence(fly), "human_scoring": mask_sequence(human),
            "fly_aligned_original": fly_aligned, "human_aligned_original": human_aligned,
            "coordinates": alignment.coordinates.tolist(), "score": float(alignment.score),
        }, handle)
        handle.write("\n")


def main():
    job_id = os.environ.get("SLURM_JOB_ID", "")
    if not re.fullmatch(r"[0-9]+", job_id):
        raise RuntimeError("Run biological alignments through SLURM")
    output = SIDE / f"alignment_{job_id}"
    if output.exists():
        raise FileExistsError(f"Refusing to overwrite {output}")
    aligner = create_aligner()
    manifest, fastas = validated_inputs()
    pairs, representatives, audit = load_pairs(MAPPING)
    output.mkdir()
    for name in ("alignments", "tables", "provenance"):
        (output / name).mkdir()
    preparation.write_rows(
        output / "provenance/input_manifest.csv", manifest, ["path", "bytes", "sha256"]
    )
    import Bio
    import numpy
    execution = {
        "started_utc": preparation.now(), "slurm_job_id": job_id,
        "python": sys.version, "biopython": Bio.__version__, "numpy": numpy.__version__,
        "python_executable": sys.executable,
        "mapping_directory": str(MAPPING.relative_to(ROOT)), "method": METHOD,
        "algorithm": aligner.algorithm, "aligner_parameters": str(aligner),
        "source_documentation": "https://biopython.org/docs/1.86/Tutorial/chapter_pairwise.html",
        "git_commit": subprocess.run(
            ["git", "-C", str(ROOT), "rev-parse", "HEAD"],
            check=True, capture_output=True, text=True,
        ).stdout.strip(),
        "packages": subprocess.run(
            [sys.executable, "-m", "pip", "--isolated", "freeze"],
            check=True, capture_output=True, text=True,
        ).stdout.splitlines(),
    }
    with (output / "provenance/execution.json").open("x") as handle:
        json.dump(execution, handle, indent=2)
        handle.write("\n")
    with (output / "provenance/BLOSUM62.txt").open("x") as handle:
        handle.write(str(aligner.substitution_matrix))
    sequences = {}
    for species in ("fly", "human"):
        genes = {
            row["FlyBaseID"] if species == "fly" else row["human_gene_id"]
            for row in pairs
        }
        sequences[species] = load_sequences(
            fastas[species], [representatives[species, gene] for gene in sorted(genes)], species
        )
    preparation.log(f"Loaded fixed representatives for {len(pairs)} pairs; starting {aligner.algorithm}")
    results = []
    with (output / "tables/gene_conservation.csv").open("x", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=METRIC_COLUMNS)
        writer.writeheader()
        for index, pair in enumerate(pairs, start=1):
            fly = sequences["fly"][pair["fly_protein_id"]]
            human = sequences["human"][pair["human_protein_id"]]
            # Take only the first traceback; the number of optimal paths can overflow.
            alignment = next(iter(aligner.align(mask_sequence(fly), mask_sequence(human))))
            fly_aligned, human_aligned = aligned_originals(alignment, fly, human)
            if fly_aligned.replace("-", "") != fly or human_aligned.replace("-", "") != human:
                raise ValueError(f"Traceback lost residues: {pair['FlyBaseID']}")
            metrics = alignment_metrics(fly_aligned, human_aligned)
            path = output / "alignments" / f"{pair['FlyBaseID']}.json"
            save_alignment(path, pair, fly, human, fly_aligned, human_aligned, alignment)
            row = {
                key: pair[key] for key in (
                    "FlyBaseID", "category", "HGNC_ID", "human_gene_id", "DIOPT_score",
                    "one_to_one", "fly_protein_id", "human_protein_id",
                )
            }
            row.update(metrics)
            row.update({
                "fly_selection_rule": representatives["fly", pair["FlyBaseID"]]["selection_rule"],
                "human_selection_rule": representatives["human", pair["human_gene_id"]]["selection_rule"],
                "fly_length_aa": len(fly), "human_length_aa": len(human),
                "fly_masked_residues": sum(aa not in STANDARD for aa in fly),
                "human_masked_residues": sum(aa not in STANDARD for aa in human),
                "fly_nonstandard_counts": json.dumps(dict(sorted(Counter(aa for aa in fly if aa not in STANDARD).items()))),
                "human_nonstandard_counts": json.dumps(dict(sorted(Counter(aa for aa in human if aa not in STANDARD).items()))),
                "alignment_score": float(alignment.score),
                "fly_sequence_sha256": representatives["fly", pair["FlyBaseID"]]["sequence_sha256"],
                "human_sequence_sha256": representatives["human", pair["human_gene_id"]]["sequence_sha256"],
                "alignment_path": str(path.relative_to(output)),
            })
            writer.writerow(row)
            results.append(row)
            if index % 100 == 0 or index == len(pairs):
                handle.flush()
                preparation.log(f"Aligned {index}/{len(pairs)} fixed pairs")
            del alignment
    preparation.write_rows(
        output / "tables/one_to_one_conservation.csv",
        [row for row in results if row["one_to_one"] == "True"], METRIC_COLUMNS,
    )
    category_counts = []
    for category in preparation.CATEGORIES:
        eligible = [row for row in audit if row["category"] == category]
        measured = [row for row in results if row["category"] == category]
        category_counts.append({
            "category": category, "eligible_genes": len(eligible),
            "excluded_before_alignment": sum(row["status"] != "ready_for_alignment" for row in eligible),
            "aligned_pairs": len(measured),
            "identity_measured": len(measured), "fly_coverage_measured": len(measured),
            "human_coverage_measured": len(measured),
            "pairs_with_masked_residues": sum(
                row["fly_masked_residues"] > 0 or row["human_masked_residues"] > 0
                for row in measured
            ),
            "one_to_one_measured": sum(row["one_to_one"] == "True" for row in measured),
        })
    preparation.write_rows(
        output / "tables/category_counts.csv", category_counts, list(category_counts[0])
    )
    for row in manifest:
        path = ROOT / row["path"]
        if path.stat().st_size != row["bytes"] or preparation.sha256(path) != row["sha256"]:
            raise RuntimeError(f"Input/script changed during alignment: {row['path']}")
    output_files = [
        {"path": str(path.relative_to(output)), "bytes": path.stat().st_size,
         "sha256": preparation.sha256(path)}
        for path in sorted(output.rglob("*")) if path.is_file()
    ]
    preparation.write_rows(
        output / "provenance/output_manifest.csv", output_files, ["path", "bytes", "sha256"]
    )
    with (output / "provenance/ALIGNMENT_COMPLETE.json").open("x") as handle:
        json.dump({
            "completed_utc": preparation.now(), "aligned_pairs": len(results),
            "category_counts": category_counts, "method": METHOD,
            "output_manifest_sha256": preparation.sha256(output / "provenance/output_manifest.csv"),
        }, handle, indent=2)
        handle.write("\n")
    preparation.log(f"Alignment complete: {len(results)} genes; no statistical tests or plotting run")


if __name__ == "__main__":
    main()
