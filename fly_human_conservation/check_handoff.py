"""Verify the Git-only local analysis inputs without HPC data or dependencies."""

import csv
import hashlib
import json
import math
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
GROUPS = (
    "SSDOs",
    "Disease orthologues excluding SSDOs",
    "Orthologues excluding disease genes",
)


def read(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t" if path.suffix == ".tsv" else ","))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def unique(rows, key):
    result = {row[key]: row for row in rows}
    require(len(result) == len(rows), f"Duplicate {key}")
    return result


def verify_hash(path, expected):
    require(path.is_file(), f"Missing required handoff file: {path}")
    require(path.stat().st_size == int(expected["bytes"]), f"Size mismatch: {path}")
    require(digest(path) == expected["sha256"], f"Checksum mismatch: {path}")


def main():
    alignment = HERE / "alignment_4935733"
    mapping = HERE / "mapping_4935338"
    completion = json.loads((alignment / "provenance/ALIGNMENT_COMPLETE.json").read_text())
    manifest_path = alignment / "provenance/output_manifest.csv"
    require(
        digest(manifest_path) == completion["output_manifest_sha256"],
        "Alignment output manifest checksum mismatch",
    )
    manifest = unique(read(manifest_path), "path")
    for relative, expected in manifest.items():
        if relative.startswith(("tables/", "provenance/")):
            verify_hash(alignment / relative, expected)

    # Historical input paths refer to the pre-move layout, not the local layout.
    prefix = "side/fly_human_conservation/"
    archived = HERE / "provenance/executed_scripts"
    for expected in read(alignment / "provenance/input_manifest.csv"):
        relative = expected["path"]
        if relative.startswith((prefix + "references/", "Data/inputs/")):
            continue
        if relative.startswith(prefix):
            suffix = relative.removeprefix(prefix)
            path = archived / suffix if "/" not in suffix else HERE / suffix
        else:
            path = ROOT / relative
        verify_hash(path, expected)
    for expected in read(mapping / "provenance/input_manifest.csv"):
        if not expected["path"].startswith("Data/inputs/"):
            verify_hash(ROOT / expected["path"], expected)

    execution = json.loads((mapping / "provenance/execution.json").read_text())
    for name, expected in execution["entry_points_sha256"].items():
        require(digest(archived / name) == expected, f"Original mapping source mismatch: {name}")

    audit = unique(read(mapping / "tables/gene_mapping_audit.csv"), "FlyBaseID")
    rows = read(alignment / "tables/gene_conservation.csv")
    measured = unique(rows, "FlyBaseID")
    sensitivity = unique(read(alignment / "tables/one_to_one_conservation.csv"), "FlyBaseID")
    require(
        sensitivity == {
            gene: row for gene, row in measured.items() if row["one_to_one"] == "True"
        },
        "One-to-one sensitivity subset changed",
    )
    require(
        set(measured) == {
            gene for gene, row in audit.items() if row["status"] == "ready_for_alignment"
        },
        "Measured and mapping-ready universes differ",
    )
    require(len(rows) == completion["aligned_pairs"] == 8435, "Aligned gene count changed")
    require(len(sensitivity) == 2299, "One-to-one gene count changed")

    def ids(name):
        return {row["FlyBaseID"] for row in read(ROOT / "Data/outputs" / name)}

    ssdo = ids("sexspecific_disease_ids.csv")
    disease = ids("all_disease_ortholog_ids.csv")
    orthologues = ids("all_ortholog_ids.csv")
    coding = {
        row["FlyBaseID"]
        for row in read(ROOT / "Data/outputs/dmel_gene_reference.csv")
        if row["gene_biotype"] == "protein_coding"
    }
    require(ssdo <= disease <= orthologues <= coding, "Canonical nesting changed")
    categories = {}
    for category, genes in zip(GROUPS, (ssdo, disease - ssdo, orthologues - disease)):
        categories.update({gene: category for gene in genes})
    require(set(audit) == set(categories), "Canonical eligible universe changed")
    for gene, row in audit.items():
        require(row["category"] == categories[gene], f"Canonical category changed: {gene}")

    for row in rows:
        gene = row["FlyBaseID"]
        require(row["category"] == categories[gene], f"Measured category changed: {gene}")
        for field in ("HGNC_ID", "human_gene_id", "fly_protein_id", "human_protein_id", "one_to_one"):
            require(row[field] == audit[gene][field], f"Mapping changed: {gene} {field}")
        fly, human, columns, paired, identical, standard_paired = (
            int(row[key]) for key in (
                "fly_length_aa", "human_length_aa", "alignment_columns",
                "residue_paired_columns", "standard_identical_columns",
                "standard_residue_paired_columns",
            )
        )
        require(
            min(fly, human) > 0 and columns == fly + human - paired
            and 0 <= identical <= standard_paired <= paired <= min(fly, human),
            f"Invalid metric denominator: {gene}",
        )
        for key, expected in (
            ("identity_percent", 100.0 * identical / columns),
            ("fly_coverage_percent", 100.0 * paired / fly),
            ("human_coverage_percent", 100.0 * paired / human),
        ):
            value = float(row[key])
            require(
                math.isfinite(value) and 0 <= value <= 100
                and math.isclose(value, expected, rel_tol=1e-12, abs_tol=1e-12),
                f"Invalid metric value: {gene} {key}",
            )
    counts = read(alignment / "tables/category_counts.csv")
    require(len(counts) == 3 and {row["category"] for row in counts} == set(GROUPS), "Category summary changed")
    for row in counts:
        category = row["category"]
        group = [item for item in rows if item["category"] == category]
        for key in ("aligned_pairs", "identity_measured", "fly_coverage_measured", "human_coverage_measured"):
            require(int(row[key]) == len(group), f"Measured count mismatch: {category} {key}")
        eligible = sum(item == category for item in categories.values())
        require(int(row["eligible_genes"]) == eligible, f"Eligible count mismatch: {category}")
        require(int(row["excluded_before_alignment"]) == eligible - len(group), f"Exclusion mismatch: {category}")
        require(
            int(row["one_to_one_measured"]) == sum(item["one_to_one"] == "True" for item in group),
            f"Sensitivity count mismatch: {category}",
        )
        require(
            int(row["pairs_with_masked_residues"]) == sum(
                int(item["fly_masked_residues"]) > 0 or int(item["human_masked_residues"]) > 0
                for item in group
            ),
            f"Masking count mismatch: {category}",
        )
    print("Handoff verified: 9,092 eligible, 8,435 measured, 2,299 one-to-one genes.")
    print("Compact outputs, archived sources, canonical inputs, and metric definitions agree.")
    print("No ignored raw inputs, reference downloads, alignments, or HPC environment required.")


if __name__ == "__main__":
    main()
