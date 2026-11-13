"""Toy alignments and compact-input checks only; no project alignments."""

import csv
import gzip
import hashlib
import json
import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import align_proteins as analysis
import prepare_references as preparation


class AlignmentTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name).resolve()
        self.addCleanup(self.tmp.cleanup)

    def test_exact_identity_and_terminal_penalties(self):
        aligner = analysis.create_aligner()
        alignment = next(iter(aligner.align("ACDE", "ACDE")))
        fly, human = analysis.aligned_originals(alignment, "ACDE", "ACDE")
        metrics = analysis.alignment_metrics(fly, human)
        self.assertEqual(metrics["identity_percent"], 100)
        self.assertEqual(metrics["fly_coverage_percent"], 100)
        self.assertEqual(metrics["human_coverage_percent"], 100)
        self.assertEqual(metrics["standard_identical_columns"], 4)
        full = next(iter(aligner.align("WWWW", "WWWWAA")))
        fly, human = analysis.aligned_originals(full, "WWWW", "WWWWAA")
        self.assertEqual(fly.replace("-", ""), "WWWW")
        self.assertEqual(human.replace("-", ""), "WWWWAA")
        metrics = analysis.alignment_metrics(fly, human)
        self.assertEqual(metrics["alignment_columns"], 6)
        self.assertEqual(metrics["standard_identical_columns"], 4)
        self.assertAlmostEqual(metrics["identity_percent"], 100 * 4 / 6)
        self.assertEqual(metrics["fly_coverage_percent"], 100)
        self.assertAlmostEqual(metrics["human_coverage_percent"], 100 * 4 / 6)
        self.assertAlmostEqual(full.score, 44 - 10 - 0.5)

    def test_internal_indel_and_masked_identity(self):
        aligner = analysis.create_aligner()
        fly, human = "WWWWACDEWWWW", "WWWWWWWW"
        alignment = next(iter(aligner.align(fly, human)))
        f, h = analysis.aligned_originals(alignment, fly, human)
        metrics = analysis.alignment_metrics(f, h)
        self.assertEqual(metrics["standard_identical_columns"], 8)
        self.assertEqual(metrics["human_gap_columns"], 4)
        self.assertEqual(metrics["alignment_columns"], 12)
        self.assertAlmostEqual(metrics["fly_coverage_percent"], 100 * 8 / 12)
        self.assertEqual(metrics["human_coverage_percent"], 100)
        self.assertEqual(analysis.mask_sequence("ABXZJUO"), "AXXXXXX")
        alignment = next(iter(aligner.align(
            analysis.mask_sequence("WUXW"), analysis.mask_sequence("WUXW"),
        )))
        f, h = analysis.aligned_originals(alignment, "WUXW", "WUXW")
        self.assertEqual(f.replace("-", ""), "WUXW")
        metrics = analysis.alignment_metrics(f, h)
        self.assertEqual(metrics["standard_identical_columns"], 2)
        self.assertEqual(metrics["residue_paired_columns"], 4)
        self.assertEqual(metrics["standard_residue_paired_columns"], 2)
        self.assertEqual(metrics["identity_percent"], 50)
        self.assertEqual(metrics["fly_coverage_percent"], 100)

    def test_invalid_sequences_and_alignment_rows_fail(self):
        for value in ("", "AA*AA", "A-A", "acde", "AA?"):
            with self.assertRaises(ValueError):
                analysis.mask_sequence(value)
        for fly, human in (("", ""), ("AA", "A"), ("A-", "A-"), ("--", "AA"), ("A?", "AA")):
            with self.assertRaises(ValueError):
                analysis.alignment_metrics(fly, human)
        aligner = analysis.create_aligner()
        aligner.mode = "local"
        alignment = next(iter(aligner.align("AAWWWWAA", "WWWW")))
        with self.assertRaisesRegex(ValueError, "not full-length global"):
            analysis.aligned_originals(alignment, "AAWWWWAA", "WWWW")

    def test_verified_fasta_reload_and_tamper_detection(self):
        sequence = "AUXW"
        representative = {
            "species": "human", "gene_id": "ENSG1",
            "protein_id": "ENSP1", "versioned_protein_id": "ENSP1.2",
            "transcript_id": "ENST1", "length_aa": "4",
            "sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest(),
            "terminal_stop_removed": "True", "nonstandard_residue_count": "2",
        }
        path = self.root / "pep.gz"
        with gzip.open(path, "wt") as handle:
            handle.write(">ENSP1.2 gene:ENSG1.2 transcript:ENST1.3\nAUXW*\n")
        self.assertEqual(
            analysis.load_sequences(path, [representative], "human"),
            {"ENSP1": sequence},
        )
        with self.assertRaisesRegex(ValueError, "differs from successful preparation"):
            analysis.load_sequences(
                path, [{**representative, "sequence_sha256": "wrong"}], "human",
            )
        with gzip.open(path, "at") as handle:
            handle.write(">ENSP1.2 gene:ENSG1.2 transcript:ENST1.3\nAUXW*\n")
        with self.assertRaisesRegex(ValueError, "Duplicate selected FASTA"):
            analysis.load_sequences(path, [representative], "human")

    def make_mapping(self):
        mapping = self.root / "mapping"
        (mapping / "tables").mkdir(parents=True)
        (mapping / "provenance").mkdir()
        pair = {
            "FlyBaseID": "FBgn1", "category": preparation.CATEGORIES[0],
            "HGNC_ID": "HGNC:1", "human_gene_id": "ENSG1", "one_to_one": "True",
            "status": "ready_for_alignment", "fly_protein_id": "FBpp1",
            "human_protein_id": "ENSP1", "human_gene_candidate_count": "1",
            "DIOPT_score": "8",
        }
        for name in ("gene_mapping_audit", "selected_pairs"):
            preparation.write_rows(mapping / f"tables/{name}.csv", [pair], list(pair))
        reps = [
            {"species": species, "gene_id": gene, "protein_id": protein,
             "sequence_sha256": "test", "selection_rule": "canonical", "length_aa": "4"}
            for species, gene, protein in (("fly", "FBgn1", "FBpp1"), ("human", "ENSG1", "ENSP1"))
        ]
        preparation.write_rows(mapping / "tables/representative_proteins.csv", reps, list(reps[0]))
        (mapping / "provenance/PREPARATION_COMPLETE.json").write_text(json.dumps({
            "category_attrition": [
                {"category": category, "ready_for_alignment": int(index == 0),
                 "one_to_one_ready": int(index == 0)}
                for index, category in enumerate(preparation.CATEGORIES)
            ],
        }))
        return mapping

    def test_mapping_universe_and_disagreement_detection(self):
        mapping = self.make_mapping()
        pairs, reps, audit = analysis.load_pairs(mapping)
        self.assertEqual(len(pairs), 1)
        self.assertEqual(pairs[0]["FlyBaseID"], "FBgn1")
        self.assertEqual(len(audit), 1)
        self.assertEqual(reps["human", "ENSG1"]["protein_id"], "ENSP1")
        path = mapping / "tables/selected_pairs.csv"
        path.write_text(path.read_text().replace("ENSP1", "ENSP2"))
        with self.assertRaisesRegex(ValueError, "disagrees with gene audit"):
            analysis.load_pairs(mapping)

    def test_input_provenance_rejects_changed_inputs_or_references(self):
        mapping = self.make_mapping()
        (mapping / "provenance/execution.json").write_text("{}")
        for name in (
            "align_proteins.py", "align_proteins.sh", "alignment_requirements.txt",
            "test_align_proteins.py", "prepare_references.py",
        ):
            (self.root / name).write_text("test script")
        canonical = self.root / "canonical.csv"
        canonical.write_text("test input")
        preparation.write_rows(
            mapping / "provenance/input_manifest.csv",
            [{"path": "canonical.csv", "bytes": canonical.stat().st_size,
              "sha256": preparation.sha256(canonical)}], ["path", "bytes", "sha256"],
        )
        entries = []
        for relative in (
            "references/fly/" + preparation.FLY_PEP, "references/human/pep.fa.bgz",
        ):
            path = self.root / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(b"synthetic reference")
            entries.append({
                "path": relative, "bytes": path.stat().st_size,
                "sha256": preparation.sha256(path), "verified": "True",
                "provider_checksum": preparation.EXPECTED[relative.removeprefix("references/")],
            })
        preparation.write_rows(
            mapping / "provenance/reference_manifest.tsv", entries, list(entries[0]),
        )
        with (
            patch.object(analysis, "ROOT", self.root),
            patch.object(analysis, "SIDE", self.root),
            patch.object(analysis, "MAPPING", mapping),
        ):
            manifest, fastas = analysis.validated_inputs()
            self.assertEqual(len(fastas), 2)
            self.assertTrue(any(row["path"] == "canonical.csv" for row in manifest))
            canonical.write_text("changed")
            with self.assertRaisesRegex(ValueError, "Canonical/orthology input changed"):
                analysis.validated_inputs()
            canonical.write_text("test input")
            fastas["human"].write_bytes(b"changed")
            with self.assertRaisesRegex(ValueError, "Peptide reference differs"):
                analysis.validated_inputs()
            fastas["human"].write_bytes(b"synthetic reference")
            output_manifest = mapping / "provenance/output_manifest.csv"
            selected = mapping / "tables/selected_pairs.csv"
            preparation.write_rows(output_manifest, [{
                "path": "tables/selected_pairs.csv", "bytes": selected.stat().st_size,
                "sha256": preparation.sha256(selected),
            }], ["path", "bytes", "sha256"])
            completion = mapping / "provenance/PREPARATION_COMPLETE.json"
            complete = json.loads(completion.read_text())
            complete["output_manifest_sha256"] = preparation.sha256(output_manifest)
            completion.write_text(json.dumps(complete))
            (self.root / "align_recovered.sh").write_text("recovered wrapper")
            with patch.object(analysis, "MAPPING", self.root / "original_mapping"):
                manifest, _ = analysis.validated_inputs(mapping)
                self.assertTrue(any(row["path"] == "align_recovered.sh" for row in manifest))
                selected.write_text(selected.read_text() + "\n")
                with self.assertRaisesRegex(ValueError, "Recovery mapping output changed"):
                    analysis.validated_inputs(mapping)

    def test_main_exports_metrics_audit_and_completion(self):
        mapping = self.make_mapping()
        pair, reps, audit = analysis.load_pairs(mapping)
        for species, gene, sequence in (("fly", "FBgn1", "WUXW"), ("human", "ENSG1", "WUXWAA")):
            row = reps[species, gene]
            row["sequence_sha256"] = hashlib.sha256(sequence.encode()).hexdigest()
        with (
            patch.object(analysis, "ROOT", self.root),
            patch.object(analysis, "SIDE", self.root),
            patch.object(analysis, "MAPPING", mapping),
        ):
            with patch.object(analysis, "validated_inputs", return_value=([], {"fly": "fly", "human": "human"})) as validate:
                with patch.object(analysis, "load_pairs", return_value=(pair, reps, audit)) as load:
                    with patch.object(analysis, "load_sequences", side_effect=[
                        {"FBpp1": "WUXW"}, {"ENSP1": "WUXWAA"},
                    ]):
                        with patch.object(analysis.platform, "platform", return_value="test"), patch.object(analysis.subprocess, "run") as process:
                            process.return_value.stdout = "test\n"
                            analysis.main(mapping, run_id="123")
            validate.assert_called_once_with(mapping)
            load.assert_called_once_with(mapping)
            with self.assertRaises(FileExistsError):
                analysis.main(run_id="123")
        output = self.root / "alignment_local_123"
        with (output / "tables/gene_conservation.csv").open(newline="") as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual(len(rows), 1)
        row = rows[0]
        self.assertEqual(row["fly_masked_residues"], "2")
        self.assertEqual(row["standard_identical_columns"], "2")
        self.assertAlmostEqual(float(row["identity_percent"]), 100 * 2 / 6)
        self.assertAlmostEqual(float(row["human_coverage_percent"]), 100 * 4 / 6)
        aligned = json.loads((output / row["alignment_path"]).read_text())
        self.assertEqual(aligned["fly_original"], "WUXW")
        self.assertEqual(aligned["fly_scoring"], "WXXW")
        self.assertEqual(aligned["human_aligned_original"].replace("-", ""), "WUXWAA")
        self.assertTrue((output / "provenance/ALIGNMENT_COMPLETE.json").is_file())
        with (output / "tables/one_to_one_conservation.csv").open(newline="") as handle:
            self.assertEqual(len(list(csv.DictReader(handle))), 1)


if __name__ == "__main__":
    unittest.main()
