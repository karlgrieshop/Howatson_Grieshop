"""Small synthetic checks; no network access or biological analysis."""

import csv
import gzip
import hashlib
import importlib.util
import os
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch


SPEC = importlib.util.spec_from_file_location(
    "prepare", Path(__file__).with_name("prepare_references.py")
)
prepare = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(prepare)


class PreparationTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.addCleanup(self.tmp.cleanup)

    def csv(self, relative, rows, columns):
        path = self.root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        prepare.write_rows(path, rows, columns)
        return path

    def gz(self, name, text):
        path = self.root / name
        with gzip.open(path, "wt") as handle:
            handle.write(text)
        return path

    def test_pair_selection_full_graph_and_missing_genes(self):
        for name, ids in (
            ("sexspecific_disease_ids", ["FBgn1", "FBgn2"]),
            ("all_disease_ortholog_ids", ["FBgn1", "FBgn2", "FBgn3"]),
            ("all_ortholog_ids", ["FBgn1", "FBgn2", "FBgn3", "FBgn4"]),
        ):
            self.csv(
                f"Data/outputs/{name}.csv",
                [{"FlyBaseID": gene} for gene in ids], ["FlyBaseID"],
            )
        self.csv(
            "Data/outputs/dmel_gene_reference.csv",
            [{"FlyBaseID": gene, "gene_biotype": "protein_coding"}
             for gene in ["FBgn1", "FBgn2", "FBgn3", "FBgn4"]],
            ["FlyBaseID", "gene_biotype"],
        )
        pairs = [
            ("FBgn1", "HGNC:2", "8"), ("FBgn1", "HGNC:10", "8"),
            ("FBgn3", "HGNC:3", "9"), ("FBgn4", "HGNC:4", "4"),
            ("FBgn99", "HGNC:3", "7"),
        ]
        rows = [
            dict(zip(prepare.PAIR_COLUMNS, (fly, human, score, "fly", "human")))
            for fly, human, score in pairs
        ]
        self.csv("Data/inputs/all_orthologs.csv", rows, prepare.PAIR_COLUMNS)
        self.csv("Data/inputs/all_disease_orthologs.csv", rows[2:3], prepare.PAIR_COLUMNS)
        categories, tables, selected = prepare.prepare_pairs(self.root)
        self.assertEqual(len(categories), 4)
        self.assertNotIn("FBgn2", selected)
        self.assertEqual(selected["FBgn1"]["HGNC_ID"], "HGNC:10")
        self.assertEqual(selected["FBgn1"]["top_score_tie_count"], 2)
        self.assertFalse(selected["FBgn3"]["one_to_one"])
        self.assertTrue(selected["FBgn4"]["one_to_one"])
        self.assertEqual(selected["FBgn1"]["DIOPT_score"], 8)
        self.assertEqual(len(tables["all_orthologs"]), 5)

    def test_xref_exact_direct_and_ambiguous_mappings(self):
        xref = self.gz(
            "xref.gz",
            "gene_stable_id\txref_id\tdb_name\tinfo_type\n"
            "ENSG1.2\tHGNC:1\tHGNC Symbol\tDIRECT\n"
            "ENSG1.2\tHGNC:1\tHGNC Symbol\tDIRECT\n"
            "ENSG2\tHGNC:2\tHGNC Symbol\tDIRECT\n"
            "ENSG3\tHGNC:2\tHGNC Symbol\tDIRECT\n"
            "ENSG9\tHGNC:1\tGeneCards\tDIRECT\n"
            "ENSG9\tHGNC:1\tHGNC Symbol\tDEPENDENT\n"
            "LRG_10\tHGNC:1\tHGNC Symbol\tDIRECT\n"
            "LRG_11\tHGNC:3\tHGNC Symbol\tDIRECT\n",
        )
        mapping, counts, excluded = prepare.hgnc_mapping(xref, {"HGNC:1", "HGNC:2", "HGNC:3"})
        self.assertEqual(mapping["HGNC:1"], {"ENSG1"})
        self.assertEqual(counts["HGNC:1"], 2)
        self.assertEqual(mapping["HGNC:2"], {"ENSG2", "ENSG3"})
        self.assertNotIn("HGNC:3", mapping)
        self.assertEqual([row["gene_stable_id"] for row in excluded], ["LRG_10", "LRG_11"])
        self.assertTrue(all(
            row["exclusion_reason"] == "LRG_reference_namespace_not_ENSG"
            for row in excluded
        ))
        categories = {f"FBgn{n}": "SSDOs" for n in range(1, 5)}
        selected = {
            f"FBgn{n}": {"HGNC_ID": f"HGNC:{n}", "DIOPT_score": 8, "one_to_one": True}
            for n in range(1, 4)
        }
        representatives = {gene: {"protein_id": gene + "p"} for gene in categories}
        rows = prepare.mapping_rows(
            categories, selected, mapping, representatives, {"ENSG1": {"protein_id": "ENSP1"}}
        )
        self.assertEqual([row["status"] for row in rows], [
            "ready_for_alignment", "HGNC_ambiguous", "HGNC_unmapped",
            "no_supplied_orthology_pair",
        ])
        self.assertEqual(rows[1]["human_gene_id"], "")

    def test_unexpected_xref_namespace_still_fails_explicitly(self):
        xref = self.gz(
            "bad_xref.gz",
            "gene_stable_id\txref_id\tdb_name\tinfo_type\n"
            "UNKNOWN_10\tHGNC:1\tHGNC Symbol\tDIRECT\n",
        )
        with self.assertRaisesRegex(ValueError, "Unexpected gene namespace"):
            prepare.hgnc_mapping(xref, {"HGNC:1"})

    def test_canonical_longest_fallback_and_invalid_proteins(self):
        gtf_lines = ["#!genome-build BDGP6.54"]
        for gene, protein, transcript, canonical in [
            ("FBgn1", "P1", "T1", True), ("FBgn1", "P2", "T2", False),
            ("FBgn2", "P3", "T3", False), ("FBgn2", "P4", "T4", False),
            ("FBgn2", "P5", "T5", False), ("FBgn3", "P6", "T6", True),
            ("FBgn3", "P7", "T7", True), ("FBgn4", "P8", "T8", True),
            ("FBgn4", "P9", "T9", False), ("FBgn5", "P10", "T10", True),
            ("FBgn6", "P11", "T11", True),
        ]:
            attr = (
                f'gene_id "{gene}"; transcript_id "{transcript}"; '
                f'protein_id "{protein}"; gene_biotype "protein_coding";'
            )
            if canonical:
                attr += ' tag "Ensembl_canonical";'
            gtf_lines.append("3R\tFlyBase\tCDS\t1\t10\t.\t+\t0\t" + attr)
        gtf = self.gz("fly.gtf.gz", "\n".join(gtf_lines) + "\n")
        fasta = self.gz(
            "pep.gz",
            ">P1 gene:FBgn1 transcript:T1\nAAA*\n"
            ">P2 gene:FBgn1 transcript:T2\nAAAAAA\n"
            ">P3 gene:FBgn2 transcript:T3\nAAAA\n"
            ">P4 gene:FBgn2 transcript:T4\nAAAA\n"
            ">P5 gene:FBgn2 transcript:T5\nAAAA*AA\n"
            ">P6 gene:FBgn3 transcript:T6\nAAA\n"
            ">P7 gene:FBgn3 transcript:T7\nAAAAA\n"
            ">P8 gene:FBgn4 transcript:T8\nAA*AA\n"
            ">P9 gene:FBgn4 transcript:T9\nAAX\n"
            ">P10.1 gene:FBgn5 transcript:T10\nAAA\n"
            ">P10.2 gene:FBgn5 transcript:T10\nAAAA\n"
            ">P11 gene:FBgn6 transcript:WRONG\nAAA\n",
        )
        proteins, canonical, coding, header = prepare.annotation(
            gtf, {f"FBgn{n}" for n in range(1, 7)}, "fly"
        )
        representatives, audit = prepare.representative_proteins(
            fasta, proteins, canonical, coding, "fly"
        )
        self.assertEqual(representatives["FBgn1"]["protein_id"], "P1")
        self.assertEqual(representatives["FBgn1"]["length_aa"], 3)
        self.assertTrue(representatives["FBgn1"]["terminal_stop_removed"])
        self.assertEqual(representatives["FBgn2"]["protein_id"], "P3")
        self.assertEqual(representatives["FBgn3"]["protein_id"], "P7")
        self.assertEqual(representatives["FBgn4"]["protein_id"], "P9")
        self.assertEqual(representatives["FBgn4"]["nonstandard_residue_count"], 1)
        self.assertNotIn("FBgn5", representatives)
        self.assertNotIn("FBgn6", representatives)
        self.assertIn("ambiguous_fasta_records", {row["status"] for row in audit})
        self.assertIn("transcript_mapping_conflict", {row["status"] for row in audit})
        with self.assertRaisesRegex(ValueError, "build header mismatch"):
            prepare.annotation(gtf, set(), "human")

    def test_checksum_verification_and_manifest_delimiter(self):
        references = self.root / "references"
        references.mkdir()
        (self.root / "provenance").mkdir()
        data = {
            "pep.fa.bgz": b"human peptide", "genes.gtf.bgz": b"human annotation",
            "xref.tsv.gz": b"human xref", prepare.FLY_PEP: b"fly peptide",
            prepare.FLY_GTF: b"fly annotation",
        }
        md5s = {name: hashlib.md5(data[name]).hexdigest()
                for name in ("pep.fa.bgz", "genes.gtf.bgz", "xref.tsv.gz")}
        sums = {}
        for name in (prepare.FLY_PEP, prepare.FLY_GTF):
            temp = self.root / name
            temp.write_bytes(data[name])
            sums[name] = " ".join(subprocess.check_output(
                ["sum", "-r", str(temp)], text=True,
            ).split()[:2])
        pinned = {
            **{"human/" + name: checksum for name, checksum in md5s.items()},
            **{"fly/" + name: checksum for name, checksum in sums.items()},
        }

        def fake_download(url, path):
            if path.name == "md5sum.txt":
                content = "".join(f"{digest}  {name}\n" for name, digest in md5s.items()).encode()
            elif path.name.endswith("CHECKSUMS"):
                name = prepare.FLY_PEP if path.name == "pep.CHECKSUMS" else prepare.FLY_GTF
                content = f"{sums[name]} {name}\n".encode()
            else:
                content = data[path.name]
            with path.open("xb") as handle:
                handle.write(content)
            return {"url": url, "retrieved_utc": "test", "last_modified": "",
                    "bytes": len(content), "sha256": prepare.sha256(path)}

        with patch.object(prepare, "download", fake_download), patch.object(prepare, "EXPECTED", pinned):
            prepare.reference_files(self.root)
        manifest = self.root / "provenance/reference_manifest.tsv"
        with manifest.open(newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(rows), 8)
        self.assertEqual(sum(row["verified"] == "True" for row in rows), 5)
        with self.assertRaises(FileExistsError):
            prepare.reference_files(self.root)
        input_path = self.root / "input.csv"
        input_path.write_bytes(b"input")
        inputs = [{
            "path": input_path.name, "bytes": input_path.stat().st_size,
            "sha256": prepare.sha256(input_path),
        }]
        prepare.write_rows(
            self.root / "provenance/input_manifest.csv",
            inputs, ["path", "bytes", "sha256"],
        )
        output = self.root / "mapping_test"
        (output / "provenance").mkdir(parents=True)
        with patch.object(prepare, "EXPECTED", pinned), patch.object(prepare, "download") as download:
            self.assertEqual(prepare.reuse_references(self.root, output, inputs), references)
            download.assert_not_called()
        reused = prepare.read_rows(
            output / "provenance/reference_manifest.tsv",
            ["path", "sha256", "retrieved_utc", "reuse_verified_utc"],
        )
        self.assertEqual(len(reused), 8)
        self.assertTrue(all(row["retrieved_utc"] == "test" for row in reused))
        with self.assertRaisesRegex(ValueError, "Input changed"):
            prepare.reuse_references(
                self.root, output, [{**inputs[0], "sha256": "changed"}],
            )
        (references / "human/pep.fa.bgz").write_bytes(b"changed")
        with patch.object(prepare, "EXPECTED", pinned):
            with self.assertRaisesRegex(ValueError, "Reference fingerprint changed"):
                prepare.reuse_references(self.root, output, inputs)
            with self.assertRaisesRegex(ValueError, "Checksum mismatch"):
                prepare.verify_provider_checksum(references, "human/pep.fa.bgz")

    def test_refuses_existing_outputs_before_downloads(self):
        for relative in prepare.INPUTS:
            path = self.root / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            path.touch()
        side = self.root / "side/fly_human_conservation"
        (side / "tables").mkdir(parents=True)
        with patch.object(prepare, "ROOT", self.root), patch.object(prepare, "SIDE", side):
            with patch.dict(os.environ, {"SLURM_JOB_ID": "test"}):
                with patch.object(prepare, "download") as download:
                    with self.assertRaisesRegex(FileExistsError, "Refusing to overwrite"):
                        prepare.main()
                    download.assert_not_called()
        (side / "mapping_test").mkdir()
        with patch.object(prepare, "ROOT", self.root), patch.object(prepare, "SIDE", side):
            with patch.dict(os.environ, {"SLURM_JOB_ID": "test"}):
                with self.assertRaisesRegex(FileExistsError, "Refusing to overwrite"):
                    prepare.main(reuse_verified_references=True)

    def test_main_exports_ready_missing_and_ambiguous_genes(self):
        for name, ids in (
            ("sexspecific_disease_ids", ["FBgn1", "FBgn2"]),
            ("all_disease_ortholog_ids", ["FBgn1", "FBgn2", "FBgn3"]),
            ("all_ortholog_ids", ["FBgn1", "FBgn2", "FBgn3", "FBgn4"]),
        ):
            self.csv(
                f"Data/outputs/{name}.csv",
                [{"FlyBaseID": gene} for gene in ids], ["FlyBaseID"],
            )
        self.csv(
            "Data/outputs/dmel_gene_reference.csv",
            [{"FlyBaseID": f"FBgn{n}", "gene_biotype": "protein_coding"} for n in range(1, 5)],
            ["FlyBaseID", "gene_biotype"],
        )
        pairs = [
            dict(zip(prepare.PAIR_COLUMNS, (f"FBgn{n}", f"HGNC:{n}", "8", "fly", "human")))
            for n in (1, 3, 4)
        ]
        self.csv("Data/inputs/all_orthologs.csv", pairs, prepare.PAIR_COLUMNS)
        self.csv("Data/inputs/all_disease_orthologs.csv", pairs[1:2], prepare.PAIR_COLUMNS)
        side = self.root / "side/fly_human_conservation"
        side.mkdir(parents=True)
        for name in ("prepare_references.py", "prepare_references.sh"):
            (side / name).write_bytes(Path(__file__).with_name(name).read_bytes())

        def fake_references(destination):
            refs = destination / "references"
            for species, members in (("fly", (1, 2, 3, 4)), ("human", (1, 3, 33))):
                (refs / species).mkdir()
                lines = (
                    ["#!genome-build BDGP6.54"] if species == "fly" else
                    ["#!genome-build GRCh38.p14", "#!genebuild-last-updated 2026-04"]
                )
                fasta = []
                for n in members:
                    gene = f"FBgn{n}" if species == "fly" else f"ENSG{n}"
                    attr = (
                        f'gene_id "{gene}"; transcript_id "T{n}"; protein_id "P{n}"; '
                        'gene_biotype "protein_coding"; tag "Ensembl_canonical";'
                    )
                    lines.append("1\ttest\tCDS\t1\t9\t.\t+\t0\t" + attr)
                    fasta.append(f">P{n} gene:{gene} transcript:T{n}\nAAA\n")
                gtf = prepare.FLY_GTF if species == "fly" else "genes.gtf.bgz"
                pep = prepare.FLY_PEP if species == "fly" else "pep.fa.bgz"
                with gzip.open(refs / species / gtf, "wt") as handle:
                    handle.write("\n".join(lines) + "\n")
                with gzip.open(refs / species / pep, "wt") as handle:
                    handle.write("".join(fasta))
            with gzip.open(refs / "human/xref.tsv.gz", "wt") as handle:
                handle.write(
                    "gene_stable_id\txref_id\tdb_name\tinfo_type\n"
                    "ENSG1\tHGNC:1\tHGNC Symbol\tDIRECT\n"
                    "ENSG3\tHGNC:3\tHGNC Symbol\tDIRECT\n"
                    "ENSG33\tHGNC:3\tHGNC Symbol\tDIRECT\n"
                    "LRG_10\tHGNC:1\tHGNC Symbol\tDIRECT\n"
                )
            return refs

        with patch.object(prepare, "ROOT", self.root), patch.object(prepare, "SIDE", side):
            with patch.dict(os.environ, {"SLURM_JOB_ID": "test"}):
                with patch.object(prepare, "reference_files", fake_references):
                    with patch.object(prepare.subprocess, "run") as git:
                        git.return_value.stdout = "synthetic-commit\n"
                        prepare.main()
        with (side / "tables/gene_mapping_audit.csv").open(newline="") as handle:
            mapped = list(csv.DictReader(handle))
        self.assertEqual([row["status"] for row in mapped], [
            "ready_for_alignment", "no_supplied_orthology_pair", "HGNC_ambiguous", "HGNC_unmapped",
        ])
        for name, n in (
            ("all_pairs", 3), ("selected_pairs", 3), ("one_to_one_pairs", 3),
            ("gene_mapping_audit", 4), ("category_attrition", 3),
        ):
            with (side / f"tables/{name}.csv").open(newline="") as handle:
                self.assertEqual(len(list(csv.DictReader(handle))), n)
        self.assertTrue((side / "provenance/PREPARATION_COMPLETE.json").is_file())
        with (side / "tables/excluded_human_xrefs.csv").open(newline="") as handle:
            excluded = list(csv.DictReader(handle))
        self.assertEqual(excluded[0]["gene_stable_id"], "LRG_10")
        self.assertEqual(len(excluded), 1)
        original_input_hash = prepare.sha256(side / "provenance/input_manifest.csv")

        def fake_reuse(source, destination, inputs):
            self.assertEqual(source, side)
            self.assertEqual(destination, side / "mapping_recovery_test")
            self.assertEqual(len(inputs), len(prepare.INPUTS))
            return source / "references"

        with patch.object(prepare, "ROOT", self.root), patch.object(prepare, "SIDE", side):
            with patch.dict(os.environ, {"SLURM_JOB_ID": "recovery_test"}):
                with patch.object(prepare, "reuse_references", fake_reuse):
                    with patch.object(prepare, "reference_files") as downloads:
                        with patch.object(prepare.subprocess, "run") as git:
                            git.return_value.stdout = "synthetic-commit\n"
                            prepare.main(reuse_verified_references=True)
                        downloads.assert_not_called()
        self.assertTrue(
            (side / "mapping_recovery_test/provenance/PREPARATION_COMPLETE.json").is_file()
        )
        self.assertEqual(
            prepare.sha256(side / "provenance/input_manifest.csv"), original_input_hash,
        )
        self.assertFalse((side / "mapping_recovery_test/references").exists())


if __name__ == "__main__":
    unittest.main()
