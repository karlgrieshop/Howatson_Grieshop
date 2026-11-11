"""Synthetic evidence checks; never run project reference scans in tests."""

import gzip
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import review_mapping_recovery as review


class RecoveryReviewTests(unittest.TestCase):
    def setUp(self):
        self.row = {
            "FlyBaseID": "FBgn1", "category": review.prep.CATEGORIES[0],
            "original_reason": "HGNC_ambiguous", "selected_HGNC_ID": "HGNC:1",
            "selected_human_symbol": "GENE", "selected_DIOPT_score": "8",
            "original_candidate_ENSGs": "ENSG1;ENSG2", "original_one_to_one": "True",
            "source_disease_records": "[]", "decision": "pending_evidence_review",
            "proposed_HGNC_ID": "", "proposed_ENSG_ID": "", "proposed_protein_ID": "",
            "evidence_source": "test", "evidence_release": review.FIXED_RELEASE,
            "evidence": "", "residual_uncertainty": "",
        }
        self.doc = {
            "hgnc_id": "HGNC:1", "symbol": "GENE", "status": "Approved",
            "ensembl_gene_id": "ENSG1", "locus_group": "protein-coding gene",
            "locus_type": "gene with protein product", "location": "1p23",
        }
        self.annotation = {
            "gene_name": "GENE", "gene_biotype": "protein_coding", "chromosome": "1",
        }

    def propose(self, **changes):
        return review.propose(
            changes.get("row", self.row), changes.get("doc", self.doc),
            changes.get("annotation", self.annotation),
            changes.get("protein", {"protein_id": "ENSP1"}),
            changes.get("direct", {"HGNC:1"}),
        )

    def test_identity_annotation_and_xref_support_a_proposal(self):
        row = self.propose()
        self.assertEqual(row["decision"], "propose_recovery_authoritative_identity")
        self.assertEqual(row["proposed_ENSG_ID"], "ENSG1")
        self.assertEqual(row["proposed_protein_ID"], "ENSP1")
        self.assertEqual(self.row["decision"], "pending_evidence_review")

    def test_protein_availability_does_not_resolve_identity(self):
        for changes in (
            {"annotation": {}},
            {"annotation": {**self.annotation, "gene_name": "PARALOGUE"}},
            {"annotation": {**self.annotation, "chromosome": "2"}},
            {"annotation": {**self.annotation, "chromosome": "KI_ALT"}},
            {"annotation": {**self.annotation, "gene_biotype": "pseudogene"}},
            {"direct": set()}, {"direct": {"HGNC:1", "HGNC:2"}},
            {"doc": {**self.doc, "ensembl_gene_id": "ENSG3"}},
        ):
            with self.subTest(changes=changes):
                self.assertEqual(self.propose(**changes)["decision"], "manual_review_required")

    def test_supported_identity_without_protein_stays_excluded(self):
        self.assertEqual(
            self.propose(protein={})["decision"], "retain_excluded_no_valid_fixed_protein",
        )
        self.assertEqual(
            self.propose(doc={**self.doc, "locus_group": "pseudogene"})["decision"],
            "retain_excluded_HGNC_nonprotein_locus",
        )
        with self.assertRaisesRegex(ValueError, "change a selected HGNC"):
            self.propose(doc={**self.doc, "hgnc_id": "HGNC:2"})

    def test_placeholder_requires_exact_partner_and_one_to_one_unchanged(self):
        row = {**self.row, "selected_HGNC_ID": "HGNC:0", "original_reason": "HGNC_unmapped"}
        proposal = self.propose(row=row)
        self.assertEqual(proposal["decision"], "propose_recovery_authoritative_identity")
        self.assertEqual(proposal["selected_HGNC_ID"], "HGNC:0")
        self.assertEqual(proposal["original_one_to_one"], "True")
        counts = review.category_counts(
            [{"FlyBaseID": "FBgn1", "category": self.row["category"],
              "status": "HGNC_unmapped", "one_to_one": "True"}],
            [proposal],
        )[0]
        self.assertEqual(counts["original_measured"], 0)
        self.assertEqual(counts["projected_mapped_if_approved"], 1)
        self.assertEqual(counts["proposed_one_to_one_additional"], 1)

    def test_full_review_with_tiny_fixed_reference(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "references/human").mkdir(parents=True)
            (root / "out/downloads").mkdir(parents=True)
            (root / "out/tables").mkdir()
            (root / "out/provenance").mkdir()
            with gzip.open(root / "references/human/genes.gtf.bgz", "wt") as handle:
                handle.write(
                    '#!genome-build GRCh38.p14\n#!genebuild-last-updated 2026-04\n'
                    '1\ttest\tgene\t1\t9\t.\t+\t.\tgene_id "ENSG1"; gene_name "GENE"; gene_biotype "protein_coding";\n'
                    '1\ttest\tCDS\t1\t9\t.\t+\t0\tgene_id "ENSG1"; gene_name "GENE"; gene_biotype "protein_coding"; '
                    'transcript_id "ENST1"; protein_id "ENSP1"; tag "Ensembl_canonical";\n'
                )
            with gzip.open(root / "references/human/pep.fa.bgz", "wt") as handle:
                handle.write(">ENSP1 gene:ENSG1 transcript:ENST1\nAAA\n")
            with gzip.open(root / "references/human/xref.tsv.gz", "wt") as handle:
                handle.write(
                    "gene_stable_id\txref_id\txref_label\tdb_name\tinfo_type\n"
                    "ENSG1\tHGNC:1\tGENE\tHGNC Symbol\tDIRECT\n"
                    "ENSG2\t98765\tOLD_LOC\tNCBI gene (formerly Entrezgene)\tDIRECT\n"
                )

            def fake_fetch(url, path):
                review.prep.write_rows(path, [self.doc], list(self.doc))
                return {"retrieved_utc": "test", "path": str(path), "url": url}

            with (
                patch.object(review, "BASE", root),
                patch.object(review, "verify_references", return_value=[
                    {"path": "references/human/genes.gtf.bgz", "sha256": "test"},
                ]),
                patch.object(review, "fetch", side_effect=fake_fetch),
            ):
                proposals, manifest = review.full_review(root / "out", [self.row])
            self.assertEqual(proposals[0]["decision"], "propose_recovery_authoritative_identity")
            self.assertEqual(proposals[0]["proposed_protein_ID"], "ENSP1")
            direct, records = review.scan_xrefs(
                root / "references/human/xref.tsv.gz", set(), {"LOC98765"},
            )
            self.assertEqual(records[0]["gene_stable_id"], "ENSG2")
            self.assertEqual(manifest[0]["retrieved_utc"], "test")


if __name__ == "__main__":
    unittest.main()
