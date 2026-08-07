from __future__ import annotations

import json
from pathlib import Path
import sys
import tempfile
import unittest

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "reanalysis/scripts"))

from common import base_accession, read_go_annotations  # noqa: E402
from audit_target_sources import lactyl_mask  # noqa: E402
from extractors import matching_accession_mappings, site_position_pairs  # noqa: E402
from run_pipeline import CATEGORIES, REGION_ORDER, attach_go, cell_type_statistics  # noqa: E402


class PipelineUnitTests(unittest.TestCase):
    def test_base_accession_removes_database_and_isoform_syntax(self) -> None:
        self.assertEqual(base_accession("sp|P49959-2|MRE11_HUMAN"), "P49959")
        self.assertEqual(base_accession("REV__CON__Q9H9Q4-3"), "Q9H9Q4")

    def test_go_reader_excludes_not_annotations(self) -> None:
        columns = [
            "GENE PRODUCT DB",
            "GENE PRODUCT ID",
            "SYMBOL",
            "QUALIFIER",
            "GO TERM",
            "GO NAME",
            "ECO ID",
            "GO EVIDENCE CODE",
            "REFERENCE",
            "WITH/FROM",
            "TAXON ID",
            "ASSIGNED BY",
            "ANNOTATION EXTENSION",
            "GO ASPECT",
        ]
        rows = [
            ["UniProtKB", "P49959", "MRE11", "NOT|involved_in", "GO:0006281", "DNA repair", "", "IDA", "PMID:1", "", "9606", "GOA", "", "P"],
            ["UniProtKB", "Q9H9Q4", "NHEJ1", "involved_in", "GO:0006303", "nonhomologous end joining", "", "IDA", "PMID:2", "", "9606", "GOA", "", "P"],
        ]
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "go.tsv"
            pd.DataFrame(rows, columns=columns).to_csv(path, sep="\t", index=False)
            retained, raw = read_go_annotations(path)
        self.assertEqual(set(retained["BaseAccession"]), {"Q9H9Q4"})
        self.assertEqual(int(raw["ExcludedNOT"].sum()), 1)

    def test_peaks_mapping_rejects_accession_mismatch(self) -> None:
        matches = [("P49959", 100), ("Q9H9Q4", 200)]
        self.assertEqual(matching_accession_mappings(matches, "P49959|MRE11_HUMAN"), [("P49959", 100)])
        self.assertEqual(matching_accession_mappings(matches, "O60934|NBN_HUMAN"), [])

    def test_site_position_pairs_reject_length_mismatch(self) -> None:
        self.assertEqual(site_position_pairs(["1", "2"], ["10", "20"]), [("1", "10"), ("2", "20")])
        self.assertEqual(site_position_pairs(["1", "2"], ["10"]), [])

    def test_go_gene_fallback_remains_human_only(self) -> None:
        proteins = pd.DataFrame(
            [{"BaseAccession": "UNKNOWN", "GeneSymbol": "SHARED", "ProteinName": "", "KlaSites": "K1"}]
        )
        go_summary = pd.DataFrame(columns=["BaseAccession", "GOSymbol", "GOTerms", "GONames", "GOEvidenceCodes", "GOReferences", "GOAnnotationCount"])
        go_raw = pd.DataFrame(
            [{
                "BaseAccession": "MOUSE1", "SYMBOL": "SHARED", "QUALIFIER": "involved_in",
                "GO TERM": "GO:0006281", "GO NAME": "DNA repair", "GO EVIDENCE CODE": "IDA",
                "REFERENCE": "PMID:1", "TAXON ID": "10090", "ExcludedNOT": False,
            }]
        )
        result = attach_go(proteins, go_summary, go_raw)
        self.assertEqual(result.at[0, "GOMatchMode"], "unmatched")

    def test_target_audit_recognizes_common_kla_encodings(self) -> None:
        table = pd.DataFrame(
            {
                "PTM": ["Lac", "", ""],
                "Modified peptide": ["PEPTIDE", "AK(La (K))R", "PEPTIDE"],
                "La (K)": ["0", "0", "1"],
            }
        )
        self.assertEqual(lactyl_mask(table).tolist(), [True, True, True])


class PipelineOutputTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.tables = PROJECT_ROOT / "reanalysis/results/tables"
        cls.intermediate = PROJECT_ROOT / "reanalysis/intermediate"
        if not (cls.tables / "dataset_analysis_summary.csv").exists():
            raise unittest.SkipTest("Run reanalysis/scripts/run_pipeline.py first")

    def test_excluded_datasets_never_enter_primary_evidence(self) -> None:
        evidence = pd.read_csv(
            self.intermediate / "kla_by_dataset/all_primary_sample_level_kla_sites.csv",
            dtype=str,
        )
        self.assertFalse(set(evidence["PXD"]) & {"PXD038880", "PXD050906"})
        self.assertEqual(
            set(evidence["PXD"]),
            {"PXD014870", "PXD028488", "PXD050470", "PXD053474", "PXD060185", "PXD078013", "PXD078736"},
        )

    def test_grouping_uses_mcf10a_as_immortalized_and_mcf7_as_tumor(self) -> None:
        grouping = pd.read_csv(PROJECT_ROOT / "reanalysis/config/grouping_schemes.csv", dtype=str).set_index("CellType")
        for column in ("teacher_requested_grouping", "biologically_conventional_grouping"):
            self.assertEqual(grouping.at["MCF10A", column], "normal_immortalized_cell_lines")
            self.assertEqual(grouping.at["MCF7", column], "tumor_cell_lines")

    def test_venn_counts_equal_exact_region_tables(self) -> None:
        for scheme in ("teacher_requested_grouping", "biologically_conventional_grouping"):
            for analysis in ("all_kla", "kla_go_ddr"):
                root = self.tables / "venn_regions" / scheme / analysis
                counts = pd.read_csv(root / "venn_region_counts.csv").set_index("Region")["ProteinCount"]
                membership = pd.read_csv(root / "venn_membership.csv", dtype=str)
                self.assertFalse(membership["BaseAccession"].duplicated().any())
                for region in REGION_ORDER:
                    frame = pd.read_csv(root / f"{region}.csv", dtype=str)
                    self.assertEqual(len(frame), int(counts[region]))
                    self.assertFalse(frame["BaseAccession"].duplicated().any())
                expected_total = 3112 if analysis == "all_kla" else 275
                self.assertEqual(int(counts.sum()), expected_total)

    def test_category_sets_reconstruct_from_regions(self) -> None:
        root = self.tables / "venn_regions/teacher_requested_grouping/all_kla"
        membership = pd.read_csv(root / "venn_membership.csv", dtype=str)
        self.assertEqual(membership["BaseAccession"].nunique(), 3112)
        for category in CATEGORIES:
            category_table = pd.read_csv(root / f"{category}_all.csv", dtype=str)
            self.assertFalse(category_table["BaseAccession"].duplicated().any())

    def test_tumor_specific_exports_match_tumor_only_regions(self) -> None:
        root = self.tables / "venn_regions/teacher_requested_grouping"
        pairs = [
            ("all_kla/tumor_only.csv", "tumor_specific_kla_proteins.csv"),
            ("kla_go_ddr/tumor_only.csv", "tumor_specific_kla_ddr_proteins.csv"),
        ]
        for region_path, export_path in pairs:
            region = set(pd.read_csv(root / region_path, dtype=str)["BaseAccession"])
            export = set(pd.read_csv(self.tables / export_path, dtype=str)["BaseAccession"])
            self.assertEqual(region, export)

    def test_three_group_combined_tables_are_unique_and_complete(self) -> None:
        expected = {
            "all_kla_three_groups_combined.csv": 3112,
            "kla_go_ddr_three_groups_combined.csv": 275,
        }
        flags = [
            "InHippocampusTissue",
            "InNormalImmortalizedCellLines",
            "InTumorCellLines",
        ]
        for filename, expected_rows in expected.items():
            frame = pd.read_csv(self.tables / filename, dtype=str).fillna("")
            self.assertEqual(len(frame), expected_rows)
            self.assertFalse(frame["BaseAccession"].duplicated().any())
            self.assertTrue(frame["VennRegion"].isin(REGION_ORDER).all())
            calculated = frame[flags].eq("Yes").sum(axis=1).astype(str)
            self.assertTrue(calculated.eq(frame["DetectedGroupCount"]).all())

    def test_non_deduplicated_combined_tables_preserve_group_rows(self) -> None:
        analyses = {
            "all_kla": "all_kla_three_groups_combined_non_deduplicated.csv",
            "kla_go_ddr": "kla_go_ddr_three_groups_combined_non_deduplicated.csv",
        }
        for analysis, filename in analyses.items():
            root = self.tables / "venn_regions/teacher_requested_grouping" / analysis
            expected_rows = sum(
                len(pd.read_csv(root / f"{category}_all.csv", dtype=str))
                for category in CATEGORIES
            )
            frame = pd.read_csv(self.tables / filename, dtype=str).fillna("")
            self.assertEqual(len(frame), expected_rows)
            self.assertEqual(set(frame["SourceCategory"]), set(CATEGORIES))
            self.assertTrue(frame["BaseAccession"].duplicated().any())

    def test_cell_type_statistics_match_primary_evidence(self) -> None:
        statistics = pd.read_csv(self.tables / "cell_type_kla_ddr_statistics.csv")
        primary = pd.read_csv(
            self.intermediate / "kla_by_dataset/all_primary_sample_level_kla_sites.csv",
            dtype=str,
        )
        go = pd.read_csv(
            self.intermediate / "go_intersection/all_pxd_kla_go_ddr_proteins.csv",
            dtype=str,
        )
        grouping = pd.read_csv(
            PROJECT_ROOT / "reanalysis/config/grouping_schemes.csv",
            dtype=str,
        )
        expected = cell_type_statistics(primary, set(go["BaseAccession"]), grouping)
        pd.testing.assert_frame_equal(statistics, expected, check_dtype=False)
        self.assertEqual(
            statistics.set_index("CellOrTissueType").loc["MCF10A", "TotalKlaProteins"],
            947,
        )
        self.assertEqual(
            statistics.set_index("CellOrTissueType").loc["MCF10A", "KlaGoDdrProteins"],
            116,
        )

    def test_pxd053474_reconciliation_rule(self) -> None:
        comparison = pd.read_csv(self.tables / "pxd053474/search_vs_supplementary_all.csv")
        counts = comparison["ComparisonCategory"].value_counts().to_dict()
        self.assertEqual(counts, {"consistent": 1275, "search_only": 826, "supplementary_only": 3})
        self.assertEqual(int(comparison["PrimaryIncluded"].sum()), 1298)

    def test_pxd014870_sensitivity_is_subset_of_main(self) -> None:
        main = pd.read_csv(self.intermediate / "kla_by_dataset/PXD014870_sample_level_kla_sites.csv", dtype=str)
        sensitivity = pd.read_csv(self.intermediate / "kla_by_dataset/PXD014870_sensitivity_localization_0.75.csv", dtype=str)
        main_keys = set(zip(main["SourceFile"], main["SourceRow"], main["BaseAccession"], main["KlaSite"]))
        sensitivity_keys = set(zip(sensitivity["SourceFile"], sensitivity["SourceRow"], sensitivity["BaseAccession"], sensitivity["KlaSite"]))
        self.assertTrue(sensitivity_keys <= main_keys)
        self.assertTrue(pd.to_numeric(sensitivity["LocalizationProb"], errors="coerce").ge(0.75).all())

    def test_target_source_audit_distinguishes_non_kla_detection(self) -> None:
        audit = pd.read_csv(self.tables / "target_protein_source_level_audit_MRE11_XLF_NBS1.csv", dtype=str).fillna("")
        original = audit[audit["PXD"].isin({"PXD014870", "PXD028488", "PXD050470", "PXD053474"})]
        self.assertTrue(original["ExtractedPrimaryRows"].astype(int).eq(0).all())
        self.assertTrue(original["SourceKlaTargetRows"].astype(int).eq(0).all())

    def test_reference_ddr_statistics_use_selected_protein_denominators(self) -> None:
        reference = pd.read_csv(self.tables / "cell_type_reference_ddr_statistics.csv")
        comparison = pd.read_csv(
            self.tables / "cell_type_kla_vs_reference_ddr_statistics.csv"
        )
        control_information = pd.read_csv(
            self.tables / "cell_type_reference_control_information.csv"
        )
        config = json.loads(
            (PROJECT_ROOT / "reanalysis/config/reference_proteome_selection.json").read_text()
        )
        expected = {
            row["cell_type"]: int(row["reference_protein_count_main"]) for row in config
        }
        self.assertEqual(len(reference), 10)
        self.assertEqual(set(reference["CellOrTissueType"]), set(expected))
        for row in reference.itertuples(index=False):
            self.assertEqual(row.ReferenceProteinCount, expected[row.CellOrTissueType])
            self.assertLessEqual(row.ReferenceDdrProteinCount, row.ReferenceProteinCount)
            self.assertAlmostEqual(
                row.ReferenceDdrFraction,
                row.ReferenceDdrProteinCount / row.ReferenceProteinCount,
            )
        self.assertEqual(len(comparison), 10)
        self.assertEqual(len(control_information), 10)
        self.assertFalse(control_information["选择理由"].fillna("").eq("").any())
        self.assertTrue(
            control_information["年份"].astype(int).between(2015, 2026).all()
        )
        self.assertTrue(
            (
                comparison["DdrFractionPercentagePointDifference"]
                - (
                    comparison["KlaGoDdrFraction"]
                    - comparison["ReferenceDdrFraction"]
                )
            )
            .abs()
            .lt(1e-12)
            .all()
        )

    def test_tall_surrogate_sensitivity_sets_are_consistent(self) -> None:
        frame = pd.read_csv(self.tables / "tall104_surrogate_ddr_sensitivity.csv")
        values = frame.set_index("TAllReferenceSet")
        self.assertEqual(values.at["TALL-1_primary", "ProteinCount"], 3383)
        self.assertEqual(values.at["Jurkat_sensitivity", "ProteinCount"], 3363)
        self.assertEqual(values.at["TALL-1_Jurkat_union", "ProteinCount"], 3835)
        self.assertEqual(values.at["TALL-1_Jurkat_intersection", "ProteinCount"], 2911)
        self.assertTrue((frame["DdrProteinCount"] <= frame["ProteinCount"]).all())


if __name__ == "__main__":
    unittest.main(verbosity=2)
