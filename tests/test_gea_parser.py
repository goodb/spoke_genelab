"""
Unit tests for GEA parser module.
"""

import os
import pytest
from pathlib import Path

# Get the path to the test data
TEST_DATA_DIR = Path(__file__).parent.parent / "E-GEOD-5305-gea"


class TestGEAParser:
    """Tests for the GEA file parser."""

    @pytest.fixture
    def experiment_dir(self):
        """Return path to test experiment directory."""
        return str(TEST_DATA_DIR)

    def test_parse_idf_file(self, experiment_dir):
        """Test parsing of IDF metadata file."""
        from scripts.gea.gea_parser import parse_idf_file

        idf_path = os.path.join(experiment_dir, "E-GEOD-5305.idf.txt")
        if not os.path.exists(idf_path):
            pytest.skip("Test data not available")

        metadata = parse_idf_file(idf_path)

        assert metadata["title"] != ""
        assert "p68" in metadata["title"].lower() or "transcription" in metadata["title"].lower()
        assert metadata["submitter_last_name"] == "Hoffman"
        assert "Eric" in metadata["submitter_first_name"]
        assert len(metadata["experimental_factors"]) > 0

    def test_parse_configuration_xml(self, experiment_dir):
        """Test parsing of configuration XML file."""
        from scripts.gea.gea_parser import parse_configuration_xml

        config_path = os.path.join(experiment_dir, "E-GEOD-5305-configuration.xml")
        if not os.path.exists(config_path):
            pytest.skip("Test data not available")

        assay_groups, contrasts = parse_configuration_xml(config_path)

        # Should have 8 assay groups (g1-g8)
        assert len(assay_groups) == 8
        assert "g1" in assay_groups
        assert "g8" in assay_groups

        # Should have 4 contrasts
        assert len(contrasts) == 4

        # Check contrast structure
        contrast_ids = {c.id for c in contrasts}
        assert "g1_g3" in contrast_ids
        assert "g2_g4" in contrast_ids

        # Check assay group has samples
        g1 = assay_groups["g1"]
        assert len(g1.samples) > 0
        assert g1.label != ""

    def test_parse_sdrf_file(self, experiment_dir):
        """Test parsing of condensed SDRF file."""
        from scripts.gea.gea_parser import parse_sdrf_file

        sdrf_path = os.path.join(experiment_dir, "E-GEOD-5305.condensed-sdrf.tsv")
        if not os.path.exists(sdrf_path):
            pytest.skip("Test data not available")

        samples_df = parse_sdrf_file(sdrf_path)

        assert not samples_df.empty
        assert "Sample" in samples_df.columns

        # Check for organism characteristic
        organism_cols = [c for c in samples_df.columns if "organism" in c.lower()]
        assert len(organism_cols) > 0

    def test_parse_analytics_file(self, experiment_dir):
        """Test parsing of analytics TSV file."""
        from scripts.gea.gea_parser import parse_analytics_file

        analytics_path = os.path.join(
            experiment_dir, "E-GEOD-5305_A-AFFY-23-analytics.tsv"
        )
        if not os.path.exists(analytics_path):
            pytest.skip("Test data not available")

        df = parse_analytics_file(analytics_path)

        assert not df.empty
        assert "gene_id" in df.columns
        assert "gene_name" in df.columns

        # Check for contrast columns
        contrast_cols = [c for c in df.columns if "g1_g3" in c]
        assert len(contrast_cols) > 0

    def test_parse_gsea_file(self, experiment_dir):
        """Test parsing of GSEA results file."""
        from scripts.gea.gea_parser import parse_gsea_file

        gsea_path = os.path.join(
            experiment_dir, "E-GEOD-5305.g1_g3.go.gsea.tsv"
        )
        if not os.path.exists(gsea_path):
            pytest.skip("Test data not available")

        df = parse_gsea_file(gsea_path)

        assert not df.empty
        # Should have p-value and effect size columns
        assert any("p_value" in c.lower() or "p-value" in c.lower() for c in df.columns)

    def test_load_gea_experiment(self, experiment_dir):
        """Test loading a complete GEA experiment."""
        from scripts.gea.gea_parser import load_gea_experiment

        if not os.path.exists(experiment_dir):
            pytest.skip("Test data not available")

        experiment = load_gea_experiment(experiment_dir)

        assert experiment.accession == "E-GEOD-5305"
        assert experiment.title != ""
        assert len(experiment.contrasts) == 4
        assert len(experiment.assay_groups) == 8
        assert "Mus musculus" in experiment.organism

    def test_find_gea_files(self, experiment_dir):
        """Test finding GEA files in a directory."""
        from scripts.gea.gea_parser import find_gea_files

        if not os.path.exists(experiment_dir):
            pytest.skip("Test data not available")

        files = find_gea_files(experiment_dir)

        assert len(files["idf"]) == 1
        assert len(files["sdrf"]) == 1
        assert len(files["config"]) == 1
        assert len(files["analytics"]) >= 1
        assert len(files["gsea_go"]) >= 1

    def test_get_gsea_enrichment_type(self):
        """Test determining enrichment type from filename."""
        from scripts.gea.gea_parser import get_gsea_enrichment_type

        assert get_gsea_enrichment_type("E-GEOD-5305.g1_g3.go.gsea.tsv") == "go"
        assert get_gsea_enrichment_type("E-GEOD-5305.g1_g3.reactome.gsea.tsv") == "reactome"
        assert get_gsea_enrichment_type("E-GEOD-5305.g1_g3.interpro.gsea.tsv") == "interpro"
        assert get_gsea_enrichment_type("other.tsv") == "unknown"

    def test_get_contrast_from_gsea_filename(self):
        """Test extracting contrast ID from GSEA filename."""
        from scripts.gea.gea_parser import get_contrast_from_gsea_filename

        assert get_contrast_from_gsea_filename("E-GEOD-5305.g1_g3.go.gsea.tsv") == "g1_g3"
        assert get_contrast_from_gsea_filename("E-GEOD-5305.g2_g4.reactome.gsea.tsv") == "g2_g4"


class TestGEAExtractors:
    """Tests for GEA extractor modules."""

    @pytest.fixture
    def experiment(self):
        """Load test experiment."""
        from scripts.gea.gea_parser import load_gea_experiment

        if not os.path.exists(TEST_DATA_DIR):
            pytest.skip("Test data not available")

        return load_gea_experiment(str(TEST_DATA_DIR))

    def test_extract_study_nodes(self, experiment):
        """Test study node extraction."""
        from scripts.gea.gea_study_extractor import extract_study_nodes

        studies = extract_study_nodes(experiment)

        assert len(studies) == 1
        assert studies.iloc[0]["identifier"] == "E-GEOD-5305"
        assert "identifier" in studies.columns
        assert "name" in studies.columns

    def test_extract_assay_nodes(self, experiment):
        """Test assay node extraction."""
        from scripts.gea.gea_assay_extractor import extract_assay_nodes

        assays = extract_assay_nodes(experiment)

        assert len(assays) == 4  # 4 contrasts
        assert "identifier" in assays.columns
        assert "contrast_id" in assays.columns

        # Check contrast IDs
        contrast_ids = set(assays["contrast_id"])
        assert "g1_g3" in contrast_ids

    def test_extract_mgene_nodes(self, experiment):
        """Test MGene node extraction."""
        from scripts.gea.gea_gene_extractor import create_mgene_nodes

        mgenes = create_mgene_nodes(experiment)

        assert not mgenes.empty
        assert "identifier" in mgenes.columns
        assert "symbol" in mgenes.columns

        # Check for Ensembl IDs
        sample_id = mgenes.iloc[0]["identifier"]
        assert sample_id.startswith("ENSMUSG")

    def test_extract_gsea_results(self, experiment):
        """Test GSEA extraction."""
        from scripts.gea.gea_gsea_extractor import extract_gsea_results

        results = extract_gsea_results(experiment, p_value_threshold=1.0)  # Get all results

        assert not results.empty
        assert "enrichment_type" in results.columns
        assert "contrast_id" in results.columns

        # Should have GO results
        go_results = results[results["enrichment_type"] == "go"]
        assert not go_results.empty


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
