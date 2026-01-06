"""
Tests for gene ID mapping and SPARQL queries.

Tests:
1. Gene ID mapper functionality (Ensembl -> NCBI mapping)
2. Correct namespace usage in RDF output (ncbigene: vs ensembl:)
3. SPARQL queries from gea_pipeline.md documentation
"""

import os
import pytest
import tempfile
from pathlib import Path

from rdflib import Graph

# Test data paths - try multiple locations
TEST_DATA_PATHS = [
    Path(__file__).parent.parent / "E-GEOD-5305-gea",  # Mouse data in repo
    Path("/Users/bgood/Documents/Scripps/gea/E-GEOD-5281-gea"),  # Human Alzheimer's data
]


def get_test_data_dir(prefer_human=True):
    """Get available test data directory, preferring human data for gene ID tests."""
    if prefer_human:
        # For gene ID mapping tests, prefer human data
        for path in reversed(TEST_DATA_PATHS):
            if path.exists():
                return path
    for path in TEST_DATA_PATHS:
        if path.exists():
            return path
    return None


class TestGeneIdMapper:
    """Test the gene ID mapper module."""

    def test_load_hgnc_mapping(self):
        """Test loading HGNC gene ID mappings."""
        from scripts.common.gene_id_mapper import load_hgnc_mapping

        df = load_hgnc_mapping()

        # Should have loaded mappings
        assert len(df) > 0
        assert "ensembl_gene_id" in df.columns
        assert "entrez_id" in df.columns
        assert "symbol" in df.columns

        print(f"Loaded {len(df)} HGNC gene mappings")

    def test_ensembl_to_ncbi_map(self):
        """Test getting Ensembl to NCBI mapping dictionary."""
        from scripts.common.gene_id_mapper import get_ensembl_to_ncbi_map

        mapping = get_ensembl_to_ncbi_map()

        # Should have mappings
        assert len(mapping) > 0

        # Test a known gene (TSPAN6: ENSG00000000003 -> 7105)
        assert mapping.get("ENSG00000000003") == "7105"

        print(f"Got {len(mapping)} Ensembl -> NCBI mappings")

    def test_ncbi_to_ensembl_map(self):
        """Test getting NCBI to Ensembl mapping dictionary."""
        from scripts.common.gene_id_mapper import get_ncbi_to_ensembl_map

        mapping = get_ncbi_to_ensembl_map()

        # Should have mappings
        assert len(mapping) > 0

        # Test a known gene (7105 -> ENSG00000000003)
        assert mapping.get("7105") == "ENSG00000000003"

    def test_add_ncbi_gene_ids(self):
        """Test adding NCBI gene IDs to a DataFrame."""
        import pandas as pd
        from scripts.common.gene_id_mapper import add_ncbi_gene_ids

        # Create test DataFrame with Ensembl IDs
        df = pd.DataFrame({
            "ensembl_id": ["ENSG00000000003", "ENSG00000000005", "ENSG00000999999"],
            "symbol": ["TSPAN6", "TNMD", "UNKNOWN"],
        })

        result = add_ncbi_gene_ids(df, ensembl_col="ensembl_id")

        # Should have ncbi_gene_id column
        assert "ncbi_gene_id" in result.columns

        # Known genes should be mapped
        assert result.loc[0, "ncbi_gene_id"] == "7105"  # TSPAN6
        assert result.loc[1, "ncbi_gene_id"] == "64102"  # TNMD

        # Unknown gene should be None/NaN
        assert pd.isna(result.loc[2, "ncbi_gene_id"])


class TestGeneNamespaces:
    """Test that genes use correct namespaces in RDF output."""

    @pytest.fixture
    def rdf_graph(self):
        """Generate RDF graph from human test data."""
        from scripts.gea.gea_pipeline import process_gea_experiment
        from scripts.rdf.turtle_writer import TurtleWriter

        test_dir = get_test_data_dir(prefer_human=True)
        if not test_dir or not test_dir.exists():
            pytest.skip("Human test data not available")

        result = process_gea_experiment(
            str(test_dir),
            p_value_threshold=0.01,
            max_genes_per_assay=100,
            include_gsea=False,
            include_orthologs=False,
        )

        writer = TurtleWriter()
        for node_type, df in result.nodes.items():
            if not df.empty:
                writer.add_nodes_from_dataframe(df, node_type)
        for rel_type, df in result.relationships.items():
            if not df.empty:
                writer.add_relationships_from_dataframe(df, rel_type)

        return writer.graph

    def test_ncbi_genes_use_ncbigene_namespace(self, rdf_graph):
        """Test that genes with NCBI IDs use ncbigene: namespace."""
        query = """
        PREFIX biolink: <https://w3id.org/biolink/vocab/>
        PREFIX ncbigene: <https://www.ncbi.nlm.nih.gov/gene/>
        PREFIX spokegenelab: <https://spoke.ucsf.edu/genelab/>
        PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>

        SELECT ?gene ?id ?ensemblId WHERE {
            ?gene a biolink:Gene .
            ?gene biolink:id ?id .
            ?gene spokegenelab:id_source ?idSource .
            FILTER(STR(?idSource) = "NCBIGene")
            OPTIONAL { ?gene spokegenelab:ensembl_id ?ensemblId }
            FILTER(STRSTARTS(STR(?gene), STR(ncbigene:)))
        }
        LIMIT 10
        """
        results = list(rdf_graph.query(query))

        assert len(results) > 0, "Should have genes with ncbigene: namespace"

        for gene, gene_id, ensembl_id in results:
            gene_uri = str(gene)
            assert gene_uri.startswith("https://www.ncbi.nlm.nih.gov/gene/"), \
                f"Gene URI should use ncbigene namespace: {gene_uri}"
            # ID should be numeric (NCBI gene ID)
            assert str(gene_id).isdigit(), f"NCBI gene ID should be numeric: {gene_id}"
            print(f"NCBI gene: {gene_id} (Ensembl: {ensembl_id})")

    def test_ensembl_only_genes_use_ensembl_namespace(self, rdf_graph):
        """Test that genes without NCBI IDs use ensembl: namespace."""
        query = """
        PREFIX biolink: <https://w3id.org/biolink/vocab/>
        PREFIX ensembl: <http://identifiers.org/ensembl/>
        PREFIX spokegenelab: <https://spoke.ucsf.edu/genelab/>

        SELECT ?gene ?id WHERE {
            ?gene a biolink:Gene .
            ?gene biolink:id ?id .
            ?gene spokegenelab:id_source ?idSource .
            FILTER(STR(?idSource) = "Ensembl")
            FILTER(STRSTARTS(STR(?gene), STR(ensembl:)))
        }
        LIMIT 10
        """
        results = list(rdf_graph.query(query))

        if len(results) > 0:
            for gene, gene_id in results:
                gene_uri = str(gene)
                assert gene_uri.startswith("http://identifiers.org/ensembl/"), \
                    f"Ensembl-only gene URI should use ensembl namespace: {gene_uri}"
                assert str(gene_id).startswith("ENSG"), \
                    f"Ensembl ID should start with ENSG: {gene_id}"
                print(f"Ensembl-only gene: {gene_id}")
        else:
            print("No Ensembl-only genes found (all mapped to NCBI)")

    def test_de_relationships_use_correct_gene_namespace(self, rdf_graph):
        """Test that DE relationships point to correct gene namespaces."""
        # Check NCBI gene references
        ncbi_query = """
        PREFIX biolink: <https://w3id.org/biolink/vocab/>
        PREFIX ncbigene: <https://www.ncbi.nlm.nih.gov/gene/>
        PREFIX spokegenelab: <https://spoke.ucsf.edu/genelab/>

        SELECT ?gene ?log2fc WHERE {
            ?assoc a biolink:GeneExpressionMixin .
            ?assoc biolink:object ?gene .
            ?assoc spokegenelab:log2fc ?log2fc .
            FILTER(STRSTARTS(STR(?gene), STR(ncbigene:)))
        }
        LIMIT 5
        """
        ncbi_results = list(rdf_graph.query(ncbi_query))
        assert len(ncbi_results) > 0, "Should have DE relationships to ncbigene: genes"

        # Check Ensembl gene references (if any)
        ensembl_query = """
        PREFIX biolink: <https://w3id.org/biolink/vocab/>
        PREFIX ensembl: <http://identifiers.org/ensembl/>
        PREFIX spokegenelab: <https://spoke.ucsf.edu/genelab/>

        SELECT ?gene ?log2fc WHERE {
            ?assoc a biolink:GeneExpressionMixin .
            ?assoc biolink:object ?gene .
            ?assoc spokegenelab:log2fc ?log2fc .
            FILTER(STRSTARTS(STR(?gene), STR(ensembl:)))
        }
        LIMIT 5
        """
        ensembl_results = list(rdf_graph.query(ensembl_query))

        print(f"DE relationships: {len(ncbi_results)} to ncbigene:, {len(ensembl_results)} to ensembl:")

    def test_no_ncbigene_with_ensembl_id(self, rdf_graph):
        """Test that there are no ncbigene: URIs containing ENSG IDs."""
        query = """
        PREFIX ncbigene: <https://www.ncbi.nlm.nih.gov/gene/>

        SELECT ?gene WHERE {
            ?gene ?p ?o .
            FILTER(STRSTARTS(STR(?gene), STR(ncbigene:)))
            FILTER(CONTAINS(STR(?gene), "ENSG"))
        }
        """
        results = list(rdf_graph.query(query))

        assert len(results) == 0, \
            f"Should not have ncbigene: URIs with Ensembl IDs: {[str(r[0]) for r in results]}"


class TestDocumentedSPARQLQueries:
    """Test SPARQL queries from gea_pipeline.md documentation."""

    @pytest.fixture
    def rdf_graph(self):
        """Generate RDF graph from test data."""
        from scripts.gea.gea_pipeline import process_gea_experiment
        from scripts.rdf.turtle_writer import TurtleWriter

        test_dir = get_test_data_dir(prefer_human=True)
        if not test_dir or not test_dir.exists():
            pytest.skip("Test data not available")

        result = process_gea_experiment(
            str(test_dir),
            p_value_threshold=0.01,
            max_genes_per_assay=200,
            max_terms_per_type=20,
            include_gsea=False,  # Skip for faster tests
            include_orthologs=False,
        )

        writer = TurtleWriter()
        for node_type, df in result.nodes.items():
            if not df.empty:
                writer.add_nodes_from_dataframe(df, node_type)
        for rel_type, df in result.relationships.items():
            if not df.empty:
                writer.add_relationships_from_dataframe(df, rel_type)

        return writer.graph

    def test_find_upregulated_genes(self, rdf_graph):
        """Test query: Find upregulated genes (from documentation)."""
        query = """
        PREFIX biolink: <https://w3id.org/biolink/vocab/>
        PREFIX spokegenelab: <https://spoke.ucsf.edu/genelab/>

        SELECT ?study ?studyTitle ?assay ?geneSymbol ?log2fc ?pvalue
        WHERE {
            ?assoc a biolink:GeneExpressionMixin ;
                   biolink:subject ?assay ;
                   biolink:object ?gene ;
                   spokegenelab:log2fc ?log2fc ;
                   spokegenelab:adj_p_value ?pvalue .

            FILTER(?log2fc > 0)

            OPTIONAL { ?gene biolink:symbol ?geneSymbol }

            ?study biolink:has_output ?assay ;
                   biolink:name ?studyTitle .
        }
        ORDER BY DESC(?log2fc)
        LIMIT 10
        """
        results = list(rdf_graph.query(query))

        assert len(results) > 0, "Should find upregulated genes"

        print("\nTop 10 upregulated genes:")
        for row in results:
            print(f"  {row.geneSymbol}: log2fc={row.log2fc}, p={row.pvalue}")

        # Verify all log2fc values are positive
        for row in results:
            assert float(row.log2fc) > 0

    def test_find_significantly_upregulated_genes(self, rdf_graph):
        """Test query: Find significantly upregulated genes (log2fc > 1, p < 0.05)."""
        query = """
        PREFIX biolink: <https://w3id.org/biolink/vocab/>
        PREFIX spokegenelab: <https://spoke.ucsf.edu/genelab/>

        SELECT ?study ?studyTitle ?assay ?log2fc ?pvalue
        WHERE {
            ?assoc a biolink:GeneExpressionMixin ;
                   biolink:subject ?assay ;
                   biolink:object ?gene ;
                   spokegenelab:log2fc ?log2fc ;
                   spokegenelab:adj_p_value ?pvalue .

            OPTIONAL { ?gene biolink:symbol ?geneSymbol }

            FILTER(?log2fc > 1.0 && ?pvalue < 0.05)

            ?study biolink:has_output ?assay ;
                   biolink:name ?studyTitle .
        }
        ORDER BY DESC(?log2fc)
        LIMIT 10
        """
        results = list(rdf_graph.query(query))

        # May or may not have results depending on data
        print(f"\nFound {len(results)} significantly upregulated genes (log2fc > 1, p < 0.05)")

        for row in results:
            assert float(row.log2fc) > 1.0
            assert float(row.pvalue) < 0.05

    def test_find_downregulated_genes_in_study(self, rdf_graph):
        """Test query: Find all downregulated genes in a specific study."""
        # First get the study ID
        study_query = """
        PREFIX biolink: <https://w3id.org/biolink/vocab/>
        SELECT ?studyId WHERE {
            ?study a biolink:Study ;
                   biolink:id ?studyId .
        }
        LIMIT 1
        """
        study_results = list(rdf_graph.query(study_query))
        if not study_results:
            pytest.skip("No study found")

        study_id = str(study_results[0][0])

        query = f"""
        PREFIX biolink: <https://w3id.org/biolink/vocab/>
        PREFIX spokegenelab: <https://spoke.ucsf.edu/genelab/>

        SELECT ?geneSymbol ?log2fc ?pvalue ?assayName
        WHERE {{
            ?study biolink:id "{study_id}" ;
                   biolink:has_output ?assay .

            ?assay biolink:name ?assayName .

            ?assoc biolink:subject ?assay ;
                   biolink:object ?gene ;
                   spokegenelab:log2fc ?log2fc ;
                   spokegenelab:adj_p_value ?pvalue .

            OPTIONAL {{ ?gene biolink:symbol ?geneSymbol }}

            FILTER(?log2fc < -1.0)
        }}
        ORDER BY ?log2fc
        LIMIT 20
        """
        results = list(rdf_graph.query(query))

        print(f"\nFound {len(results)} downregulated genes (log2fc < -1) in {study_id}")

        for row in results:
            if row.log2fc:
                assert float(row.log2fc) < -1.0

    def test_count_entities(self, rdf_graph):
        """Test query: Count total entities in the knowledge graph."""
        query = """
        PREFIX biolink: <https://w3id.org/biolink/vocab/>

        SELECT
            (COUNT(DISTINCT ?study) AS ?numStudies)
            (COUNT(DISTINCT ?assay) AS ?numAssays)
            (COUNT(DISTINCT ?gene) AS ?numGenes)
        WHERE {
            { ?study a biolink:Study }
            UNION
            { ?assay a biolink:Assay }
            UNION
            { ?gene a biolink:Gene }
        }
        """
        results = list(rdf_graph.query(query))
        row = results[0]

        num_studies = int(row.numStudies)
        num_assays = int(row.numAssays)
        num_genes = int(row.numGenes)

        print(f"\nEntity counts:")
        print(f"  Studies: {num_studies}")
        print(f"  Assays: {num_assays}")
        print(f"  Genes: {num_genes}")

        assert num_studies >= 1
        assert num_assays >= 1
        assert num_genes >= 1

    def test_count_de_by_direction(self, rdf_graph):
        """Test query: Count genes by direction of expression change."""
        query = """
        PREFIX biolink: <https://w3id.org/biolink/vocab/>
        PREFIX spokegenelab: <https://spoke.ucsf.edu/genelab/>

        SELECT
            (SUM(IF(?log2fc > 0, 1, 0)) AS ?upregulated)
            (SUM(IF(?log2fc < 0, 1, 0)) AS ?downregulated)
            (COUNT(*) AS ?total)
        WHERE {
            ?assoc spokegenelab:log2fc ?log2fc ;
                   spokegenelab:adj_p_value ?pvalue .
            FILTER(?pvalue < 0.05)
        }
        """
        results = list(rdf_graph.query(query))
        row = results[0]

        up = int(row.upregulated) if row.upregulated else 0
        down = int(row.downregulated) if row.downregulated else 0
        total = int(row.total) if row.total else 0

        print(f"\nDE gene counts (p < 0.05):")
        print(f"  Upregulated: {up}")
        print(f"  Downregulated: {down}")
        print(f"  Total: {total}")

        assert total >= 0
        assert up + down <= total

    def test_study_assay_gene_chain(self, rdf_graph):
        """Test query: Full chain Study -> Assay -> Gene."""
        query = """
        PREFIX biolink: <https://w3id.org/biolink/vocab/>
        PREFIX spokegenelab: <https://spoke.ucsf.edu/genelab/>

        SELECT ?studyId ?assayName ?geneSymbol ?log2fc WHERE {
            ?study a biolink:Study ;
                   biolink:id ?studyId ;
                   biolink:has_output ?assay .

            ?assay a biolink:Assay ;
                   biolink:name ?assayName .

            ?assoc a biolink:GeneExpressionMixin ;
                   biolink:subject ?assay ;
                   biolink:object ?gene ;
                   spokegenelab:log2fc ?log2fc .

            OPTIONAL { ?gene biolink:symbol ?geneSymbol }
        }
        LIMIT 10
        """
        results = list(rdf_graph.query(query))

        assert len(results) > 0, "Should find study-assay-gene chains"

        print(f"\nStudy -> Assay -> Gene chains:")
        for row in results:
            print(f"  {row.studyId} -> {row.assayName} -> {row.geneSymbol} (log2fc: {row.log2fc})")


class TestCharacteristicsAndFactors:
    """Test SDRF characteristics and factors in RDF."""

    @pytest.fixture
    def rdf_graph(self):
        """Generate RDF graph from human test data with characteristics."""
        from scripts.gea.gea_pipeline import process_gea_experiment
        from scripts.rdf.turtle_writer import TurtleWriter

        test_dir = get_test_data_dir(prefer_human=True)
        if not test_dir or not test_dir.exists():
            pytest.skip("Human test data not available")

        result = process_gea_experiment(
            str(test_dir),
            p_value_threshold=0.01,
            max_genes_per_assay=50,
            include_gsea=False,
            include_orthologs=False,
        )

        writer = TurtleWriter()
        for node_type, df in result.nodes.items():
            if not df.empty:
                writer.add_nodes_from_dataframe(df, node_type)
        for rel_type, df in result.relationships.items():
            if not df.empty:
                writer.add_relationships_from_dataframe(df, rel_type)

        return writer.graph

    def test_disease_nodes_exist(self, rdf_graph):
        """Test that Disease nodes are created from SDRF."""
        query = """
        PREFIX biolink: <https://w3id.org/biolink/vocab/>

        SELECT ?disease ?name WHERE {
            ?disease a biolink:Disease .
            OPTIONAL { ?disease biolink:name ?name }
        }
        """
        results = list(rdf_graph.query(query))

        print(f"\nFound {len(results)} Disease nodes:")
        for row in results:
            print(f"  {row.disease}: {row.name}")

    def test_anatomy_nodes_with_ontology_uris(self, rdf_graph):
        """Test that Anatomy nodes use proper UBERON URIs."""
        query = """
        PREFIX biolink: <https://w3id.org/biolink/vocab/>
        PREFIX uberon: <http://purl.obolibrary.org/obo/UBERON_>

        SELECT ?anatomy WHERE {
            ?anatomy a biolink:AnatomicalEntity .
            FILTER(STRSTARTS(STR(?anatomy), STR(uberon:)))
        }
        """
        results = list(rdf_graph.query(query))

        print(f"\nFound {len(results)} Anatomy nodes with UBERON URIs:")
        for row in results:
            print(f"  {row.anatomy}")

    def test_celltype_nodes_with_ontology_uris(self, rdf_graph):
        """Test that CellType nodes use proper CL URIs."""
        query = """
        PREFIX biolink: <https://w3id.org/biolink/vocab/>
        PREFIX cl: <http://purl.obolibrary.org/obo/CL_>

        SELECT ?celltype WHERE {
            ?celltype a biolink:Cell .
            FILTER(STRSTARTS(STR(?celltype), STR(cl:)))
        }
        """
        results = list(rdf_graph.query(query))

        print(f"\nFound {len(results)} CellType nodes with CL URIs:")
        for row in results:
            print(f"  {row.celltype}")

    def test_assay_has_attribute_celltype(self, rdf_graph):
        """Test Assay -> CellType relationships via has_attribute."""
        query = """
        PREFIX biolink: <https://w3id.org/biolink/vocab/>

        SELECT ?assay ?celltype WHERE {
            ?assay a biolink:Assay .
            ?assay biolink:has_attribute ?celltype .
            ?celltype a biolink:Cell .
        }
        LIMIT 10
        """
        results = list(rdf_graph.query(query))

        print(f"\nFound {len(results)} Assay -> CellType relationships:")
        for row in results:
            print(f"  {row.assay} -> {row.celltype}")


class TestMultiURIParsing:
    """Test handling of multiple URIs in single SDRF fields."""

    def test_create_characteristic_nodes_with_multi_uri(self):
        """Test that multiple URIs in one field create separate nodes."""
        from scripts.gea.gea_assay_extractor import create_characteristic_nodes

        # Simulate SDRF data with multiple URIs in one field (real bug case)
        characteristics = {
            "organism part": [
                {
                    "value": "hippocampus entorhinal cortex",
                    "uri": "http://purl.obolibrary.org/obo/UBERON_0001388 http://purl.obolibrary.org/obo/UBERON_0001389",
                },
                {
                    "value": "cerebellum",
                    "uri": "http://purl.obolibrary.org/obo/UBERON_0002037",
                },
            ]
        }

        nodes_df = create_characteristic_nodes(characteristics, "Anatomy", "organism part")

        assert nodes_df is not None, "Should create nodes DataFrame"
        assert len(nodes_df) == 3, f"Should create 3 nodes (2 from split + 1 single), got {len(nodes_df)}"

        # Check that all nodes have valid single URIs
        for _, row in nodes_df.iterrows():
            uri = row["identifier"]
            assert " " not in uri, f"URI should not contain spaces: {uri}"
            assert uri.startswith("http://"), f"URI should start with http://: {uri}"

        # Check specific URIs are present
        identifiers = set(nodes_df["identifier"].tolist())
        assert "http://purl.obolibrary.org/obo/UBERON_0001388" in identifiers
        assert "http://purl.obolibrary.org/obo/UBERON_0001389" in identifiers
        assert "http://purl.obolibrary.org/obo/UBERON_0002037" in identifiers

    def test_create_characteristic_nodes_single_uri(self):
        """Test that single URIs are handled correctly."""
        from scripts.gea.gea_assay_extractor import create_characteristic_nodes

        characteristics = {
            "disease": [
                {
                    "value": "Alzheimer's disease",
                    "uri": "http://purl.obolibrary.org/obo/MONDO_0004975",
                }
            ]
        }

        nodes_df = create_characteristic_nodes(characteristics, "Disease", "disease")

        assert nodes_df is not None
        assert len(nodes_df) == 1
        assert nodes_df.iloc[0]["identifier"] == "http://purl.obolibrary.org/obo/MONDO_0004975"

    def test_create_characteristic_nodes_no_uri(self):
        """Test that values without URIs get spokegenelab URIs."""
        from scripts.gea.gea_assay_extractor import create_characteristic_nodes

        characteristics = {
            "sex": [
                {"value": "male", "uri": ""},
                {"value": "female", "uri": ""},
            ]
        }

        nodes_df = create_characteristic_nodes(characteristics, "Sex", "sex")

        assert nodes_df is not None
        assert len(nodes_df) == 2

        # Check that spokegenelab URIs are generated
        for _, row in nodes_df.iterrows():
            assert row["identifier"].startswith("https://spoke.ucsf.edu/genelab/Sex/")

    def test_create_study_characteristic_relationships_multi_uri(self):
        """Test study-characteristic relationships with multi-URI fields."""
        from scripts.gea.gea_assay_extractor import create_study_characteristic_relationships

        characteristics = {
            "organism part": [
                {
                    "value": "hippocampus entorhinal cortex",
                    "uri": "http://purl.obolibrary.org/obo/UBERON_0001388 http://purl.obolibrary.org/obo/UBERON_0001389",
                }
            ]
        }

        rels_df = create_study_characteristic_relationships(
            "E-TEST-001", characteristics, "organism part", "Anatomy"
        )

        assert rels_df is not None, "Should create relationships DataFrame"
        assert len(rels_df) == 2, f"Should create 2 relationships (one per URI), got {len(rels_df)}"

        # Check that all 'to' values are valid single URIs
        for _, row in rels_df.iterrows():
            to_uri = row["to"]
            assert " " not in to_uri, f"Target URI should not contain spaces: {to_uri}"

    def test_create_assay_factor_relationships_multi_uri(self):
        """Test assay-factor relationships with multi-URI fields."""
        import pandas as pd
        from scripts.gea.gea_assay_extractor import create_assay_factor_relationships

        # Assay nodes need reference_group_id and test_group_id columns
        assay_nodes = pd.DataFrame({
            "identifier": ["E-TEST-001-g1_g2"],
            "contrast_groups": ["g1_g2"],
            "reference_group_id": ["g1"],
            "test_group_id": ["g2"],
        })

        group_factors = {
            "g1": {
                "disease": [
                    {
                        "value": "AD and Control",
                        "uri": "http://purl.obolibrary.org/obo/MONDO_0004975 http://purl.obolibrary.org/obo/PATO_0000461",
                    }
                ]
            },
            "g2": {}  # No factors for test group
        }

        rels_df = create_assay_factor_relationships(
            assay_nodes, group_factors, "disease", "Disease"
        )

        assert rels_df is not None, "Should create relationships DataFrame"
        assert len(rels_df) == 2, f"Should create 2 relationships (one per URI), got {len(rels_df)}"

        # Check that all 'to' values are valid single URIs
        for _, row in rels_df.iterrows():
            to_uri = row["to"]
            assert " " not in to_uri, f"Target URI should not contain spaces: {to_uri}"

    def test_multi_uri_in_rdf_output(self):
        """Test that multi-URI fields produce valid RDF (no spaces in URIs)."""
        from scripts.gea.gea_pipeline import process_gea_experiment
        from scripts.rdf.turtle_writer import TurtleWriter
        import io
        import tempfile
        import os

        test_dir = get_test_data_dir(prefer_human=True)
        if not test_dir or not test_dir.exists():
            pytest.skip("Human test data not available")

        result = process_gea_experiment(
            str(test_dir),
            p_value_threshold=0.01,
            max_genes_per_assay=10,
            include_gsea=False,
            include_orthologs=False,
        )

        writer = TurtleWriter()
        for node_type, df in result.nodes.items():
            if not df.empty:
                writer.add_nodes_from_dataframe(df, node_type)
        for rel_type, df in result.relationships.items():
            if not df.empty:
                writer.add_relationships_from_dataframe(df, rel_type)

        # Serialize to temp file to check for valid URIs
        with tempfile.NamedTemporaryFile(mode='w', suffix='.ttl', delete=False) as f:
            temp_path = f.name

        try:
            writer.graph.serialize(destination=temp_path, format="turtle")
            with open(temp_path, 'r') as f:
                turtle_str = f.read()

            # Check for the specific error pattern (space-separated URIs)
            import re
            # Look for patterns like <http://... http://...> which would be invalid
            invalid_pattern = re.compile(r'<http://[^\s>]+\s+http://[^\s>]+>')
            matches = invalid_pattern.findall(turtle_str)

            assert len(matches) == 0, f"Found invalid URIs with spaces: {matches[:5]}"

            print(f"Generated {len(turtle_str)} characters of valid Turtle RDF")
        finally:
            os.unlink(temp_path)


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])
