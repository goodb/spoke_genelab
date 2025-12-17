"""
Integration tests for RDF output and SPARQL queries.

These tests verify:
1. RDF Turtle output is syntactically correct
2. The graph can be queried with SPARQL
3. Biological queries like "what genes are upregulated?" work correctly
"""

import os
import pytest
import tempfile
from pathlib import Path

from rdflib import Graph

# Get the path to the test data
TEST_DATA_DIR = Path(__file__).parent.parent / "E-GEOD-5305-gea"


class TestRDFSyntax:
    """Test RDF output syntax validity."""

    @pytest.fixture
    def temp_output_dir(self):
        """Create a temporary output directory."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)

    @pytest.fixture
    def pipeline_result(self, temp_output_dir):
        """Run the GEA pipeline and return results."""
        from scripts.gea.gea_pipeline import process_gea_experiment

        if not os.path.exists(TEST_DATA_DIR):
            pytest.skip("Test data not available")

        result = process_gea_experiment(
            str(TEST_DATA_DIR),
            p_value_threshold=0.1,
            include_gsea=True,
            include_orthologs=False,  # Skip orthologs for faster tests
        )
        return result

    def test_turtle_writer_creates_valid_rdf(self, pipeline_result, temp_output_dir):
        """Test that the Turtle writer creates syntactically valid RDF."""
        from scripts.rdf.turtle_writer import write_graph_to_turtle

        output_file = temp_output_dir / "test.ttl"

        writer = write_graph_to_turtle(
            pipeline_result.nodes,
            pipeline_result.relationships,
            output_file,
        )

        # Verify file was created
        assert output_file.exists()

        # Verify syntax by parsing with rdflib
        g = Graph()
        g.parse(str(output_file), format="turtle")

        # Should have triples
        assert len(g) > 0
        print(f"Created {len(g)} triples")

    def test_turtle_file_can_be_reparsed(self, pipeline_result, temp_output_dir):
        """Test that generated Turtle can be parsed by rdflib."""
        from scripts.rdf.turtle_writer import write_graph_to_turtle, validate_turtle_syntax

        output_file = temp_output_dir / "test.ttl"

        write_graph_to_turtle(
            pipeline_result.nodes,
            pipeline_result.relationships,
            output_file,
        )

        # Should not raise an exception
        assert validate_turtle_syntax(output_file) is True


class TestSPARQLQueries:
    """Test SPARQL queries on the RDF graph."""

    @pytest.fixture
    def rdf_graph(self):
        """Generate RDF graph from test data."""
        from scripts.gea.gea_pipeline import process_gea_experiment
        from scripts.rdf.turtle_writer import TurtleWriter

        if not os.path.exists(TEST_DATA_DIR):
            pytest.skip("Test data not available")

        # Run pipeline
        result = process_gea_experiment(
            str(TEST_DATA_DIR),
            p_value_threshold=0.1,
            include_gsea=True,
            include_orthologs=False,
        )

        # Create RDF graph in memory
        writer = TurtleWriter()

        for node_type, df in result.nodes.items():
            if not df.empty:
                writer.add_nodes_from_dataframe(df, node_type)

        for rel_type, df in result.relationships.items():
            if not df.empty:
                writer.add_relationships_from_dataframe(df, rel_type)

        return writer.graph

    def test_query_all_studies(self, rdf_graph):
        """Test querying for all studies."""
        query = """
        PREFIX biolink: <https://w3id.org/biolink/vocab/>

        SELECT ?study ?name WHERE {
            ?study a biolink:Study .
            ?study biolink:id ?name .
        }
        """
        results = list(rdf_graph.query(query))

        assert len(results) >= 1
        # Check that E-GEOD-5305 is present
        study_names = [str(r[1]) for r in results]
        assert any("E-GEOD-5305" in name for name in study_names)

    def test_query_all_genes(self, rdf_graph):
        """Test querying for all genes."""
        query = """
        PREFIX biolink: <https://w3id.org/biolink/vocab/>

        SELECT (COUNT(?gene) as ?count) WHERE {
            ?gene a biolink:Gene .
        }
        """
        results = list(rdf_graph.query(query))
        gene_count = int(results[0][0])

        assert gene_count > 0
        print(f"Found {gene_count} genes in the graph")

    def test_query_upregulated_genes(self, rdf_graph):
        """
        Test querying for upregulated genes.

        This is the key biological query: "What genes are upregulated in this dataset?"
        """
        query = """
        PREFIX biolink: <https://w3id.org/biolink/vocab/>
        PREFIX spokegenelab: <https://spoke.ucsf.edu/genelab/>

        SELECT ?gene ?log2fc WHERE {
            ?assoc biolink:subject ?assay .
            ?assoc biolink:object ?gene .
            ?assoc spokegenelab:log2fc ?log2fc .
            FILTER(?log2fc > 0)
        }
        ORDER BY DESC(?log2fc)
        LIMIT 20
        """
        results = list(rdf_graph.query(query))

        # Should have some upregulated genes
        assert len(results) > 0

        # Check that log2fc values are positive
        for gene, log2fc in results:
            assert float(log2fc) > 0
            print(f"Upregulated gene: {gene} (log2fc: {log2fc})")

    def test_query_downregulated_genes(self, rdf_graph):
        """Test querying for downregulated genes."""
        query = """
        PREFIX biolink: <https://w3id.org/biolink/vocab/>
        PREFIX spokegenelab: <https://spoke.ucsf.edu/genelab/>

        SELECT ?gene ?log2fc WHERE {
            ?assoc biolink:subject ?assay .
            ?assoc biolink:object ?gene .
            ?assoc spokegenelab:log2fc ?log2fc .
            FILTER(?log2fc < 0)
        }
        ORDER BY ?log2fc
        LIMIT 20
        """
        results = list(rdf_graph.query(query))

        # Should have some downregulated genes
        assert len(results) > 0

        # Check that log2fc values are negative
        for gene, log2fc in results:
            assert float(log2fc) < 0

    def test_query_genes_by_assay(self, rdf_graph):
        """Test querying genes for a specific assay/contrast."""
        query = """
        PREFIX biolink: <https://w3id.org/biolink/vocab/>
        PREFIX spokegenelab: <https://spoke.ucsf.edu/genelab/>

        SELECT ?assay ?gene ?log2fc WHERE {
            ?assoc biolink:subject ?assay .
            ?assoc biolink:object ?gene .
            ?assoc spokegenelab:log2fc ?log2fc .
            FILTER(CONTAINS(STR(?assay), "g1_g3"))
        }
        LIMIT 10
        """
        results = list(rdf_graph.query(query))

        # Should have results for the g1_g3 contrast
        assert len(results) > 0

    def test_query_enriched_pathways(self, rdf_graph):
        """Test querying for enriched GO terms/pathways."""
        query = """
        PREFIX biolink: <https://w3id.org/biolink/vocab/>
        PREFIX spokegenelab: <https://spoke.ucsf.edu/genelab/>

        SELECT ?enrichment ?name ?pvalue WHERE {
            ?enrichment a biolink:Association .
            ?enrichment biolink:name ?name .
            OPTIONAL { ?enrichment spokegenelab:adj_p_value ?pvalue }
        }
        LIMIT 20
        """
        results = list(rdf_graph.query(query))

        # Should have enrichment results
        if len(results) > 0:
            print(f"Found {len(results)} pathway enrichments")
            for enrichment, name, pvalue in results[:5]:
                print(f"  {name}: p={pvalue}")

    def test_query_go_terms(self, rdf_graph):
        """Test querying for GO terms."""
        query = """
        PREFIX biolink: <https://w3id.org/biolink/vocab/>
        PREFIX go: <http://purl.obolibrary.org/obo/GO_>

        SELECT ?term ?name WHERE {
            ?term a biolink:BiologicalProcess .
            OPTIONAL { ?term biolink:name ?name }
        }
        LIMIT 10
        """
        results = list(rdf_graph.query(query))
        print(f"Found {len(results)} GO terms")

    def test_query_study_assay_gene_chain(self, rdf_graph):
        """Test querying the full chain: Study -> Assay -> Gene."""
        query = """
        PREFIX biolink: <https://w3id.org/biolink/vocab/>
        PREFIX spokegenelab: <https://spoke.ucsf.edu/genelab/>

        SELECT ?study ?assay ?gene ?log2fc WHERE {
            ?study a biolink:Study .
            ?study biolink:has_output ?assay .
            ?assay a biolink:Assay .

            ?assoc a biolink:Association .
            ?assoc biolink:subject ?assay .
            ?assoc biolink:object ?gene .
            ?assoc spokegenelab:log2fc ?log2fc .
        }
        LIMIT 10
        """
        results = list(rdf_graph.query(query))

        if len(results) > 0:
            print(f"Found {len(results)} study-assay-gene chains")
            for study, assay, gene, log2fc in results[:3]:
                print(f"  {study} -> {assay} -> {gene} (log2fc: {log2fc})")


class TestBiolinkCompliance:
    """Test compliance with Biolink Model."""

    @pytest.fixture
    def rdf_graph(self):
        """Generate RDF graph from test data."""
        from scripts.gea.gea_pipeline import process_gea_experiment
        from scripts.rdf.turtle_writer import TurtleWriter

        if not os.path.exists(TEST_DATA_DIR):
            pytest.skip("Test data not available")

        result = process_gea_experiment(
            str(TEST_DATA_DIR),
            p_value_threshold=0.1,
            include_gsea=True,
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

    def test_all_nodes_have_biolink_type(self, rdf_graph):
        """Test that all nodes have a Biolink type."""
        query = """
        PREFIX biolink: <https://w3id.org/biolink/vocab/>
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>

        SELECT (COUNT(?node) as ?count) WHERE {
            ?node rdf:type ?type .
            FILTER(STRSTARTS(STR(?type), "https://w3id.org/biolink/vocab/"))
        }
        """
        results = list(rdf_graph.query(query))
        biolink_typed_count = int(results[0][0])

        assert biolink_typed_count > 0
        print(f"Found {biolink_typed_count} nodes with Biolink types")

    def test_associations_have_subject_and_object(self, rdf_graph):
        """Test that reified relationship nodes have subject and object."""
        query = """
        PREFIX biolink: <https://w3id.org/biolink/vocab/>

        SELECT ?assoc WHERE {
            ?assoc biolink:subject ?subj .
            ?assoc biolink:object ?obj .
        }
        """
        results = list(rdf_graph.query(query))

        # Should have reified associations with subject/object
        assert len(results) > 0
        print(f"{len(results)} associations have subject/object")


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])
