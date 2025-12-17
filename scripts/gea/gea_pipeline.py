"""
GEA Pipeline orchestrator.

Main orchestration module for processing GEA experiment directories
and producing graph outputs (Neo4j CSV and/or RDF).
"""

import os
from pathlib import Path
from typing import Dict, List, Optional, Union

import pandas as pd

from .gea_parser import GEAExperiment, load_gea_experiment, find_gea_files
from .gea_study_extractor import extract_study_nodes, get_study_summary
from .gea_gene_extractor import (
    create_mgene_nodes,
    create_gene_and_ortholog_data,
    extract_differential_expression,
)
from .gea_assay_extractor import (
    extract_assay_nodes,
    extract_materials_from_samples,
    map_materials_to_ontology,
    assign_materials_to_assays,
    create_anatomy_nodes,
    create_celltype_nodes,
    create_study_assay_relationships,
    create_assay_anatomy_relationships,
    create_assay_celltype_relationships,
)
from .gea_gsea_extractor import extract_all_gsea_data

from ..common.csv_writer import save_graph_to_csv
from ..common.config import Config


class GEAPipelineResult:
    """Container for GEA pipeline results."""

    def __init__(self):
        self.nodes: Dict[str, pd.DataFrame] = {}
        self.relationships: Dict[str, pd.DataFrame] = {}
        self.experiment: Optional[GEAExperiment] = None
        self.errors: List[str] = []

    def add_nodes(self, node_type: str, df: pd.DataFrame):
        """Add or merge nodes of a given type."""
        if df.empty:
            return
        if node_type in self.nodes:
            self.nodes[node_type] = pd.concat(
                [self.nodes[node_type], df], ignore_index=True
            ).drop_duplicates(subset="identifier" if "identifier" in df.columns else None)
        else:
            self.nodes[node_type] = df

    def add_relationships(self, rel_type: str, df: pd.DataFrame):
        """Add or merge relationships of a given type."""
        if df.empty:
            return
        if rel_type in self.relationships:
            self.relationships[rel_type] = pd.concat(
                [self.relationships[rel_type], df], ignore_index=True
            ).drop_duplicates()
        else:
            self.relationships[rel_type] = df

    def get_stats(self) -> Dict[str, int]:
        """Get statistics about the pipeline results."""
        stats = {}
        for node_type, df in self.nodes.items():
            stats[f"nodes_{node_type}"] = len(df)
        for rel_type, df in self.relationships.items():
            stats[f"relationships_{rel_type}"] = len(df)
        return stats


def process_gea_experiment(
    experiment_dir: Union[str, Path],
    config: Optional[Config] = None,
    bioportal_apikey: Optional[str] = None,
    p_value_threshold: float = 0.1,
    include_gsea: bool = True,
    include_orthologs: bool = True,
) -> GEAPipelineResult:
    """
    Process a single GEA experiment directory.

    Args:
        experiment_dir: Path to experiment directory
        config: Optional Config object
        bioportal_apikey: BioPortal API key for ontology mapping
        p_value_threshold: Adjusted p-value threshold for filtering
        include_gsea: Whether to include GSEA/pathway enrichment data
        include_orthologs: Whether to map orthologs to human genes

    Returns:
        GEAPipelineResult with all nodes and relationships
    """
    result = GEAPipelineResult()

    # Load experiment
    print(f"Loading experiment from {experiment_dir}")
    try:
        experiment = load_gea_experiment(str(experiment_dir))
        result.experiment = experiment
    except Exception as e:
        result.errors.append(f"Failed to load experiment: {e}")
        return result

    # Print summary
    summary = get_study_summary(experiment)
    print(f"  Accession: {summary['accession']}")
    print(f"  Title: {summary['title']}")
    print(f"  Organism: {summary['organism']}")
    print(f"  Contrasts: {summary['num_contrasts']}")

    # Extract Study nodes
    print("Extracting Study nodes...")
    study_nodes = extract_study_nodes(experiment)
    result.add_nodes("Study", study_nodes)

    # Extract Assay nodes
    print("Extracting Assay nodes...")
    assay_nodes = extract_assay_nodes(experiment)

    # Map materials to ontology if API key provided
    if bioportal_apikey:
        print("Mapping materials to ontology...")
        materials = extract_materials_from_samples(experiment)
        if not materials.empty:
            mapped_materials = map_materials_to_ontology(materials, bioportal_apikey)
            assay_nodes = assign_materials_to_assays(assay_nodes, experiment, mapped_materials)

    result.add_nodes("Assay", assay_nodes)

    # Create Study-Assay relationships
    study_assay_rels = create_study_assay_relationships(assay_nodes)
    result.add_relationships("Study-PERFORMED_SpAS-Assay", study_assay_rels)

    # Create Anatomy and CellType nodes and relationships
    anatomy_nodes = create_anatomy_nodes(assay_nodes)
    result.add_nodes("Anatomy", anatomy_nodes)

    celltype_nodes = create_celltype_nodes(assay_nodes)
    result.add_nodes("CellType", celltype_nodes)

    assay_anatomy_rels = create_assay_anatomy_relationships(assay_nodes)
    result.add_relationships("Assay-INVESTIGATED_ASiA-Anatomy", assay_anatomy_rels)

    assay_celltype_rels = create_assay_celltype_relationships(assay_nodes)
    result.add_relationships("Assay-INVESTIGATED_ASiCT-CellType", assay_celltype_rels)

    # Extract MGene nodes
    print("Extracting MGene nodes...")
    mgene_nodes = create_mgene_nodes(experiment)
    result.add_nodes("MGene", mgene_nodes)

    # Map to human orthologs
    if include_orthologs and not mgene_nodes.empty:
        print("Mapping orthologs...")
        ortholog_data = create_gene_and_ortholog_data(mgene_nodes)
        result.add_nodes("Gene", ortholog_data.get("Gene", pd.DataFrame()))
        result.add_relationships(
            "MGene-IS_ORTHOLOG_MGiG-Gene",
            ortholog_data.get("MGene-IS_ORTHOLOG_MGiG-Gene", pd.DataFrame())
        )

    # Extract differential expression relationships
    print("Extracting differential expression data...")
    de_rels = extract_differential_expression(experiment, p_value_threshold)
    result.add_relationships("Assay-MEASURED_DIFFERENTIAL_EXPRESSION_ASmMG-MGene", de_rels)

    # Extract GSEA data
    if include_gsea:
        print("Extracting GSEA/pathway enrichment data...")
        gsea_data = extract_all_gsea_data(experiment, p_value_threshold)

        for key, df in gsea_data.items():
            if "-" in key:  # Relationship
                result.add_relationships(key, df)
            else:  # Node
                result.add_nodes(key, df)

    # Print stats
    stats = result.get_stats()
    print("\nPipeline Results:")
    for key, count in sorted(stats.items()):
        print(f"  {key}: {count}")

    return result


def process_gea_batch(
    experiments_dir: Union[str, Path],
    config: Optional[Config] = None,
    **kwargs,
) -> GEAPipelineResult:
    """
    Process all GEA experiments in a directory.

    Args:
        experiments_dir: Directory containing experiment subdirectories
        config: Optional Config object
        **kwargs: Additional arguments passed to process_gea_experiment

    Returns:
        GEAPipelineResult with combined data from all experiments
    """
    combined_result = GEAPipelineResult()
    experiments_dir = Path(experiments_dir)

    # Find experiment directories
    exp_dirs = [d for d in experiments_dir.iterdir() if d.is_dir() and "-gea" in d.name]

    print(f"Found {len(exp_dirs)} experiment directories")

    for exp_dir in exp_dirs:
        try:
            result = process_gea_experiment(exp_dir, config, **kwargs)

            # Merge results
            for node_type, df in result.nodes.items():
                combined_result.add_nodes(node_type, df)
            for rel_type, df in result.relationships.items():
                combined_result.add_relationships(rel_type, df)

            combined_result.errors.extend(result.errors)

        except Exception as e:
            combined_result.errors.append(f"{exp_dir.name}: {e}")

    return combined_result


def run_gea_pipeline(
    input_dir: Union[str, Path],
    output_dir: Union[str, Path],
    config: Optional[Config] = None,
    output_csv: bool = True,
    output_rdf: bool = False,
    **kwargs,
) -> GEAPipelineResult:
    """
    Main entry point for GEA pipeline.

    Args:
        input_dir: Input directory (single experiment or batch)
        output_dir: Output directory for CSV/RDF files
        config: Optional Config object
        output_csv: Whether to output Neo4j CSV files
        output_rdf: Whether to output RDF Turtle files
        **kwargs: Additional arguments

    Returns:
        GEAPipelineResult with all processed data
    """
    input_path = Path(input_dir)
    output_path = Path(output_dir)

    # Determine if single experiment or batch
    if "-gea" in input_path.name:
        # Single experiment
        result = process_gea_experiment(input_path, config, **kwargs)
    else:
        # Batch processing
        result = process_gea_batch(input_path, config, **kwargs)

    # Output CSV files
    if output_csv:
        node_dir = output_path / "nodes"
        rel_dir = output_path / "relationships"
        node_dir.mkdir(parents=True, exist_ok=True)
        rel_dir.mkdir(parents=True, exist_ok=True)

        save_graph_to_csv(
            result.nodes,
            result.relationships,
            node_dir,
            rel_dir,
        )

    # Output RDF files
    if output_rdf:
        from ..rdf.turtle_writer import write_graph_to_turtle

        rdf_dir = output_path / "rdf"
        rdf_dir.mkdir(parents=True, exist_ok=True)

        write_graph_to_turtle(
            result.nodes,
            result.relationships,
            rdf_dir / "gxa_rdf.ttl",
        )

    return result


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="GEA Pipeline")
    parser.add_argument("input_dir", help="Input directory (experiment or batch)")
    parser.add_argument("--output-dir", "-o", default="./output", help="Output directory")
    parser.add_argument("--csv", action="store_true", help="Output Neo4j CSV files")
    parser.add_argument("--rdf", action="store_true", help="Output RDF Turtle files")
    parser.add_argument("--p-value", type=float, default=0.1, help="P-value threshold")
    parser.add_argument("--bioportal-key", help="BioPortal API key")
    parser.add_argument("--no-gsea", action="store_true", help="Skip GSEA extraction")
    parser.add_argument("--no-orthologs", action="store_true", help="Skip ortholog mapping")

    args = parser.parse_args()

    result = run_gea_pipeline(
        args.input_dir,
        args.output_dir,
        output_csv=args.csv or not args.rdf,  # Default to CSV if neither specified
        output_rdf=args.rdf,
        bioportal_apikey=args.bioportal_key,
        p_value_threshold=args.p_value,
        include_gsea=not args.no_gsea,
        include_orthologs=not args.no_orthologs,
    )

    if result.errors:
        print("\nErrors:")
        for error in result.errors:
            print(f"  {error}")
