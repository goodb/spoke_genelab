"""
Main CLI orchestrator for SPOKE-GeneLab pipeline.

Provides a unified command-line interface for running both OSDR and GEA
data processing pipelines with options for Neo4j CSV and RDF output.
"""

import argparse
import sys
from pathlib import Path

from .common.config import load_config, Config


def run_osdr_pipeline(args, config: Config):
    """Run the OSDR (NASA Open Science Data Repository) pipeline."""
    print("OSDR pipeline not yet implemented as scripts.")
    print("Please use the Jupyter notebooks in notebooks/ directory.")
    return 1


def run_gea_pipeline(args, config: Config):
    """Run the GEA (Gene Expression Atlas) pipeline."""
    from .gea.gea_pipeline import run_gea_pipeline as gea_run

    result = gea_run(
        input_dir=args.input_dir,
        output_dir=args.output_dir,
        config=config,
        output_csv=args.csv or not args.rdf,
        output_rdf=args.rdf,
        bioportal_apikey=args.bioportal_key or config.bioportal_api_key,
        p_value_threshold=args.p_value,
        include_gsea=not args.no_gsea,
        include_orthologs=not args.no_orthologs,
    )

    if result.errors:
        print("\nErrors encountered:")
        for error in result.errors:
            print(f"  - {error}")
        return 1

    return 0


def main():
    """Main entry point for the CLI."""
    parser = argparse.ArgumentParser(
        description="SPOKE-GeneLab Pipeline CLI",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Process a GEA experiment and output CSV files
  python -m scripts.main --source gea --input-dir ./E-GEOD-5305-gea --csv

  # Process and generate RDF output
  python -m scripts.main --source gea --input-dir ./E-GEOD-5305-gea --rdf

  # Process with both CSV and RDF output
  python -m scripts.main --source gea --input-dir ./E-GEOD-5305-gea --csv --rdf

  # Process with custom p-value threshold
  python -m scripts.main --source gea --input-dir ./E-GEOD-5305-gea --p-value 0.05 --rdf
        """,
    )

    # Required arguments
    parser.add_argument(
        "--source",
        choices=["osdr", "gea"],
        required=True,
        help="Data source: 'osdr' for NASA Open Science Data Repository, 'gea' for Gene Expression Atlas",
    )

    # Input/output arguments
    parser.add_argument(
        "--input-dir", "-i",
        help="Input directory (required for GEA, ignored for OSDR)",
    )
    parser.add_argument(
        "--output-dir", "-o",
        default="./output",
        help="Output directory (default: ./output)",
    )
    parser.add_argument(
        "--config",
        help="Path to configuration file (YAML)",
    )

    # Output format arguments
    parser.add_argument(
        "--csv",
        action="store_true",
        help="Output Neo4j CSV files",
    )
    parser.add_argument(
        "--rdf",
        action="store_true",
        help="Output RDF Turtle files",
    )

    # Processing arguments
    parser.add_argument(
        "--p-value",
        type=float,
        default=0.1,
        help="Adjusted p-value threshold for significance (default: 0.1)",
    )
    parser.add_argument(
        "--bioportal-key",
        help="BioPortal API key for ontology mapping",
    )
    parser.add_argument(
        "--no-gsea",
        action="store_true",
        help="Skip GSEA/pathway enrichment extraction",
    )
    parser.add_argument(
        "--no-orthologs",
        action="store_true",
        help="Skip ortholog mapping to human genes",
    )

    # Parse arguments
    args = parser.parse_args()

    # Validate arguments
    if args.source == "gea" and not args.input_dir:
        parser.error("--input-dir is required for GEA source")

    # Load configuration
    config = load_config(config_file=args.config)

    # Ensure output directory exists
    Path(args.output_dir).mkdir(parents=True, exist_ok=True)

    # Run appropriate pipeline
    if args.source == "osdr":
        return run_osdr_pipeline(args, config)
    elif args.source == "gea":
        return run_gea_pipeline(args, config)

    return 0


if __name__ == "__main__":
    sys.exit(main())
