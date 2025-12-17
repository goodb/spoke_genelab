"""
SPOKE-GeneLab Pipeline Scripts

This package contains modular scripts for processing gene expression data
from NASA OSDR and Gene Expression Atlas (GEA) into a Neo4j knowledge graph
and RDF/Biolink output.

Subpackages:
    - pipeline: OSDR data processing pipeline (refactored from notebooks)
    - gea: Gene Expression Atlas local file import pipeline
    - rdf: RDF/Biolink model output generation
    - common: Shared utilities and configuration
"""

__version__ = "0.2.0"
