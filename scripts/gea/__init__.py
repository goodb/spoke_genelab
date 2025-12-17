"""
Gene Expression Atlas (GEA) local file import pipeline.

This package processes local GEA experiment directories (e.g., E-GEOD-5305-gea/)
and converts them to the SPOKE-GeneLab knowledge graph format.

Modules:
    - gea_parser: Parse GEA file formats (IDF, SDRF, configuration.xml, analytics)
    - gea_study_extractor: Extract study/experiment metadata
    - gea_gene_extractor: Extract gene and differential expression data
    - gea_assay_extractor: Extract assay and contrast information
    - gea_gsea_extractor: Extract GSEA pathway enrichment results
    - gea_pipeline: Main GEA pipeline orchestrator
"""
