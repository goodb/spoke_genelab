# NASA SPOKE-GeneLab Knowledge Graph

This repository contains the code and metadata needed to build a **Knowledge Graph (KG)** for [NASA GeneLab](https://www.nasa.gov/osdr-genelab-about/) omics datasets hosted on the [Open Science Data Repository (OSDR)](https://osdr.nasa.gov/bio/repo/search?q=&data_source=cgene,alsda&data_type=study).

---

## 🚀 Features

- **Automated graph construction** from datasets in the OSDR
- **Incremental update** for new datasets
- **Statistical filtering** of results for significance
- **Species selection** via a configurable whitelist
- **Versioned metadata** for reproducibility (v0.1.0)
- **Federated query** using Neo4j Fabric with the [Scalable Precision Medicine Open Knowledge Engine (SPOKE) KG](https://spoke.ucsf.edu/)

---

## 🧪 Supported Data Types

| Measurement                  | Technology                                              | Property         | Selection Criteria |
| ---------------------------- | ------------------------------------------------------- | -----------------|-----------------|
| Transcription profiling      | RNA Sequencing (RNA‑Seq)                                | Log2 fold change | Adjusted p-value <= 0.1 |
| Transcription profiling      | DNA microarray                                          | Log2 fold change | Adjusted p-value <= 0.1 |
| DNA methylation profiling    | Whole Genome Bisulfite Sequencing                       | Methylation difference % | q-value <= 0.1 |
| DNA methylation profiling    | Reduced‑Representation Bisulfite Sequencing (RRBS)      | Methylation difference % | q-value <= 0.1 |

---

## ⚙️ How It Works

1. **Fetch** omics study records using the OSDR API  
2. **Filter** datasets by statistical thresholds and target species  
3. **Map** model organism genes to human genes
4. **Map** cell and tissue types to the [Cell (CL)](https://bioportal.bioontology.org/ontologies/CL) and [Uber Anatomy Ontology (UBERON)](https://bioportal.bioontology.org/ontologies/UBERON) ontology, respectively 
5. **Export** CSV files for graph database upload
6. **Import** CSV files into a Neo4j Graph database

---

## 🕸️ Graph Schema

![Simplified Graph Schema](docs/spoke-genelab-v0.0.5-simplified.png)

**Figure**: Schematic overview of the GeneLab knowledge graph structure, highlighting key node types (circles) and relationships (arrows).

The `MEASURED_DIFFERENTIAL_EXPRESSION` relationship encodes Log₂ fold changes derived from transcription profiling assays, while the `MEASURED_DIFFERENTIAL_METHYLATION` relationship captures methylation differences identified through DNA methylation assays. The `METHYLATED_IN` relationship links model organism genes (`MGene`) to 1,000 base pair genomic regions (`MethylationRegion`) exhibiting differential methylation.

Proxy nodes (shown in gray) represent standardized identifiers for human genes (ENTREZ ID), anatomical structures (UBERON ID), and cell types (CL ID), enabling integration with external Neo4j databases and supporting composite graph database construction.

Diagram generated using [arrows.app](https://arrows.app).

---

## 📁 Metadata Directory Structure

The following node and relationship metadata files define the graph schema.

- **Nodes**  
  [kg/v0.1.0/metadata/nodes/](kg/v0.1.0/metadata/nodes/)

- **Relationships**   
  [kg/v0.1.0/metadata/relationships/](kg/v0.1.0/metadata/relationships/)

The organization and conventions for defining the metadata and data are described in the [kg-import](https://github.com/sbl-sdsc/kg-import) Git repository.

---

## 🔗 SPOKE - GeneLab Composite Database

![](docs/spoke-genlab-v0.0.3-composite.png)

**Figure**: Integration of the SPOKE and GeneLab knowledge graphs using proxy nodes.  
The **GeneLab** graph (right), a knowledge graph representing spaceflight omics datasets, depicts key experimental entities: `Assay`, `Study`, `Mission`, `MGene`, and `MethylationRegion`, along with their relationships. 
**Proxy nodes** (gray) represent external identifiers (ENTREZ, UBERON, CL) and enable linkage to the **[SPOKE](https://spoke.ucsf.edu/)** graph (left), a rich biomedical knowledge graph comprising biological processes, molecular functions, diseases, compounds, and more. The dashed lines indicate mappings to enable the construction of a [composite Neo4j graph database](https://neo4j.com/docs/operations-manual/current/tutorial/tutorial-composite-database/). The composite graph enables federated queries across multiple KGs.

---

## ⚙️ Data Import Into Neo4j Knowledge Graph

### Setup Neo4j Desktop

1. Download the Neo4j Desktop application from the [Neo4j Download Center](https://neo4j.com/download-center/#desktop) and follow the installation instructions.

2. When the installation is complete, Neo4j Desktop will launch. Click the `New` button to create a new project.

![](docs/new_project.png)

3. Hover the cursor over the created project, click the edit button, and change the project name from `Project` to `spoke-genelab`.

![](docs/rename_project.png)

4. Click the `ADD` button and select `Local DBMS`. **Select Neo4j version 5.23.0.**

![](docs/add_graph_dbms.png)

5. Enter the password `neo4jdemo` and click `Create`.
    
![](docs/create_dbms.png)
    
6. Select `Terminal` to open a terminal window.
    
![](docs/open_terminal.png)

7. Type `pwd` in the terminal window to show the path to the `NEO4J_INSTALL_PATH` directory. This path is required in the `.env` file, see the next section.
 
![](docs/get_path.png)


------

### Setup the Environment

Prerequisites: Miniconda3 (light-weight, preferred) or Anaconda3 and Mamba (faster than Conda)

* Install [Miniconda3](https://docs.conda.io/en/latest/miniconda.html)
* Update an existing miniconda3 installation: ```conda update conda```
* Install Mamba: ```conda install mamba -n base -c conda-forge```
* Install Git (if not installed): ```conda install git -n base -c anaconda```
------

1. Clone this Repository

```
git clone https://github.com/BaranziniLab/spoke_genelab.git
cd spoke_genelab
```

2. Create a Conda environment

The file `environment.yml` specifies the Python version and all required dependencies.

```
mamba env create -f environment.yml
```

3. Create an account in [BioPortal](https://bioportal.bioontology.org/) and copy the API key. BioPortal is used to map terms to ontologies.

   
4. Copy the file `env_template` to `.env`

5. Edit the file `.env` and set the following variables

KG version number

`KG_VERSION=v0.1.0`

Path to the cloned git repository

`KG_GIT=/Users/.../spoke_genelab/`

Path to the Neo4J instance in Neo4j Desktop (in quotes). Make sure to enclose the path in quotes.

`NEO4J_INSTALL_PATH="/Users/.../Library/Application Support/Neo4j Desktop/Application/relate-data/dbmss/dbms-3d4b95d1-0219-480b-a3c4-ee5a409cc383"`

BioPortal API Key

`BIOPORTAL_API_KEY=<bioportal api key>`

------

### Download and Process Datasets and upload to Neo4J Graph Database

1. Start the spoke-genelab Graph DBMS

![](docs/start_dbms.png)

2. Activate the conda environment

```
conda activate spoke-genelab
```

3. Launch Jupyter Lab

```
jupyter lab
```

4. Navigate to the `notebooks` directory and run the following notebooks

| Notebook                   |   Description           |
|----------------------------|-------------------------|
| 1_download_datasets.ipynb  | Downloads datasets     |
| 2_create_study_mission_nodes.ipynb | Creates Study and Mission nodes and their relationships |
| 3_create_gene_nodes.ipynb  | Creates MGene (model organism) and mapped Gene (human) gene nodes |
| 4_create_assay_nodes.ipynb | Creates Assay nodes and their relationships |
| 5_import_to_neo4j.ipynb    | Imports the formatted data into a Neo4j KG |
| 6_query_examples.ipynb     | Runs example queries (optional) |

5. When the import is completed, click the `Refresh` button in Neo4j Desktop. The newly created database `spoke-genelab-v0.1.0` will be listed.

![](docs/db_imported.png)

6. Click the `Open` button to launch the database.

![](docs/open_dbms.png)

7. Click on the database icon on the left.

![](docs/select_db_icon.png)

8. Use the pull-down menu to select a version of `spoke-genelab-v0.1.0` database. Wait for about 30+ seconds until the database is loaded and the nodes are listed as shown below.
   
![](docs/db_ready.png)

9. Set the Graph Stylesheet

Drag the file kg/v0.1.0/style.grass onto the Neo4j Browser window to set the node colors, sizes, and labels.

10. Now you are ready to run Cypher queries on the selected database.

11. When you are finished, stop the database in the Neo4j Desktop.

To stop the conda environment, type

```conda deactivate```

------

### Dump Neo4J Graph Database
1. Stop the database

2. Hover the cursor over the `spoke-genelab-v0.1.0` database and select `Dump` from the menu.

![](docs/dump_db.png)

3. When the dump is complete, click the `Reveal files in Finder` button to open the directory that contains the `spoke-genelab-v0.1.0.dump` file.

![](docs/dump_location.png)

This database dump will be used to create the SPOKE-GeneLab composite database.

------

## 🧬 Gene Expression Atlas (GEA) Import & RDF Export

In addition to the OSDR pipeline, this repository includes a **command-line pipeline** for importing local [Gene Expression Atlas (GEA)](https://www.ebi.ac.uk/gxa/) experiment data and exporting to **RDF (Turtle format)** following the [Biolink Model](https://biolink.github.io/biolink-model/).

### Supported GEA Data

The GEA pipeline processes experiment directories containing:
- **IDF files** (`.idf.txt`) - Experiment metadata
- **SDRF files** (`.condensed-sdrf.tsv`) - Sample metadata with ontology URIs
- **Configuration XML** (`-configuration.xml`) - Assay groups and contrasts
- **Analytics TSV** (`-analytics.tsv`) - Differential expression results
- **GSEA files** (`.gsea.tsv`) - Gene set enrichment analysis for GO, Reactome, InterPro

### Installation

```bash
# Install additional dependencies
pip install rdflib pytest
```

### CLI Usage

```bash
# Process a GEA experiment directory and output Neo4j CSV files
python -m scripts.main --source gea --input-dir ./E-GEOD-5305-gea --csv

# Process and generate RDF Turtle output
python -m scripts.main --source gea --input-dir ./E-GEOD-5305-gea --rdf

# Both CSV and RDF output
python -m scripts.main --source gea --input-dir ./E-GEOD-5305-gea --csv --rdf

# Custom p-value threshold
python -m scripts.main --source gea --input-dir ./E-GEOD-5305-gea --rdf --p-value 0.05

# Skip GSEA extraction (faster processing)
python -m scripts.main --source gea --input-dir ./E-GEOD-5305-gea --rdf --no-gsea
```

### Output Structure

```
output/
├── nodes/                    # Neo4j node CSV files
│   ├── Study_*.csv
│   ├── Assay_*.csv
│   ├── MGene_*.csv
│   ├── PathwayEnrichment_*.csv
│   ├── GOTerm_*.csv
│   ├── ReactomePathway_*.csv
│   └── InterProDomain_*.csv
├── relationships/            # Neo4j relationship CSV files
│   ├── Study-PERFORMED_SpAS-Assay_*.csv
│   ├── Assay-MEASURED_DIFFERENTIAL_EXPRESSION_ASmMG-MGene_*.csv
│   └── ...
└── rdf/
    └── gxa_rdf.ttl           # RDF Turtle file (Biolink compliant)
```

### RDF Biolink Model Mapping

The RDF output follows the [Biolink Model](https://biolink.github.io/biolink-model/) ontology:

| Graph Node | Biolink Class | URI Pattern |
|------------|---------------|-------------|
| Study | `biolink:Study` | `spokegenelab:Study/{id}` |
| Assay | `biolink:Assay` | `spokegenelab:Assay/{id}` |
| MGene | `biolink:Gene` | `ncbigene:{ensembl_id}` |
| Gene | `biolink:Gene` | `ncbigene:{entrez_id}` |
| PathwayEnrichment | `biolink:Association` | `spokegenelab:Enrichment/{id}` |
| GOTerm | `biolink:BiologicalProcess` | `go:{GO_id}` |
| ReactomePathway | `biolink:Pathway` | `reactome:{id}` |

Differential expression relationships are **reified** as Biolink Associations with properties:
- `spokegenelab:log2fc` - Log2 fold change
- `spokegenelab:adj_p_value` - Adjusted p-value

### Querying RDF with SPARQL

The RDF output can be loaded into any SPARQL-compatible triple store or queried directly with Python:

```python
from rdflib import Graph

# Load the RDF graph
g = Graph()
g.parse("output/rdf/gxa_rdf.ttl", format="turtle")

# Query for upregulated genes
query = """
PREFIX biolink: <https://w3id.org/biolink/vocab/>
PREFIX spokegenelab: <https://spoke.ucsf.edu/genelab/>

SELECT ?gene ?log2fc WHERE {
    ?assoc biolink:subject ?assay .
    ?assoc biolink:object ?gene .
    ?assoc spokegenelab:log2fc ?log2fc .
    FILTER(?log2fc > 1.0)
}
ORDER BY DESC(?log2fc)
LIMIT 20
"""
results = g.query(query)
for row in results:
    print(f"Gene: {row.gene}, Log2FC: {row.log2fc}")
```

#### Example SPARQL Queries

**Find all upregulated genes (log2fc > 0):**
```sparql
PREFIX biolink: <https://w3id.org/biolink/vocab/>
PREFIX spokegenelab: <https://spoke.ucsf.edu/genelab/>

SELECT ?gene ?log2fc ?pvalue WHERE {
    ?assoc biolink:object ?gene .
    ?assoc spokegenelab:log2fc ?log2fc .
    ?assoc spokegenelab:adj_p_value ?pvalue .
    FILTER(?log2fc > 0)
}
ORDER BY DESC(?log2fc)
```

**Find downregulated genes (log2fc < 0):**
```sparql
PREFIX biolink: <https://w3id.org/biolink/vocab/>
PREFIX spokegenelab: <https://spoke.ucsf.edu/genelab/>

SELECT ?gene ?log2fc WHERE {
    ?assoc biolink:object ?gene .
    ?assoc spokegenelab:log2fc ?log2fc .
    FILTER(?log2fc < -1.0)
}
ORDER BY ?log2fc
```

**Find genes for a specific assay/contrast:**
```sparql
PREFIX biolink: <https://w3id.org/biolink/vocab/>
PREFIX spokegenelab: <https://spoke.ucsf.edu/genelab/>

SELECT ?gene ?log2fc WHERE {
    ?assoc biolink:subject ?assay .
    ?assoc biolink:object ?gene .
    ?assoc spokegenelab:log2fc ?log2fc .
    FILTER(CONTAINS(STR(?assay), "g1_g3"))
}
```

**List enriched GO terms:**
```sparql
PREFIX biolink: <https://w3id.org/biolink/vocab/>
PREFIX spokegenelab: <https://spoke.ucsf.edu/genelab/>

SELECT ?term ?name ?pvalue WHERE {
    ?enrichment a biolink:Association .
    ?enrichment biolink:name ?name .
    ?enrichment spokegenelab:adj_p_value ?pvalue .
    ?enrichment biolink:object ?term .
    FILTER(CONTAINS(STR(?term), "GO_"))
}
ORDER BY ?pvalue
```

### Running Tests

```bash
# Run all tests
python -m pytest tests/ -v

# Run only GEA parser tests
python -m pytest tests/test_gea_parser.py -v

# Run RDF integration tests
python -m pytest tests/test_rdf_integration.py -v
```

### Scripts Directory Structure

```
scripts/
├── __init__.py
├── main.py                    # CLI entry point
├── common/
│   ├── config.py             # Configuration management
│   ├── graph_builder.py      # Node/relationship utilities
│   └── csv_writer.py         # Neo4j CSV output
├── gea/
│   ├── gea_parser.py         # GEA file parsing
│   ├── gea_study_extractor.py
│   ├── gea_gene_extractor.py
│   ├── gea_assay_extractor.py
│   ├── gea_gsea_extractor.py # Pathway enrichment
│   └── gea_pipeline.py       # GEA orchestrator
└── rdf/
    ├── rdf_config.py         # Namespace definitions
    ├── biolink_mapper.py     # Biolink Model mappings
    └── turtle_writer.py      # RDF Turtle generation
```

------

## 📚 Citation

PW Rose, CA Nelson, SG Gebre, AM Saravia-Butler, K Soman, KA Grigorev, LM Sanders, SV Costes, SE Baranzini, NASA SPOKE-GeneLab Knowledge Graph. Available online: https://github.com/BaranziniLab/spoke_genelab (2025)

CA Nelson, PW Rose, K Soman, LM Sanders, SG Gebre, SV Costes, SE Baranzini, Nasa Genelab-Knowledge Graph Fabric Enables Deep Biomedical Analysis of Multi-Omics Datasets, https://ntrs.nasa.gov/citations/20250000723 (2025)

------

## 💰 Funding
NSF Award number [2333819](https://www.nsf.gov/awardsearch/showAward?AWD_ID=2333819), Proto-OKN Theme 1: Connecting Biomedical information on Earth and in Space via the SPOKE knowledge graph.

