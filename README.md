# ChemGraphBuilder [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.13148036.svg)](https://doi.org/10.5281/zenodo.13148036)

`chemgraphbuilder` is a Python package designed for transforming chemical data into knowledge graphs.  
It leverages **PubChem** for data extraction and **Neo4j** for building graph databases, enabling researchers to efficiently extract, process, and visualize complex chemical relationships with precision.  
The package is designed for easy extension to include other data sources in future releases.

---

## Table of Contents
- [ChemGraphBuilder ](#chemgraphbuilder-)
  - [Table of Contents](#table-of-contents)
  - [Neo4j Requirements](#neo4j-requirements)
  - [🚀 Quick Start](#-quick-start)
    - [1️⃣ Install and Run Neo4j](#1️⃣-install-and-run-neo4j)
    - [2️⃣ Install `chemgraphbuilder`](#2️⃣-install-chemgraphbuilder)
    - [3️⃣ Connect to Neo4j in Python](#3️⃣-connect-to-neo4j-in-python)
    - [4️⃣ First Example](#4️⃣-first-example)
  - [Usage](#usage)
    - [From Python](#from-python)
    - [From Command Line](#from-command-line)
  - [Features](#features)
  - [Documentation](#documentation)
  - [Contributing](#contributing)
  - [License](#license)
  - [Contact](#contact)
  - [Acknowledgments](#acknowledgments)

---

## Neo4j Requirements

`chemgraphbuilder` requires a running **Neo4j** database that is accessible via **Bolt URI**, **username**, and **password**.

You can run Neo4j:
- **Locally** (Neo4j Desktop or Docker)
- **Remotely** (Neo4j Aura Cloud)

**Default Bolt port:** `7687`  
**Default Web UI port:** `7474`

---

## 🚀 Quick Start

Follow these steps to get up and running with `chemgraphbuilder` and Neo4j in under 5 minutes.

### 1️⃣ Install and Run Neo4j

**Option A – Docker (fastest)**
```bash
docker run \
  --name neo4j \
  -p 7474:7474 -p 7687:7687 \
  -e NEO4J_AUTH=neo4j/testpassword \
  neo4j:5.14
````

* Bolt URI: `bolt://localhost:7687`
* Username: `neo4j`
* Password: `testpassword`
* Web UI: [http://localhost:7474](http://localhost:7474)

**Option B – Neo4j Desktop**

1. Download from: [https://neo4j.com/download/](https://neo4j.com/download/)
2. Create a **new project** and database.
3. Note the **Bolt URI**, **username**, and **password**.

---

### 2️⃣ Install `chemgraphbuilder`

```bash
pip install chemgraphbuilder
```
Or visit the [PyPI Project Page](https://pypi.org/project/chemgraphbuilder) for the latest release.

---

### 3️⃣ Connect to Neo4j in Python

```python
from chemgraphbuilder import Neo4jBase

# Connection details
NEO4J_URI = "bolt://localhost:7687"
NEO4J_USER = "neo4j"
NEO4J_PASSWORD = "testpassword"

# Connect and test
db = Neo4jBase(uri=NEO4J_URI, user=NEO4J_USER, password=NEO4J_PASSWORD)
db.test_connection()
```

---

### 4️⃣ First Example

```python
# Create a simple test node
db.run_query("CREATE (:Test {name: 'Hello Neo4j'})")
print("Node created!")
```

Check in the Neo4j Browser:

```cypher
MATCH (n) RETURN n;
```

---

## Usage

### From Python

```python
from chemgraphbuilder.setup_data_folder import SetupDataFolder
from chemgraphbuilder.node_collector_processor import NodesCollectorProcessor
from chemgraphbuilder.relationship_collector_processor import RelationshipsCollectorProcessor
from chemgraphbuilder.graph_nodes_loader import GraphNodesLoader
from chemgraphbuilder.graph_relationships_loader import GraphRelationshipsLoader

# Setup data folder
setup_folder = SetupDataFolder()
setup_folder.setup()

# Collect nodes
collector = NodesCollectorProcessor(node_type=node_type, enzyme_list=enzyme_list, start_chunk=0)
collector.collect_and_process_data()

# Collect relationships
collector = RelationshipsCollectorProcessor(relationship_type=relationship_type, start_chunk=0)
collector.collect_relationship_data()

# Load nodes into Neo4j
graph_nodes_loader = GraphNodesLoader(uri, username, password)
graph_nodes_loader.load_data_for_node_type(label)
graph_nodes_loader.close()

# Load relationships into Neo4j
graph_relationships_loader = GraphRelationshipsLoader(uri, username, password)
graph_relationships_loader.add_relationships(relationship_type)
graph_relationships_loader.close()
```

### From Command Line

```bash
setup-data-folder
collect-process-nodes --node_type Compound --enzyme_list gene1,gene2 --start_chunk 0
collect-process-relationships --relationship_type Assay_Compound --start_chunk 0
load-graph-nodes --uri bolt://localhost:7687 --username neo4j --password password --label Compound
load-graph_relationships --uri bolt://localhost:7687 --username neo4j --password password --relationship_type Assay_Gene
```

More examples: [Usage Examples](https://asmaa-a-abdelwahab.github.io/ChemGraphBuilder/1.node_collector_processor/).

---

## Features

* **Node Representation:** Compounds, genes, proteins, bioassays.
* **Comprehensive Relationships:** Includes assay-compound, assay-gene, compound similarity, co-occurrence, inhibitor/activator/ligand, etc.
* **Data Integration:** Schema supports adding new sources.
* **Flexible Access:** Command line & Python API.

---

## Documentation

Full docs: [ChemGraphBuilder Documentation](https://asmaa-a-abdelwahab.github.io/ChemGraphBuilder/)

---

## Contributing

Issues: [GitHub Issues](https://github.com/asmaa-a-abdelwahab/chemgraphbuilder/issues)
Pull requests welcome.

---

## License

GPL-3.0 – see [LICENSE](https://github.com/asmaa-a-abdelwahab/chemgraphbuilder/blob/main/LICENSE).

---

## Contact

**Asmaa A. Abdelwahab** – [asmaa.a.abdelwahab@gmail.com](mailto:asmaa.a.abdelwahab@gmail.com)

---

## Acknowledgments

* **PubChem** – for chemical and bioassay data.
* **Neo4j** – for graph database capabilities.
