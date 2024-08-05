# ChemGraphBuilder
chemgraphbuilder is a Python package for transforming chemical data into knowledge graphs. Leveraging PubChem for data extraction and Neo4j for graph databases, it enables researchers to easily extract, process, and visualize complex chemical relationships with precision.

## Table of Contents

1. [Installation](#installation)
2. [Usage](#usage)
3. [Features](#features)
4. [Documentation](#documentation)
5. [Contributing](#contributing)
6. [License](#license)
7. [Contact](#contact)
8. [Acknowledgments](#acknowledgments)

## Installation

To install ChemGraphBuilder, use pip:

```bash
pip install chemgraphbuilder
```

## Usage

### From Python

```python
from chemgraphbuilder.setup_data_folder import SetupDataFolder
from chemgraphbuilder.node_collector_processor import NodesCollectorProcessor
from chemgraphbuilder.relationship_collector_processor import RelationshipsCollectorProcessor
from chemgraphbuilder.graph_nodes_loader import GraphNodesLoader
from chemgraphbuilder.graph_relationships_loader import GraphRelationshipsLoader

# Initialize and setup the data directory before collecting any data
setup_folder = SetupDataFolder()
setup_folder.setup()

# Initialize the collector & Collect and process the data
collector = NodesCollectorProcessor(node_type=node_type, enzyme_list=enzyme_list, start_chunk=0)
collector.collect_and_process_data()

# Initialize the collector & Collect and process the relationship data
collector = RelationshipsCollectorProcessor(relationship_type=relationship_type, start_chunk=0)
collector.collect_relationship_data()

# Initialize the loader & load nodes into neo4j database
graph_nodes_loader = GraphNodesLoader(uri, username, password)
graph_nodes_loader.load_data_for_node_type(label)
graph_nodes_loader.close()

# Initialize the loader & load relationships into neo4j database
graph_relationships_loader = GraphRelationshipsLoader(uri, username, password)
graph_relationships_loader.add_relationships(relationship_type)
graph_relationships_loader.close()
```

### From Command Line

```bash
setup-data-folder
collect-process-nodes --node_type Compound --enzyme_list gene1,gene2 --start_chunk 0 # the default start-chunk is 0
collect-process-relationships --relationship_type Assay_Compound --start_chunk 0
load-graph-nodes --uri bolt://localhost:7687 --username neo4j --password password --label Compound
load-graph_relationships --uri bolt://localhost:7687 --username neo4j --password password --relationship_type Assay_Gene
```

For more detailed examples, visit the [Usage Examples](https://asmaa-a-abdelwahab.github.io/ChemGraphBuilder/1.node_collector_processor/).

## Features

- **Node Representation:** Incorporates diverse nodes such as compounds, genes, proteins, and bioassays.
- **Comprehensive Relationships:** Maps out various interactions, including gene-protein relationships, bioassay-gene relationships, bioassay-compound relationships, compound similarities, compound co-occurrences in literature, and more nuanced interactions like inhibitor, activator, ligand, and other roles between compounds and genes.
- **Data Integration:** The knowledge graph schema is designed to support the integration of additional data sources, enhancing the depth and accuracy of the knowledge graph.
- **Command Line and Programmatic Access:** Provides flexibility in usage, allowing for integration into larger workflows or standalone analyses.

## Documentation

Full documentation is available at [ChemGraphBuilder Documentation](https://asmaa-a-abdelwahab.github.io/ChemGraphBuilder/).

## Contributing

Contributions are welcome! If any issues are found or suggestions for improvements arise, they can be reported via the [GitHub Issues](https://github.com/asmaa-a-abdelwahab/chemgraphbuilder/issues) page. Contributions to the codebase through pull requests are also encouraged.

## License

This project is licensed under the GPL-3.0 License - see the [LICENSE](https://github.com/asmaa-a-abdelwahab/ChemGraphBuilder#GPL-3.0-1-ov-file) file for details.

## Contact

For questions or support, please contact [Asmaa A. Abdelwahab](mailto:asmaa.a.abdelwahab@gmail.com).

## Acknowledgments

This project utilizes the PubChem Database and its API for accessing chemical and bioassay data. We acknowledge the efforts of the PubChem team for maintaining such a valuable resource.
