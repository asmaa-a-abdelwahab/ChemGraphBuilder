# The `GraphNodesLoader` class is usable from Python code, and the command line interface.

### Usage in Different Contexts

#### 1. **Python Code**

```python
from chemgraphbuilder.graph_nodes_loader import GraphNodesLoader

# User-provided connection details and label
uri = "****"
username = "****"
password = "****"
label = "Compound"  # you can add other labels such as (BioAssay, Gene, Protein)

graph_nodes_loader = GraphNodesLoader(uri, username, password)
graph_nodes_loader.load_data_for_node_type(label)
```

#### 2. **Command Line**

After installing the package using `pip install .` or `pip install chemgraphbuilder`, you can use the command-line interface:

```sh
load-graph-nodes --uri **** --username **** --password **** --label Compound
```
