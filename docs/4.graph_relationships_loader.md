# The `GraphRelationshipsLoader` class is usable from Python code, and the command line interface.

### Usage in Different Contexts

#### **Available Relationships**
- `Assay_Gene`
- `Assay_Compound`
- `Compound_Gene`
- `Compound_Transformation`
- `Gene_Enzyme`
- `Compound_Similarities`
- `Cpd_Cpd_CoOccurence`
- `Cpd_Gene_CoOccurence`

#### 1. **Python Code**

```python
from chemgraphbuilder.graph_relationships_loader import GraphRelationshipsLoader

uri = "bolt://localhost:7687"
username = "neo4j"
password = "password"

loader = GraphRelationshipsLoader(uri, username, password)
loader.add_relationships("Assay_Gene")
loader.close()
```

#### 2. **Command Line**

After installing the package using `pip install .` or `pip install chemgraphbuilder`, you can use the command-line interface:

```sh
load-graph_relationships --uri bolt://localhost:7687 --username neo4j --password password --relationship_type Assay_Gene

```
