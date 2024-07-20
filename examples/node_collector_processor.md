You can use the `NodesCollectorProcessor` class in both Python and the command line:

### Python Script Example

You can use the `NodesCollectorProcessor` class within a Python script as follows:

```python
from chemgraphbuilder.nodes_collector_processor import NodesCollectorProcessor
from chemgraphbuilder.setup_data_folder import SetupDataFolder

# Set the connection details and parameters
uri = "bolt://localhost:7689"
username = "neo4j"
password = "your_password"
node_type = "Compound"  # Change to "BioAssay", "Gene", or "Protein" as needed
enzyme_list = ['CYP2D6', 'CYP3A4']

# Initialize and setup the data directory
setup_folder = SetupDataFolder()
setup_folder.setup()

# Initialize the collector
collector = NodesCollectorProcessor(uri=uri, username=username, password=password, node_type=node_type, enzyme_list=enzyme_list)

# Collect and process the data
collector.collect_and_process_data()

# Close the connection
collector.close()
```


### Command Line Interface (CLI) Example

You can use the `NodesCollectorProcessor` class from the command line by executing the script with the necessary arguments:

```sh
setup-data-folder 
collect-process-nodes --uri bolt://localhost:7689 --username neo4j --password your_password --node_type Compound --enzyme_list CYP2D6,CYP3A4
```

Make sure to replace `"your_password"` with the actual password for your Neo4j database. The `node_type` argument can be one of `"Compound"`, `"BioAssay"`, `"Gene"`, or `"Protein"`, depending on the type of data you want to collect. The `enzyme_list` should be a comma-separated string of enzyme names.
