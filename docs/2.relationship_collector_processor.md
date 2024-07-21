# Relationship Collector and Processor

This module collects and processes relationship data for different types of relationships using the `RelationshipPropertiesExtractor` and `RelationshipDataProcessor` classes.

## Available Relationship Types

- Assay_Compound
- Assay_Enzyme
- Gene_Enzyme
- Compound_Gene
- Compound_Similarity
- Compound_Cooccurrence
- Compound_Transformation

## Usage Examples

### Python

To use the `RelationshipCollectorProcessor` class in Python, follow the example below:

```python
from relationship_collector_processor import RelationshipCollectorProcessor

# Set the relationship type
relationship_type = "Assay_Compound"  # Change to the desired relationship type

# Initialize the collector
collector = RelationshipCollectorProcessor(relationship_type=relationship_type)

# Collect and process the relationship data
collector.collect_relationship_data()
```

### Command Line

To collect and process relationship data from the command line, use the following command:

```sh
# Collect and process the relationship data for the specified relationship type
collect-process-relationships --relationship_type Assay_Compound
```

Replace `Assay_Compound` with any of the available relationship types to collect different data.
