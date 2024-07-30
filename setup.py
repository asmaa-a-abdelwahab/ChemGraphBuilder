from setuptools import setup, find_packages

setup(
    name='chemgraphbuilder',
    version='1.9.0',
    packages=find_packages(),
    install_requires=[
        "beautifulsoup4",
        "matplotlib",
        "numpy",
        "pandas",
        "neo4j",
        "requests",
        "xmltodict",
        "setuptools",
        "dask",
    ],
    author='Asmaa Ali Abdelwahab',
    author_email='asmaa.a.abdelwahab@gmail.com',
    description='A package to build chemical knowledge graphs using data from PubChem and Neo4j',
    url='https://github.com/asmaa-a-abdelwahab/ChemGraphBuilder',
    entry_points={
        'console_scripts': [
            'setup-data-folder=chemgraphbuilder.setup_data_folder:main',
            'collect-process-nodes=chemgraphbuilder.node_collector_processor:main',
            'collect-process-relationships=chemgraphbuilder.relationship_collector_processor:main',
            'load-graph-nodes=chemgraphbuilder.graph_nodes_loader:main',
            'load-graph-relationships=chemgraphbuilder.graph_relationships_loader:main',
        ],
    },
)
