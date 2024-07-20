from setuptools import setup, find_packages

setup(
    name='chemgraphbuilder',
    version='0.1',
    packages=find_packages(),
    install_requires=[
        'neo4j',
        'pandas'
    ],
    entry_points={
        'console_scripts': [
            'load-graph-nodes=chemgraphbuilder.graph_nodes_loader:main',
            'load-graph-relationships=chemgraphbuilder.graph_relationships_loader:main',
        ],
    },
)
