
import logging
import getpass
import os
from neo4j import GraphDatabase

class Neo4jConnectionError(Exception):
    """Custom exception for Neo4j connection errors."""
    pass

class Neo4jBase:
    """
    Base class to manage connections with the Neo4j database.

    Attributes:
    - uri: The connection URI for the Neo4j database.
    - user: The username to use for authentication.
    - driver: The driver object used to interact with the Neo4j database.

    Methods:
    - connect_to_neo4j: Establish a connection to the Neo4j database.
    - close: Close the connection to the Neo4j database.
    """

    def __init__(self, logger=None, uri="neo4j+s://f8875ed1.databases.neo4j.io:7687", user="neo4j"):
        self.uri = uri
        self.user = user
        self.driver = None

        # Set up logging configuration
        logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
        self.logger = logging.getLogger(__name__)

    def connect_to_neo4j(self):
        """Establish a connection to the Neo4j database using provided URI and username."""
        password = os.getenv("NEO4J_PASSWORD")  # Check if password is set in environment variables
        if not password:
            password = getpass.getpass(prompt="Enter Neo4j password: ")

        try:
            self.driver = GraphDatabase.driver(self.uri, auth=(self.user, password))
            self.logger.info("Successfully connected to the Neo4j database.")
        except Exception as e:
            self.logger.error(f"Failed to connect to the Neo4j database: {e}")
            raise Neo4jConnectionError("Failed to connect to the Neo4j database.") from e

    def close(self):
        """Close the connection to the Neo4j database."""
        if self.driver:
            self.driver.close()
            self.logger.info("Neo4j connection closed successfully.")
