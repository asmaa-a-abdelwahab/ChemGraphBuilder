"""
Tests for the AddGraphNodes class from the add_graph_nodes module.
"""

import unittest
from unittest.mock import patch, MagicMock
import os
import sys
import pandas as pd

# Add the parent directory to the system path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from chemgraphbuilder.add_graph_nodes import AddGraphNodes

class TestAddGraphNodes(unittest.TestCase):
    """
    Test case for the AddGraphNodes class.
    """

    @patch('chemgraphbuilder.add_graph_nodes.Neo4jBase.__init__')
    @patch('chemgraphbuilder.add_graph_nodes.GraphDatabase.driver')
    def setUp(self, mock_driver, _mock_neo4j_base):
        """
        Set up the test case with a mocked Neo4j driver.
        """
        self.mock_driver = mock_driver
        self.mock_driver.session = MagicMock()
        self.add_graph_nodes = AddGraphNodes(driver=self.mock_driver)

    @patch('chemgraphbuilder.add_graph_nodes.GraphDatabase.driver')
    def test_create_uniqueness_constraint(self, mock_driver):
        """
        Test the create_uniqueness_constraint method.
        """
        mock_session = MagicMock()
        mock_driver.session.return_value.__enter__.return_value = mock_session

        AddGraphNodes.create_uniqueness_constraint(
            mock_driver, 'TestLabel', 'test_property')

        mock_session.run.assert_called_once_with(
            'CREATE CONSTRAINT IF NOT EXISTS FOR (n:TestLabel) '
            'REQUIRE n.test_property IS UNIQUE'
        )

    def test_generate_property_string(self):
        """
        Test the _generate_property_string method.
        """
        self.assertEqual(AddGraphNodes._generate_property_string(123), 123)
        self.assertEqual(AddGraphNodes._generate_property_string(123.456), 123.456)
        self.assertEqual(AddGraphNodes._generate_property_string('text'), "'text'")
        self.assertEqual(
            AddGraphNodes._generate_property_string("text'with\nspecial"),
            "'text\\'with\\nspecial'"
        )

    @patch('chemgraphbuilder.add_graph_nodes.pd.read_csv')
    def test_read_csv_file(self, mock_read_csv):
        """
        Test the read_csv_file method.
        """
        mock_df = pd.DataFrame({
            'unique_id': [1, 2],
            'property1': ['value1', 'value2'],
            'property2': [10.0, 20.0]
        })
        mock_read_csv.return_value = mock_df

        result = self.add_graph_nodes.read_csv_file('dummy_path.csv', 'unique_id')
        expected = {
            1: {'property1': 'value1', 'property2': 10.0},
            2: {'property1': 'value2', 'property2': 20.0}
        }
        self.assertEqual(result, expected)

    @patch('chemgraphbuilder.add_graph_nodes.pd.read_csv')
    def test_combine_csv_files(self, mock_read_csv):
        """
        Test the combine_csv_files method.
        """
        mock_df1 = pd.DataFrame({
            'unique_id': [1],
            'property1': ['value1'],
            'property2': [10.0]
        })
        mock_df2 = pd.DataFrame({
            'unique_id': [2],
            'property1': ['value2'],
            'property2': [20.0]
        })
        mock_read_csv.side_effect = [mock_df1, mock_df2]

        with patch('os.listdir', return_value=['file1.csv', 'file2.csv']):
            result = self.add_graph_nodes.combine_csv_files('dummy_directory')

        expected = pd.concat([mock_df1, mock_df2], ignore_index=True)
        pd.testing.assert_frame_equal(result, expected)

    @patch.object(AddGraphNodes, '_generate_property_string', side_effect=lambda x: f'processed_{x}')
    def test_generate_cypher_queries(self, _mock_generate_property_string):
        """
        Test the generate_cypher_queries method.
        """
        node_dict = {
            1: {'property1': 'value1', 'property2': 10.0},
            2: {'property1': 'value2', 'property2': 20.0}
        }
        queries = list(self.add_graph_nodes.generate_cypher_queries(
            node_dict, 'TestLabel', 'unique_id'))

        expected_queries = [
            'MERGE (n:TestLabel {unique_id: 1}) SET n.property1 = processed_value1, '
            'n.property2 = processed_10.0',
            'MERGE (n:TestLabel {unique_id: 2}) SET n.property1 = processed_value2, '
            'n.property2 = processed_20.0'
        ]

        self.assertEqual(queries, expected_queries)

    @patch.object(AddGraphNodes, 'execute_queries')
    @patch.object(AddGraphNodes, 'generate_cypher_queries')
    @patch.object(AddGraphNodes, 'read_csv_file')
    def test_process_and_add_nodes(self, mock_read_csv_file, mock_generate_cypher_queries, mock_execute_queries):
        """
        Test the process_and_add_nodes method.
        """
        mock_read_csv_file.return_value = {1: {'property1': 'value1', 'property2': 10.0}}
        mock_generate_cypher_queries.return_value = ['dummy_query']

        self.add_graph_nodes.process_and_add_nodes(
            'dummy_path.csv', 'TestLabel', 'unique_id')

        mock_read_csv_file.assert_called_once_with('dummy_path.csv', 'unique_id')
        mock_generate_cypher_queries.assert_called_once_with(
            {1: {'property1': 'value1', 'property2': 10.0}}, 'TestLabel', 'unique_id'
        )
        mock_execute_queries.assert_called_once_with(['dummy_query'])

    @patch.object(AddGraphNodes, 'process_and_add_nodes')
    @patch.object(AddGraphNodes, 'combine_csv_files')
    @patch('chemgraphbuilder.add_graph_nodes.os.remove')
    @patch('chemgraphbuilder.add_graph_nodes.pd.DataFrame.to_csv')
    def test_process_and_add_nodes_from_directory(
        self, mock_to_csv, mock_os_remove, mock_combine_csv_files, mock_process_and_add_nodes):
        """
        Test the process_and_add_nodes_from_directory method.
        """
        mock_df = pd.DataFrame({'unique_id': [1], 'property1': ['value1'], 'property2': [10.0]})
        mock_combine_csv_files.return_value = mock_df

        self.add_graph_nodes.process_and_add_nodes_from_directory(
            'dummy_directory', 'TestLabel', 'unique_id')

        mock_combine_csv_files.assert_called_once_with('dummy_directory')
        mock_to_csv.assert_called_once_with(
            os.path.join('dummy_directory', 'combined_temp.csv'), index=False)
        mock_process_and_add_nodes.assert_called_once_with(
            os.path.join('dummy_directory', 'combined_temp.csv'), 'TestLabel', 'unique_id')
        mock_os_remove.assert_called_once_with(
            os.path.join('dummy_directory', 'combined_temp.csv'))

if __name__ == '__main__':
    unittest.main()
