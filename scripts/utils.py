"""
Utility functions for parsing and data processing.
"""


def parse_nodes_elements(file_path):
    """
    Parse nodes and elements from a file.

    :param file_path: Path to the file containing node and element data
    :return: Tuple of (nodes dict, elements dict)
    """
    nodes = {}
    elements = {}
    with open(file_path, 'r') as f:
        lines = f.readlines()

    node_section = False
    element_section = False
    for line in lines:
        line = line.strip()
        if line.startswith('*Node'):
            node_section = True
            element_section = False
            continue
        if line.startswith('*Element'):
            node_section = False
            element_section = True
            continue
        if node_section and line:
            parts = line.split(',')
            node_id = int(parts[0])
            x, y = map(float, parts[1:3])
            nodes[node_id] = (x, y)
        if element_section and line:
            parts = line.split(',')
            element_id = int(parts[0])
            node_ids = list(map(int, parts[1:]))
            elements[element_id] = node_ids
    return nodes, elements
