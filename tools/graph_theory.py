import networkx
from networkx.generators.atlas import graph_atlas_g
from networkx.generators import cycle_graph, path_graph
import shelve
import os
import json

from settings import global_settings


_two_core_path = os.path.join(global_settings['package_path'], 'tools', 'two_core_names.json')
if os.path.exists(_two_core_path):
    with open(_two_core_path, 'r') as inf:
        two_core_names = json.loads(inf.read())
else:
    two_core_names = None

atlas_graph_list = graph_atlas_g()


def get_unknown_graph_shelf():
    # unknown_graphs contains details of graphs (including the graph objects themselves)
    #  in atlas table (or to be added) that are not in list_of_graphs()
    unknown_graph_shelf = os.path.join(global_settings['package_path'], 'tools', 'unknown_graphs')
    unknown_graph_shelf_filepath = '{0}.db'.format(unknown_graph_shelf)
    if not os.path.exists(unknown_graph_shelf_filepath):
        message = 'File not found: {0}\n'.format(unknown_graph_shelf_filepath)
        message += 'File should have unknown_graph names as keys, and dictionaries of graph information as values'
        raise IOError(message)
    else:
        return unknown_graph_shelf


def list_of_graphs(unknown_graphs=False):
    """

    Parameters
    ----------
    unknown_graphs : bool
        If True, appends the graphs from the unknown_graphs.db shelf.

    Returns
    -------
    graph_list : list(networkx.Graph)
        List of graph objects with .name attributes corresponding to their names in the Atlas of Graphs.
        These names are used in the coeus database table AtlasDB.

    """
    # atlas graphs
    graph_list = atlas_graph_list
    # cycle and path graphs
    max_n = 100
    for n in range(8, max_n + 1):
        c = cycle_graph(n)
        p = path_graph(n)
        c.name = 'C{0}'.format(n)
        p.name = 'P{0}'.format(n)
        graph_list.append(c)
        graph_list.append(p)
    # Add other custom graphs here.
    # If unknown_graphs, append them to the list from shelf.
    if unknown_graphs:
        graph_list += get_unknown_graph_list()
    return graph_list


def graph_to_plain_graph(g):
    # construct h fully in case of unorderable/unsortable edges.
    h = networkx.Graph()
    h.add_nodes_from(range(len(g.nodes())))
    edges = [(g.nodes().index(e1), g.nodes().index(e2)) for e1, e2 in g.edges()]
    h.add_edges_from(edges)
    return h


def get_unknown_graph_list():
    unknown_graph_shelf = get_unknown_graph_shelf()
    s = shelve.open(unknown_graph_shelf)
    unknown_graph_list = list(s.values())
    s.close()
    return unknown_graph_list


def isomorphism_checker(g, graph_list=None):
    """ Finds the name of the graph by checking for isomorphism with the graphs in the graph_list.

    Notes
    -----
    If g is a MultiGraph, a DiGraph or a MultiDiGraph, the isomorphism is run against the underlying Graph.
    In other words, all the edge annotations (including directions) are ignored.

    Parameters
    ----------
    g : A networkx Graph.
    graph_list :
        A list of networkx.Graph objects against which to check for isomorphism.
        If the graphs in the list do not have names, then their index will be returned as a string

    Returns
    -------
    iso_name : str, or None
        The name of the graph from the graph_list that is isomorphic to g.
        If the graph does not have a 'name' attribute, its index in the list will be returned as a string.
        None if no such isomorph is found in the list_of_graphs.

    """
    # if no graph_list is supplied, run list_of_graphs for the Graph Atlas list.
    if graph_list is None:
        graph_list = list_of_graphs()
    # run isomorphism check on the Graph of g (i.e. not the DiGraph or MultiGraph).
    h = graph_to_plain_graph(g)
    isomorph = next(filter(lambda x: networkx.is_isomorphic(h, x), graph_list), None)
    if isomorph is None:
        return
    else:
        iso_name = isomorph.name
        if iso_name is None:
            iso_name = str(graph_list.index(isomorph))
    return iso_name


def sorted_connected_components(g, include_trivials=False):
    """ List of connected component subgraphs of graph g, ordered in decreasing number of nodes.

    Parameters
    ----------
    g : networkx.Graph
    include_trivials : bool
        If True, trivial connected components (i.e. singular nodes) will be included.

    Returns
    -------
    [networkx.Graph]
        List of connected component subgraphs.
    """
    h = graph_to_plain_graph(g)
    components = sorted(networkx.connected_component_subgraphs(h, copy=False),
                        key=lambda x: len(x.nodes()), reverse=True)
    if not include_trivials:
        components = [x for x in components if len(x.nodes()) > 1]
    return components


def get_graph_name(g):
    # construct h fully in case of unorderable/unsortable edges.
    h = graph_to_plain_graph(g)
    name = isomorphism_checker(h)
    if name is None:
        unknown_graph_list = get_unknown_graph_list()
        name = isomorphism_checker(h, graph_list=unknown_graph_list)
        if name is None:
            name = 'U{0}'.format(len(unknown_graph_list) + 1)
    return name


def store_graph(g):
    """ Good practice to use this function over _add_graph_to_shelf as it deals with two_cores properly. """
    storage_changed = False
    h = graph_to_plain_graph(g)
    two_core = networkx.k_core(h, k=2)
    # add_two_core first
    added_twocore = _add_graph_to_shelf(two_core)
    added_graph = _add_graph_to_shelf(h)
    if added_graph:
        storage_changed = True
        two_core_name = get_graph_name(two_core)
        graph_name = get_graph_name(h)
        add_two_core_name_to_json(graph_name=graph_name, two_core_name=two_core_name)
        if added_twocore:
            add_two_core_name_to_json(graph_name=two_core_name, two_core_name=two_core_name)
    return storage_changed


def _add_graph_to_shelf(g):
    """ Runs isomorphism checker against stored dictionary of non-Atlas graphs.

    Notes
    -----
    The dictionary of non-Atlas graphs is stored as a .db (output from shelve) file.
        Keys are unknown graph names.
        Values are dictionaries storing the unknown graphs themselves as well as some basic properties of the graphs.
    If graph is not in Atlas, nor in this dictionary of non-Atlas graphs, then it is added to this dictionary.
    Unknown graphs are named 'Un', where n is the order in which they have been encountered.
    The number n does not tell you anything about the properties of the graph.

    Parameters
    ----------
    g : A networkx Graph.

    Returns
    -------
    added : bool
        True if graph was added to shelf
    """
    added = False
    h = graph_to_plain_graph(g)
    name = get_graph_name(h)
    if name[0] == 'U':
        unknown_graph_shelf = get_unknown_graph_shelf()
        s = shelve.open(unknown_graph_shelf)
        if name not in s.keys():
            h.name = name
            s[name] = h
            added = True
        s.close()
    return added


def add_two_core_name_to_json(graph_name, two_core_name, force_add=False):
    """ Add an amino acid to the two_core_names.json file used to populate the atlas table.

    Parameters
    ----------
    graph_name : str
        Name of graph from graph Atlas, e.g. 'G163'.
    two_core_name : str
        Name of the two core of that corresponds to to graph_name.
    force_add : bool, optional
        If True, will over-write existing dictionary value for graph_name if already in two_core_names.json.
        If False, then an IOError is raised if graph_name is already in two_core_names.json.

    Raises
    ------
    IOError
        If code is already in two_core_names.json and force_add is False.

    Returns
    -------
    None

    """
    if not two_core_names:
        return
    # If code is already in the dictionary, raise an error
    if (not force_add) and graph_name in two_core_names.keys():
        raise IOError("{0} is already in the two_core_names dictionary, with two core name: {1}"
                      .format(graph_name, two_core_names[graph_name]))
    two_core_names[graph_name] = two_core_name
    # Write over json file with updated dictionary.
    with open(_two_core_path, 'w') as foo:
        foo.write(json.dumps(two_core_names))
    return


__author__ = 'Jack W. Heal'
