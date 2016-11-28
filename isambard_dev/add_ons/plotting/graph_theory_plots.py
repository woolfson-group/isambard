import networkx
import matplotlib.pyplot as plt


# TODO Completely re-write this function for dealing with protein.tags['graph'] objects.
def draw_graph(G):
    """ Plot a graph object returned from get_graph

    :param G: A graph, returned from get_graph
    :return:
    """
    e_p = []
    e_ap = []
    e_perp = []
    for e in G.edges(data=True):
        a = e[2]['angle']
        if a < 60:
            e_p.append(e)
        elif a > 120:
            e_ap.append(e)
        else:
            e_perp.append(e)
    plt.clf()
    fig, ax = plt.subplots()
    ax.set_xticklabels([])
    ax.set_xticks([])
    ax.set_yticklabels([])
    ax.set_yticks([])

    if G.number_of_nodes() > 10:
        pos = networkx.spring_layout(G)
    else:
        pos = networkx.circular_layout(G)
    labels = {}
    for i in G.nodes():
        labels[i] = i
    networkx.draw_networkx_nodes(G, pos, node_color='c')
    networkx.draw_networkx_labels(G, pos, labels)
    if len(e_p) > 0:
        networkx.draw_networkx_edges(G, pos, edgelist=e_p, width=2.0)
    if len(e_ap) > 0:
        networkx.draw_networkx_edges(G, pos, edgelist=e_ap, width=2.0, style='dotted')
    if len(e_perp) > 0:
        networkx.draw_networkx_edges(G, pos, edgelist=e_perp, width=2.0, style='dashed')

    line_labels = dict([((u, v), d['weight']) for u, v, d in G.edges(data=True)])
    networkx.draw_networkx_edge_labels(G, pos, edge_labels=line_labels, label_pos=0.3)

    return fig

__author__ = 'Jack W. Heal'
__status__ = 'Development'

