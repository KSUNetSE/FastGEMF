import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from sortedcontainers import SortedList

def plot_results(T, StateCount, compartments, line_styles=['-', '--', ':'], 
                font_size=12, font_family='serif', grid=True):
    
    plt.figure(figsize=(10, 6))
    num_compartments = len(compartments)
    colors = plt.cm.rainbow(np.linspace(0, 1, num_compartments))
    
    for i, compartment in enumerate(compartments):
        plt.plot(T, StateCount[i, :], color=colors[i], 
                 label=compartment)
    
    plt.xlabel('Time', fontsize=font_size, fontfamily=font_family)
    plt.ylabel('Number of Nodes', fontsize=font_size, fontfamily=font_family)
    plt.title('State Count Over Time', fontsize=font_size+2, fontfamily=font_family)
    plt.legend(fontsize=font_size-2)
    plt.grid(grid)
    plt.tight_layout()
    plt.show()

def extend_to_max_length(data_dict):
    sorted_times = SortedList((-len(data['T']), key) for key, data in data_dict.items())
    max_length,sim_no = -1*sorted_times[0][0],sorted_times[0][1]

    for sim, data in data_dict.items():
        T = data['T']
        statecount = data['statecount']


        if len(T) < max_length:
            T_extended = np.pad(T, (0, max_length - len(T)), mode='edge')
            statecount_extended = np.pad(statecount, ((0,0),(0, max_length - statecount.shape[1])), mode='edge')

            data_dict[sim]['T'] = T_extended
            data_dict[sim]['statecount'] = statecount_extended

    return data_dict, sim_no

def plot_multiple_results(results,  compartments, font_size=12, font_family='serif', grid=True):
    results, sim_no=extend_to_max_length(results)

    T = list(results.values())[sim_no]['T']
    state_counts = np.array([result['statecount'] for result in results.values()])
    num_compartments = len(compartments)
    colors = plt.cm.rainbow(np.linspace(0, 1, num_compartments))

    mean_statecount = np.mean(state_counts, axis=0)
    min_statecount  = np.min(state_counts, axis=0)
    max_statecount  = np.max(state_counts, axis=0)

    plt.figure(figsize=(10, 6))

    for i, compartment in enumerate(compartments):
        plt.plot(T, mean_statecount[i, :], color=colors[i], 
                 label=compartment)
        plt.fill_between(T, min_statecount[i, :], max_statecount[i, :],  color=colors[i], alpha=0.3)
        


    plt.xlabel('Time', fontsize=font_size, fontfamily=font_family)
    plt.ylabel('Population of States', fontsize=font_size, fontfamily=font_family)
    plt.legend(fontsize=font_size-2)
    plt.grid(grid)
    plt.tight_layout()
    plt.show()

    

def draw_model_graph(model):
    N = len(model.compartments)
    angle = 2*np.pi / N
    pos = {compartment: (np.cos(i * angle), np.sin(i * angle)) 
           for i, compartment in enumerate(model.compartments)}


    plt.figure(figsize=(5, 5))
    G_node = nx.MultiDiGraph()
    for compartment in model.compartments:
        G_node.add_node(compartment)
    colors = plt.cm.rainbow(np.linspace(0, 1, N))
    color_map = {compartment: color for compartment, color in zip(model.compartments, colors)}
    
    edge_curves = {}
    for nt in model.node_transitions:
        if (nt.from_state, nt.to_state) in edge_curves:
            edge_curves[(nt.from_state, nt.to_state)] += 0.1
        else:
            edge_curves[(nt.from_state, nt.to_state)] = 0.1
        G_node.add_edge(nt.from_state, nt.to_state, style='dashed', 
                        label=f"{nt.name}, rate: ({nt.rate})", 
                        connectionstyle=f'arc3,rad={edge_curves[(nt.from_state, nt.to_state)]}')
    
    nx.draw_networkx_nodes(G_node, pos, node_color=[color_map[node] for node in G_node.nodes()], node_size=3000, alpha=0.6)
    nx.draw_networkx_labels(G_node, pos, font_size=18, font_family="sans-serif")
    
    edge_labels = {(u, v): d['label'] for u, v, d in G_node.edges(data=True) if d['label']}
    for (u, v, d) in G_node.edges(data=True):
        nx.draw_networkx_edges(G_node, pos, edgelist=[(u, v)], width=3, alpha=0.7, 
                               edge_color='grey', style='dashed', 
                               connectionstyle=d['connectionstyle'])
    nx.draw_networkx_edge_labels(G_node, pos, edge_labels=edge_labels, font_color='blue')
    
    plt.title("Node-based Transitions")
    plt.axis('off')
    plt.margins(0.2)
    plt.tight_layout()
    plt.show()
    layers = set(et.network_layer for et in model.edge_transitions)
    for layer in layers:
        plt.figure(figsize=(6, 5))
        G_edge = nx.MultiDiGraph()
        for compartment in model.compartments:
            G_edge.add_node(compartment)
        
        edge_curves = {}
        for et in model.edge_transitions:
            if et.network_layer == layer:
                if (et.from_state, et.to_state) in edge_curves:
                    edge_curves[(et.from_state, et.to_state)] += 0.1
                else:
                    edge_curves[(et.from_state, et.to_state)] = 0.1
                G_edge.add_edge(et.from_state, et.to_state, style='solid', 
                                label=f"{et.name} (Inf.: {et.inducer}, rate: {et.rate})",
                                connectionstyle=f'arc3,rad={edge_curves[(et.from_state, et.to_state)]}')
        
        nx.draw_networkx_nodes(G_edge, pos, node_color=[color_map[node] for node in G_edge.nodes()], node_size=3000, alpha=0.6)
        nx.draw_networkx_labels(G_edge, pos, font_size=18, font_family="sans-serif")
        
        edge_labels = {(u, v): d['label'] for u, v, d in G_edge.edges(data=True) if d['label']}
        for (u, v, d) in G_edge.edges(data=True):
            nx.draw_networkx_edges(G_edge, pos, edgelist=[(u, v)], width=3, alpha=0.7, 
                                   edge_color='black', 
                                   connectionstyle=d['connectionstyle'])
        nx.draw_networkx_edge_labels(G_edge, pos, edge_labels=edge_labels, font_color='blue')
        
        plt.title(f"Edge-based Transitions - Layer: {layer}")
        plt.axis('off')
        plt.margins(0.2)
        plt.tight_layout()
        plt.show()