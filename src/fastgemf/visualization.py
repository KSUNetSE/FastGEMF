import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from sortedcontainers import SortedList
import matplotlib
import sys
import time
import datetime
import matplotlib.pyplot as plt
import matplotlib
from scipy.interpolate import interp1d
from os.path import join as ospj
import os
from typing import Literal  
matplotlib.rcParams['backend'] = 'Agg'  # Default to Agg
if 'matplotlib.backends' in sys.modules:
    import importlib
    importlib.reload(matplotlib.backends)
matplotlib.use('Agg', force=True)  # Force Agg at the start
def choose_matplotlib_backend():
    """
    Choose backend: Agg for non-interactive (e.g., exec()), TkAgg only for confirmed GUI.
    """
    is_interactive = hasattr(sys, 'ps1') or 'IPYTHON' in sys.argv
    
    is_script = sys.argv[0] != ''  # empty argv[0] means interactive shell
    has_display = 'DISPLAY' in os.environ and os.environ['DISPLAY'] != ''
    
    #print(f"Interactive: {is_interactive}, Script: {is_script}, Display: {has_display}")
    
    if not is_interactive and is_script:
       #print("Detected script-like context (e.g., exec()), using Agg")
        return 'Agg'
    
    # TkAgg if GUI interactive and GUI seems available
    if is_interactive and has_display:
        try:
            import tkinter
            tk_root = tkinter.Tk()
            tk_root.withdraw()
            tk_root.update()
            tk_root.destroy()
            return 'TkAgg'
        except (ImportError, tkinter.TclError) as e:
            return 'Agg'
    
    return 'Agg'
# Set backend dynamically, but only switch to TkAgg if GUI is confirmed
backend = choose_matplotlib_backend()
matplotlib.use(backend, force=True)


def plot_results(T, StateCount, compartments, line_styles=['-', '--', ':'], 
                font_size=12,
                font_family='serif',
                show_figure=True,
                save_figure=True,
                title="",
                save_path= ospj(os.getcwd(),"results_figure.pdf"),
                y_axis="fraction",
                grid=True):
    title= title if title is not "" else "State Count Over Time"
    plt.figure(figsize=(10, 6))
    num_compartments = len(compartments)
    colors = plt.cm.rainbow(np.linspace(0, 1, num_compartments))
    N=int( np.sum(StateCount[:,0],axis=0))
    if y_axis=="fraction":
        StateCount=StateCount/N        
    for i, compartment in enumerate(compartments):
        plt.plot(T, StateCount[i, :], color=colors[i], 
                 label=compartment)
    
    plt.xlabel('Time', fontsize=font_size, fontfamily=font_family)
    plt.ylabel(f'{y_axis } of states (Number of Nodes: {N})', fontsize=font_size, fontfamily=font_family)
    plt.title(title, fontsize=font_size+2, fontfamily=font_family)
    plt.legend(fontsize=font_size-2)
    plt.grid(grid)
    plt.tight_layout()
    if save_figure:
        plt.savefig(save_path, dpi=600, bbox_inches='tight')
    if show_figure and matplotlib.get_backend() == 'TkAgg':
        plt.show(block=True)
    elif show_figure:
        print("Warning: Interactive display not available with current backend. Figure saved instead.")
    
    plt.close('all')  # Always clean up

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

def plot_minmax_results(results,  compartments, font_size=12, font_family='serif',show_figure=False,save_figure=True, title="",
                          save_path= ospj(os.getcwd(),"simulation_figure.pdf"),grid=True,y_axis="fraction"):#y_axis="population"
    
    
    T = [np.array(result['T']) for result in results.values()]
    state_counts = [np.array(result['statecount']) for result in results.values()]
    num_compartments = len(compartments)
    title= title+"\n(Solid = Mean, Shaded = Min-Max Shading)" if title is not "" else "State Count Over Time\n(Solid = Mean, Shaded = Min-Max Shading)"
    colors = plt.cm.rainbow(np.linspace(0, 1, num_compartments))

    # determine the common time grid based on the longest simulation
    max_length = max(len(t) for t in T)  
    final_time = max(t[-1] for t in T)   
    common_time = np.linspace(0, final_time, max_length)

    # interpolate all simulations to the common time grid
    interpolated_state_counts = []
    for i in range(len(results)):
        interp_funcs = []
        for j in range(num_compartments):
            # Interpolate each compartment's state count
            interp_func = interp1d(T[i], state_counts[i][j, :], kind='previous', fill_value="extrapolate")
            interp_funcs.append(interp_func(common_time))
        interpolated_state_counts.append(np.array(interp_funcs))  # Shape: (num_compartments, max_length)

    interpolated_state_counts = np.array(interpolated_state_counts)  # Shape: (num_sims, num_comparts, max_length)
    if y_axis=="fraction":
        N=int( np.sum(state_counts[0],axis=0)[0])
        interpolated_state_counts=interpolated_state_counts/N
    mean_state_counts = np.mean(interpolated_state_counts, axis=0)  # Shape: (num_comparts, max_length)
    min_state_counts=np.min(interpolated_state_counts, axis=0)
    max_state_counts = np.max(interpolated_state_counts, axis=0)    # Shape: (num_comparts, max_length)

    plt.figure(figsize=(10, 6), dpi=300)
        
    #  individual simulation trajectories (faint lines)
    """for i in range(len(results)):
        for j in range(num_compartments):
            plt.plot(
                common_time,
                interpolated_state_counts[i][j, :],
                color=colors[j],
                alpha=0.3,  # Lower alpha for individual runs
                label=compartments[j] if i == 0 else ""  # Label only on first iteration
            )
    """
    #  mean as solid lines and shaded region for mean-to-max
    for j in range(num_compartments):
        plt.plot(
            common_time,
            mean_state_counts[j, :],
            color=colors[j],
            lw=.5,  # thick line for mean
            alpha=1,    
            
        )
        plt.fill_between(
            common_time,
            min_state_counts[j, :],
            max_state_counts[j, :],
            color=colors[j],
            label=f"{compartments[j]}",
            alpha=0.3,
            #=f"{compartments[j]} (Mean-Max)"
        )
    grid and plt.grid(True)
    plt.title(title, fontsize=font_size+2, fontfamily= font_family)
    plt.xlabel('Time', fontsize=font_size, fontfamily=font_family)
    plt.ylabel(f'{y_axis } of States (Number of Nodes: {N})', fontsize=font_size, fontfamily=font_family)
    plt.legend(fontsize=font_size)
    if save_figure:
        plt.savefig(save_path, dpi=600, bbox_inches='tight')
    if show_figure and matplotlib.get_backend() == 'TkAgg':
        plt.show(block=True)
    elif show_figure and  not save_figure:
        print("Warning: Interactive display not available with current backend. you can save the figure insted, just pass kwarg: save_figure:True, and save_path:your\path\here ")
    elif show_figure and save_figure:
        print("Warning: Interactive display not available with current backend. Figure saved instead at directory: ",save_path)
        

def plot_shaded_results(
        times,
        mean_statecounts,
        statecounts_variations,
        compartments,
        variation_type="range",
        font_size=12,
        font_family='serif',
        show_figure=False,
        save_figure=True,
        save_path=None,
        title="",
        grid=True,
        y_axis="fraction"
    ):


    plt.figure(figsize=(10, 6), dpi=300)

    if variation_type == "iqr":
        variation_label = "IQR (25–75%)"
    elif variation_type == "90ci":
        variation_label = "90% CI (5–95%)"
    elif variation_type == "std":
        variation_label = "Mean ± 1 STD"
    elif variation_type == "range":
        variation_label = "Range (Min–Max)"
    else:
        variation_label = ""

    if title.strip() == "":
        title = f"State Count Over Time\n(Solid = Mean, Shaded = {variation_label})"
    else:
        title = f"{title}\n(Solid = Mean, Shaded = {variation_label})"

    if save_path is None:
        save_path = ospj(os.getcwd(), "simulation_figure.pdf")

    num_compartments = len(compartments)
    N = int(np.sum(mean_statecounts, axis=0)[0])

    colors = plt.cm.rainbow(np.linspace(0, 1, num_compartments))

    for i, compartment in enumerate(compartments):

        mean = mean_statecounts[i, :]
        lower = statecounts_variations[0, i, :]
        upper = statecounts_variations[1, i, :]

        if y_axis == "fraction":
            mean = mean / N
            lower = lower / N
            upper = upper / N

        plt.plot(times, mean, color=colors[i], label=f"{compartment}")

        plt.fill_between(times, lower, upper, color=colors[i], alpha=0.3)

    if grid:
        plt.grid(True)

    plt.xlabel('Time', fontsize=font_size, fontfamily=font_family)
    plt.ylabel(f'{y_axis} of States (Number of Nodes: {N})',
               fontsize=font_size, fontfamily=font_family)

    plt.title(title, fontsize=font_size + 2, fontfamily=font_family)
    plt.legend(fontsize=font_size - 2)
    
    
    if save_figure:
        plt.savefig(save_path, dpi=600, bbox_inches='tight')
    if show_figure and matplotlib.get_backend() == 'TkAgg':
        plt.show(block=True)
    elif show_figure and  not save_figure:
        print("Warning: Interactive display not available with current backend. you can save the figure insted, just pass kwarg: save_figure:True, and save_path:your\path\here ")
    elif show_figure and save_figure:
        print("Warning: Interactive display not available with current backend. Figure saved instead at directory: ",save_path)
    
    
    
    
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