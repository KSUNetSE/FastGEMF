
import scipy.sparse as sp
import networkx as nx
import time
from .error import StopLoopException 
import numpy as np
from sortedcontainers import SortedList
from .times_structure import TimeSorted,TimeNp
import numpy as np

def sample(rates):
    """
    Draw a random variable from a discrete distribution defined by 'probability'

    Parameters:
        probability (list or np.array): A list or array of probabilities associated 
                                        with each state. 

    """
    sum_rates=np.sum(rates)
    if sum_rates>1e-6:
        probability=rates/sum_rates
        return np.searchsorted(np.cumsum(probability), np.random.random())
    else: 
        raise StopLoopException("Infection Is Over! All infected Are Removed!")

def find_key_for_value(dictionary, value):
        for key, val in dictionary.items():
            if val == value:
                return key
        return None

def invert_influence_dict(influence_dict):
    inverted_dict = {}
    for influencer, layers in influence_dict.items():
        for layer in layers:
            if layer not in inverted_dict:
                inverted_dict[layer] = []
            inverted_dict[layer].append(influencer)
    return inverted_dict


def setup_simulation_matrices(cfg):
    (cfg.network_layers).sort()
    layer_indices = {layer: index for index, layer in enumerate(cfg.network_layers)}

    M = len(cfg.compartments)
    L = len(cfg.network_layers)
    q=[]
    #compartments=sorted(cfg.compartments)
    compartments=cfg.compartments
    infuenced_layers={compartments.index(t.inducer):[] for t  in cfg.edge_transitions_valued}
    layers_influencer={}
    for t in cfg.edge_transitions_valued:
        l=[]
        if t.inducer in compartments :
            l=layer_indices[t.network_layer]
            if l not in infuenced_layers[compartments.index(t.inducer)]:
                infuenced_layers[compartments.index(t.inducer)].append(l)
            if compartments.index(t.inducer) not in q:
                 q.append(compartments.index(t.inducer))

    layers_influencer=invert_influence_dict(infuenced_layers)
    A_delta = np.zeros((M, M), dtype=object)  # Using dtype=object to store rates as strings
    for t in cfg.node_transitions_valued:
        i = compartments.index(t.from_state)
        j = compartments.index(t.to_state)
        A_delta[i, j] = t.rate
    
    A_beta = np.zeros((L, M, M), dtype=object)  # 3D array for edge transitions
    for t in cfg.edge_transitions_valued:
        l = layer_indices[t.network_layer]
        i = compartments.index(t.from_state)
        j = compartments.index(t.to_state)
        A_beta[l, i, j] = t.rate

    return {
        #'q': np.array(q),
        'q':infuenced_layers,
        'A_delta': A_delta,
        'A_beta': A_beta,
        'M': M,
        'L': L,
        'layers_influencer':layers_influencer
    }

def extract_transition_parameters(arr):
    result = {}
    dim=arr.ndim
    if dim==2:
        for i in range(arr.shape[0]):
            for j in range(arr.shape[1]):
                if isinstance(arr[i, j], str):
                    result[arr[i, j]] = (i, j)
        return result
    
    elif dim==3:
        for l in range(arr.shape[0]):
            for i in range(arr.shape[1]):
                for j in range(arr.shape[1]):
                    if isinstance(arr[l,i, j], str):
                        result[arr[l,i, j]] = (l,i, j)
        return result
    

def pick_fastest_method(nodes, edges):


    result=1/(1+np.exp(-(nodes-120_000)))
    
    if result <=.5:
        return TimeNp.generate_times
    elif result > .5:
        return TimeSorted.generate_times

