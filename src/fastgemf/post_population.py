import numpy as np
import copy
import pandas as pd
from typing import Dict, Tuple
from scipy.interpolate import interp1d
from typing import Literal
def post_population(x0, model_matrices,event_data, nodes ):
    """
    The post_population caclulcates the evolution of state through the total time of simulation.

    Parameters:
    -----------
    x0 : np.ndarray
        The initial state of network

    model_matrices : ModelMatrices
        A_beta, A_delta, q, b_l, b_il, etc

    event_data : EventData
        Contains the information about event times(min_time) and state transitions (states_k and states_k_plus_1) and sampled nodes in each iteration 

    nodes : int
        Total number of nodes in the network.

    Returns:
    --------
    time : np.array
        absolute times

    StateCount : np.ndarray
        An M x T array (M compartments, T time steps), where each element represents the number of
        nodes in each compartment at each time step.
    ts : np.ndarray
       interarrival times

    states_k : np.ndarray
        The compartment states before each event (states_k).

    states_k_plus_1 : np.ndarray
        The compartment states after each event (states_k_plus_1).

    """
    M=model_matrices.M
    N=nodes
    time=copy.deepcopy((event_data.min_times))
    ts=np.array(time)
    time.insert(0,0)
    states_k=event_data.states_k
    states_k_plus_1=event_data.states_k_plus_1

    ts[1:]=ts[1:]-ts[:-1]
    X0 = np.zeros((M, N))
    x0=x0.astype(int)
    col_i = np.arange(N) 
    np.add.at(X0,(x0,col_i),1)
    StateCount = np.zeros((M, len(ts) + 1))
    StateCount[:, 0] = X0.sum(axis=1)
    DX = np.zeros((M, len(ts)))
    np.subtract.at(DX, (states_k, np.arange(len(ts))), 1)
    np.add.at(DX, (states_k_plus_1, np.arange(len(ts))), 1)
    StateCount[:, 1:] = np.cumsum(DX, axis=1) + StateCount[:, 0][:, np.newaxis]
    return time, StateCount[:,:],  ts,event_data.sampled_nodes, states_k, states_k_plus_1



def final_results_stats(times: Dict[int, np.ndarray],
                        state_counts: Dict[int, np.ndarray]
                        ,number_of_compartments: int
                        ,variation_type: Literal["iqr", "90ci", "std", "range"] = "90ci"
                       ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, pd.DataFrame]:
    """
    Process multiple stochastic simulations with different event-based time lengths.

    Parameters
    times : dict[int -> np.ndarray]
        Dictionary where each entry is the time array of a simulation.
        Example: times[i] = array([t0, t1, ..., tK_i])

    state_counts : dict[int -> np.ndarray]
        Dictionary of state count matrices.
        state_counts[i] is shape (M, K_i+1)
        where M = number_of_compartments
        
    variation_type : str
        Type of variation to compute across simulations:    "iqr", "90ci", "std", "range"

    Returns

    unified_time : np.ndarray
        The union of all time points across simulations (sorted, unique).

    state_mean : np.ndarray
        Array of shape (M, len(unified_time)) containing the mean across simulations.

    state_std : np.ndarray
        Array of shape (M, len(unified_time)) containing the std deviation across simulations.

    full_results_df : pd.DataFrame
        A DataFrame where each row corresponds to a simulation,
        containing: simulation id, full time array, full state matrix.

    Notes: Different simulations produce different number of events.
    To average them correctly, we  interpolate all state counts
    to a common time grid.
    """

    sim_ids = sorted(times.keys())
    num_sims = len(sim_ids)
    all_times = np.unique(np.concatenate([times[i] for i in sim_ids]))
    unified_time = np.sort(all_times)
    M = number_of_compartments
    state_matrix_all = np.zeros((num_sims, M, len(unified_time)))

    for idx, sim in enumerate(sim_ids):
        t = times[sim]
        X = state_counts[sim]  # shape (M, len(t))

        for m in range(M):
            interp_func = interp1d(t, X[m], kind='previous', fill_value="extrapolate")
            state_matrix_all[idx, m] = interp_func(unified_time)
    statecounts_mean = np.mean(state_matrix_all, axis=0)        
    if variation_type == "iqr": # interquartile range 25th to 75th percentile
        lower_percentile = 25
        upper_percentile = 75       
        statecounts_variations = np.percentile(state_matrix_all, [lower_percentile, upper_percentile], axis=0)
        
    elif variation_type == "90ci": # 90% confidence interval
        lower_percentile = 5
        upper_percentile = 95       
        statecounts_variations = np.percentile(state_matrix_all, [lower_percentile, upper_percentile], axis=0)  
         
    elif variation_type == "std":
        state_std = np.std(state_matrix_all, axis=0)
        statecounts_variations = np.array([
            statecounts_mean - state_std,
            statecounts_mean + state_std
        ])
    else:  # range
        min_states = np.min(state_matrix_all, axis=0)
        max_states = np.max(state_matrix_all, axis=0)
        statecounts_variations = np.array([min_states, max_states])  
        
    

    rows = []
    for sim in sim_ids:
        rows.append({
            "simulation": sim,
            "time": times[sim],
            "state_count": state_counts[sim]
        })
    full_results_df = pd.DataFrame(rows)

    return unified_time, statecounts_mean, statecounts_variations, full_results_df