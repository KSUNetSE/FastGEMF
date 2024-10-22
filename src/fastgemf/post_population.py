import numpy as np
import copy

def post_population(x0, model_matrices,event_data, nodes  ):
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
    return time, StateCount[:,:],  ts, states_k, states_k_plus_1