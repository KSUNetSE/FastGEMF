
from .utils import sample
import numpy as np
from .error import StopLoopException 


def sample_event(times, event_data,model_matrices, rate_arrays, X, iteration,Tf):
    """
    Samples the next event by finding the node with the minimum absolute time, updating its state, and adjusting its rate.

    Parameters:
    -----------
    times : TimeNP or TimeSorted
        An object responsible for handling absolute times.

    event_data : EventData
        An object that tracks node states before and after events.

    model_matrices : ModelMatrices
        Contains the model matrices used in the GEMF framework.

    rate_arrays : RateArray
        Arrays containing the transition rates for each node.

    X : np.ndarray
        A NumPy array representing the current states of the nodes in the network.

    iteration : np.ndarray
        A NumPy array indicating the current iteration number of the simulation.

    Tf : np.ndarray
        A NumPy array representing the current simulation time.

    Description:
    ------------
    This function samples the next event from the time structure, updates the corresponding node's state, and adjusts the transition rates accordingly. It first identifies the node scheduled for the next event, updates its state using the information from `model_matrices`, and then modifies the rates in `rate_arrays` to reflect the state change.
    
    """

    sampled_node, min_time=times.get_min()
    state_k=X[sampled_node].astype(int)
    state_k_plus_1=sample_new_state(model_matrices, rate_arrays, state_k, sampled_node )
    event_data.add_event(state_k, state_k_plus_1, sampled_node, min_time )
    X[sampled_node]=state_k_plus_1
    rate_arrays.R -=rate_arrays.Rn[sampled_node]
    iteration +=1
    rate_arrays.Rn[sampled_node]=update_node_rate(model_matrices, rate_arrays,
                                                   state_k_plus_1, sampled_node )
    
    times.pop_and_push(sampled_node, min_time, rate_arrays.Rn[sampled_node]) 
    
    rate_arrays.R += rate_arrays.Rn[sampled_node]        
    Tf[0]=min_time
    


def update_network(times, event_data,model_matrices, rate_arrays, networks, X):
    """
    Cautious update approach: Updates affected neighbors after a node transitions.

    Parameters:
    -----------
    times : TimeNP or TimeSorted
        An object responsible for handling absolute times.

    event_data : EventData
        An object that tracks node states before and after transitions.

    model_matrices : ModelMatrices
        Contains the model matrices used in the GEMF framework.

    rate_arrays : RateArray
        Arrays containing the transition rates for each node.

    X : np.ndarray
        A NumPy array representing the current states of the nodes in the network.

    iteration : np.ndarray
        A NumPy array indicating the current iteration number of the simulation.

    Tf : np.ndarray
        A NumPy array representing the current simulation time.

    Description:
    ------------
    This function cautiously updates the network by first processing a node transition, then updating the rates of neighboring nodes that are affected by this transition. The method ensures that any changes in a node's state are reflected in its neighbors' transition rates, using the provided model matrices and rate arrays.
    
    """
    
    layers=model_matrices.q.get(event_data.states_k_plus_1[-1])
    if layers:
        for layer in layers: 
            out_neighbors=networks.get_out_neighbors(layer,event_data.sampled_nodes[-1])
            neighbor_weights=networks.get_out_weight(layer,event_data.sampled_nodes[-1])
            rate_arrays.Nq[layer][out_neighbors] +=neighbor_weights
            delta_rates=model_matrices.bil[X[out_neighbors],layer]*neighbor_weights
            nonzero_indices=delta_rates.nonzero()[0]
            rate_arrays.Rn[out_neighbors[nonzero_indices]] +=delta_rates[nonzero_indices] 
            rate_arrays.R += np.sum(delta_rates)
            times.update(out_neighbors[nonzero_indices],event_data.min_times[-1],
                        rate_arrays.Rn[out_neighbors[nonzero_indices]] )
        
            
    layers=model_matrices.q.get(event_data.states_k[-1])
    if layers:    
        for layer in layers: 
            out_neighbors=networks.get_out_neighbors(layer,event_data.sampled_nodes[-1])
            neighbor_weights=networks.get_out_weight(layer,event_data.sampled_nodes[-1])
            rate_arrays.Nq[layer][out_neighbors] -=neighbor_weights
            delta_rates=model_matrices.bil[X[out_neighbors],layer]*neighbor_weights
            nonzero_indices=delta_rates.nonzero()[0]
            rate_arrays.Rn[out_neighbors[nonzero_indices]] -=delta_rates[nonzero_indices]
            rate_arrays.R -= np.sum(delta_rates) 
            times.update(out_neighbors[nonzero_indices],event_data.min_times[-1],
                        rate_arrays.Rn[out_neighbors[nonzero_indices]] )        


def sample_new_state(model_matrices, rate_arrays, state_k,sampled_node):
    """
    Samples the new state for a node based on current rates.
    -----------
    Returns:
    -----------
    int
    an integer in  {0,1,...,M} 
    
    """

    return (sample( np.ravel(model_matrices.A_delta[state_k,:].T
                    +np.dot(model_matrices.bi[state_k],
                            rate_arrays.Nq[:,sampled_node] )) ) )

def update_node_rate(model_matrices, rate_arrays, state_k_plus_1, sampled_node):
    """
    Updates the rate of transitions for a node based on its new state.

    Returns:
    -----------
    float
        The updated transition rate for the node.
    
    """

    return ( model_matrices.di[state_k_plus_1]
                +np.dot(model_matrices.bil[state_k_plus_1,:],
                        rate_arrays.Nq[:,sampled_node]) )   
