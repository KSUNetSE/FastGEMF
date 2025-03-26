import numpy as np
import random


def randomly_assign_states(N, states_number,hubs=None):
    """
    Assigns states to nodes in the network randomly based on the given state numbers.

    Parameters:
    -----------
    N : int
        The total number of nodes in the network.

    states_number : dict
        A dictionary with states as keys and the number of nodes in each state as values.( e.g., {0: 10, 1: 10} , where 0 state is by default susceptible.)

    
    hubs : list of int, optional, default is None
        A list of node with highest degree .

    Returns:
    -------
    np.ndarray
        An array of length `N` with node states assigned randomly.

    Raises:
    ------
    ValueError
        If the sum of the number of nodes in each state exceeds the total number of nodes in the network.

    """
    state_number={state:number for state, number in states_number.items() if number>1e-6 and state>1e-6}
    #print(state_number)
    x0=np.zeros(N)
    states=list(state_number.keys())
    number=list(state_number.values())
    #print(states_number,number)
    if sum(number)>N:
        raise ValueError("The sum of the number of nodes in each state should not exceed the total number of nodes in the network!")
    if type(hubs) is not np.ndarray and type(hubs) is not list or hubs is None:
        nodes=np.arange(0,N,dtype=int)
        random.shuffle(nodes)
    else:
        nodes=hubs
    index=np.cumsum([0]+number) 
    index=np.array(index)
    for  i ,state  in enumerate(states):
        x0[nodes[index[i]:index[i+1]]]=state 
    return np.array(x0,dtype=int)
def initial_condition_func(networks, inst, initial_condition,counter):
    """
    Generates the initial conditions for nodes in the network based on the specified method in `initial_condition`.

    Parameters:
    -----------
    networks : Network
        The network object that defines the structure of the network layers.

    inst : ModelConfiguration
        The ModelConfiguration instance that defines the compartments and transitions.

    initial_condition : dict
        A dictionary specifying the method for initializing node states. Methods for the SIS model include:
        - 'default_percentage': {} 
        Randomly assigns 10 percent of the population to the inducing state (e.g., 10 percent infected) and the remaining 90 percent equally and randomly to other states (e.g., 90 percent susceptible).
        - 'percentage': {'I': 5, 'S': 95} 
        User-defined percentages for each compartment, assigning nodes randomly to states based on these percentages (e.g., 5 percent infected, 95 percent susceptible).
        - 'number': {'I': 10, 'S': 90} 
        User-defined number of nodes to be assigned to specific states (e.g., 10 nodes percent infected, 90 nodes percent susceptible).
        - 'hubs_number': {'I': 10} 
        User-defined number of hubs (nodes with the most connections) to be assigned to specific states.
        - 'exact': x0
        Directly sets the initial state from an existing array `x0`.

    counter : int
        A counter to track how many times the function has been called, used to prevent infinite recursion.

    Returns:
    --------
    x0 : np.ndarray
        An N x M array representing the initial state of each node in the network, where each element is the 
        compartment/state index for that node.
    
    Description:
    ------------
    This function assigns initial states to nodes based on the selected method:
    - 'default_percentage': Randomly assigns 10% of the nodes to inducer compartments and 90% to susceptible( 0 ).
    - 'percentage': Distributes nodes into compartments according to user-defined percentages.
    - 'hubs_number': Places a fixed number of nodes into specific compartments, targeting nodes with the highest degree.
    - 'exact': Uses an existing array `x0` to directly set the initial state of each node.

    If an invalid method is provided, the function defaults to the 'default_percentage' method.

    Raises:
    -------
    ValueError:
        Raised if the percentages provided in the 'percentage' method do not sum to 1 or 100.
    """

    method = list(initial_condition.keys())[0]
    if method == 'default_percentage': 
        N = networks.nodes
        J = list(inst.q.keys())
        x0 = np.zeros(N)
        #others = [i for i, _ in enumerate(inst.compartments) if i not in J]  
        percentage_inducer = 0.1
        inducers_numbers={inducer:int(number) for inducer, number in zip(J,([percentage_inducer*N]*len(J)))}
        x0=randomly_assign_states(N, inducers_numbers)
        """
        NJ = int(N * percentage_inducer)
        NO = int(N * percentage_susceptible)
        x0 = np.zeros(N)
        
        if NJ > N:
            return 'Oops! Initial infection is more than the total population'
        else:
            temp = np.random.permutation(N)
            nj = temp[:NJ]
            no = temp[NJ:NJ + NO]  # Remaining nodes for 'others'
            index = np.arange(0, NJ, NJ / len(J)).astype(int)
            index = np.append(index, NJ)
            for i in range(len(J)):
                x0[nj[index[i]:index[i + 1]]] = J[i]

            if others:  
                nodes_per_other = NO // len(others) 
                leftover = NO % len(others) 
                index_others = np.arange(0, NO, nodes_per_other).astype(int)
                index_others = np.append(index_others, NO)
                for i, state in enumerate(others):
                    x0[no[index_others[i]:index_others[i+1]]] = state

                if leftover:
                    extra_indices = np.random.choice(others, size=leftover, replace=True)
                    for i, state in enumerate(extra_indices):
                        x0[no[-leftover + i]] = state
        """
        return x0.astype(int)
    if method == 'number':
        x0=np.zeros(networks.nodes)
        N = networks.nodes
        state_numbers={i: int(initial_condition['number'].get(compartment, 0)) for i, compartment in enumerate(inst.compartments)}
        x0=randomly_assign_states(N, state_numbers)
        return x0

        
    if method == 'percentage': 
        N = networks.nodes
        state_percentage={i: initial_condition['percentage'].get(compartment, 0) for i, compartment in enumerate(inst.compartments)}
        J = list(inst.q.keys())
        if sum( list( state_percentage.values() ) ) not in [1,100]:
            raise ValueError ("The percentage for states must add-up to 1 or 100")
        if sum( list( state_percentage.values() ) ) == 1:
            normalizer=1
        else:
            normalizer=100 
        state_numbers={state: int(percent/normalizer*N) for state, percent in state_percentage.items()}   
        x0=randomly_assign_states(N, state_numbers)
        return x0
    
    if method == 'hubs_number':
        N = networks.nodes
        inst.compartments
        state_numbers={i: initial_condition['hubs_number'].get(compartment, 0) for i, compartment in enumerate(inst.compartments)}
        x0=np.zeros(N)
        tot_no_nodes=sum(list(state_numbers.values()))
        states=list(state_numbers.keys())
        numbers=list(state_numbers.values())
        nodes=networks.get_highest_degree_nodes(0,tot_no_nodes)
        #x0[nodes]=random.choices(states, weights=numbers, k=tot_no_nodes)
        x0=randomly_assign_states(N, state_numbers,nodes)
        return x0.astype(int)

    if method == 'exact':
        x0=initial_condition['exact']    
        return x0.astype(int)
    else :
        if counter<1:
            print("Please enter a valid option for initial condition! \n ") 
            print ('exact' ',''percentage'',''random' ',''hubs_number' '  '  'are only valid options')
            print("\n For unvalid options FastGEMF automatically put 10 percent at inducing states \n and distributes 90 percent equally  between other states")
        return initial_condition_func(networks, inst, 
                                    initial_condition={'default_percentage':[]}, counter=counter)

