import numpy  as np
import scipy.sparse as sp
import copy as copy
from .utils import  pick_fastest_method
from .network import Network
from .initial_condition import initial_condition_func
from .times_structure import TimeSorted,TimeNp
from .modelconfiguration import ModelConfiguration
from  dataclasses import dataclass, field
from  typing import List, Any, Dict
from .stop_conditions import stop_cond

@dataclass
class EventData:
    """
    This class  stores event data during the simulation.

    Attributes:
    -----------
    states_k : List[int]
        the state of node sampled node before transition.

    states_k_plus_1 : List[int]
        the state of node sampled node after transition.

    sampled_nodes : List[int]
        a list sampled nodes in each iteration.

    min_times : List[int]
        List of absolute times at each event (minimum times when events occur).

    Methods:
    --------
    add_event(current_state, next_state, sampled_node, min_time):
        Adds an event to the event data by recording the current and next states, the sampled node, and the event time.
    """    
    states_k: List[int] = field(default_factory=list)
    states_k_plus_1: List[int] = field(default_factory=list)
    sampled_nodes: List[int] = field(default_factory=list)
    min_times: List[int] = field(default_factory=list)

    def add_event(self, current_state, next_state, sampled_node, min_time):
        self.states_k.append(current_state)
        self.states_k_plus_1.append(next_state)
        self.sampled_nodes.append(sampled_node)
        self.min_times.append(min_time)


@dataclass
class ModelMatrices:
        """
    this class store the GEMF matrices used in the simulation.

    Attributes:
    -----------
    A_beta : np.ndarray
        Edge-based transititon adjacency matrix
    bi : np.ndarray
        M x M x L
    bil : np.ndarray
        M x L array, the sum of A_beta rows for each state.
    A_delta : np.ndarray
        Node-based adjaceny matrix
    di : np.ndarray
        The sum of the node-based rates
    q , layers_infuencer : Dict
        mapping influencers to layers and vice versa
    L : int
        Number of layers.
    M : int
        Number of compartments.

    Methods:
    --------
    from_instant(cls, cfg, get_model_arrays_func):
        Factory method to create ModelMatrices instance using ModelConfiguration instance that user has defined for it mechansitic model.
        """        
        A_beta:np.ndarray
        bi: np.ndarray # M x M x L
        bil: np.ndarray # M x L array, each col. is sum of A_beta rows for that layyer
        A_delta: np.ndarray
        di: np.ndarray # M    sum of the A_dela rows
        q: Dict                # key: influencer, Value: layer
        layers_infuencer: Dict # key: layer, Value: influencer
        L: int
        M: int

        @classmethod
        def from_instant(cls, cfg: Any, get_model_arrays_func):
                A_beta, bi, bil, A_delta, di, q,layers_infuencer, L, M = get_model_arrays_func(cfg)
                return cls(A_beta=A_beta, bi=bi,bil=bil, A_delta=A_delta,di=di,
                q=q,layers_infuencer=layers_infuencer, L=L, M=M )
        
@dataclass
class RateArray:
        """
    A class to store and manage the rate arrays used in the simulation.

    Attributes:
        -----------
    Nq : np.ndarray
        N x L array, storing sum of link weights of influencing neighbors.

    Rn : np.ndarray
        1D array of node rates, representing the transition rates for each node.

    R : float
        The total rate of events in the system, calculated as the sum of the node rates.

    Methods:
        --------
    from_instant(cls, cfg, X, networks, model_matrices, get_rates_func):
        A factory method to create a `RateArray` instance by calculating rates using a given configuration,
        node states, network structure, and model matrices.

        """
        Nq: np.ndarray  # N x L, Sum of  (Weights*Edge_interaction_rate) for eah node at each layer
        Rn: np.ndarray  # Node rates
        R: float        # Total rate

        @classmethod
        def from_instant(cls, cfg: Any, X: np.ndarray, networks: Network,model_matrices: ModelMatrices,  get_rates_func):
                Nq, Rn, R = get_rates_func(cfg, X, networks,model_matrices)
                return cls(Nq=Nq, Rn=Rn, R=R)
    

class Initialize():
        """
    A class to initializing the simulation environment and preparing all necessary components 
    such as event data, network structure, model matrices, and rate arrays.

    Attributes:
    -----------
    event_data : EventData
        Object that stores events occurring during the simulation.

    networks : Network
        The network structure using CSR and CSC data structure for effecient search through network.

    model_matrices : ModelMatrices
        Contains the GEMF matrices used in the simulation.

    X0 : np.ndarray 
        N x M array, Initial state of all nodes, representing which compartment each node is in at the start.

    X : np.ndarray
        A copy of the initial state `X0` used to keep track of the current states during the simulation.

    rate_arrays : RateArray
        Stores node transition rates related values.

    times : np.ndarray
       absolute times of nodes.

    iteration : np.ndarray
        counting the number of events.

    stop_condition : dict
        Condition that defines when the simulation should stop, such as the number of events or a time limit.

    Tf : np.ndarray
        Final simulation time.

    Methods:
    --------
    __init__(model_inst, sim_cfg, generate_times_func, counter):
        Initializes the simulation environment by setting up event data, networks, model matrices,
        initial conditions, and rate arrays. Time steps are generated based on the chosen method.
    
    _get_networks(model_inst):
        Retrieves network data either from file paths or directly from pre-defined environment objects

    _get_rates(model_inst, X, net, model_matrices, initial_inf_nodes=None):
        Computes transition rates for nodes and interactions based on the current state of the network 
        and model matrices.

    _get_modified_rates(cfg):
        Modifies and returns the model matrices used for transitions and rates based on the configuration.

    _initial_condition(model_inst, sim_cfg, counter):
        Sets up the initial condition for the node states based on the network and model configuration.
                """
        def __init__(self, model_inst:ModelConfiguration, sim_cfg, generate_times_func, counter) :

                self.event_data=EventData()
                self.networks=Network(self._get_networks(model_inst))
                self.model_matrices=ModelMatrices.from_instant(model_inst, self._get_modified_rates)
                self.X0=self._initial_condtion(model_inst, sim_cfg, counter)
                self.X=copy.deepcopy(self.X0)
                self.rate_arrays=RateArray.from_instant(model_inst, self.X, self.networks,
                                                         self.model_matrices, self._get_rates)
        
                if generate_times_func: 
                        self.times=generate_times_func(self.rate_arrays.Rn)
                else:
                        generate_times_func=pick_fastest_method(self.networks.nodes, 
                                                                self.networks.edges[0] )  
                        self.times=generate_times_func(self.rate_arrays.Rn)

                self.iteration=np.array(-1)
                self.stop_condition= sim_cfg['stop_condition']
                self.Tf=np.array([0.0]).astype(float)
                
                
 

        def _get_networks(self,model_inst):
               if  isinstance(model_inst.networks.values(),str):
                        return [sp.load_npz(dir) for dir in  model_inst.network_directories.values()]
               else:
                        return [network for network in model_inst.networks.values()]


        def _get_rates(self,model_inst,X, net: Network, model_matrices: ModelMatrices, initial_inf_nodes=None ):
                L, q= model_matrices.L, model_matrices.q
                N=net.nodes
                di, bil=model_matrices.di, model_matrices.bil
                Nq = np.zeros((L,N))       
                for n in range(N):
                        for l in range(L):
                                influencers=model_matrices.layers_infuencer.get(l)
                                for influencer in influencers: 
                                        Nln=net.get_in_neighbors(l,n)
                                        weights=net.get_in_weight(l,n)
                                        Nq[l][n]+=sum((np.array(X[Nln]==influencer))* weights)
                Rn = np.zeros(N)
                Rn = np.array([di[X[n]] + np.dot(bil[X[n], :], Nq[:, n]) for n in range(N)])
                R=Rn.sum(axis=0)        
                return  Nq, Rn, R       


        def _get_modified_rates (self,cfg):
                cfg._get_GEMF_matrices()
                di = cfg.A_delta.sum(axis=1) 
                layers_infuencer=cfg.layers_infuencer
                bil = np.array([np.sum(mat, axis=1) for mat in cfg.A_beta]).T
                bi = np.zeros(( cfg.M, cfg.M, cfg.L ))
                for i in range(cfg.M):
                        for l in range(cfg.L):
                                bi[i, :, l] = cfg.A_beta[l][i,:]
                return cfg.A_beta, bi, bil, cfg.A_delta,di, cfg.q,layers_infuencer, cfg.L, cfg.M 


        def _initial_condtion(self,model_inst, sim_cfg, counter):
                X0= initial_condition_func(self.networks, model_inst, 
                                              sim_cfg['initial_condition'], counter)
                counter+=1
                return X0
                