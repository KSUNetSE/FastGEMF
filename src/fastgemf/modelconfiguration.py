from  .utils import setup_simulation_matrices, extract_transition_parameters
import numpy as np
import yaml
from dataclasses import dataclass
from typing import List, Dict
import scipy.sparse as sp

@dataclass
class NodeTransitionValued:
    """
    A class representing a node-based transition which has  a value for the transiton.

    Attributes:
    -----------
    name : str
        The name of the transition (e.g., "recovery" ).

    from_state : str
        The state the node is transitioning from (e.g., "infected" or "S").

    to_state : str
        The state the node is transitioning to (e.g., "susceptible").

    rate : str
        The rate expression or identifier for the transition (e.g., "beta").

    value : float
        The numerical value associated with the rate (e.g., 0.01).

    Methods:
    --------
    __str__():
        Returns a formatted string for printing  the node-based transition.
    """
    name: str
    from_state: str
    to_state: str
    rate: str
    value: float

    def __str__(self):
        """
        Returns a string representation of the edge-based transition.
        for example:
        "recovery: infected --> susceptible  (delta: 0.01)".
        """
        return f"{self.name}: {self.from_state} -> {self.to_state} ({self.rate}: {self.value})"

@dataclass
class EdgeTransitionValued:
    """
    A class representing an edge-based transition of nodes with an associated value,

    Attributes:
    -----------
    name : str
        The name of the edge transition (e.g., "infection").

    from_state : str
        The state the source node is transitioning from (e.g., "susceptible").
    to_state : str
        The state the target node is transitioning to (e.g., "Infected").

    inducer : str
        The node or compartment that induces the transition (e.g., "Infected").

    network_layer : str
        The layer of the network where the interaction occurs (e.g., "social_contacts").

    rate : str
        The rate expression or identifier for the transition (e.g., "beta").
    value : float
        The numerical value associated with the rate (e.g., 0.05).

    Methods:
    --------
    __str__():
        Returns a formatted string representation of the edge-based transition.
    """
    name: str
    from_state: str
    to_state: str
    inducer: str
    network_layer: str
    rate: str
    value: float

    def __str__(self):
        """
        Returns a string  of the edg-based transition.
        "infection: susceptible --(Infected)-> Infected on social_contacts (beta:  0.05)".
        """
        return f"{self.name}: {self.from_state} --({self.inducer})-> {self.to_state} on {self.network_layer} ({self.rate}: {self.value})"
    
        
class ModelConfiguration:
    """
    A class for speciyng the values for the parametric mechanistic model deifned using ModelSchema. User can convert to/from YAML format for easier reproducibility of the results and 

    Attributes:
    -----------
    model : Any, optional
        An instance of the model that contains predefined settings, such as name, compartments, and network layers.

    name : str
        The name of the model.

    compartments : List[str]
        List of compartments in the model, representing the possible states of nodes.

    network_layers : List[str]
        The different layers of the network representing different modes of interaction (e.g., social, physical,sexual, etc.).

    node_transitions_parametric : List[dict]
        A list of parametric node transitions before parameter values are assigned.

    edge_transitions_parametric : List[dict]
        A list of parametric edge transitions before parameter values are assigned.

    parameters : Dict[str, float]
        A dictionary of parameters and their values that are used in the model's transitions.

    node_transitions_valued : List[NodeTransitionValued]
        A list of node transitions with assigned parameter values.

    edge_transitions_valued : List[EdgeTransitionValued]
        A list of edge transitions with assigned parameter values.

    network_directories : Dict[str, str]
        A dictionary mapping network layers to their file directories for loading network data.

    networks : Dict[str, Any]
        A dictionary mapping network layers to their network graph objects 

    Methods:
    --------
    get_networks(\**graphs):
        Updates the network graphs for the specified layers by receiving the CSR fomrat or Adjacency matrix.

    load_network_directory(\**direct):
        Loads networks from file directories for specified layers.

    add_parameter(\**kwargs):
        Adds model parameters and updates transitions based on these parameters.

    _update_transitions():
        Updates node and edge transitions with the current parameter values.

    to_dict() -> dict:
        Converts the model configuration into a dictionary format.

    to_yaml(file_path: str):
        Saves the model configuration as a YAML file.

    from_yaml(yaml_file: str) -> 'ModelConfiguration':
        Loads a model configuration from a YAML file.

    __str__() -> str:
        Provides a string representation of the model, including compartments, network layers, and transitions.

    _get_transition_matrices():
        Generates simulation matrices based on the model's transitions and parameters.

    _get_GEMF_matrices():
        Constructs GEMF  transition matrices using the current parameters.

    _print_GEMF_matrices():
        Prints the GEMF matrices and the model's configuration for debugging and verification purposes.
        Comment: Sine we map states to numerical values, susceptible is always mapped to zero.
    """
    def __init__(self, model=None, name=None, compartments=None, network_layers=None):
        """
        Initializes the ModelConfiguration object. 
        Parameters:
        -----------
        model : ModelSchema, optional if YAML of ModelSchema instance has passed. (by using from_yaml() function).
            An instance of a model that defines name, compartments, and network layers.

        name : str, optional.
            Name of the model (used if a model instance is not provided).

        compartments : List[str], optional.
            A list of compartments representing the states (used if a model instance is not provided).

        network_layers : List[str], optional.
            A list of network layers defining the interaction structures (used if a model instance is not provided).
        """
        if model: # when model is passed
            self.model = model
            self.name = model.name
            self.compartments = model.compartments
            self.network_layers = model.network_layers
            self.node_transitions_parametric: list[dict]=model.node_transitions
            self.edge_transitions_parametric: list[dict]=model.edge_transitions
        else: # when yaml is passed
            self.model = None
            self.name = name
            self.compartments : List[str]=compartments
            self.network_layers : List[str]=network_layers
            self.node_transitions_parametric: list[dict]=[]
            self.edge_transitions_parametric: list[dict]=[]

        self.network_layers.sort()
        self.parameters: Dict[str, float] = {}
        self.node_transitions_valued: List[NodeTransitionValued] = []
        self.edge_transitions_valued: List[EdgeTransitionValued] = []
        self.network_directories = {layer: 'Add the directory for graph here!' for layer in self.network_layers}
        self.networks = {layer: 'Insert the graph Compressed Sparse Row object here!' for layer in self.network_layers}


    def get_networks(self, **graphs):
        """

        update network graphs by recieving network objects for each user specified layeryers.
        for example: 
        SIR=ModeConfiguration.from_yaml("SIR.yaml") # creating the ModelConfiguration instance
        SIR.get_networks(layer_1=network_1) # network_1 should be either adjaceny matrix or CSR format. of not, transform before passing.

        """

        for layer_name, graph in graphs.items(): 
            for name in  self.network_layers:
                if layer_name==name:
                    self.networks.update({name:graph})
        return self
    
    def load_network_directory(self, **direct):
        """

        Update network graphs by loading network from their directories.
        for example: 
        SIR=ModeConfiguration.from_yaml("SIR.yaml") # creating the ModelConfiguration instance
        SIR.load_network_directory(layer_1=usr\defined\path\network_1.npz) # network_1 should be either adjaceny matrix or CSR format save as .npz format.

        """        
        for layer_name, network_dir in direct.items():
            for name in  self.network_layers:
                if layer_name==name:
                    self.network_directories.update({name:network_dir})
                    graph = sp.load_npz(network_dir)
                    self.get_networks(**{name: graph})
        return self
    
    def add_parameter(self, **kwargs):
        """

        updates node and edge transitions based on the new values.
        for example: 
        SIR=ModeConfiguration.from_yaml("SIR.yaml") # creating the ModelConfiguration instance
        SIR.add_parameter(beta=.05, delta=.01) # network_1 should be either adjaceny matrix or CSR format save as .npz format.

        """   
        self.parameters.update(kwargs)
        self._update_transitions()
        return self

    def _update_transitions(self):
        self.node_transitions_valued.clear()
        self.edge_transitions_valued.clear()

        for transition in self.node_transitions_parametric:
            if transition.rate in self.parameters.keys():
                valued_transition = NodeTransitionValued(
                    transition.name,
                    transition.from_state,
                    transition.to_state,
                    transition.rate,
                    self.parameters[transition.rate]
                )
                self.node_transitions_valued.append(valued_transition)

        for transition in self.edge_transitions_parametric:
            if transition.rate in self.parameters.keys():
                valued_transition = EdgeTransitionValued(
                    transition.name,
                    transition.from_state,
                    transition.to_state,
                    transition.inducer,
                    transition.network_layer,
                    transition.rate,
                    self.parameters[transition.rate]
                )
                self.edge_transitions_valued.append(valued_transition)
 
    def to_dict(self) -> dict:
        return {
            'name': self.name,
            'compartments': self.compartments,
            'parameters':self.parameters,
            'network_layers': self.network_directories,
            'node_transitions': [{ 
                'name': t.name,
                'from_state': t.from_state, 
                'to_state': t.to_state,
                'rate': t.rate,
                'value': t.value
            } for t in self.node_transitions_valued],
            'edge_transitions': [{
                'name': t.name,
                'from_state': t.from_state,
                'to_state': t.to_state,
                'inducer': t.inducer,
                'network_layer': t.network_layer,
                'rate': t.rate,
                'value': t.value
            } for t in self.edge_transitions_valued],

        }

    def to_yaml(self, file_path: str):
        """Saves the model configuration as a YAML file."""
        with open(file_path, 'w') as file:
            yaml.dump(self.to_dict(),   )



            
    @classmethod
    def from_yaml(cls, yaml_file: str) -> 'ModelConfiguration':
        """Loads a model configuration from a YAML file."""

        with open(yaml_file, 'r') as file:
            data = yaml.safe_load(file)
        
        instance = cls(
            name=data.get('name', 'Unnamed Model'),
            compartments=data.get('compartments', []),
            network_layers=list(data.get('network_layers', {}).keys())
        )
        instance.parameters=data.get('parameters', {})
        network_layers = data.get('network_layers', {})
        instance.load_network_directory(**network_layers)

        for trans in data.get('node_transitions', []):
            instance.node_transitions_valued.append(NodeTransitionValued(**trans))
            instance.node_transitions_parametric.append(NodeTransitionValued(**trans))

        for trans in data.get('edge_transitions', []):
            instance.edge_transitions_valued.append(EdgeTransitionValued(**trans))
            instance.edge_transitions_parametric.append(EdgeTransitionValued(**trans))
            
        instance.add_parameter(**instance.parameters)

        return instance



    def __str__(self):
        node_trans_str = '\n'.join(str(t) for t in self.node_transitions_valued)
        edge_trans_str = '\n'.join(str(t) for t in self.edge_transitions_valued)
        networks = []
        for name, graph in self.networks.items():
            networks.append(f"{name}: Network With {graph.shape[0]} Nodes and {graph.nnz} Edges")
        networks_str = '\n'.join(networks)
        return (f"Model: {self.name}\n"
                f"Compartments: {', '.join(self.compartments)}\n"
                f"Network Layers:\n{networks_str}\n"
                f"Node-Based Transitions:\n{node_trans_str}\n"
                f"Edge-Based Transitions:\n{edge_trans_str}\n"

                )

                                        ############################################### 
    
    def _get_transition_matrices(self):
            self.q=setup_simulation_matrices(self)["q"]
            self.layers_infuencer=setup_simulation_matrices(self)['layers_influencer']
            self.L=setup_simulation_matrices(self)["L"]
            self.M=setup_simulation_matrices(self)["M"]
            self.A_beta=setup_simulation_matrices(self)["A_beta"]
            self.A_delta=setup_simulation_matrices(self)["A_delta"]
        
    def _get_GEMF_matrices(self):
        self._get_transition_matrices()
        self.map_edge_parameters =extract_transition_parameters(self.A_beta)
        self.map_node_parameters =extract_transition_parameters(self.A_delta)
        self.A_beta, self.A_delta=np.zeros((self.L,self.M,self.M)),np.zeros((self.M,self.M))

        for key, value in self.parameters.items():
            if key in self.map_edge_parameters.keys():
                self.A_beta[self.map_edge_parameters[key]]=value

            elif key in self.map_node_parameters.keys() :
                self.A_delta[self.map_node_parameters[key]]=value

        return self
    
    def _print_GEMF_matrices(self):
        self._get_GEMF_matrices()
        print(f"Configuration for Model: {self.name}\n"
            f"influencer_of_layers:\n {self.q}\n"  
            f"Edge_based_matrix:\n {self.A_beta}\n"
            f"Node_based_matrix:\n {self.A_delta}\n"
            f"Networks: {self.network_layers}\n")