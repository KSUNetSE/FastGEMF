from typing import List, Dict
import yaml
from pathlib import Path
from typing import List, Dict, Any
from dataclasses import dataclass, field
from .visualization import draw_model_graph

@dataclass
class NodeTransition:
    """
    A class representing a parametric node-based transition .

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

    Methods:
    --------
    __str__():
        Returns a formatted string for printing  the node-based transition.
    """
    name: str
    from_state: str
    to_state: str
    rate: str
    def __str__(self):
        """
        Returns a string representation of the edge-based transition.
        for example:
        "recovery: infected --> susceptible  (delta)".
        """
        return f"{self.name}: {self.from_state} --> {self.to_state} (rate: {self.rate})"

@dataclass
class EdgeTransition:
    """
    A class representing parametric edge-based transition of nodes.

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

    def __str__(self):
        """
        Returns a string  of the edg-based transition.
        "infection: susceptible --(Infected)-> Infected on social_contacts (beta)".
        """
        return f"{self.name}: {self.from_state} --({self.inducer})-> {self.to_state} on {self.network_layer} (rate: {self.rate})"
    
class ModelSchema:
    """
    A class representing the schema of a custom parametric network-based mechanistic model, including compartments,
    network layers, node transitions, and edge-based interactions. This class is basically a schema or template of the model, allowing users
    to define general mechanistic model, which can later configure using ModelConfiguration class.

    Attributes:
    -----------
    name : str
        The name of the model.
    compartments : List[str]
        A list of compartments representing the states in the model (e.g., Susceptible, Infected).

    network_layers : List[str]
        A list of network layers' names, each representing a different mode of interaction (e.g., social, physical).

    node_transitions : List[NodeTransition]
        A list of node-based transitions object 

    edge_transitions : List[EdgeTransition]
        A list of edge-based transitions object.

    Methods:
    --------
    define_compartment(compartments: List[str]) -> 'ModelSchema':
        Defines compartments in the model by appending them if not already present.
    
    add_node_transition(name: str, from_state: str, to_state: str, rate: str) -> 'ModelSchema':
        Adds a node-based transition between compartments with an associated parameteric rate(e.g., 'beta', 'delta').
    
    add_edge_interaction(name: str, from_state: str, to_state: str, inducer: str, network_layer: str, rate: str) -> 'ModelSchema':
        Adds an edge-based interaction between state across speicif network layer with a specified inducer and parameteric rate.
    
    add_network_layer(\*layers: str) -> 'ModelSchema':
        Adds one or more network layers to the model.

    to_dict() -> dict:
        Converts the model schema to a dictionary format for serialization.

    to_yaml(file_path: str):
        Serializes the model schema to a YAML file.

    from_yaml(cls, yaml_file: str) -> 'ModelSchema':
        Class method to load a model schema from a YAML file.

    draw_model_graph():
        Visualizes the model's transitions using an external visualization library.

    __str__() -> str:
        Returns a string representation of the model schema, including compartments, network layers, and transitions.
    """
    def __init__(self, name: str = "Custom Model"):
        self.name = name
        self.compartments: List[str] = []
        self.network_layers: List[str] = []
        self.node_transitions: List[NodeTransition] = []
        self.edge_transitions: List[EdgeTransition] = []

    def define_compartment(self, compartments: List[str]) -> 'ModelSchema':
        for compartment in compartments:
            if compartment not in self.compartments:
                self.compartments.append(compartment)
        return self
    
    def add_node_transition(self, name: str, from_state: str, to_state: str, rate: str) -> 'ModelSchema':
        transition = NodeTransition(name, from_state, to_state, rate)
        self.node_transitions.append(transition)
        return self

    def add_edge_interaction(self, name: str, from_state: str, to_state: str, inducer: str, network_layer: str, rate: str) -> 'ModelSchema':
        transition = EdgeTransition(name, from_state, to_state, inducer, network_layer, rate)
        self.edge_transitions.append(transition)
        return self
    
    def add_network_layer(self, *layers: str) -> 'ModelSchema':
        for layer in layers:
            if layer not in self.network_layers:
                    self.network_layers.append(layer)
        return self
    

    
    def to_dict(self) -> dict:
        return {
            'name': self.name,
            'compartments': self.compartments,
            'network_layers': self.network_layers,
            'node_transitions': [{ 
                'name': t.name,
                'from_state': t.from_state, 
                'to_state': t.to_state,
                'rate': t.rate
            } for t in self.node_transitions],
            'edge_transitions': [{
                'name': t.name,
                'from_state': t.from_state,
                'to_state': t.to_state,
                'inducer': t.inducer,
                'network_layer': t.network_layer,
                'rate': t.rate
            } for t in self.edge_transitions]
        }

    def to_yaml(self, file_path: str):
        # Serialize the model to a YAML file
        with open(file_path, 'w') as file:
            yaml.dump(self.to_dict(), file, default_flow_style=False)

    @classmethod
    def from_yaml(cls, yaml_file: str) -> 'ModelSchema':
        # Deserialize a YAML file back into a Model instance
        with open(yaml_file, 'r') as file:
            data = yaml.safe_load(file)
        model = cls(name=data.get('name', 'Unnamed Model'))
        model.define_compartment(data['compartments'])
        model.add_network_layer(data['network_layers'])
        for trans in data.get('node_transitions', []):
            model.node_transitions.append(NodeTransition(**trans))
        for trans in data.get('edge_transitions', []):
            model.edge_transitions.append(EdgeTransition(**trans))
        return model
    
    def draw_model_graph(self):
        draw_model_graph(self)
        
    def __str__(self):
        node_trans_str = '\n'.join(str(t) for t in self.node_transitions)
        edge_trans_str = '\n'.join(str(t) for t in self.edge_transitions)
        return (f"Model: {self.name}\n"
                f"Compartments: {', '.join(self.compartments)}\n"
                f"Network Layers:\n {self.network_layers}\n"
                f"Node-Based Transitions:\n{node_trans_str}\n"
                f"Edge-Based Transitions:\n{edge_trans_str}")
    
