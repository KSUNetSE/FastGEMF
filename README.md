
[![License](https://img.shields.io/github/license/KSUNetSE/FastGEMF?style=plastic)](https://github.com/KSUNetSE/FastGEMF/blob/main/LICENSE)
![GitHub last commit](https://img.shields.io/github/last-commit/KSUNetSE/FastGEMF?style=plastic&cacheSeconds=60)
[![PyPI Downloads](https://static.pepy.tech/badge/fastgemf)](https://pepy.tech/projects/fastgemf)

<div align="center">
  <h1> FastGEMF</h1>
  <h1> Scalable High-Speed Simulation of Stochastic Spreading Processes over Complex Multilayer Networks</h1>
 <p><img src="images/logo.png" alt="Alt text" width="450"/> </p>
</div>


FastGEMF is scalable  spread process simulator  for small to large-scale multi-layer complex networks. You can define mechanistic model with mutiple infleuncer through simple steps and simualte the difussion process over networks.



## Table of Contents

- [Features](#features)
- [Installation](#installation)
- [Quick Implementation](#quick-implementation)
- [Usage](#usage)
- [Repository Structure](#repository-structure)
- [How To Cite](#how-to-cite)

## Features 
- Fast and efficient event based simulator( logarithmic time complexity w.r.t. number of nodes in networks).
- Multi-layer networks are supported.
- Networks can be weighted/unweighted or directed/undirected.
- Mechanistic models can have different influencing agents competing over same or different layers.

## Installation

FastGEMF can be installed on Linux, macOS, and Windows operating systems. Ensure you have Python and pip installed on your system by following the official pip documentation.

Using PyPI
To install FastGEMF from the Python Package Index (PyPI), execute the following command in your terminal:
```sh
python3 -m pip install fastgemf --user --upgrade
```
Using the Source Code Repository: 
If you prefer to install the latest version directly from the source, you can use pip to install FastGEMF from GitHub:

```sh
python3 -m pip install https://github.com/KSUNetse/FastGEMF/archive/main.zip --user --upgrade
```
Alternatively, you can clone the repository and install FastGEMF locally:

```sh
git clone https://github.com/KSUNetse/FastGEMF.git
cd FastGEMF
python3 -m pip install . --user --upgrade
```
After installation, you can verify that FastGEMF is correctly installed by running:
```sh
python -c "import fastgemf; print(fastgemf.__version__)"
```

## Quick implementation
In the following example, we define a simple SIR model over a single layer contact network.

```python
import fastgemf as fg
import networkx as nx
# to suppress warnings
import warnings
warnings.filterwarnings('ignore', category=RuntimeWarning)

# Define the SIR model schema
sir_model = (
    fg.ModelSchema("SIR")
    .define_compartment(['S', 'I', 'R'])
    .add_network_layer('contact_network')
    .add_node_transition(
        name='recovery',
        from_state='I',
        to_state='R',
        rate='delta'
    )
    .add_edge_interaction(
        name='infection',
        from_state='S',
        to_state='I',
        inducer='I',
        network_layer='contact_network',
        rate='beta'
    )
)

# Print the model schema
print(sir_model)

# Visualize the model graph
sir_model.draw_model_graph()

# Generate a Barabási-Albert graph as an example network using popular module NetworkX
num_nodes = 10000
num_edges = 6
contact_network = nx.barabasi_albert_graph(num_nodes, num_edges)

# Convert the network to a sparse matrix format
contact_network_csr= nx.to_scipy_sparse_array(contact_network)

# Create a model configuration instance
sir_instance = (
    fg.ModelConfiguration(sir_model)
    .add_parameter(beta=0.3, delta=0.1)
    .get_networks(contact_network=contact_network_csr)
)

# Print the configured model
print(sir_instance)


# Creating  the Simulation object
sim = fg.Simulation(sir_instance, initial_condition={'percentage': {'I': 10, 'S': 90}}, stop_condition={'time': 30}, nsim=1)
sim.run()
sim.plot_results()
```

## Usage

By using Two modules `ModelSchema` and `ModelConfiguration` we can define any competing-like mechanistic model over different or same layers. `ModelSchema` module enabes user to create any multi-agent mechanistic model. The module can be from `FastGEMF` as:
```pyhton
import fastgemf as fg 

#Here we define a simple, yet popular `SIR` model.
SIR_model = (fg.ModelSchema(name='SIR') # model name
    .define_compartment(compartments=["S", "I", "R"]) # the compartments of the model
    .add_network_layer("contact_network") # Defining the layers' 
    .add_node_transition(name="recovery", from_state="I", to_state="R", rate="delta") # Defining node-based transition
    .add_edge_interaction(name="infection", from_state="S", to_state="I", inducer="I", network_layer="contact_network", rate="beta") # Defining edge-based transition
)
```
The user can now save the parametric mechanistic model defined in a `yaml` file, as a template for easy future use.
```pyhton
SIR_model.to_yaml("your_directory\SIR.yaml") # to save the defined model in `your_directory`
```
The user can draw  or print the model they defined as:
```python
SIR_model.draw_model_graph()
```
Below are two visual representations of the SIR model we defined:
<p align="center">
  <img src="images/SIR_node_based.png" alt="SIR Node Transition" width="60%" />
<\p>
<p align="center">
  <img src="images/SIR_edge_based.png" alt="SIR Edge Interaction" width="60%" />
</p>
or user can easily print the details  as: 
  
```python
print(SIR_model)
```
The printed result is as:
```
Model: SIR
Compartments: S, I, R
Network Layers: ['contact_network']
Node-Based Transitions: recovery: I --> R (rate: delta)
Edge-Based Transitions: infection: S --(I)-> I on contact_network (rate: beta)
```
In the next step, we use the `ModelConfiguration` module to assign values to the parameters of the model instane we defined(If you want to save to YAML make sure you have loaded the directories for networks):
```python
SIR_instnace=(fg.ModelConfiguration(SIR_model)
.add_parameter(beta=.05,delta=.05) # set the defined parameters to numeric values
.get_networks(contact_network=network1 )) # passing the network object to model instnace
#.load_network_directory(contact_network="directory/for/contact_network") # If you want to save to YAML make sure you have loaded the directories for network
SIR_instnace.to_yaml("your_directory\SIR_instance.yaml") # saving the instance of the model as a yaml file
```
and the instance of `ModelConfiguration` by printing it:
```python
print(SIR_instance)
```
The printed result is as:
```
Model: SIR
Compartments: S, I, R
Network Layers:
contact_network: Network With N Nodes and L Edges # N and L will replace with real size of the passed network
Node-Based Transitions:
recovery: I -> R (delta: 0.05)
Edge-Based Transitions:
infection: S --(I)-> I on contact_network (beta: 0.05)

```
Finally to run the simulation and plotting the results, we need to create an object of the `Simulation` class and create a YAML file of the object for reprocibility of the results as follows:

```python
# Creating  the Simulation object
sim = fg.Simulation(SIR_instance, initial_condition={'percentage': {'I': 10, 'S': 90}}, stop_condition={'time': 30}, nsim=1)
sim.run()
sim.plot_results()
sim.to_yaml('sim.yaml')
```
`initial_condition={'percentage': {'I': 10, 'S': 90}` specifies the percentage of nodes to be initially and randomly at different states , `stop_condition={'time': 30}` can be number of `events` as number of transitions(each transition is counted as one event) or `time` as in this example. `nsim` determines the number of simulation to run.

Also, we can pass  two YAML files for `ModelCofiguration` and `Simulation` instances as:
```python
sim = fg.Simulation.from_yaml('SIR_instance.yaml','sim.yaml')
sim.run()
sim.plot_results()
```

Now Similiarly now we can define the `SI_1I_2` model, when two influencers are cometing over two different layers, while each of the influencer have different rate for transmiting the disease and different curing rates. We can define this model as:
```python
# Define the SI_1I_2 model schema
si1i2_model = (
    fg.ModelSchema("SI_1I_2")
    .define_compartment(['S', 'I_1', 'I_2'])
    .add_network_layer('primary_layer')
    .add_network_layer('secondary_layer')
    .add_node_transition(name='recovery1', from_state='I_1', to_state='S', rate='delta1')
    .add_node_transition(name='recovery2', from_state='I_2', to_state='S', rate='delta2')
    .add_edge_interaction(
        name='infection1', from_state='S', to_state='I_1', inducer='I_1',
        network_layer='primary_layer', rate='beta1'
    )
    .add_edge_interaction(
        name='infection2', from_state='S', to_state='I_2', inducer='I_2',
        network_layer='secondary_layer', rate='beta2'
    )
)

print(si1i2_model)
si1i2_model.draw_model_graph()

# Generate two random geometric graphs
primary_network = nx.random_geometric_graph(1000, 0.02)
secondary_network = nx.random_geometric_graph(1000, 0.05)

# Convert networks to sparse matrix format
primary_network_csr = nx.to_scipy_sparse_array(primary_network)
secondary_network_csr = nx.to_scipy_sparse_array(secondary_network)

# Create an isntance of ModelConfiguration, as SI_1I_2 model which also
si1i2_instance = (
    fg.ModelConfiguration(si1i2_model)
    .add_parameter(beta1=0.1, delta1=0.2, beta2=0.05, delta2=0.1)
    .get_networks(primary_layer=primary_network_csr, secondary_layer=secondary_network_csr)
)

print(si1i2_instance)

# Creating  the Simulation object
sim = fg.Simulation(si1i2_instance, initial_condition={'percentage': {'I_1': 5, 'I_2': 5, 'S':90}}, stop_condition={'time': 100}, nsim=1)
sim.run()
sim.plot_results()
```

## Repository Structure

```
FastGEMF/                  
├── docs/                      
├── predefined_instants/       
├── src/                       
│   └── fastgemf/              
│       ├── __init__.py        
│       ├── error.py          
│       ├── GEMFCore.py       
│       ├── GEMFSimulation.py  
│       ├── initial_condition.py
│       ├── initializer.py     
│       ├── modelconfiguration.py
│       ├── modelschema.py    
│       ├── network.py         
│       ├── post_population.py 
│       ├── stop_conditions.py 
│       ├── times_structure.py 
│       ├── utils.py           
│       ├── visualization.py   
├── tests/                     
├── pyproject.toml             
├── README.md                  
├── requirements.txt          
├── setup.cfg                  
└── setup.py
```
## How to cite
```
@ARTICLE{10876117,
  author={Hossein Samaei, Mohammad and Darabi Sahneh, Faryad and Scoglio, Caterina},
  journal={IEEE Access}, 
  title={FastGEMF: Scalable High-Speed Simulation of Stochastic Spreading Processes Over Complex Multilayer Networks}, 
  year={2025},
  volume={13},
  number={},
  pages={27112-27125},
  keywords={Stochastic processes;Epidemics;Nonhomogeneous media;Computational modeling;Analytical models;Scalability;Python;Complexity theory;Accuracy;Solid modeling;Complex networks;Markov process;epidemic spreading;mechanistic models;simulation},
  doi={10.1109/ACCESS.2025.3539345}}

  ```
