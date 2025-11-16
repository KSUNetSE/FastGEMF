import fastgemf as fg
import networkx as nx
import os
from os.path import join as ospj    
# to suppress warnings
import warnings
warnings.filterwarnings('ignore', category=RuntimeWarning)

# Define the SIR model schema
sir_model = (
    fg.ModelSchema("SIS")
    .define_compartment(['S', 'I'])
    .add_network_layer('contact_network')
    .add_node_transition(
        name='recovery',
        from_state='I',
        to_state='S',
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
num_nodes = 1000
num_edges = 6
contact_network = nx.barabasi_albert_graph(num_nodes, num_edges)

# Convert the network to a sparse matrix format
contact_network_csr= nx.to_scipy_sparse_array(contact_network)

# Create a model configuration instance
sir_instance = (
    fg.ModelConfiguration(sir_model)
    .add_parameter(beta=0.01, delta=0.1)
    .get_networks(contact_network=contact_network_csr)
)

# Print the configured model
print(sir_instance)


# Creating  the Simulation object
sim = fg.Simulation(sir_instance, initial_condition={'percentage': {'I': 1, 'S': 99}}, stop_condition={'time': 200}, nsim=10)
sim.run()
#if nsim>1, you get the mean and variation of the results, instead of single run results, modify variation_type to get different types of variations
variation_type=["iqr", "90ci", "std", "range"
                ] # if nsim>1, specify the type of variation to get the results and for save the plots
for var_type in variation_type:
    print(f"Plotting with variation type: {var_type}")
    time, statecounts_mean, *_ = sim.get_results(variation_type=var_type)
    statecounts_variations= None
    sim.plot_results(time, statecounts_mean, statecounts_variations, variation_type=var_type,
                            title=f"SIR Model Simulation - Variation: {var_type}",
                            y_axis="fraction",
                            grid=True,
                            save_path= ospj(os.getcwd(),f"simulation_figure-{var_type}.png")
                           )
