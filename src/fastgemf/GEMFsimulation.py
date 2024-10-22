import numpy  as np
import yaml 
from .error import StopLoopException 
from .times_structure import TimeSorted,TimeNp
from .post_population import post_population
from .initializer import Initialize
from .modelconfiguration import ModelConfiguration
from .GEMFcore import sample_event, update_network
from .stop_conditions import stop_cond
import scipy.sparse as sp
import copy as copy
from tqdm import tqdm
from .visualization import  plot_multiple_results,plot_results
from dataclasses import dataclass, field
from typing import List, Any, Dict


    

class Simulation:
    """
    The Simulation class represents a framework for epidemic modeling using GEMFcore.

    This class handles the initialization, configuration, and execution of simulations
    with customizable time structures, initial conditions, and stopping criteria. It allows
    running multiple simulations, extracting results, and plotting outcomes.

    Attributes:
    -----------
    inst : ModelConfiguration, optional
        Instance of the model configuration that defines the epidemic model's structure (compartments, transitions, etc.).
        If passing a YAML file, this parameter should not be provided.
    time_structure_option : str, optional, default is 'Auto'
        Specifies the data structure for managing absolute times ('Auto', 'Array', 'SortedList').
    initial_condition : dict, optional
        Defines the initial state distribution, such as the percentage of inducers and others.
        There are four options for specifying the initial condition:
        
        - 'percentage': dict, Defines the percentage for each compartment. For example, in the SIS model: {'I': 5, 'S': 95}, which means 5% infected and 95% susceptible, randomly distributed.
        - 'exact': np.ndarray or list, Specifies the exact state of the system. Each element represents a node with a value between 0 and M.
        - 'hubs_number': dict, Randomly sets a specific number of hubs to a certain state. For example: {'I': 10} means 10 hubs are infected, and all others are susceptible.
        
        If passing a YAML file, this parameter should not be provided.
    stop_condition : dict, optional
        Defines the condition that halts the simulation. The simulation stops if there are no more transitions (i.e., the sum of rates equals zero) or if user-defined conditions are met:
        
        - "events": int, Number of events to simulate, stopping when this threshold is reached. For example: 10000 events.
        - "time": float, Specifies the end time of the simulation. For example: 5.3 units of time.
        
        If passing a simulation YAML file, this parameter should not be provided.
    nsim : int, optional, default is 1
        Number of simulation runs to perform.

    Methods:
    --------
    __init__(inst=None, time_structure_option='Auto', initial_condition=None, stop_condition=None, nsim=1):
        Initializes the simulation with a model instance and configurations, 
        based on either YAML files containing model and simulation information or the passed parameters (inst, initial_condition, and stop_condition).

    initialize(inst=None, sim_cfg=None, time_generator_function=None, yaml=None):
        Sets up the simulation environment using either a provided model instance or YAML configuration.

    reset():
        Resets the simulation by re-initializing it with the original user-defined settings.

    _choose_method(time_structure):
        Selects the method for generating time steps based on the chosen time structure option.
        Time steps can be handled using Numpy arrays or SortedList data structures.

    forward():
        Simulates one step ahead, updating the network state based on event sampling.

    run_single_time():
        Runs a single simulation until the stop condition is met.

    run():
        Runs multiple simulations (if specified) and stores the results.

    from_yaml(cls, model_yaml, simulation_yaml):
        Class method to create a Simulation instance from YAML configuration files.

    get_results():
        Extracts and returns the simulation results, including:
        - Absolute times (np.ndarray)
        - State population history (np.ndarray)
        - Interarrival times (np.ndarray)
        - The state from which a node transitioned at each event (np.ndarray)
        - The state to which a node transitioned at each event (np.ndarray)

    stop_condition():
        Checks if the simulation should terminate based on transition rates or custom stop conditions.

    plot_results():
        Plots the simulation results for either a single run or multiple runs, if applicable.

    to_yaml(file_path):
        Saves the simulation configuration to a YAML file for reproducibility and future use.
    """

    def __init__(self, inst=None , time_structure_option='Auto',
                  initial_condition={'default_percentage':{'inducers':10, 'other':90 }},
                    stop_condition={'events':1000},
                    nsim=1):
        
        self.time_structure_option=time_structure_option
        self.inst=inst
        self.sim_cfg={
                        'time_structure_option':time_structure_option,
                        'initial_condition':initial_condition,
                        'stop_condition':stop_condition,
                        'nsim':nsim        
                    }
        self.time_generator_func=self._choose_method(self.time_structure_option)
        self.setup=None
        self.counter=np.array([0])
        self.initialize(self.inst,self.sim_cfg, self.time_generator_func)
        

    def initialize(self,inst=None, sim_cfg=None, time_generator_function=None, yaml=None):
        if not yaml and inst:
            self.setup=Initialize(inst,sim_cfg, time_generator_function, self.counter)

        else:
            print('Please pass an instant of the model to simulator!')


    def reset(self):  
        self.initialize(self.inst, self.sim_cfg, self.time_generator_func)


    def _choose_method(self,time_structure):
        if   time_structure=='SortedList':
            return TimeSorted.generate_times
        
        elif time_structure=='Array': #numpy array
            return TimeNp.generate_times
        
        elif time_structure=='Auto': 
            return None
        else:
            print("Please enter a valid option for time structure:")
            print("Array, SortedList, and Auto are only accepted methods")


    
    def forward(self):
        sample_event(self.setup.times, self.setup.event_data, self.setup.model_matrices, self.setup.rate_arrays, 
                          self.setup.X, self.setup.iteration,self.setup.Tf)
        
        update_network(self.setup.times, self.setup.event_data,
                            self.setup.model_matrices, self.setup.rate_arrays,self.setup.networks, self.setup.X)
        
    def run_signle_time(self):
        self.reset()
        try:
            while not self.stop_condition():
                self.forward()
        except  StopLoopException as e:
            print(f"Loop terminated: {e}")


    def run(self):
        results = {}
        if self.sim_cfg['nsim'] > 1:
            for i in tqdm(range(self.sim_cfg['nsim'])):
                self.run_signle_time()  
                T, statecount, *_=self.get_results()
                results[i]={'T':T, 'statecount': statecount }
        else:
            self.run_signle_time()
        self.results=results     

    @classmethod
    def from_yaml(cls, model_yaml: str, simulation_yaml: str):
        model_inst=ModelConfiguration.from_yaml(model_yaml)

        with open(simulation_yaml, 'r') as file:
            sim_cfg=yaml.safe_load(file)

        return cls( inst=model_inst,
                    time_structure_option=sim_cfg['time_structure_option'],
                    initial_condition=sim_cfg['initial_condition'],
                    stop_condition=sim_cfg['stop_condition'],
                    nsim=sim_cfg['nsim']   )
    
    
    def get_results(self):
        return post_population(self.setup.X0, self.setup.model_matrices, self.setup.event_data,
                                self.setup.networks.nodes )
        
    def stop_condition(self):
        return (self.setup.rate_arrays.R < 1e-6 or 
                stop_cond(self.setup)
                )
    
    def plot_results(self):
        if  self.sim_cfg['nsim']<2:
            T, StateCounts,*_ =self.get_results( )
            plot_results(T,StateCounts,self.inst.compartments)
        elif self.sim_cfg['nsim']>=2:
            plot_multiple_results(self.results, self.inst.compartments)
    

    def to_yaml(self, file_path: str):
        with open(file_path, 'w') as file:
            yaml.dump(self.sim_cfg, file, default_flow_style=False)

    