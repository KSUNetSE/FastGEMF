
__version__ = '0.1.4'

r"""
FastGEMF: Scalable High-Speed Simulation of Stochastic Spreading Processes over Complex Multilayer Networks
"""
from .GEMFsimulation import Simulation
from .network import edgelist_to_csr, Network
from .modelschema import ModelSchema
from .modelconfiguration import ModelConfiguration
