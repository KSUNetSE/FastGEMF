import scipy.sparse as sp
import networkx as nx
import numpy as np

def edgelist_to_csr(file_path, directed=False):
    """
    Converts an edge list from a file to a CSR matrix using NetworkX library functions.
    
    Parameters:
    file_path (str): Path to the file containing the edge list.
    directed (bool): Whether the graph is directed or not.
    
    Returns:
    scipy.sparse.csr_matrix: CSR matrix representation of the graph.
    """
    G = nx.DiGraph() if directed else nx.Graph()
    with open(file_path, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) == 3:
                u, v, w = map(int, parts)
                G.add_edge(u, v, weight=w)

    G_csr = nx.to_scipy_sparse_array(G, format='csr', weight='weight')
    
    return G_csr


    

class Network:

    def __init__(self, mats):
        """
        Initializes the Network with adjacency matrices or CSR matracies for each layer.

        Args:
            mats (list): A list of  adjacency matrices or CSR matrices for different layers.
        """
        self.csrMats = [sp.csr_matrix(mat) for mat in mats]
        self.cscMats = [mat.tocsc() for mat in self.csrMats]
        self.nodes = self.csrMats[0].shape[0] if mats else 0
        self.edges = [network.nnz for network in self.csrMats ]
        self.weights = self.csrMats[0].data

    def get_out_neighbors(self, l, n):
        """
        Get the outgoing neighbors of a node in layer l.

        Args:
            l (int): Layer index.
            n (int): Node index.

        Returns:
            np.ndarray: Indices of nodes connected by outgoing edges from node n.
        """        
        start_pointer = self.csrMats[l].indptr[n]
        end_pointer = self.csrMats[l].indptr[n + 1]
        return self.csrMats[l].indices[start_pointer:end_pointer]
    
    def get_in_neighbors(self, l, n):
        """
        Get the incoming neighbors of a node in layer l.

        Args:
            l (int): Layer index.
            n (int): Node index.

        Returns:
            np.ndarray: Indices of nodes connected by incoming edges to node n.
        """
        start_pointer, end_pointer = self.cscMats[l].indptr[n], self.cscMats[l].indptr[n + 1]
        return self.cscMats[l].indices[start_pointer:end_pointer]

    def get_out_weight(self, l, n):
        """
        Get the weights of outgoing edges from a node in layer l.

        Args:
            l (int): Layer index.
            n (int): Node index.

        Returns:
            np.ndarray: Weights of the outgoing edges from node n.
        """
        start_pointer = self.csrMats[l].indptr[n]
        end_pointer = self.csrMats[l].indptr[n + 1]
        return self.csrMats[l].data[start_pointer:end_pointer]

    def get_in_weight(self, l, n):
        """
        Get the weights of incoming edges to a node in layer l.

        Args:
            l (int): Layer index.
            n (int): Node index.

        Returns:
            np.ndarray: Weights of the incoming edges to node n.
        """
        start_pointer = self.cscMats[l].indptr[n]
        end_pointer = self.cscMats[l].indptr[n + 1]
        return self.cscMats[l].data[start_pointer:end_pointer]
    
    def get_highest_degree_nodes(self, l, m):
        """
            Returns the `m` highest degree nodes in layer `l` of the network.
            
            Args:
                l (int): The index of the layer (0-indexed).
                m (int): The number of top degree nodes to return.
            
            Returns:
                numpy.ndarray: An array of the `m` highest degree nodes in layer `l`.

        Warning: In large networks, this inilization can significantly increase the running time. It requires sorting all node degrees.

        """
        degrees = self.csrMats[l].sum(axis=1).A1
        
        hubs = np.argsort(degrees)[-m:]
        
        return hubs



    