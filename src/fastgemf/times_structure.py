import numpy as np
from sortedcontainers import SortedList

class TimeNp:
    """
    Using NumPy Library for sampling the nodes(storing absolute times and finding the node with min time)

    Attributes:
    -----------
    times_array : np.ndarray
        Array storing the absolute times for nodes.

    Methods:
    --------
    generate_times(rates):
        Class method that initializes a TimeNp instance by generating times for each node based on the provided rates.

    get_min():
        Returns the node with the minimum event time and its corresponding absolute time.

    pop_and_push(node, tmin, Rn_node):
        Updates the time for a node after an event and pushes the new time back into the array.

    update(affected_nodes, tmin, Rn_affected):
        Updates the abosolute times for affected nodes.

    __len__():
        Returns the number of nodes being tracked in the times array.
    """

    def __init__(self, rates):
        self.times_array = -np.log(np.random.rand(len(rates))) / np.absolute(rates)
    
    @classmethod
    def generate_times(cls,rates):
        return cls(rates)

    def get_min(self):
        sampled_node=np.argmin(self.times_array) 
        return sampled_node, self.times_array[sampled_node]
    
    def pop_and_push(self,node, tmin, Rn_node):
        self.times_array[node]=tmin - np.log(np.random.rand()) / np.absolute(Rn_node)

    def update(self, affected_nodes: np.ndarray, tmin:float, Rn_affected: np.ndarray ):
        if  isinstance(affected_nodes,(float, int, list, tuple)):
            affected_nodes=np.array(affected_nodes)
            self.times_array[affected_nodes] = tmin - np.log(np.random.rand(len(affected_nodes))) / np.absolute(Rn_affected)
        else:
            self.times_array[affected_nodes] = tmin - np.log(np.random.rand(len(affected_nodes))) / np.absolute(Rn_affected)            

    def __len__(self):
        return len(self.times_array)
    
class TimeSorted:
    """
    Using SortedList Library for sampling the nodes(finding the node with min abosulte time)

    Attributes:
    -----------
    times_array : np.ndarray
        Array storing the abosulte times for each node. for keeping the track of absolute times and their corresponding nodes.

    times_sorted : SortedList
        Sorted list of nodes and their abosulte times.

    minArg : int
        The node with the current minimum abosulte time.

    minValue : float
        The minimum abosulte time.

    Methods:
    --------
    generate_times(Rn):
        Class method that initializes a TimeSorted instance by generating abosulte times for each node based on the provided rates.

    get_min():
        Returns the node with the minimum abosulte time and the corresponding time.

    pop_and_push(node, tmin, Rn_node):
        Updates the abosulte time for a node and reinserts it into the sorted list.

    update(affected_nodes, tmin, Rn_affected):
        Updates the abosulte times for affected nodes in the sorted list.

    __len__():
        Returns the number of nodes being tracked in the sorted list.
    """
    def __init__(self, Rn):

        self.times_array = -np.log(np.random.rand(len(Rn))) / np.absolute(Rn)
        self.times_sorted = SortedList(enumerate([]), key=lambda x: x[1])
        
        for node in Rn.nonzero()[0]:
            self.times_sorted.add((node, self.times_array[node]))
        
        self.minArg = self.times_sorted[0][0]
        self.minValue = self.times_sorted[0][1]

    @classmethod
    def generate_times(cls, Rn):
        return cls(Rn)

    def get_min(self):
        return self.minArg, self.minValue
    
    def pop_and_push(self,node, tmin, Rn_node):
        new_time=tmin - np.log(np.random.rand()) / np.absolute(Rn_node)

        if new_time==float('inf') and  self.times_sorted:
            self._remove_node(node,new_time)

        elif  self.times_sorted:
            self._update_node(node, new_time)
            
        if not self.times_sorted:
            print('\nInfection has ended. All nodes are removed!')
        else:
            self.minValue = self.times_sorted[0][1]
            self.minArg = self.times_sorted[0][0] 

    def update(self, affected_nodes: np.ndarray, tmin:float, Rn_affected: np.ndarray):

        new_times = tmin - np.log(np.random.rand(len(affected_nodes))) / np.absolute(Rn_affected)
    
        for k, n in enumerate(affected_nodes):

            if new_times[k] == float('inf') and  self.times_sorted:
                self._remove_node(n,new_times[k])

            elif self.times_array[n] == float('inf') and  self.times_sorted:
                self._add_node(n, new_times[k])

            elif  self.times_sorted:
                self._update_node(n, new_times[k])
        if not self.times_sorted:

            print('\nInfection has ended. All nodes are removed!')
        else:                        
            self.minValue = self.times_sorted[0][1]
            self.minArg = self.times_sorted[0][0]

    def _remove_node(self, node, new_time):
        self.times_sorted.discard((node, self.times_array[node]))
        self.times_array[node] = new_time

    def _update_node(self, node, new_time):
        self.times_sorted.discard((node, self.times_array[node]))
        self.times_sorted.add((node, new_time))
        self.times_array[node] = new_time

    def _add_node(self, node, new_time):
        self.times_sorted.add((node, new_time))
        self.times_array[node] = new_time

    def __len__(self):
        return len(self.times_sorted)
        
