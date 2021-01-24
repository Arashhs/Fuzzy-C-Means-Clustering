import random

class Point:
    def __init__(self, values):
        self.values = values # values for each point
        self.unit_intervals = [] # unit intervals for this point in regards to each center


class KCM:
    def __init__(self, points, min_clusters_num, max_clusters_num, m):
        self.points = points
        self.min_clusters_num = min_clusters_num
        self.max_clusters_num = max_clusters_num
        self.centers = []
        self.m = m

    
    # simply generating random centers
    def init_centers(self, c, dim):
        ret_centers = []
        for _ in range(c):
            rnd = [random.random() for i in range(dim)]
            ret_centers.append(rnd)
        return ret_centers

    
    # resetting unit intervals
    def clear_unit_intervals(self, c):
        for point in self.points:
            point.unit_intervals = [0] * c

    
    # initializing our environment for running KCM with the new number of centers
    def init_single_kcm(self, c):
        self.centers = self.init_centers(c, len(self.points[0].values))
        self.clear_unit_intervals(c)
        

    
    # running the KCM algorithm with the given c and finding the centers
    def run_cluster(self, c):
        self.init_single_kcm(c)


    
    # running the KCM algorithm for different number of centers and finding the appropriate one
    def kcm_cluster(self):
        self.run_cluster(3)




