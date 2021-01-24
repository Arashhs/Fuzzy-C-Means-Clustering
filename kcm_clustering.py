import random, math

class Point:
    def __init__(self, values):
        self.values = values # values for each point
        self.unit_intervals = [] # unit intervals for this point in regards to each center


class KCM:
    def __init__(self, points, min_clusters_num, max_clusters_num, m, convergence_limit=0.01):
        self.points = points
        self.min_clusters_num = min_clusters_num
        self.max_clusters_num = max_clusters_num
        self.centers = []
        self.m = m
        self.convergence_limit = convergence_limit

    
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


    # calculating the distance between two points
    def calculate_distance(self, a, b):
        sum_square = 0
        for i in range(len(a)):
            sum_square += (a[i] - b[i])**2
        return math.sqrt(sum_square)


    # calculating the u_ik for X_k's unit interval in regards to the center V_i
    def update_uik(self, i, k, c):
        vi = self.centers[i]
        xk = self.points[k].values
        uik = 0
        term1 = self.calculate_distance(xk, vi)
        for j in range(c):
            vj = self.centers[j]
            term2 = self.calculate_distance(xk, vj)
            uik += (term1 / term2) ** (2/(self.m - 1))
        uik = 1/uik
        return uik


    # updating center vi
    def update_vi(self, i):
        vi = self.centers[i]
        n = len(self.points)
        m = self.m

        sigma1 = [0] * len(vi)        
        sigma2 = 0
        # calculating the denominator of vi calculation formula
        for k in range(n):
            uik = self.points[k].unit_intervals[i]
            sigma2 += uik**m

        # calculating the numerator of vi calculation formula
        for k in range(n):
            uik = self.points[k].unit_intervals[i]
            xk = self.points[k].values
            for ind in range(len(vi)):
                sigma1[ind] += (uik**m)*xk[ind]
                
        # calculating vi
        for i in range(len(vi)):
            vi[i] = sigma1[i] / sigma2

        return vi

    
    # checking if all the coordinates of the new center have converged
    def is_converged(self, old_centers, new_centers):
        limit = self.convergence_limit
        for i in range(len(new_centers)):
            for j in range(len(new_centers[0])):
                if abs(new_centers[i][j] - old_centers[i][j]) > limit:
                    return False
        return True


    
    # running the KCM algorithm with the given c and finding the centers
    def run_cluster(self, c):
        self.init_single_kcm(c)

        # while centers are not converged, run the algorithm
        while True:
            old_centers = self.centers.copy()

            # updating u_ik values
            for k in range(len(self.points)):
                for i in range(c):
                    self.points[k].unit_intervals[i] = self.update_uik(i, k, c)

            # updating v_i values (centers)
            for i in range(c):
                self.centers[i] = self.update_vi(i)

            if self.is_converged(old_centers, self.centers):
                break
        
        return self.centers

    
    # calculating the Entropy for the given points and centers
    def calculate_entropy(self):
        entropy = 0
        for i in range(len(self.centers)):
            for k in range(len(self.points)):
                uik = self.points[k].unit_intervals[i]
                entropy = entropy - uik*math.log(uik)
        return entropy



    
    # running the KCM algorithm for different number of centers and finding the appropriate one
    def kcm_cluster(self):
        entropies = []
        for i in range(self.min_clusters_num, self.max_clusters_num + 1):
            self.run_cluster(i)
            entropy = self.calculate_entropy()
            entropies.append(entropy)
            print("Number of centers: {}, Entropy: {}".format(i, entropy))
        print(self.centers)




