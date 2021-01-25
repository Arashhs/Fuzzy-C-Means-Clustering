import random, math, copy
import matplotlib.pyplot as plt


class Point:
    def __init__(self, values):
        self.values = values # values for each point
        self.unit_intervals = [] # unit intervals for this point in regards to each center
        self.label = None # label that shows each point belongs to which cluster


class KCM:
    def __init__(self, points, min_clusters_num, max_clusters_num, m, convergence_limit=0.005):
        self.points = points
        self.min_clusters_num = min_clusters_num
        self.max_clusters_num = max_clusters_num
        self.centers = []
        self.m = m
        self.convergence_limit = convergence_limit
        self.clusters = []

    
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
            old_centers = copy.deepcopy(self.centers)

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
        c = len(self.centers)
        for i in range(len(self.centers)):
            for k in range(len(self.points)):
                uik = self.points[k].unit_intervals[i]
                entropy = entropy - uik*math.log(uik)
        return entropy/math.log(c)


    # label each point to show the cluster it belongs to
    def label_points(self):
        for point in self.points:
            max_index = 0
            for i in range(len(point.unit_intervals)):
                if point.unit_intervals[i] > point.unit_intervals[max_index]:
                    max_index = i
            point.label = max_index


    def build_clusters(self):
        self.label_points()
        self.clusters = [[] for i in range(len(self.centers))]
        for point in self.points:
            self.clusters[point.label].append(point)

    
    # printing the centers
    def print_centers(self):
        centers = self.centers
        for i in range(len(centers)):
            pr = "Center {}: ".format(i+1)
            for j in range(len(centers[i])):
                pr += str(centers[i][j]) + ', '
            pr = pr[:-1]
            print(pr)



    
    # running the KCM algorithm for different number of centers and finding the appropriate one
    def kcm_cluster(self):
        entropies = []
        all_points = []
        all_centers = []
        all_clusters = []
        for i in range(self.min_clusters_num, self.max_clusters_num + 1):
            self.run_cluster(i)
            entropy = self.calculate_entropy()
            entropies.append(entropy)
            print("Number of centers: {}, Entropy: {}".format(i, entropy))
            self.build_clusters()
            all_points.append(copy.deepcopy(self.points))
            all_centers.append(copy.deepcopy(self.centers))
            all_clusters.append(copy.deepcopy(self.clusters))
        min_index = 0
        for i in range(len(entropies)):
            if entropies[i] < entropies[min_index]:
                min_index = i
        self.points = all_points[min_index]
        self.centers = all_centers[min_index]
        self.clusters = all_clusters[min_index]
        print("\nBest number of centers:", min_index + self.min_clusters_num)
        self.print_centers()

    
    # plot the result
    def kcm_plot(self):
        colors = ["green","blue","yellow","pink","black","orange","purple","beige","brown","gray","cyan","magenta"]
        centers_x = [c[0] for c in self.centers]
        centers_y = [c[1] for c in self.centers]
        for i in range(len(self.clusters)):
            x_values = [y.values[0] for y in [x for x in self.clusters[i]]]
            y_values = [y.values[1] for y in [x for x in self.clusters[i]]]
            plt.scatter(x_values, y_values, c=colors[i])
        plt.scatter(centers_x, centers_y, c="red", s=40)
        plt.show()







