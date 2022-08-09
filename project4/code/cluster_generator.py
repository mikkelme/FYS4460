import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import measurements

class TwoDim_cluster:
    def __init__(self, L, p):
        if p < 0:
            print("Input error: p most be a positive number")
            exit()
        self. L = L
        self.p = p
        self.generate_cluster()

    def generate_cluster(self):
        # np.random.seed(974412678)
        self.z = np.random.rand(self.L,self.L)
        self.m = self.z < self.p
        self.lw, self.num = measurements.label(self.m)
        self.area = measurements.sum(self.m, self.lw, index = np.arange(self.lw.max()+1))

    def show(self):
        areaImg = self.area[self.lw]
        plt.imshow(areaImg, origin = "upper")
        plt.colorbar()
        plt.show()

    def find_percolating_cluster(self):
        """ find percolating cluster
            and returns ...             """
        #Find percolating cluster labels
        horisontal_connect  = np.intersect1d(self.lw[:,0], self.lw[:,-1])
        vertical_connect = np.intersect1d(self.lw[0,:], self.lw[-1,:])
        union_connect = np.union1d(horisontal_connect, vertical_connect)
        perc_labels = union_connect.ravel()[np.flatnonzero(union_connect)]
        return perc_labels

    def is_spanning(self):
        perc_labels = self.find_percolating_cluster()
        perculation = len(perc_labels) > 0
        return perculation

    def rm_percolating_cluster(self):
        """ removes percolating cluster
            from m, lw & area arrays    """
        perc_labels = self.find_percolating_cluster()
        perculation = len(perc_labels) > 0

        if perculation: #Remove possible perculation cluster
            perc_idx = np.argwhere(np.isin(self.lw, perc_labels))
            self.m[perc_idx[:,0], perc_idx[:,1]] = False
            self.lw, self.num = measurements.label(self.m)
            self.area = measurements.sum(self.m, self.lw, index = np.arange(self.lw.max()+1))


    def percolating_cluster_mass(self, MC_cycles):
        #Collect perculating cluster areas (M)
        M = np.zeros(MC_cycles)
        for i in range(MC_cycles):
            perc_labels = self.find_percolating_cluster()
            #Handles more than one perc_cluster by taking mean
            M[i] = np.mean(self.area[perc_labels])
            self.generate_cluster()
        return np.mean(M[np.isfinite(M)])


    def maximal_logbin_size(self, a):
        """ Calculate maximal base exponent for
            a^x logaritmic binning for L X L system
            to be used in self.logarithmic_bins()      """
        sum = 0; x = 0
        while sum < self.L**2:
            x += 1
            sum += a**x
        self.logbin_max = x


    def logarithmic_bins(self, s, a):
        """ Return logaritmic bins
            based on s-array         """
        #Maximal bin array estimate using s_max <= L^2
        tmp_bins = np.cumsum(a**np.arange(self.logbin_max))

        #Find number of bins corresponding to s_max
        num_bin_edges = len(tmp_bins[tmp_bins < s.max()]) + 1

        #Calculate bin proporties
        bin_edges = np.cumsum(a**np.arange(num_bin_edges))
        bin_width = a**np.arange(num_bin_edges)[1:]
        bin_center = bin_edges[0:-1] + 1/2*bin_width
        return bin_edges, bin_width, bin_center


    def density(self, a, MC_cycles):
        """ Return cluster density n """
        self.rm_percolating_cluster()

        #Collect cluster areas (s)
        s = self.area.tolist()
        for i in range(1, MC_cycles):
            self.generate_cluster()
            self.rm_percolating_cluster()
            s += self.area.tolist()
        s = np.array(s)[np.nonzero(s)] #remove zeros

        #Define histogram binning
        self.maximal_logbin_size(a)
        bin_edges, bin_width, bin_center = self.logarithmic_bins(s, a)

        #Calculate histogram and cluster density (n)
        hist, hist_bins = np.histogram(s, bin_edges)
        n = hist/bin_width/MC_cycles/self.L**2

        return bin_center, n

    def exwalk(self, plot = False):
        """ Walks the exterior paths of the
            spanning cluster (from top to bottom)
            with a right and a left turning walker. """
        from walk import walk

        #--- Ensure top-bottom spanning cluster ---#
        ncount = 0
        perc_labels = []
        while (len(perc_labels) == 0):
            self.generate_cluster()
            ncount += 1
            if (ncount > 1000):
                print("Could not make percolation cluster...")
                break
            vertical_connect = np.intersect1d(self.lw[0,:], self.lw[-1,:])
            perc_labels = vertical_connect[np.where(vertical_connect > 0)]

        if len(perc_labels) > 0: #Perform walk on spanning cluster
            zz = (self.lw == perc_labels[0]) #Chose first possible spanning cluster

            #--- Walk ---#
            left_walker,right_walker = walk(zz) # Perform walk
            SC = left_walker * right_walker     # Single Connected sites (where l*r is nonzero)
            M_SC = np.sum(SC > 0)               # Mass of Single Connected sites
            zz_SC= zz + SC                      # Spanning cluster + Single Connected sites
            BB = left_walker + right_walker

            if plot:
                plt.figure(num=5, dpi=80, facecolor='w', edgecolor='k')
                plt.title("Right walker")
                plt.imshow(right_walker, interpolation='nearest', origin='upper') # Display spanning cluster
                plt.colorbar()

                plt.figure(num=4, dpi=80, facecolor='w', edgecolor='k')
                plt.title("Left walker")
                plt.imshow(left_walker, interpolation='nearest', origin='upper') # Display spanning cluster
                plt.colorbar()

                plt.figure(num=3, dpi=80, facecolor='w', edgecolor='k')
                plt.title("Backbone (left walker + right walker)")
                plt.imshow(BB, interpolation='nearest', origin='upper') # Display spanning cluster
                plt.colorbar()

                plt.figure(num=2, dpi=80, facecolor='w', edgecolor='k')
                plt.title(r"Single connected bonds ($M_{SC} =$ " + f"{M_SC})")
                plt.imshow(SC, interpolation='nearest', origin='upper') # Display spanning cluster
                plt.colorbar()

                plt.figure(num=1, dpi=80, facecolor='w', edgecolor='k')
                plt.title("Spanning Cluster + Single connected bonds")
                plt.imshow(zz_SC, interpolation='nearest', origin='upper') # Display spanning cluster
                plt.colorbar()

                plt.show()

            return M_SC


if __name__ == "__main__":
    L = 100
    p_c = 0.59275
    p = 0.59275


    cluster = TwoDim_cluster(L,p)
    cluster.exwalk(True)








# b = np.arange(1, lw.max()+1)
# np.random.shuffle(b)
# b = np.concatenate((np.array([0]), b))
# shuffledLw = b[lw]
