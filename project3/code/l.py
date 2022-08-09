from cluster_generator import *
from plot_set import *
from scipy.interpolate import interp1d



def PI(L, p_try, MC_cycles):
    if isinstance(p_try, (float,int)): # single value
        # print("CLOSED SECTION")
        # exit()
        cluster = TwoDim_cluster(L, p_try)
        perculation = np.zeros(MC_cycles)
        for j in range(MC_cycles):
            perc_labels = cluster.find_percolating_cluster()
            perculation[j] = len(perc_labels) > 0
            cluster.generate_cluster()
        pi = np.mean(perculation)
        return pi

    else: # array, list, tuple
        pi = np.zeros(len(p_try))
        for i in range(len(p_try)):
            cluster = TwoDim_cluster(L, p_try[i])
            perculation = np.zeros(MC_cycles)
            for j in range(MC_cycles):
                perc_labels = cluster.find_percolating_cluster()
                perculation[j] = len(perc_labels) > 0
                cluster.generate_cluster()
            pi[i] = np.mean(perculation)
    return pi

def PI_equal_x(x, L_try, p_span, MC_cycles):
    pi_eq_x = np.zeros(len(L_try))
    for i in range(len(L_try)):
        print(f"\r {i}/{len(L_try)}, x = {x}, L = {L_try[i]}", end = "")
        pi = PI(L_try[i], p_span, MC_cycles)
        ipc = np.argmax(pi > x) #Find first value where pi > x

        #Interpolate from icp-1 to ipc to find intersections
        slope = (pi[ipc]-pi[ipc-1])/(p_span[ipc] - p_span[ipc-1])
        gain = x - pi[ipc-1]
        pi_eq_x[i] = p_span[ipc-1] + gain/slope

        plot = False
        if plot:
            curve = plt.plot(p_span, pi, label = f"L = {L_try[i]}")
            plt.plot(pi_eq_x[i], x, "o", color = curve[0].get_color())
    print(f"\r {i+1}/{len(L_try)}, x = {x}, L = {L_try[i]}")
    return pi_eq_x

def plot(L_try, p_span, X, PI_MC_cycles, save_data):
    if save_data:
        np.save("l_data/L_try.npy", L_try)

    plt.figure(num=0, dpi=80, facecolor='w', edgecolor='k')
    for x in X:
        pi_eq_x = PI_equal_x(x, L_try, p_span, PI_MC_cycles)
        plt.plot(L_try, pi_eq_x, "o", label = f"x = {x}")
        if save_data:
            np.save("l_data/PI_eq_" + str(x) + ".npy", pi_eq_x)

    plt.xlabel(r"$L$", fontsize=14)
    plt.ylabel(r"$p_{\Pi=x}$", fontsize=14)
    plt.legend(fontsize = 13)
    plt.tight_layout(pad=1.1, w_pad=0.7, h_pad=0.2)
    # plt.savefig("../article/figures/l.pdf", bbox_inches="tight")






if __name__ == "__main__":
    L_try = [25, 50, 100, 200, 400, 800]
    p_span = np.linspace(0.4,0.75, 50)
    PI_MC_cycles = 200
    X = [0.3, 0.8]

    plot(L_try, p_span, X, PI_MC_cycles, save_data = True)
    plt.savefig("../article/figures/l.pdf", bbox_inches="tight")
    plt.show()


    # load test
    # L = np.load("l_data/L_try.npy")
    # pi_eq_1  = np.load("l_data/PI_eq_0.3.npy")
    # pi_eq_2  = np.load("l_data/PI_eq_0.8.npy")
    # print(L, pi_eq_1, pi_eq_2)



#
