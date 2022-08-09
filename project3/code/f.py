from cluster_generator import *
from plot_set import *



def plot(L, a, MC_cycles, p_try):
    """ Plots (n,s) for multiple p values
        in p_try averaged over some MC_cycles """


    plt.figure(num=0, dpi=80, facecolor='w', edgecolor='k')

    if isinstance(p_try, (float,int)): # single value
        cluster = TwoDim_cluster(L,p_try)
        s, n = cluster.density(a, MC_cycles)
        plt.plot(s, n, "--", marker = "o", label = f"p = {p_try:g}")

    else: # array, list, tuple
        p_num = 0
        print(f"\r {p_num}/{len(p_try)}", end = "")
        for p in p_try:
            p_num += 1
            cluster = TwoDim_cluster(L,p)
            s, n = cluster.density(a, MC_cycles)
            label = f"p = {p:g}"
            if p == p_c:
                label = r"$p_c$" + f" = {p:g}"
            print(f"\r {p_num}/{len(p_try)}", end = "")
            plt.plot(s, n, "--", marker = "o", label = label)
        print()

    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel(r"$s$", fontsize=14)
    plt.ylabel(r"$n(s,p)$", fontsize=14)
    plt.legend(fontsize = 13)
    plt.tight_layout(pad=1.1, w_pad=0.7, h_pad=0.2)



if __name__ == "__main__":
    L = 1000
    a = 1.6
    MC_cycles = 100

    p_c = 0.59275
    diff = 0.01
    span = 0.05
    p_below = np.linspace(p_c - (diff + span), p_c - diff, 5)
    p_above = np.linspace(p_c + (diff + span), p_c + diff, 5)


    plot(L, a, MC_cycles, p_below)
    #plt.savefig("../article/figures/f_p_below.pdf", bbox_inches="tight")
    plt.show()
    exit()
    print("Next one")

    plot(L, a, MC_cycles, p_above)
    #plt.savefig("../article/figures/f_p_above.pdf", bbox_inches="tight")
    plt.show()
