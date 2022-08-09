from cluster_generator import *
from plot_set import *


def plot(L_try, p, MC_cycles):

    M = np.zeros(len(L_try))
    for i in range(len(L_try)):
        print(f"\r {i+1}/{len(L_try)}, L = {L_try[i]}", end = "")
        cluster = TwoDim_cluster(L_try[i], p)
        M[i] = cluster.percolating_cluster_mass(MC_cycles)
    print()

    a, b, a_err, b_err = lin_fit(np.log10(L_try),np.log10(M))
    decimals_a = int(np.ceil(-np.log10(a_err)))

    plt.figure(num=0, dpi=80, facecolor='w', edgecolor='k')
    plt.plot(L_try, M, "o", label = "Datapoints")
    plt.plot(L_try, L_try**a * 10**b, "--",  label = f"Linear fit (log-plot)\nSlope = {a:.{decimals_a}f} " r"$\pm$" f" {a_err:.1g}")
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel(r"$L$", fontsize=14)
    plt.ylabel(r"$M$", fontsize=14)
    plt.legend(fontsize = 13)
    plt.tight_layout(pad=1.1, w_pad=0.7, h_pad=0.2)

if __name__ == "__main__":
    p_c = 0.59275
    L_try = 2**np.arange(4,11+1)
    MC_cycles = 500

    plot(L_try, p_c, MC_cycles)
    plt.savefig("../article/figures/i.pdf", bbox_inches="tight")
    plt.show()
