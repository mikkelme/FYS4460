from cluster_generator import *
from plot_set import *

def Find_M(L, p, MC_cycles):
    M = np.zeros(len(L))

    for i in range(len(L)):
        print(f"\r {i}/{len(L)}, L = {L[i]}", end = "")
        cluster = TwoDim_cluster(L[i],p)
        M_SC = cluster.exwalk()
        for j in range(1, MC_cycles):
            cluster.generate_cluster()
            M_SC += cluster.exwalk()
        M[i] = M_SC/MC_cycles # Average
    print(f"\r {i+1}/{len(L)}, L = {L[i]}")

    return M

def Find_D_SC(L, M):
    x, y = np.log10(L), np.log10(M)
    a, b, a_err, b_err = lin_fit(x,y)

    decimals_a = int(np.ceil(-np.log10(a_err)))
    print(f"D_SC = {a:.{decimals_a}f} +- {a_err:.1g}")

    x_fit = np.linspace(np.min(x), np.max(x), int(1e3))
    plt.figure(num=0, dpi=80, facecolor='w', edgecolor='k')
    plt.plot(L, M, "o", label = "Datapoints")
    plt.plot(10**x_fit, 10**(a*x_fit + b), label = f"Linear fit (log-plot)\nSlope = {a:.{decimals_a}f} " r"$\pm$" f" {a_err:.1g}")
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel(r"L", fontsize=14)
    plt.ylabel(r"$\langle M_{SC} \rangle (L)$", fontsize=14)
    plt.legend(fontsize = 13)
    plt.tight_layout(pad=1.1, w_pad=0.7, h_pad=0.2)








if __name__ == "__main__":
    L = [25, 50, 100, 200, 400, 800]
    p = 0.59275
    MC_cycles = 200

    M = Find_M(L, p, MC_cycles)
    Find_D_SC(L, M)

    plt.savefig("../article/figures/q.pdf", bbox_inches="tight")
    plt.show()
