from conductivity_measurements import *


def mu(G, pVals, p_c):
    x, y = np.log10(pVals-p_c), np.log10(G)
    a, b, a_err, b_err = lin_fit(x,y)
    decimals_a = decimals(a_err)
    x_fit = np.array([x[0], x[-1]])

    plt.figure(num=0, dpi=80, facecolor='w', edgecolor='k')
    plt.plot(pVals-p_c, G,'o', label='$G$')
    plt.plot(10**x_fit, 10**(a*x_fit + b), label = f"Linear fit (log-plot)\nSlope = {a:.{decimals_a}f} " r"$\pm$" f" {a_err:.1g}")
    plt.xlabel(r"$p-p_c$", fontsize=14)
    plt.ylabel(r"$G$", fontsize=14)
    plt.xscale("log")
    plt.yscale("log")
    plt.legend(fontsize = 13)
    plt.tight_layout(pad=1.1, w_pad=0.7, h_pad=0.2)



if __name__ == "__main__":
    L = 400
    nSamples = 100

    p_c = 0.59275
    p_start = 0.61; p_end = 0.85
    p_num = 20
    pVals = p_c + np.logspace(np.log10(p_start-p_c), np.log10(p_end-p_c), p_num)

    G, P = conductivity_density(L, pVals, nSamples, plot = False)
    mu(G, pVals, p_c)

    #plt.savefig("../article/figures/mu.pdf", bbox_inches="tight")
    plt.show()
