from conductivity_measurements import *


def G_L(LVals, p_c, nSamples):

    G = np.zeros(len(LVals))
    for i in range(len(LVals)):
        print(f"L = {LVals[i]}")
        G[i] = conductivity_density(LVals[i], [p_c], nSamples, plot = False)[0][0]

    x, y = np.log10(LVals), np.log10(G)
    a, b, a_err, b_err = lin_fit(x,y)
    decimals_a = decimals(a_err)
    x_fit = np.array([x[0], x[-1]])

    zeta = -4/3*a
    zeta_err = 4/3*a_err
    decimals_zeta = decimals(zeta_err)

    print(f"zeta_R = {zeta:.{decimals_zeta}f} +- {zeta_err:.1g}")


    plt.figure(num=0, dpi=80, facecolor='w', edgecolor='k')
    plt.plot(LVals, G,'o', label='$G$')
    plt.plot(10**x_fit, 10**(a*x_fit + b), label = f"Linear fit (log-plot)\nSlope = {a:.{decimals_a}f} " r"$\pm$" f" {a_err:.1g}")
    plt.xlabel(r"$L$", fontsize=14)
    plt.ylabel(r"$G(p=p_c, L)$", fontsize=14)
    plt.xscale("log")
    plt.yscale("log")
    plt.legend(fontsize = 13)
    plt.tight_layout(pad=1.1, w_pad=0.7, h_pad=0.2)


    #
    #
    # x, y = np.log10(pVals-p_c), np.log10(G)
    # a, b, a_err, b_err = lin_fit(x,y)
    # decimals_a = decimals(a_err)
    # x_fit = np.array([x[0], x[-1]])
    #
    # plt.figure(num=0, dpi=80, facecolor='w', edgecolor='k')
    # plt.plot(pVals-p_c, G,'o', label='$G$')
    # plt.plot(10**x_fit, 10**(a*x_fit + b), label = f"Linear fit (log-plot)\nSlope = {a:.{decimals_a}f} " r"$\pm$" f" {a_err:.1g}")
    # plt.xlabel(r"$p-p_c$", fontsize=14)
    # plt.ylabel(r"$G$", fontsize=14)
    # plt.xscale("log")
    # plt.yscale("log")
    # plt.legend(fontsize = 13)
    # plt.tight_layout(pad=1.1, w_pad=0.7, h_pad=0.2)



if __name__ == "__main__":
    LVals = 2**np.arange(5, 11)
    nSamples = 100
    p_c = 0.59275

    G_L(LVals, p_c, nSamples)
    # plt.savefig("../article/figures/G_L.pdf", bbox_inches="tight")
    plt.savefig("../article/figures/G_L_alt.pdf", bbox_inches="tight")

    plt.show()
