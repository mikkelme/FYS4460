import numpy as np
import matplotlib
from plot_set import *






if __name__ == "__main__":
    L = np.load("l_data/L_try.npy")
    pi_eq_1  = np.load("l_data/PI_eq_0.3.npy")
    pi_eq_2  = np.load("l_data/PI_eq_0.8.npy")

    x = np.log10(L)
    y = np.log10(pi_eq_2 - pi_eq_1)
    a, b, a_err, b_err = lin_fit(x,y)

    decimals_a = int(np.ceil(-np.log10(a_err)))
    nu = -1/a
    nu_std = a_err/a**2
    decimals_nu = int(np.ceil(-np.log10(nu_std)))
    print(f"Sigma = {nu:.{decimals_nu}f} +- {nu_std:.1g}")
    x_fit = np.linspace(x.min(), x.max(), int(1e4))

    plt.figure(num=0, dpi=80, facecolor='w', edgecolor='k')
    plt.plot(L, pi_eq_2 - pi_eq_1, "o")
    plt.plot(10**x_fit, 10**(a*x_fit + b), label = f"Linear fit (log-plot)\nSlope = {a:.{decimals_a}f} " r"$\pm$" f" {a_err:.1g}")

    plt.xscale("log")
    plt.yscale("log")

    plt.xlabel(r"$L$", fontsize=14)
    plt.ylabel(r"$p_{\Pi = 0.8} - p_{\Pi = 0.3}$", fontsize=14)
    plt.legend(fontsize = 13)
    plt.tight_layout(pad=1.1, w_pad=0.7, h_pad=0.2)
    plt.savefig("../article/figures/m.pdf", bbox_inches="tight")
    plt.show()
