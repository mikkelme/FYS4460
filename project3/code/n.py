import numpy as np
import matplotlib
from plot_set import *


if __name__ == "__main__":
    L = np.load("l_data/L_try.npy")
    nu = 4/3
    pi_eq_1  = np.load("l_data/PI_eq_0.3.npy")
    pi_eq_2  = np.load("l_data/PI_eq_0.8.npy")
    label_x = [0.3, 0.8]

    x = L**(-1/nu)
    y = [pi_eq_1, pi_eq_2]
    b = np.zeros(len(y))
    b_err = np.zeros(len(y))

    plt.figure(num=0, dpi=80, facecolor='w', edgecolor='k')
    for i in range(len(y)):
        a, b[i], a_err, b_err[i] = lin_fit(x, y[i])
        decimals_b = int(np.ceil(-np.log10(b_err[i])))
        # label = "$\Pi_{x=" + str(label_x[i]) + "}$"
        curve = plt.plot(x, y[i], "o", label = f"x = {label_x[i]}")
        plt.plot(x, x*a + b[i], color = curve[0].get_color(), label = f"Linear fit\nIntersection = {b[i]:.{decimals_b}f} " r"$\pm$" f" {b_err[i]:.1g}")
    plt.xlabel(r"$L^{-1/\nu}$", fontsize=14)
    plt.ylabel(r"$p_{\Pi=x}$", fontsize=14)
    plt.legend(fontsize = 13)
    plt.tight_layout(pad=1.1, w_pad=0.7, h_pad=0.2)

    pc = np.mean(b)
    pc_err = np.sqrt(np.sum(b_err**2))
    decimals_pc = int(np.ceil(-np.log10(pc_err)))
    print(f"Pc = {pc:.{decimals_pc}f} +- {pc_err:.1g}")




    plt.savefig("../article/figures/n.pdf", bbox_inches="tight")
    plt.show()
