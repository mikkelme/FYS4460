from cluster_generator import *
from plot_set import *
from scipy.interpolate import interp1d






def s_xi(p, L, a, s_pc, n_pc, plot = True):
    cluster = TwoDim_cluster(L, p)
    s, n = cluster.density(a, MC_cycles)
    end_idx = np.min((len(s),len(s_pc)))

    n_over_npc =  n[:end_idx]/n_pc[:end_idx]
    isfinite = np.isfinite(np.log10(n_over_npc))
    x, y = np.log10(s[:end_idx][isfinite]), np.log10(n_over_npc[isfinite])
    f = interp1d(x,y , kind = "linear")
    s_log = np.linspace(np.log10(s[:end_idx][isfinite][0]), np.log10(s[:end_idx][isfinite][-1]), int(1e3))
    xi_candidates = s_log[np.where(10**f(s_log) < 0.5)]


    if plot:
        points = plt.plot(s[isfinite], n_over_npc[isfinite], "o", label = f"p = {p:g}")
        plt.plot(10**s_log, 10**f(s_log), "--", color = points[0].get_color())

    if len(xi_candidates) < 3:
        if len(s) < len(s_pc):
            #plt.plot(np.abs(p-0.59275), s[isfinite][-1], "o")
            return s[:end_idx][isfinite][-1]
        else:
            return np.nan
    else:
        return 10**np.min(xi_candidates)


def s_xi_p(p_try, L, a, s_pc, n_pc, plot = True):
    p_c = 0.59275
    xi = np.zeros(len(p_try))
    for i in range(len(p_try)):
        print(f"\r {i+1}/{len(p_try)}, p = {p_try[i]:g}", end = "")
        xi_element = s_xi(p_try[i], L, a, s_pc, n_pc, plot)
        if np.isfinite(xi_element):
            xi[i] = xi_element
        else:
            print(f"break at: p =  {p_try[i]:g}")
            break
    print()

    if plot:
        plt.hlines(y = 0.5, xmin = np.min(s_pc), xmax = np.max(s_pc), label = r"$\frac{n(s,p)}{n(s,p_c)}= 0.5$")
        plt.plot(xi, xi*0 + 0.5, color = "black", linestyle = "None", marker = "o", mfc = "None", label = r"$s_{\xi}$")
        plt.xscale("log")
        plt.yscale("log")
        plt.xlabel(r"$s$", fontsize=14)
        plt.ylabel(r"$n(s,p) \ / \ n(s,p_c)$", fontsize=14)
        plt.legend(fontsize = 13)
        plt.tight_layout(pad=1.1, w_pad=0.7, h_pad=0.2)
        plt.savefig("../article/figures/h_n_over_nc.pdf", bbox_inches="tight")
        plt.show()

    x, y = np.log10(np.abs(p_try-p_c)), np.log10(xi)
    isfinite = np.isfinite(y)
    x,y = x[isfinite],y[isfinite]
    a, b, a_err, b_err = lin_fit(x,y)

    decimals_a = int(np.ceil(-np.log10(a_err)))
    sigma = -1/a
    sigma_std = a_err/a**2
    decimals_sigma = int(np.ceil(-np.log10(sigma_std)))
    print(f"Sigma = {sigma:.{decimals_sigma}f} +- {sigma_std:.1g}")

    x_fit = np.linspace(np.min(x), np.max(x), 50)
    plt.figure(num=0, dpi=80, facecolor='w', edgecolor='k')
    plt.plot(10**x, 10**y, "o", alpha = 1, label = "Datapoints")
    plt.plot(10**x_fit, 10**(a*x_fit + b), label = f"Linear fit (log-plot)\nSlope = {a:.{decimals_a}f} " r"$\pm$" f" {a_err:.1g}")
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel(r"$|p - p_c|$", fontsize=14)
    plt.ylabel(r"$s_{\xi}$", fontsize=14)
    plt.legend(fontsize = 13)
    plt.tight_layout(pad=1.1, w_pad=0.7, h_pad=0.2)




    return sigma



if __name__ == "__main__":
    p_c = 0.59275

    L = 1000
    a = 1.6
    MC_cycles = 500


    # p_below = np.linspace(0.56, 0.58, 20)
    p_start = 0.56
    p_end = 0.58
    num = 20
    p_below = p_c -np.logspace(np.log10(p_c-p_start), np.log10(p_c-p_end), num)


    cluster = TwoDim_cluster(L, p_c)
    s_pc, n_pc = cluster.density(a, MC_cycles)

    # s_xi(p, L, a, s_pc, n_pc)

    sigma = s_xi_p(p_below, L, a, s_pc, n_pc, plot = False)
    print(sigma)

    #plt.savefig("../article/figures/h_sigma.pdf", bbox_inches="tight")
    #plt.show()



    #
