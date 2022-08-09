from cluster_generator import *
from plot_set import *


def plot(p, a, MC_cycles, L_try):

    plt.figure(num=0, dpi=80, facecolor='w', edgecolor='k')

    if isinstance(L_try, (float,int)): # single value
        cluster = TwoDim_cluster(L_try,p_c)
        s, n = cluster.density(a, MC_cycles)
        plt.plot(s,n,  "--", marker = "o", label = f"p = {p_try:g}")

    else: # array, list, tuple
        L_num = 0
        print(f"\r {L_num}/{len(L_try)}", end = "")
        for L in L_try:
            L_num += 1
            cluster = TwoDim_cluster(L,p_c)
            s, n = cluster.density(a, MC_cycles)
            print(f"\r {L_num}/{len(L_try)}", end = "")

            plot = plt.plot(s, n,  "--", marker = "o", label = f"L = {L:g}")
            if L_num == len(L_try): #Last L
                s_min = 1
                s_max = 4
                s_idx = np.where(np.logical_and(s_min < np.log10(s), np.log10(s) < s_max, np.log10(n) != -np.inf))
                x, y = np.log10(s[s_idx]), np.log10(n[s_idx])
                a, b, a_err, b_err = lin_fit(x,y)
                decimals_a = decimals(a_err)

                hex_color = plot[0].get_color().lstrip('#')
                rgb_color = np.array(tuple(int(hex_color[i:i+2], 16) for i in (0, 2, 4)))/255
                rgb_color *= 0.6 #give shade
                plt.plot(s[s_idx], s[s_idx]**(a)*10**(b), color = rgb_color, label = f"Linfit for L = {L_try[-1]}\nSlope = {a:.{decimals_a}f} " r"$\pm$" f" {a_err:.1g}")
        print()


    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel(r"$s$", fontsize=14)
    plt.ylabel(r"$n(s,p_c)$", fontsize=14)
    plt.legend(fontsize = 13)
    plt.tight_layout(pad=1.1, w_pad=0.7, h_pad=0.2)






if __name__ == "__main__":
    p_c = 0.59275
    a = 1.3
    MC_cycles = 1000


    L_try = 2**np.arange(4,9+1)

    plot(p_c, a, MC_cycles, L_try)
    plt.savefig("../article/figures/g.pdf", bbox_inches="tight")
    plt.show()
