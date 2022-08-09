from l import *



if __name__ == "__main__":
    L_try = [25, 50, 100, 200, 400, 800]
    # L_try = [25, 50, 100, 200]

    # p_c = 0.59275
    nu = 4/3
    p_c = 0.594
    MC_cycles = 200

    span = 0.05
    p_span = np.linspace(p_c - span, p_c + span, 20)
    PI_MC_cycles = 200

    plt.figure(num=0, dpi=80, facecolor='w', edgecolor='k')
    for i in range(len(L_try)):
        print(f"\r {i}/{len(L_try)}, L = {L_try[i]}", end = "")
        pi = PI(L_try[i], p_span, MC_cycles)
        plt.plot(L_try[i]**(1/nu)*(p_span-p_c), pi, "o", label = f"L = {L_try[i]}")
    print(f"\r {i}/{len(L_try)}, L = {L_try[i]}")
    plt.xlabel(r"$L^{1/\nu}(p-p_c)$", fontsize=14)
    plt.ylabel(r"$\Pi(p,L)$", fontsize=14)
    plt.legend(fontsize = 13)
    plt.tight_layout(pad=1.1, w_pad=0.7, h_pad=0.2)
    plt.savefig("../article/figures/n2.pdf", bbox_inches="tight")
    plt.show()
