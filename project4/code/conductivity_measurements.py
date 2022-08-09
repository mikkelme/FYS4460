from pylab import *
from scipy.ndimage import measurements
from perculation_flow import *
from plot_set import *

def find_perculation_cluster(L, p):
    ncount = 0
    perc = []
    while (ncount < 1000):
        ncount += 1
        z=rand(L, L)<p
        lw,num = measurements.label(z)
        perc_x = intersect1d(lw[0,:],lw[-1,:])
        perc = perc_x[where(perc_x > 0)]
        if len(perc) > 0: # Found spanning cluster
            return z, lw, perc
    print("Couldn’t make percolation cluster...")
    exit()

def conductivity_density(L, pVals, nSamples, plot = True):

    G = zeros(len(pVals)) # Condutivity
    P = zeros(len(pVals)) # Density

    for p_idx in range(len(pVals)):
        p = pVals[p_idx]
        for j in range(nSamples):
            print(f"\r p: {p_idx+1}/{len(pVals)}, MC Cycles: {j+1}/{nSamples} ", end= "")
            z, lw, perc = find_perculation_cluster(L, p)
            area = measurements.sum(z, lw, perc[0])
            P[p_idx] += area # Find P(p,L)
            zz = asarray((lw == perc[0])) # zz=spanning cluster
            zzz = zz.T
            g = sitetobond (zzz) # Generate bond lattice
            Pvec, c_eff = FIND_COND(g, L, L) # Find conducance
            G[p_idx] += c_eff
    print("\nDone")
    G /= nSamples
    P /= (nSamples*L*L)

    if plot:
        plt.figure(num=0, dpi=80, facecolor='w', edgecolor='k')
        plot(pVals,G,'-o', label='$G$')
        plot(pVals,P,'-o', label='$P$')

        plt.xlabel(r"$p$", fontsize=14)
        plt.ylabel(r"$G,P$", fontsize=14)
        plt.legend(fontsize = 13)
        plt.tight_layout(pad=1.1, w_pad=0.7, h_pad=0.2)
    return G, P


if __name__ == "__main__":
    L = 400
    pVals = logspace(log10(0.58), log10(0.85), 20)
    nSamples = 100
    print(pVals)
    print(np.linspace(0.58, 0.85,20))
    exit()
    conductivity_density(L, pVals, nSamples)

    plt.savefig("../article/figures/GP(p).pdf", bbox_inches="tight")
    plt.show()
