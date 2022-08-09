import numpy as np
import matplotlib.pyplot as plt
import plot_set
import statsmodels.api as sm



def read_fix_ave_time(filename):
    lines = []
    info = []
    with open(filename, "r") as infile:
        for line in infile:
            if line[0] == "#":
                info.append(line)
            else:
                lines.append(np.array(line.split(), dtype=float))
                break
        for line in infile:
            lines.append(np.array(line.split(), dtype=float))
    return np.array(lines)

def lin_fit(x,y):
    x = sm.add_constant(x)
    model = sm.OLS(y, x)
    res = model.fit()
    b, a = res.params
    b_err, a_err = res.bse
    return a, b, a_err, b_err


def plot_msd(timesteps, msd, dt):
    t = (timesteps - timesteps[0])*dt
    fit_start = int(0.05*len(t))
    a, b, a_err, b_err = lin_fit(t[fit_start:],msd[fit_start:])
    D = a/6
    print(f"Estimated diffusion constant: D = {D:.5f}")

    plt.figure(num=0, dpi=80, facecolor='w', edgecolor='k')
    plt.plot(t, msd, label = f"Data")
    plt.plot(t[fit_start:], t[fit_start:]*a + b, linestyle ="--", alpha = 0.8, label = f"Linear fit: a = {a:.5f}" + r" $\pm$ " + f"{a_err:1.0e}\nb = {b:.2f}" + r" $\pm$ " + f"{b_err:1.0e}")
    plt.xlabel(r"$t/\tau$", fontsize = 14)
    plt.ylabel(r"$\langle r^2 (t)\rangle/\sigma^2$", fontsize = 14)
    plt.tight_layout(pad=1.1, w_pad=0.7, h_pad=0.2)
    plt.legend(fontsize = 13)
    plt.savefig("../article/figures/msd.pdf", bbox_inches="tight")
    plt.show()





if __name__ == "__main__":
    filename = "msd.txt"
    dt = 0.005
    data = read_fix_ave_time(filename)
    plot_msd(data[:,0], data[:,1], dt)
