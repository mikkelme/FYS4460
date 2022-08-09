import numpy as np
import matplotlib.pyplot as plt
import plot_set
import statsmodels.api as sm



def read_radial_bins(filename):
    with open(filename, "r") as infile:
        for line in infile:
            if line[0] != "#":
                step, num_chunks, N = np.array(line.split(), dtype = int)
                break
        data = []
        while True:
            data_step_i = np.zeros((num_chunks, 2))
            for i in range(num_chunks):
                #Each line: Chunk Coord1 Coord2 Ncount vx
                data_step_i[i] = np.array(infile.readline().split(), dtype=float)[[1,4]]
            data.append(data_step_i)
            info = infile.readline().split()
            if len(info) < 3:
                break
        return np.array(data)

def lin_fit(x,y):
    x = sm.add_constant(x)
    model = sm.OLS(y, x)
    res = model.fit()
    b, a = res.params
    b_err, a_err = res.bse
    return a, b, a_err, b_err


def plot_radial_distribution(data, avg_last_steps):
    R = data[-1,-1,0] + data[-1,0,0]
    x, y = (R**2 - data[-1:,:,0]**2)[0], np.mean(data[-avg_last_steps:,:,1], axis = 0)

    a, b, a_err, b_err = lin_fit(x,y) # y = a(R^2 - r^2)
    r_fit = np.linspace(0, R, int(1e4))

    plt.figure(num=0, dpi=80, facecolor='w', edgecolor='k')
    plt.plot(x,y, "o", label = "datapoints")
    plt.plot(x, x*a + b, label = "linear fitr: $v(r) = a(R^2 - r^2)$" + f"\na = {a:.4f}" + " $\pm$ " + f"{a_err:1.0e} " + r"$\sigma \tau$")
    plt.tight_layout(pad=1.1, w_pad=0.7, h_pad=0.2)
    plt.legend(fontsize = 13)
    plt.xlabel(r"$(a^2 - r^2)/\sigma^2$", fontsize = 14)
    plt.ylabel(r"$v_x(r)/\sigma \tau^{-1}$", fontsize = 14)
    plt.savefig("../article/figures/flow_profile_linfit.pdf", bbox_inches="tight")

    plt.figure(num=1, dpi=80, facecolor='w', edgecolor='k')
    plt.plot(np.mean(data[-avg_last_steps:,:,0], axis = 0), np.mean(data[-avg_last_steps:,:,1], axis = 0), label = f"datapoints (average over last {100*avg_last_steps:.0e} timesteps)")
    plt.plot(r_fit, a*(R**2 - r_fit**2), label = r"Profile from fit: $v(r) = a(R^2 - r^2)$")
    plt.xlabel(r"$r/\sigma$", fontsize = 14)
    plt.ylabel(r"$v(r)/\sigma \tau^{-1}$", fontsize = 14)
    plt.legend(fontsize = 13)
    plt.tight_layout(pad=1.1, w_pad=0.7, h_pad=0.2)
    plt.savefig("../article/figures/flow_profile.pdf", bbox_inches="tight")

    return a, b, a_err, b_err



def cal_viscosity(a, a_err):
    n =  1/2*4/(5.72/3.4)**3
    F_x = 0.1
    mu = n*F_x/(4*a)
    Delta_mu = n*F_x/(4*a**2)*a_err
    decimals = int(np.ceil(-np.log10(Delta_mu)))
    print(f"mu = {mu:.{decimals}f} +- {Delta_mu:.1g}")


def plot_inner_product(data, eff_dt):
    vv = np.sum(data[-1]**2)
    t = np.linspace(0, len(data)*eff_dt, len(data))
    inner_product = np.zeros(len(data))
    for i in range(len(data)):
        inner_product[i] = np.sum(data[i,:]*data[-1,:])/vv

    plt.figure(num=3, dpi=80, facecolor='w', edgecolor='k')
    plt.plot(t, inner_product)
    plt.xlabel(r"$t/\tau$", fontsize = 14)
    plt.ylabel("Inner product", fontsize = 14)
    plt.tight_layout(pad=1.1, w_pad=0.7, h_pad=0.2)
    plt.savefig("../article/figures/flow_inner_product.pdf", bbox_inches="tight")







if __name__ == "__main__":
    filename = "radial_vel.txt"
    eff_dt = 100*0.005
    data = read_radial_bins(filename)[:] #Array[steps, chunks, (coord, vx)]
    #plot_inner_product(data, eff_dt)
    a, b, a_err, b_err = plot_radial_distribution(data, avg_last_steps = 1000)
    cal_viscosity(a, a_err)
    #plt.show()
