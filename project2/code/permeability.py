import numpy as np
import matplotlib.pyplot as plt
from plot_set import *
from msd import read_fix_ave_time
import os


def find_perm(data_folder):
    avg_last_steps = 10
    phi = np.zeros(len(filenames))
    U = np.zeros(len(filenames))
    for i in range(len(filenames)):
        data = read_fix_ave_time(filenames[i]) #step, U
        phi[i] = float(filenames[i].split("_")[-1].strip(".txt"))
        U[i] = np.mean(data[-avg_last_steps:, 1])

        # print(phi[i])
        # plt.plot(data[:,0], data[:,1])
        # plt.show()

    const = 1
    mu = 1.19
    k = U*mu/const

    x, y = np.log10(phi), np.log10(k)
    a, b, a_err, b_err = lin_fit(x,y)
    decimals_a = int(np.ceil(-np.log10(a_err)))

    
    x_fit = np.linspace(np.min(x), np.max(x), int(1e3))
    plt.figure(num=0, dpi=80, facecolor='w', edgecolor='k')
    plt.plot(phi, k, "o", label = "Datapoints")
    plt.plot(10**x_fit, 10**(a*x_fit + b), label = f"Linear fit (log-plot)\nSlope = {a:.{decimals_a}f} " r"$\pm$" f" {a_err:.1g}")
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel(r"$\phi$", fontsize=14)
    plt.ylabel(r"k", fontsize=14)
    plt.legend(fontsize = 13)
    plt.tight_layout(pad=1.1, w_pad=0.7, h_pad=0.2)

    plt.show()

def directories_from_folder(data_folder):
    filenames = []
    for filename in os.listdir(data_folder):
         if filename.endswith(".txt"):
             filenames.append(data_folder + "/" + filename)
    return filenames



if __name__ == "__main__":
    data_folder = "U_data"
    filenames = directories_from_folder(data_folder)
    find_perm(filenames)



    # plt.plot(data[:,0], data[:,1])
    # plt.show()






#
