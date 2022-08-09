import numpy as np
import matplotlib.pyplot as plt
from plot_set import *
import statsmodels.api as sm
import subprocess
import os




def read_fix_ave_time(filename):
    lines = []
    with open(filename, "r") as infile:
        info = infile.readline() + infile.readline()
        for line in infile:
            lines.append(np.array(line.split(), dtype=float))
    return np.array(lines)

def get_D(data_folder, dt):
    T = []; D = []
    start_cut_off = 0.5 # [pct]

    for filename in os.listdir(data_folder):
         if filename.endswith(".txt"):
            #data
            data = read_fix_ave_time(data_folder + "/" + filename)
            start_idx = int(start_cut_off*len(data))
            t = (data[start_idx:,0] - data[0,0])*dt
            g = data[start_idx:,2]
            mean_T = np.mean(data[start_idx:,1])

            #linear fit
            x,y = t, g
            print(filename)
            plt.plot(t,g)
            plt.show()
            x = sm.add_constant(x)
            model = sm.OLS(y, x)
            res = model.fit()
            b, a = res.params
            b_err, a_err = res.bse

            #output
            T.append(mean_T)
            D.append(a/6)

    return np.array(T), np.array(D)

def run_simulations(temp):
    script_file = "l.in"
    setting_file = "l_run_settings.in"
    T_str = str(temp).replace(".","")
    with open(setting_file, "w") as outfile:
        outfile.write(f"variable myTemp equal {temp}")
    subprocess.run(["lmp_serial < " +  script_file], shell = True)




if __name__ == "__main__":
    RUN = False
    READ = True


    run_temp = np.linspace(1500, 2000, 5)
    # print(run_temp)
    # exit()

    data_folder = "l_data"
    dt = 1.0e-3

    if RUN:
        for temp in run_temp:
            run_simulations(temp)

    if READ:
        T, D = get_D(data_folder, dt)
        sort_idx = np.argsort(T)
        T = T[sort_idx]
        D = D[sort_idx]


        #find T_m
        min_T = 1400
        # min_T = 2200
        start_idx = np.argwhere(T >= min_T)[0][0]


        x,y = T[start_idx:], D[start_idx:]
        a, b, a_err, b_err = lin_fit(x,y)
        T_m = - b/a


        #plot
        plt.figure(num=0, dpi=80, facecolor='w', edgecolor='k')
        plt.plot(T,D, "o", label = "Datapoints")
        # plt.plot([T[start_idx], T[-1]], [a*T[start_idx] + b, a*T[-1] + b], "--", label = r"Linear fit for $T > $" + f"{min_T} K")
        # plt.plot(T_m, 0, "x", label = "Estimated melting temperature\n" + r"$T_m = $" + f"{T_m:.2f} at linear fit x-axis intercept")
        plt.xlabel("$T$ $[K]$", fontsize = 14)
        plt.ylabel("$D$ $[Å^2/ps]$")
        plt.tight_layout(pad=1.1, w_pad=0.7, h_pad=0.2)
        plt.legend(fontsize = 13)
        # plt.savefig("../article/figures/D(T)_alternative.pdf", bbox_inches="tight")
        plt.show()
