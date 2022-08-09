import numpy as np
import matplotlib.pyplot as plt
import plot_set
import statsmodels.api as sm
import subprocess
import os

from l import read_fix_ave_time, get_D
from g import read_file, color_cycle, smooth_curve
from scipy import signal
from scipy.interpolate import interp1d
from scipy.signal import argrelextrema


def run_simulations(temp, script_file):
    # script_file = "spce-water-system.in"
    setting_file = "m_run_settings.in"
    T_str = str(temp).replace(".","")
    with open(setting_file, "w") as outfile:
        outfile.write(f"variable myTemp equal {temp}")
    subprocess.run(["lmp_serial < " +  script_file], shell = True)

def get_rdf(data_folder):
    data_list = []
    T = []
    for filename in os.listdir(data_folder):
         if filename.endswith(".txt"):
             data_list.append(read_file(data_folder + "/" + filename))
             T.append(float(filename.strip("m_rdf_.txt")))
    return T, np.array(data_list)




if __name__ == "__main__":
        RUN = False
        READ = True


        # run_temp = np.linspace(50, 250, 5)
        run_temp = np.array([300.0, 500.0])
        # print(run_temp)
        # exit()

        dt = 1.0
        data_folder_diffusion = "m_data_diffusion"
        data_folder_radial = "m_data_radial"

        # script_file = "spce-water-system.in"
        script_file = "m_radial.in"
        if RUN:
            for temp in run_temp:
                run_simulations(temp, script_file)

        if READ:
            T, D = get_D(data_folder_diffusion, dt)
            sort_idx = np.argsort(T)
            T = T[sort_idx]
            D = D[sort_idx]

            #plot
            plt.figure(num=0, dpi=80, facecolor='w', edgecolor='k')
            plt.plot(T,D, "o", label = "Datapoints")
            #plt.plot([T[start_idx], T[-1]], [a*T[start_idx] + b, a*T[-1] + b], "--", label = r"Linear fit for $T >= $" + f"{min_T}")
            # plt.plot(T_m, 0, "o", label = "Estimated melting temperature\n" + r"$T_m = $" + f"{T_m:.2f}")
            plt.xlabel("$T$ $[K]$", fontsize = 14)
            plt.ylabel("$D$ $[Å^2/ps]$")
            plt.tight_layout(pad=1.1, w_pad=0.7, h_pad=0.2)
            plt.legend(fontsize = 13)
            plt.savefig("../article/figures/water_D(T).pdf", bbox_inches="tight")
            #plt.show()


            window_length = [1001, 1001, 1001]
            polyorder = [5, 5, 5]
            order = [500, 500, 500]

            T, data = get_rdf(data_folder_radial)
            data_length = len(data)
            plt.figure(num=1,  figsize = (8,8), dpi=80, facecolor='w', edgecolor='k')

            for i in range(data_length):
                #Data area
                x = data[i,-200:,1]
                y = data[i,-200:,2]

                #Interpolation
                f = interp1d(x, y, kind = "cubic")
                x_p1d = np.linspace(data[i,-200,1],data[i,-1,1], int(1e4))
                y_p1d = f(x_p1d)

                #Savgol filter (smooth curve)
                smooth_start = 2500
                smooth_end = 4000
                x_smooth = x_p1d[smooth_start:smooth_end]
                y_smooth = smooth_curve(y_p1d[smooth_start:smooth_end], window_length[i], polyorder[i])

                #Find local maxima
                local_max = argrelextrema(y_smooth, np.greater, order = order[i])
                x_max = x_smooth[local_max]
                y_max = y_smooth[local_max]


                #Plot curves
                plt.subplot(data_length,1,i+1)
                plt.plot(data[i,-200:,1], data[i,-200:,2], "o", markersize = 4, color =  color_cycle(i), label = f"Datapoints, T = {T[i]} K")
                #plt.plot(x_p1d, y_p1d, color =  color_cycle(i), label = "Interpolation curve")
                plt.plot(x_smooth, y_smooth, "--", color =  color_cycle(i),  label = "Savgol filer")

                #Plot text
                x_offset = 0.2
                y_offset = 0.05
                for j in range(len(x_max)):
                    plt.text(x_max[j] + x_offset, y_max[j] + y_offset, f"{x_max[j]:.2f}")

                plt.plot(x_max, y_max, "o", color =  color_cycle(3), label = "Local $r_{max}$ (savgol filter)" )


                #Labels
                if i == 1:
                    plt.ylabel(r"$g(r)$ (normalized)", fontsize = 14)
                plt.legend(fontsize = 13)
            plt.tight_layout(pad=1.1, w_pad=0.7, h_pad=0.2)
            plt.xlabel(r"$r [Å]$", fontsize = 14)
            plt.savefig("../article/figures/water_rdf.pdf", bbox_inches="tight")

            plt.show()
