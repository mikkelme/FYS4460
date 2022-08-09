import numpy as np
import matplotlib.pyplot as plt
from read_dump import read_dump_ovito
import subprocess
from plot_set import *
from scipy import signal




def read_energy(energy_file):
    lines = []
    with open(energy_file, "r") as infile:
        info = infile.readline() + infile.readline()
        for line in infile:
            lines.append(np.array(line.split(), dtype=float))
    return np.array(lines)

def color_cycle(num_color):
    """ get color from matplotlib
        color cycle """
    color = plt.rcParams['axes.prop_cycle'].by_key()['color']
    return color[num_color]


def run_simulations(runtime, dt):
    script_file = "b.in"
    setting_file = "b_run_settings.in"
    dt_str = str(dt).split(".")[-1]
    with open(setting_file, "w") as outfile:
        outfile.write(f"variable run equal {int(runtime/dt)}\
        \nvariable timestep equal {dt}\
        \nvariable outfile string etotal{dt_str}.txt")
    subprocess.run(["lmp_serial < " +  script_file], shell = True)

def smooth_curve(input, window_length = 801, polyorder = 5):
    output = signal.savgol_filter(input, window_length, polyorder)
    return output

if __name__ == "__main__":
    Run = False
    Read = False
    Read2 = False
    Read3 = True
    runtime = 20
    # dt_run = [0.01, 0.008, 0.002]
    dt_run = np.logspace(np.log10(0.01), np.log10(0.02), 10)
    dt_read = dt_run


    if Run:
        for dt in dt_run:
            run_simulations(runtime, dt)
    if Read:
        energy_files = []
        for dt in dt_read:
            dt_str = str(dt).split(".")[-1]
            energy_files.append("etotal" + dt_str + ".txt")

        # Plot total energy over time
        fig_Etot = plt.figure(num = 1, dpi=80, facecolor='w', edgecolor='k')
        for i in range(len(energy_files)):
            energy_file = energy_files[i]
            data = read_energy(energy_file)
            dt = float("0." + energy_file.strip("etotal.txt"))
            plt.plot(data[:,0]*dt, data[:,3], color = color_cycle(i+1), label = f"dt = {dt}")

        plt.xlabel(r"$t/\tau$", fontsize = 14)
        plt.ylabel("$E_{tot}/\epsilon$", fontsize = 14)
        plt.tight_layout(pad=1.1, w_pad=0.7, h_pad=0.2)
        plt.legend(fontsize = 13)
        fig_Etot.savefig("../article/figures/Etot2.pdf", bbox_inches="tight")
        plt.show()

    if Read2:
        energy_files = []
        for dt in dt_read[0:]:
            dt_str = str(dt).split(".")[-1]
            energy_files.append("etotal" + dt_str + ".txt")

        window_length = [1001, 1201, 1901]
        polyorder = [5, 5, 3]
        # Plot total energy over time
        fig_SE = plt.figure(num = 1, dpi=80, facecolor='w', edgecolor='k')
        for i in range(len(energy_files)):
            data = read_energy(energy_files[i])
            dt = float("0." + energy_files[i].strip("etotal.txt"))
            line = plt.plot(data[:,0]*dt, data[:,3], alpha = 0.7, color = color_cycle(i+1), label = f"dt = {dt}")


            hex_color = line[0].get_color().lstrip('#')
            rgb_color = np.array(tuple(int(hex_color[i:i+2], 16) for i in (0, 2, 4)))/255
            rgb_color *= 0.8 #give shade

            fit_start = 0
            start_idx = np.argmin(np.abs(data[:,0]*dt - fit_start))

            x = data[start_idx:,0]*dt
            y = data[start_idx:,3]
            a, b, a_err, b_err = lin_fit(x,y)
            line = plt.plot(x, a*x+b, "--", color = rgb_color, label = f"Linear fit.: SE = {a_err:.2e}")


            # y_smooth = smooth_curve(data[:,3], window_length[i], polyorder[i])
            # plt.plot(data[:,0]*dt, y_smooth, color = rgb_color, label =  f"smooth (dt ={dt})")


        plt.xlabel(r"$t/\tau$", fontsize = 14)
        plt.ylabel("$E_{tot}/\epsilon$", fontsize = 14)
        plt.tight_layout(pad=1.1, w_pad=0.7, h_pad=0.2)
        plt.legend(fontsize = 13)

        fig_SE.savefig("../article/figures/Etot_SE2.pdf", bbox_inches="tight")
        plt.show()

    if Read3:
        energy_files = []
        for dt in dt_read[0:]:
            dt_str = str(dt).split(".")[-1]
            energy_files.append("etotal" + dt_str + ".txt")

        dt_array = np.zeros(len(energy_files))
        SE_array = np.zeros(len(energy_files))
        # Plot total energy over time
        fig_SE_trend = plt.figure(num = 1, dpi=80, facecolor='w', edgecolor='k')
        for i in range(len(energy_files)):
            data = read_energy(energy_files[i])
            dt = float("0." + energy_files[i].strip("etotal.txt"))
            fit_start = 0
            start_idx = np.argmin(np.abs(data[:,0]*dt - fit_start))

            x = data[start_idx:,0]*dt
            y = data[start_idx:,3]
            a, b, a_err, b_err = lin_fit(x,y)

            dt_array[i] = dt
            SE_array[i] = a_err

        x = np.log10(dt_array)
        y = np.log10(SE_array)
        a, b, a_err, b_err = lin_fit(x,y)
        decimals_a = decimals(a_err)

        plt.plot(dt_array, SE_array, "o", label = "Datapoints")
        plt.plot(10**x, 10**(a*x + b), label = f"Linear fit (log-plot)\nSlope = {a:.{decimals_a}f} " r"$\pm$" f" {a_err:.1g}")

        plt.xlabel(r"$dt/\tau$", fontsize = 14)
        plt.ylabel("$SE(dt)$", fontsize = 14)
        plt.xscale("log")
        plt.yscale("log")

        plt.tight_layout(pad=1.1, w_pad=0.7, h_pad=0.2)
        plt.legend(fontsize = 13)

        fig_SE_trend.savefig("../article/figures/Etot_SE_trend.pdf", bbox_inches="tight")
        plt.show()
