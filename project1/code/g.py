import numpy as np
import matplotlib.pyplot as plt
from plot_set import *
from scipy import signal
from scipy.interpolate import interp1d
from scipy.signal import argrelextrema

def read_file(filename):
    lines = []
    info = []
    exp_len = 4
    with open(filename, "r") as infile:
        while True: #Read lines beginning with "#"
            line = infile.readline()
            if line[0] == "#":
                info.append(line)
            else:
                #Notice that we itentionally skip a line
                break
        for line in infile:
            words = line.split()
            if len(words) == exp_len:
                lines.append(np.array(words, dtype=float))
    return np.array(lines)

def smooth_curve(input, window_length = 1001, polyorder = 5):
    output = signal.savgol_filter(input, window_length, polyorder)
    return output

def color_cycle(num_color):
    """ get coloer from matplotlib
        color cycle """
    color = plt.rcParams['axes.prop_cycle'].by_key()['color']
    return color[num_color]

def plot_rdf(filenames):
    names = ["liquid", "solid"]
    plt.figure(num=0,  figsize = (8,8), dpi=80, facecolor='w', edgecolor='k')
    window_length = [1001, 301]
    polyorder = [5, 5]
    order = [500, 100]
    for num in range(2):
        data = read_file(filenames[num])
        x = data[-200:,1]
        y = data[-200:,2]
        f = interp1d(x, y, kind = "cubic")
        x_p1d = np.linspace(data[-200,1],data[-1,1], int(1e4))
        y_1pd = f(x_p1d)

        smooth_start = 1970
        x_smooth = x_p1d[smooth_start:]
        y_smooth = smooth_curve(y_1pd[smooth_start:], window_length[num], polyorder[num])

        local_max = argrelextrema(y_smooth, np.greater, order = order[num])
        x_max = x_p1d[smooth_start:][local_max]
        y_max = y_1pd[smooth_start:][local_max]

        plt.subplot(2,1,num+1)
        plt.plot(data[-200:,1], data[-200:,2], "o", markersize = 4, label = "Datapoints")
        plt.plot(x_p1d, y_1pd, label = "Interpolation curve")
        #plt.plot(x_smooth, y_smooth, "--", label = "Smooth curve")

        x_offset = 0.05
        y_offset = 0.05
        for i in range(len(x_max)):
            plt.text(x_max[i] + x_offset, y_max[i] + y_offset, f"{x_max[i]:.2f}")

        avg_spacing = np.mean(x_max[1:]-x_max[:-1])
        SE = np.std(x_max[1:]-x_max[:-1])
        decimal_SE = decimals(SE)

        plt.plot(x_max, y_max, "o", color =  color_cycle(3), label = "Local maximum (savgol filter)\n" + r"$\langle \Delta r_{max} \rangle$ = " + f"{avg_spacing:.{decimal_SE}f} " + r"$\pm$" + f" {SE:.{decimal_SE}f}"  )
        print(f"Average distance {avg_spacing:.4f}" )


        plt.title(names[num])
        plt.xlabel(r"$r / \sigma$", fontsize = 14)
        plt.ylabel(r"$g(r)$ (normalized)", fontsize = 14)
        plt.tight_layout(pad=1.1, w_pad=0.7, h_pad=0.2)
        plt.legend(fontsize = 13)
    plt.subplots_adjust(hspace = 0.30)
    plt.savefig("../article/figures/rdf.pdf", bbox_inches="tight")





if __name__ == "__main__":
    filename_liquid = "rdf_liquid.txt"
    filename_solid = "rdf_solid.txt"
    filenames = [filename_liquid, filename_solid]
    plot_rdf(filenames)
    plt.show()
