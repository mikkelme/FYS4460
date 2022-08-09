import numpy as np
import matplotlib.pyplot as plt
import subprocess
import plot_set
import statsmodels.api as sm

def run_simulations(temp):
    script_file = "d.in"
    setting_file = "d_run_settings.in"
    T_str = str(temp).replace(".","")


    with open(setting_file, "w") as outfile:
        outfile.write(f"variable temp equal {temp} \
        \nvariable outfile string d_etot_T{T_str}.txt")
    subprocess.run(["lmp_serial < " +  script_file], shell = True)



if __name__ == "__main__":
    Run = False
    Read = True

    temp_run = [2.5, 3.0, 3.5, 4.0, 4.5, 5.0]
    dt = 0.002
    size = 10
    from b import read_energy

    if Run:
        for temp in temp_run:
            run_simulations(temp)



    if Read:
        energy_files = ["d_etot_T" + str(temp).replace(".","") + ".txt" for temp in temp_run]

        avg_temp = []
        avg_pressure = []
        for energy_file in energy_files:
            data = read_energy(energy_file)
            start_idx = len(data)//2
            N = 4*size**3
            Temp = 2/3*data[start_idx:,1]/(N)
            avg_temp.append(np.mean(Temp))
            avg_pressure.append(np.mean(data[start_idx:,4]))

        x,y = avg_temp, avg_pressure
        x = sm.add_constant(x)
        model = sm.OLS(y, x)
        res = model.fit()
        b, a = res.params
        b_err, a_err = res.bse


        temp_space = np.linspace(avg_temp[0], avg_temp[-1], int(1e3))
        plt.plot(avg_temp, avg_pressure, "o", label = "datapoints (simulation average)")
        plt.plot(temp_space, temp_space*a + b, linestyle ="--", label = f"Linear fit:\na = {a:.6f}" + r" $\pm$ " + f"{a_err:.1g}\nb = {b:.5f}" + r" $\pm$ " + f"{b_err:.1g}")
        plt.xlabel(r"$T/\epsilon k_B^{-1}$ ", fontsize = 14)
        plt.ylabel(r"$P/\epsilon \sigma^{-3}$", fontsize = 14)
        plt.tight_layout(pad=1.1, w_pad=0.7, h_pad=0.2)
        plt.legend(fontsize = 13)
        plt.savefig("../article/figures/P(T).pdf", bbox_inches="tight")

        plt.show()
