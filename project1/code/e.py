import numpy as np
import matplotlib.pyplot as plt
import subprocess
import plot_set
import statsmodels.api as sm
import os

def run_simulations(temp, density ):
    script_file = "e.in"
    setting_file = "e_run_settings.in"
    T_str = str(temp).replace(".","")
    d_str = str(density).replace(".","")
    with open(setting_file, "w") as outfile:
        outfile.write(f"variable temp equal {temp}\
        \nvariable density equal {density}\
        \nvariable outfile string e_datafile_T{T_str}_d{d_str}.txt")
    subprocess.run(["lmp_serial < " +  script_file], shell = True)


def van_der_waals_parameters(N, T, P, V):
    a_min = 0; a_max = 10
    b_min = 0; b_max = 10
    a = np.nan; b = np.nan
    x = T
    slope = N
    MSE = 1e12

    BigRounds = 10; trynum = 1e2;
    a_try = np.linspace(a_min, a_max, int(trynum))
    b_try = np.linspace(b_min, b_max, int(trynum))
    for i in range(BigRounds):
        for ai in a_try:
            for bi in b_try:
                y = (P + ai*(N/V)**2)*(V - N*bi)
                intercept = np.mean(y - slope*x)
                y_fit = slope*x + intercept
                MSE_try = np.sum((y - y_fit)**2)/len(y)
                if MSE_try < MSE:
                    MSE = MSE_try
                    a = ai
                    b = bi
        if i == 0 and a == a_min or a == a_max or b == b_min or b == b_max:
            print("Searh interval insufficient")
            print(f"a = {a:.3f}, b = {b:.3f}")
            print(f"a_min = {a_min}, a_max = {a_max}")
            print(f"b_min = {b_min}, b_max = {b_max}")
            exit()
        eps = 2*np.max([a_try[1]-a_try[0], b_try[1]-b_try[0] ])
        a_try = np.linspace(a - eps, a + eps, int(trynum))
        b_try = np.linspace(b - eps, b + eps, int(trynum))


    # print(f"van der Walls parameters (best fit):\n a = {a:.{BigRounds}f}, b = {b:.{BigRounds}f}")
    return a, b





if __name__ == "__main__":
    Run = False
    Read = True

    # temp_run = [2.5, 3.0, 3.5, 4.0, 4.5, 5.0]
    # density_run = [0.01, 0.02, 0.03, 0.04, 0.05]
    temp_run = np.linspace(1, 10, 10)
    density_run = np.linspace(0.02, 0.38, 10)



    data_folder = "e_data"
    dt = 0.005
    size = 10
    N = 4*size**3



    from b import read_energy
    if Run:
        for temp in temp_run:
            for density in density_run:
                #print(density)
                run_simulations(temp, density)

    if Read:
        avg_temp = []
        density = []
        avg_pressure = []
        for filename in os.listdir(data_folder):
             if filename.endswith(".txt"):
                data = read_energy(data_folder + "/" + filename)
                avg_temp.append(2/3*data[-1,1]/N)
                # print(filename.split("T")[-1].split("_")[0], avg_temp[-1])
                avg_pressure.append(data[-1,4])
                density.append(float( "0." + filename.split("d")[-1].strip(".txt")[1:]))
                # print(filename, density[-1])
                #print(filename, )
                # print(f"{avg_temp[-1]:.2f}, {density[-1]:.2f}, {avg_pressure[-1]*N/density[-1]:.2f}")
        # exit()
        T = np.array(avg_temp)
        rho = np.array(density)
        P = np.array(avg_pressure)
        V = N/rho

        figIdeal = plt.figure(num=1, dpi=80, facecolor='w', edgecolor='k')
        x,y = T, P*V
        slope = N
        intercept = np.mean(y - slope*x)
        y_fit = x*slope + intercept
        RMSE = np.sqrt(np.sum((y - y_fit)**2)/len(y))
        plt.plot(x, y, "o", label = f"Ideal gas datapoints")
        plt.plot(x, y_fit, linestyle ="-", label = f"Best linear fit with slope = N\nRMSE = {RMSE:.2g}")
        plt.xlabel(r"$T/\epsilon k_B^{-1}$", fontsize = 14)
        plt.ylabel(r"$PV/\epsilon$", fontsize = 14)
        plt.tight_layout(pad=1.1, w_pad=0.7, h_pad=0.2)
        plt.legend(fontsize = 13)
        plt.savefig("../article/figures/TRP_ideal.pdf", bbox_inches="tight")



        figWaal = plt.figure(num=2, dpi=80, facecolor='w', edgecolor='k')
        a, b = van_der_waals_parameters(N, T, P, V)
        slope = N
        x,y = T, (P + a*(N/V)**2)*(V - N*b)
        intercept = np.mean(y - slope*x)
        y_fit = x*slope + intercept
        RMSE = np.sqrt(np.sum((y - y_fit)**2)/len(y))

        plt.plot(x, y, "o", label = f"Van der Waals datapoints, a = {a:.3f}, b = {b:.3f}")
        plt.plot(x, y_fit, linestyle = "-", label = f"Best linear fit with slope = N\nRMSE = {RMSE:.2g}")
        plt.xlabel(r"$T/\epsilon k_B^{-1}$", fontsize = 14)
        plt.ylabel(r"$(P + a(N/V)^2)(V - Nb)/\epsilon$", fontsize = 14)
        plt.tight_layout(pad=1.1, w_pad=0.7, h_pad=0.2)
        plt.legend(fontsize = 13)
        plt.savefig("../article/figures/TRP_waal.pdf", bbox_inches="tight")



        plt.show()
