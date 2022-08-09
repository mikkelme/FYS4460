import numpy as np
import matplotlib.pyplot as plt
import subprocess
import plot_set
import statsmodels.api as sm

def run_simulations(temp):
    script_file = "f.in"
    setting_file = "f_run_settings.in"
    T_str = str(temp).replace(".","")
    with open(setting_file, "w") as outfile:
        outfile.write(f"variable temp equal {temp}\
        \nvariable outfile string f_datafile_T{T_str}.txt")
    subprocess.run(["lmp_serial < " +  script_file], shell = True)


#run_simulations(2.5)
if __name__ == "__main__":
    Run = False
    Read = True

    temp_run = [0.9, 1.0, 1.1, 1.2]
    dt = 0.005

    from b import read_energy

    if Run:
        for temp in temp_run:
            run_simulations(temp)


    if Read:
        data_files = []
        for temp in temp_run:
            T_str = str(temp).replace(".","")
            data_files.append(f"f_datafile_T{T_str}.txt")


        avg_r2 = []
        for data_file in data_files:
            data = read_energy(data_file)
            avg_r2.append([data[:,0],data[:,1]])
        avg_r2 = np.array(avg_r2)

        fig1 = plt.figure(num=1, dpi=80, facecolor='w', edgecolor='k')
        fit_start = 20# [t]
        slopes = []
        for i in range(len(avg_r2)):
            t = (avg_r2[i,0,:] - avg_r2[i,0,0])*dt
            plot = plt.plot(t, avg_r2[i,1,:], label = f"T = {temp_run[i]}")

            x,y = t[int(fit_start/dt):], avg_r2[i,1,int(fit_start/dt):]
            x = sm.add_constant(x)
            model = sm.OLS(y, x)
            res = model.fit()
            b, a = res.params
            b_err, a_err = res.bse
            t_space = np.linspace(t[int(fit_start/dt)], t[-1], int(1e4))
            hex_color = plot[0].get_color().lstrip('#')
            rgb_color = np.array(tuple(int(hex_color[i:i+2], 16) for i in (0, 2, 4)))/255
            rgb_color *= 0.85 #give shade
            err_precision = int(np.ceil(-np.log10(a_err)))
            plt.plot(t_space, t_space*a + b, color = rgb_color, linestyle ="--", label = f"Linear fit: a = {a:.{err_precision}f}" + r" $\pm$ " + f"{a_err:.1g}")
            slopes.append([temp_run[i], a/6, a_err/6])
        plt.xlabel(r"$t/\tau$", fontsize = 14)
        plt.ylabel(r"$\langle r^2(t) \rangle / \sigma^2$", fontsize = 14)
        plt.tight_layout(pad=1.1, w_pad=0.7, h_pad=0.2)
        plt.legend(fontsize = 13)
        plt.savefig("../article/figures/diffusion.pdf", bbox_inches="tight")


        fig2 = plt.figure(num=2, dpi=80, facecolor='w', edgecolor='k')
        slopes = np.array(slopes)
        plt.plot(slopes[:,0], slopes[:,1], "o", label = r"Estimated diffusion constant $D$")


        x,y = slopes[:,0], slopes[:,1]
        x = sm.add_constant(x)
        model = sm.OLS(y, x)
        res = model.fit()
        b, a = res.params
        b_err, a_err = res.bse
        T_space = np.linspace(slopes[0,0], slopes[-1,0], int(1e3))
        err_precision = int(np.ceil(-np.log10(a_err)))
        plt.plot(T_space, T_space*a + b, color = rgb_color, linestyle ="--", label = f"Linear fit: a = {a:.{err_precision}f}" + r" $\pm$ " + f"{a_err:.1g}")




        plt.xlabel(r"$T/\epsilon k_B^{-1}$", fontsize = 14)
        plt.ylabel(r"$D / \tau \sigma{^-2} $", fontsize = 14)
        plt.tight_layout(pad=1.1, w_pad=0.7, h_pad=0.2)
        plt.legend(fontsize = 13)
        plt.savefig("../article/figures/diffusion_relation.pdf", bbox_inches="tight")






        fig3 = plt.figure(num=3, dpi=80, facecolor='w', edgecolor='k')
        t_end = 0.5
        end_idx = int(t_end/dt)
        for i in range(len(avg_r2)):
            t = (avg_r2[i,0,:] - avg_r2[i,0,0])*dt
            plot = plt.plot(t[:end_idx], avg_r2[i,1,:end_idx], label = f"T = {temp_run[i]}")

        ax = plt.gca()
        ylim = ax.get_ylim()
        v_0 = 1.3
        qua = v_0**2*t**2
        plt.plot(t[:end_idx], qua[:end_idx], "--", label =  "Quadratic")
        ax.set_ylim(ylim)

        plt.xlabel(r"$t/\tau$", fontsize = 14)
        plt.ylabel(r"$\langle r^2(t) \rangle / \sigma^2$", fontsize = 14)
        plt.tight_layout(pad=1.1, w_pad=0.7, h_pad=0.2)
        plt.legend(fontsize = 13)
        plt.savefig("../article/figures/diffusion_start.pdf", bbox_inches="tight")
        #


        plt.show()
