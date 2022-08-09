import numpy as np
import matplotlib.pyplot as plt
# from read_dump import read_dump_ovito
# import subprocess
import plot_set


from b import read_energy



dt_read = [0.02, 0.01, 0.008, 0.002]
energy_files = []
for dt in dt_read:
    dt_str = str(dt).split(".")[-1]
    energy_files.append("etotal" + dt_str + ".txt")

# Plot temperature over time
fig_Temp1 = plt.figure(num = 0, dpi=80, facecolor='w', edgecolor='k')
k = 1; N = 4000

temp_exact = 2.5
Avg_temp = np.zeros(len(energy_files))
std = np.zeros(len(energy_files))
for i in range(len(energy_files)):
    data = read_energy(energy_files[i])
    start_idx = len(data)//4
    dt = float("0." + energy_files[i].strip("etotal.txt"))
    Temp = 2/3*data[:,1]/(N*k)
    Avg_temp[i] = np.mean(Temp[start_idx:])
    std[i] = np.std(Temp[start_idx:])
    # label = f"dt = {dt:.3f}" + r", $\langle T \rangle$ = " + f"{Avg_temp[i]:.4f}, " + r"$\sigma$ = " + f"{std[i]:.4f}, " + r"$\Delta T = $ " + f"{Avg_temp[i]-temp_exact:.2f}"
    label = f"dt = {dt:.3f}"
    # print(f"dt = {dt}, avg T = {Avg_temp[i]:.4f}, SE = {std[i]:.4f}, offset = {Avg_temp[i]-temp_exact:.4f}")
    plt.plot(data[start_idx:,0]*dt, Temp[start_idx:], label = label)
plt.xlabel(r"$t/\tau$", fontsize = 14)
plt.ylabel(r"$T/\epsilon k_B^{-1}$ ", fontsize = 14)
plt.tight_layout(pad=1.1, w_pad=0.7, h_pad=0.2)
plt.legend(loc = "best", fontsize = 13)


#fig_Temp1.savefig("../article/figures/temp_pp.pdf", bbox_inches="tight")

# fig_Temp1_avg = plt.figure(num = 1, dpi=80, facecolor='w', edgecolor='k')
# for i in range(len(energy_files)):
#     label = f"dt = {dt}" + r", $\langle T \rangle$ = " + f"{Avg_temp[i]:.4f}, " + r"SE = " + f"{std[i]:.4f}"
#     plt.plot(data[start_idx:,0]*dt, Avg_temp[i] + 0*data[start_idx:,0]*dt, label = label)
#     plt.xlabel(r"$t/\tau$", fontsize = 14)
#     plt.ylabel(r"$T/\epsilon k_B^{-1}$ ", fontsize = 14)
#     plt.tight_layout(pad=1.1, w_pad=0.7, h_pad=0.2)
#     plt.legend(loc = "best", fontsize = 13)


#Plot temperature over time depending on system size
energy_files = ["etot_size10.txt", "etot_size12.txt", "etot_size14.txt", "etot_size16.txt"]
system_size = [10, 12, 14, 16]
k = 1
fig_Temp2 = plt.figure(num = 2, dpi=80, facecolor='w', edgecolor='k')
for i in range(len(energy_files)):
    energy_file = energy_files[i]
    data = read_energy(energy_file)
    start_idx = len(data)//2
    size = system_size[i]
    N = 4*size**3
    Temp = 2/3*data[:,1]/(N*k)
    Avg_temp = np.mean(Temp[start_idx:])
    std = np.std(Temp[start_idx:])
    # label = f"Size: {size}x{size}x{size}, N = {N}, " + r"$\sigma$ = " + f"{std:.4f}, " + r"$\Delta T = $ " + f"{Avg_temp[i]-temp_exact:.2f}"
    label = f"Size: {size}x{size}x{size}, N = {N}"
    print(f"Size = {size}, avg T = {Avg_temp:.4f}, SE = {std:.4f}, offset = {Avg_temp-temp_exact:.4f}")
    plt.plot(data[start_idx:,0]*dt, Temp[start_idx:], label = label)
plt.xlabel(r"$t/\tau$", fontsize = 14)
plt.ylabel(r"$T/\epsilon k_B^{-1}$ ", fontsize = 14)
plt.tight_layout(pad=1.1, w_pad=0.7, h_pad=0.2)
plt.legend(fontsize = 13)
fig_Temp2.savefig("../article/figures/temp_size_pp.pdf", bbox_inches="tight")
plt.show()







#
