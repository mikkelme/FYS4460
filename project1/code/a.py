import numpy as np
import matplotlib.pyplot as plt
from read_dump import read_dump_ovito
import plot_set


filename = "1a.data"
dump_file = read_dump_ovito(filename)
vel = dump_file.get_velocities()
numTimesteps = dump_file.numTimesteps

vel_len = np.linalg.norm(vel[:,:,-3:], axis = -1)
bins = 50

def Maxwell_Boltzmann(v):
    m = 1
    T = 2.5
    k = 1
    return np.sqrt((m/(2*np.pi*k*T))**3)*4*np.pi*v**2*np.exp(-(m*v**2)/(2*k*T))


#Histogram array
hist_array = np.zeros((numTimesteps,bins))
for i in range(numTimesteps):
    hist = np.histogram(vel_len[i], bins=bins, range=None, normed=None)
    hist_array[i] = hist[0]

#Calculate inner product
HH = np.sum(hist_array[-1]**2)
inner_product = np.zeros(numTimesteps)
for i in range(numTimesteps):
    inner_product[i] = np.sum(hist_array[i,:]*hist_array[-1,:])/HH



# Plotting
fig1 = plt.figure(num = 1, dpi=80, facecolor='w', edgecolor='k')
dt = 0.005 #tau
t = np.linspace(0, numTimesteps-1, numTimesteps)*10*dt
plt.plot(t, inner_product)
plt.xlabel(r"$t/\tau$", fontsize = 14)
plt.ylabel("Inner product", fontsize = 14)
plt.tight_layout(pad=1.1, w_pad=0.7, h_pad=0.2)


fig2 = plt.figure(num = 2, dpi=80, facecolor='w', edgecolor='k')
v = np.linspace(0, 10, int(1e3))
plt.hist(vel_len[-1], bins = bins, density = True, label = "Numerical distribution")
plt.plot(v, Maxwell_Boltzmann(v), label = "Theoretical distribution")
plt.xlabel(r"$v/\sigma \tau^{-1}$", fontsize = 14)
plt.ylabel("Normalized final distribution", fontsize = 14)
plt.legend(fontsize = 13)
plt.tight_layout(pad=1.1, w_pad=0.7, h_pad=0.2)


fig3 = plt.figure(num = 3, dpi=80, facecolor='w', edgecolor='k')
v = np.linspace(0, 10, int(1e3))
plt.hist(vel_len[0], bins = bins, density = True, label = "Numerical distribution")
plt.plot(v, Maxwell_Boltzmann(v), label = "Theoretical distribution")
plt.xlabel(r"$v/\sigma \tau^{-1}$", fontsize = 14)
plt.ylabel("Normalized initial distribution", fontsize = 14)
plt.legend(fontsize = 13)
plt.tight_layout(pad=1.1, w_pad=0.7, h_pad=0.2)

fig1.savefig("../article/figures/inner_product.pdf", bbox_inches="tight")
fig2.savefig("../article/figures/MB_dist.pdf", bbox_inches="tight")
fig3.savefig("../article/figures/MB_initial_dist.pdf", bbox_inches="tight")

# print(dump_file)
plt.show()
