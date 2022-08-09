from cluster_generator import *


if __name__ == "__main__":
    big_L = 1000
    big_cluster = TwoDim_cluster(big_L,0)



    L = 250
    p_c = 0.59275

    n = big_L/L

    if n%int(n) == 0:
        n = int(n)
    else:
        print(f"{big_L}/{L} = {n} must be integer")


    for i in range(n):
        for j in range(n):
            spanning = False
            print(f"\r i = {i}, j = {j}", end = "")
            while spanning == False:
                sub_cluster = TwoDim_cluster(L,p_c)
                spanning = sub_cluster.is_spanning()


            big_cluster.m[i*L:(i+1)*L, j*L:(j+1)*L] = sub_cluster.m
            big_cluster.lw, big_cluster.num = measurements.label(big_cluster.m)
            big_cluster.area = measurements.sum(big_cluster.m, big_cluster.lw, index = np.arange(big_cluster.lw.max()+1))


            plt.figure(num=0, dpi=100, facecolor='w', edgecolor='k')
            areaImg = big_cluster.area[big_cluster.lw]
            plt.imshow(areaImg, origin = "upper")
            # plt.colorbar()
            plt.tight_layout(pad=1.1, w_pad=0.7, h_pad=0.2)
            plt.savefig(f"../article/figures/subset_gif/subset_{i}{j}.png", bbox_inches="tight")
            plt.clf()
            # big_cluster.show()
    print("\n done")

# plt.figure(num=0, dpi=80, facecolor='w', edgecolor='k')
# plt.xlabel(r"$x$", fontsize=14)
# plt.ylabel(r"$y$", fontsize=14)
# plt.legend(fontsize = 13)



big_cluster.lw, big_cluster.num = measurements.label(big_cluster.m)
big_cluster.show()
