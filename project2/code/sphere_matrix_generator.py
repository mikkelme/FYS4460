import numpy as np
import matplotlib.pyplot as plt
import random

class matrix_generator:
    def __init__(self, N_pores, R0_min, R0_max, L_min, L_max, convert = False):
        self.N_pores = N_pores
        self.R0_min = R0_min
        self.R0_max = R0_max
        self.L_min = L_min
        self.L_max = L_max

        #Scale by sigma
        if convert != False:
            self.R0_min /= convert
            self.R0_max /= convert
            self.L_min /= convert
            self.L_max /= convert

    def generate(self):
        """ Assumes cubic box """
        self.pos = self.L_min + np.random.rand(self.N_pores,3)*(self.L_max-self.L_min)
        self.R0 = self.R0_min + np.random.rand(self.N_pores)*(self.R0_max-self.R0_min)
        return self.pos, self.R0

    def show_matrix(self):
        axis = np.array([[0,1], [0,2], [1,2]])
        axis_text = [r"x [$\sigma$]", r"y [$\sigma$]", r"z [$\sigma$]"]
        for i in range(3):
            fig = plt.figure(num = i); ax = fig.gca()
            plt.plot(self.pos[:,axis[i,0]], self.pos[:,axis[i,1]], "o", label = ("sphere center"))
            for j in range(self.N_pores):
                circle = plt.Circle((self.pos[j,0], self.pos[j,1]), self.R0[j], alpha = 0.2)
                ax.add_patch(circle)
            plt.axis("equal")
            plt.xlim([self.L_min, self.L_max])
            plt.ylim([self.L_min, self.L_max])
            plt.xlabel(axis_text[axis[i,0]])
            plt.ylabel(axis_text[axis[i,1]])
            plt.legend()
        plt.show()

    def inside_sphere(self, x,y, z):
        for i in range(len(self.R0)):
            if (x - self.pos[i,0])**2 + (y - self.pos[i,1])**2 + (z - self.pos[i,2])**2  < self.R0[i]**2:
                return 1
        return 0


    def calculate_porosity(self, num):
        """ numerical calculation of phi """
        #define grid
        N_atoms = num**3
        L = self.L_max - self.L_min
        dL = L/num

        x_ = np.linspace(1/2*dL, L -1/2*dL, num)
        y_ = np.linspace(1/2*dL, L -1/2*dL, num)
        z_ = np.linspace(1/2*dL, L -1/2*dL, num)

        count = 0
        print("\n#--- Calculating phi ---#")
        print(f"Number of pseudo atoms: {N_atoms}")
        for x in x_:
            print(f"\rProgress: {x/x_[-1]*100:.2f}%", end = "")

            for y in y_:
                for z in z_:
                    count += self.inside_sphere(x,y,z)
        self.phi = count * dL**3 / L**3
        phi_ana = self.N_pores*4/3*np.pi*np.mean(self.R0)**3/L**3

        print(f"\nNumerical phi = {self.phi:g}")
        print(f"Analytical phi (no overlap or cut off) = {phi_ana:g}")
        return self.phi

    def write_for_lammps(self, filename = "pore_matrix.in"):
        region_free = f"region r_free intersect {self.N_pores}"
        region_locked = f"region r_locked union {self.N_pores}"
        with open(filename, "w") as infile:
            for i in range(self.N_pores):
                infile.write(f"region r_in{i} sphere {self.pos[i,0]} {self.pos[i,1]} {self.pos[i,2]} {self.R0[i]} side in units box\n")
                infile.write(f"region r_out{i} sphere {self.pos[i,0]} {self.pos[i,1]} {self.pos[i,2]} {self.R0[i]} side out units box\n")
                region_free += f" r_out{i}"
                region_locked += f" r_in{i}"

            infile.write("\n" + region_free)
            infile.write("\n" + region_locked)
            infile.write("\n" + "group g_free region r_free" )
            infile.write("\n" + "group g_locked region r_locked")


            # infile.write(group_free)




if __name__ == "__main__":
    N_pores = 20    # Number of pores
    b = 5.72        # Unit cell length [Å]
    R0_min = 20     # min radius [Å]
    R0_max = 30     # max radius [Å]
    L_min = 0       # min pos [Å]
    L_max = 20*b    # max pos [Å]
    sigma = 3.4     # scaling factor [Å]

    generator = matrix_generator(N_pores, R0_min, R0_max, L_min, L_max, convert = sigma)
    generator.generate()
    # generator.show_matrix()
    generator.calculate_porosity(50)
