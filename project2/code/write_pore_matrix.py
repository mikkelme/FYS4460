from sphere_matrix_generator import matrix_generator
import numpy as np




if __name__ == "__main__":
    N_pores = 2    # Number of pores
    b = 5.72        # Unit cell length [Å]
    R0_min = 20     # min radius [Å]
    R0_max = 30     # max radius [Å]
    L_min = 0       # min pos [Å]
    L_max = 20*b    # max pos [Å]
    sigma = 3.4     # scaling factor [Å]

    generator = matrix_generator(N_pores, R0_min, R0_max, L_min, L_max, convert = sigma)
    generator.generate()
    generator.write_for_lammps()
