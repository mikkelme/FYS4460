import numpy as np
import ovito

class read_dump_ovito:
    """ Class for reading lammps file using ovito module
        -   Each instance should correspond
            to a single dump file import
    """
    def __init__(self, filename):
        self.filename = filename
        self.pipeline = ovito.io.import_file(filename, multiple_frames = True)
        self.data0 = self.pipeline.compute(0)
        self.properties = list(self.data0.particles.keys())
        self.numTimesteps = self.pipeline.source.num_frames
        self.numParticles = self.data0.particles.count

    def get_positions(self):
        self.position = np.zeros((self.numTimesteps, self.numParticles, 3))
        for frame in range(self.numTimesteps):
            frame_data = self.pipeline.compute(frame)
            self.position[frame] = np.array(frame_data.particles['Position'])
        return self.position

    def get_velocities(self):
        self.velocity = np.zeros((self.numTimesteps, self.numParticles, 3))
        for frame in range(self.numTimesteps):
            frame_data = self.pipeline.compute(frame)
            self.velocity[frame] = np.array(frame_data.particles['Velocity'])
        return self.velocity

    def __str__(self): #Print data info
        string = f"#--- Data info ---#\
        \nFilename: {self.filename}\
        \nnumber of timesteps: {self.numTimesteps}\
        \nnumber of particles: {self.numParticles}\
        \nProperties: {self.properties}\n"
        return string


def read_dump_manually(filename, sort_id = True):
    """ Reads lammps dump file
        by manually indexing """

    # Open file
    infile = open(filename, "r")

    # Find nmber of atoms and number of lines to skip between timesteps
    numAtoms = 0
    skip = 0
    for line in infile:
        skip += 1
        if line == "ITEM: NUMBER OF ATOMS\n":
            line = infile.readline()
            skip += 1
            numAtoms = int(line)
        if line[:11] == "ITEM: ATOMS":
            break

    # Read data
    data = [] #timestep, atom, property
    while True:
        placeholder = []
        for i in range(numAtoms): # Append line to data list
            placeholder.append(np.array(infile.readline().split()).astype(float))
        data.append(placeholder)
        try: #Skip info lines
            for i in range(skip):
                next(infile)
        except: #Break when at the bottom of file
            print("Done reading")
            break
    data = np.array(data) # list to array
    numTimesteps = len(data)


    # Sort by id
    idx = np.argsort(data[0,:,0])
    if sort_id:
        for i in range(numTimesteps):
            data[i] = data[i, idx]
    return data




def read_binary_dump(filename, colnames=("id", "type", "x", "y", "z", "vx", "vy", "vz")):
    """ Not done yet """
    import ase.io
    frames = ase.io.read(filename, index=":", format="lammps-dump-binary", colnames=colnames)


    vel = frames[0].get_velocities()
    for i in range(1, len(frames)):
        vel = np.append(vel, frames[i].get_velocities(), axis = 0)

    pos = frames[0].get_positions()
    for i in range(1, len(frames)):
        pos = np.append(vel, frames[i].get_positions(), axis = 0)

    return pos, vel






if __name__ == "__main__":
    print("in read_dump")
    filename = "dump.data"

    # ovito_file = read_dump_ovito(filename)
    # vel = ovito_file.get_velocities()
    # pos = ovito_file.get_positions()
