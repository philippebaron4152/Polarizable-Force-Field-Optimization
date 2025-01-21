import numpy as np
import lammps_logfile

density = np.mean(np.array(lammps_logfile.File("log_run.lammps").get("Density"))[1000:])*1000
with open("density.txt", "w+") as f:
    f.write(str(density) + "\n")
