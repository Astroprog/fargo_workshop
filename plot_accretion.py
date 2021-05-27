import os
import numpy as np
import matplotlib.pyplot as plt


path = "./fargo3d/outputs/accretion/"
field_name = "gasdens"
frames = 11

r = 5.2 * np.loadtxt(os.path.join(path, "domain_y.dat"))[3:-3] # ignore ghost cells
r_c = 0.5*(r[1:] + r[:-1])


plt.figure(figsize=(12, 7))

colors = plt.cm.viridis(np.linspace(0, 1, frames))

for frame in range(frames):
    field = np.fromfile(os.path.join(path, field_name + str(frame) + ".dat"))
    plt.loglog(r_c, field, color=colors[frame])

plt.xlabel("r [au]")
plt.ylabel("Surface Density")
plt.show()
