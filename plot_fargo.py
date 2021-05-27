import os

import numpy as np
import matplotlib.pyplot as plt



def plotField(path, field_name="gasdens", output_number=0, cartesian=True, vmin=None, vmax=None, logscale=True):
    """
    Parameters:
        path:           path to output folder (e.g. fargo3d/outputs/fargo)
        cartesian:      plot in cartesian frame (circular plot)
        vmin, vmax:     lower and upper variable limits of the field to plot 

    """

    # load cell interface locations
    r = np.loadtxt(os.path.join(path, "domain_y.dat"))[3:-3] # ignore ghost cells
    phi = np.loadtxt(os.path.join(path, "domain_x.dat"))

    # load field data
    field = np.fromfile(os.path.join(path, field_name+str(output_number)+".dat")).reshape(len(r)-1, len(phi)-1)
    # print(field.shape, field)

    r_m, phi_m = np.meshgrid(r, phi)

    if cartesian:
        x_m = r_m * np.cos(phi_m)
        y_m = r_m * np.sin(phi_m)
    else:
        x_m = phi_m
        y_m = r_m

    if logscale:
        plt.pcolormesh(x_m, y_m, np.log10(field.T), cmap="inferno", vmin=vmin, vmax=vmax)
    else:
        plt.pcolormesh(x_m, y_m, field.T, cmap="inferno", vmin=vmin, vmax=vmax)
    # limit = 150
    # plt.xlim(-limit, limit)
    # plt.ylim(-limit, limit)
    plt.colorbar(pad=0)


plt.figure(figsize=(8, 7))
plotField("./fargo3d/outputs/fargo_multifluid/", field_name="dust1dens", output_number=1, cartesian=True, logscale=True)
plt.show()
# plt.savefig("jupiter_0.05.png")
