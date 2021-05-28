import os

import numpy as np
import cmasher as cmr
import matplotlib.pyplot as plt




def plotField(path, field_name="gasdens", output_number=0, vmin=None, vmax=None, logscale=True, **kwargs):
    """
    Parameters:
        path:           path to output folder (e.g. fargo3d/outputs/fargo)
        cartesian:      plot in cartesian frame (circular plot)
        vmin, vmax:     lower and upper variable limits of the field to plot 

    """

    # load cell interface locations
    r = np.loadtxt(os.path.join(path, "domain_y.dat"))[3:-3] # ignore ghost cells
    th = np.loadtxt(os.path.join(path, "domain_z.dat"))[3:-3]

    # load field data
    field = np.fromfile(os.path.join(path, field_name+str(output_number)+".dat")).reshape(len(th)-1, len(r)-1)

    r_m, th_m = np.meshgrid(r, th)

    x_m = r_m * np.sin(th_m)
    y_m = r_m * np.cos(th_m)

    if logscale:
        plt.pcolormesh(x_m, y_m, np.log10(field), vmin=vmin, vmax=vmax, **kwargs)
    else:
        plt.pcolormesh(x_m, y_m, field, vmin=vmin, vmax=vmax, **kwargs)

    plt.colorbar(pad=0)


frame = 1
plt.style.use("dark_background")
plt.figure(figsize=(12, 7))
plotField("./fargo3d/outputs/vsi/", 
        field_name="gasvz", 
        output_number=frame, 
        vmin=-1e-2, 
        vmax=1e-2, 
        logscale=False,
        cmap="cmr.iceburn")
plt.show()

plt.style.use("dark_background")
plt.figure(figsize=(8, 7))
plotField("./fargo3d/outputs/vsi/", 
        field_name="gasdens", 
        output_number=frame, 
        vmin=None, 
        vmax=None, 
        logscale=True,
        cmap="cmr.amber")
plt.show()
