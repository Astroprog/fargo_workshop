from makedustopac import *
import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt
import os

def createOpacityFile(grain_size, wavelengths):
    logawidth = 0.05          # Smear out the grain size by 5% in both directions
    na        = 20            # Use 20 grain size samples
    chop      = 5.            # Remove forward scattering within an angle of 5 degrees
    optconst  = "pyrmg70"     # The optical constants name
    optconstfile= optconst+'.lnk'
    matdens   = 3.0           # The material density in gram / cm^3
    extrapol  = True          # Extrapolate optical constants beyond its wavelength grid, if necessary
    verbose   = True          # If True, then write out status information
    ntheta    = 181           # Number of scattering angle sampling points
    theta     = np.linspace(0.,180.,ntheta)

    print("Running the code. Please wait...")
    opac       = compute_opac_mie(optconstfile, matdens, grain_size, wavelengths, theta=theta,
                              extrapolate=extrapol, logawidth=logawidth, na=na,
                              chopforward=chop, verbose=verbose)
    print("Writing the opacity to kappa file")
    write_radmc3d_kappa_file(opac, "{:.1e}_micron".format(grain_size*1e4))


wavelengths = np.logspace(-1,3,200)*1e-4
grain_size = 1e-2

createOpacityFile(grain_size, wavelengths)

