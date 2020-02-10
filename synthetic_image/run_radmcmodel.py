import numpy as np
from astropy import constants as const
from astropy import units as u

from radmcmodel import RadmcModel

parameters = dict()
parameters["nphot"] = 1e7
parameters["nphot_scat"] = 1e6
parameters["aspectratio"] = 0.05
parameters["flaringindex"] = 0.25
parameters["alpha"] = 1e-3
parameters["stellarmass"] = 1.0
parameters["stellartemperature"] = 1e4
parameters["meanmolecularweight"] = 2.4
parameters["inclination"] = 20.0
parameters["posang"] = 180.0
parameters["frequency"] = 230 * u.GHz
parameters["wavelength"] = 1300
parameters["beam_fwhm"] = [40.0, 40.0] 
parameters["beam_posang"] = 0.0 
parameters["distance"] = 100.0
parameters["rhosolid"] = 3.0
parameters["fluids"] = ["dust3"]
parameters["stokes_numbers"] = np.array([0.1])
sigma0 = 1e-4 * const.M_sun.cgs / (5.2 * const.au.cgs)**2
a = 2.0 * parameters["stokes_numbers"] / np.pi / parameters["rhosolid"] * sigma0.value
print("Using grain sizes: {}".format(a))

parameters["grain_sizes"] = a

simpath = "/Users/peter/projects/phd/docs/talks/transdisk_bormio/exercises/repo/fargo3d/outputs/fargo_multifluid"
model = RadmcModel(parameters, simpath)
model.create_setup(300, dustSettling=True)
model.make_image()

