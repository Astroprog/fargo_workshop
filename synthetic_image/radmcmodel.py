import os
import numpy as np
import multiprocessing
from astropy import constants as const
from astropy import units as u
from scipy.interpolate import griddata
from radmc3dPy import image
from makedustopac import *



class RadmcModel:
    def __init__(self, parameters, simpath):
        self.parameters = parameters
        self.simpath = simpath
        self.model_dir = "./"

    def create_setup(self, output_number, dustSettling=True):
        self.dustSettling = dustSettling
        self.create_grid(nthetac=32)
        densities_2d = self.get_densities_2d(output_number)
        densities_3d, temperature = self.puff_up_2d_densities_temperature(densities_2d)
        self.write_opacities()
        self.write_radmc_input_files(densities_3d, temperature, data_dtype=float, binary=True)

    def make_image(self, image_name="output.fits", npix=600, phi=0.0, convolve=False):

        inclination = self.parameters["inclination"]
        posang = self.parameters["posang"]
        wavelength = self.parameters["wavelength"]
        print(self.domain_size())
        image.makeImage(npix=npix, incl=inclination, posang=posang,
                        phi=phi, wav=wavelength, binary=True,
                        sizeau=self.domain_size())   # This calls radmc3d


        im = image.readImage(binary=True)
        self.parameters["sizepix_x"] = im.sizepix_x
        self.parameters["sizepix_y"] = im.sizepix_y
        if len(im.x) == 0:
            raise RuntimeError("Image creation failed, image is empty")

        if convolve:
            im = im.imConv(fwhm=1e-3*np.array(self.parameters["beam_fwhm"]), 
                           pa=1e-3*self.parameters["beam_posang"], 
                           dpc=self.parameters["distance"])

        im.raw_input = lambda s: "y"
        image.plotImage(im, log=True, maxlog=2, dpc=101.0, cmap='inferno')
        # im.writeFits(image_name, dpc=self.parameters["distance"], 
        #              coord='17h56m21.2882188601s -21d57m21.872343282s')

    def get_densities_2d(self, output_number):
        densities_2d = dict()
        for key in self.parameters["fluids"]:
            densities_2d[key] = np.fromfile(os.path.join(self.simpath, key+"dens"+str(output_number)+".dat")).reshape(len(self.rc), len(self.phic)) * const.M_sun.cgs / (5.2 * const.au.cgs)**2
        return densities_2d

    def get_temperature(self):
        aspectRatio = self.parameters["aspectratio"]
        flaringIndex = self.parameters["flaringindex"]
        stellarMass = self.parameters["stellarmass"]
        mean_mol_weight = self.parameters["meanmolecularweight"]

        cs = aspectRatio * (self.rc / (5.2 * const.au.cgs))**flaringIndex * np.sqrt(const.GM_sun * stellarMass / self.rc) #Warning: use (Mstar+Mplanet)!
        T = (mean_mol_weight * const.m_p / const.k_B * cs**2)
        T = np.tile(T, (self.nphic, self.nthetac, 1))
        return T.cgs.value

    def puff_up_2d_densities_temperature(self, densities_2d):
        """Given a planar profile (for now {r,φ} only), returns an {r,θ,φ} profile.
        This is done by inflating the disk along the meridian direction (θ) using
        the midplane scale height H:

        ρ(r,θ,φ) = Σ(r,φ) * exp(-(z(r,θ)/H(r,φ))^2/2),

        where z = r*cos(θ)

        The scale height needs to be calculated:

        H = c_s/Ω_K = sqrt(R*T*r^3/(μ*G*Msys))

        Therefore, the function will attempt to read a temperature file.

        If it can't, it will roll back to using a constant aspect ratio profile.

        Input:
            densities_2d: the 2D data structure containing the Σ(r,φ) profiles for each dust species. 
            Msys: (stellar+planet) mass in solar masses. Needed for the scale height.
            mu: mean molecular weight. Needed for the scale height.
        Output:
            densities_3d: the corresponding 3D density profiles.
            temperature_3d: the corresponding 3D temperature profile.
        """

        alpha = self.parameters["alpha"]
        mu = self.parameters["meanmolecularweight"]
        Msys = self.parameters["stellarmass"]

        aspectratio = self.parameters["aspectratio"]
        flaringindex = self.parameters["flaringindex"]
        H = self.rc * (aspectratio * self.rc / (5.2*const.au))**flaringindex     #dims: Nr*1 -> Nr

        temperature_3d = self.get_temperature()
        # temperature_3d = np.rollaxis(temperature_3d, 2) #dims: (Nθ,Nφ,Nr) -> (Nr,Nθ,Nφ)
        temperature_3d = np.swapaxes(temperature_3d, 0, 2) #dims: (Nθ,Nφ,Nr) -> (Nr,Nθ,Nφ)

        # density_gas_2d_T = densities_2d["gas"].data.T
        densities_3d = dict()
        for fluid, key in enumerate(densities_2d):
            if key != "gas":
                density_3d = np.zeros((self.nthetac, self.nphic, self.nrc))            #dims: Nθ,Nφ,Nr
                density_2d_T = densities_2d[key].T

                # a = grain_sizes[fluid-1]
                # St = 0.5 * np.pi * a * u.cm * self.parameters["rhosolid"] * u.g / u.cm**3 / density_gas_2d_T
                St = self.parameters["stokes_numbers"][fluid]
                density_2d_normalized = density_2d_T / (np.sqrt(2*np.pi)*H)            #dims: (Nr,Nφ).Τ/Nr -> (Nφ,Nr)/Nr -> Nφ,Nr

                if self.dustSettling is True:
                    dustSettlingCorrection = np.sqrt(alpha / (St * (1.0 + St**2)))
                    H *= dustSettlingCorrection

                for k in range(self.nthetac):
                    zs = self.rc * np.cos(self.thetac[k])                #dims: Nr/1 -> Nr
                    puff_factor = np.exp(-0.5*(zs/H)**2)                 #dims: Nr/(Nφ,Nr) -> Nφ,Nr OR Nr/Nr -> Nr
                    density_3d[k] = density_2d_normalized * puff_factor  #dims: (Nφ,Nr) * ((Nφ,Nr) OR (Nr)) -> Nφ,Nr
            
                density_3d = np.rollaxis(density_3d, 2)         #dims: (Nθ,Nφ,Nr) -> (Nr,Nθ,Nφ)
                densities_3d[key] = density_3d

        return densities_3d, temperature_3d


    def create_grid(self, nthetac=1):
        self.load_coordinates()
        self.grid_rc, self.grid_phic = np.meshgrid(self.rc, self.phic)
        self.grid_rc = np.ravel(self.grid_rc)
        self.grid_phic = np.ravel(self.grid_phic)

        self.nrc = len(self.rc)
        self.nphic = len(self.phic)
        self.nr = self.nrc + 1
        self.nphi = self.nphic + 1
        self.ntheta = nthetac + 1
        self.nthetac = nthetac

        theta_m  = np.pi * 0.5 - 0.1
        theta_p  = np.pi * 0.5 + 0.1
        self.thetai   = np.linspace(theta_m, theta_p, nthetac + 1)
        self.thetac   = 0.5 * (self.thetai[0:nthetac] + self.thetai[1:nthetac + 1])

    def load_coordinates(self):
        self.ri = np.loadtxt(os.path.join(self.simpath, "domain_y.dat"))[3:-3] * 5.2 * const.au.cgs# ignore ghost cells
        self.phii = np.loadtxt(os.path.join(self.simpath, "domain_x.dat"))
        self.rc = 0.5 * (self.ri[1:] + self.ri[:-1])
        self.phic = 0.5 * (self.phii[1:] + self.phii[:-1])
        self.map_phi()

    def map_phi(self):
        # RADMC3D wants to have phi in [0, 2pi] not in [-pi,pi]
        #self.phic += 2*np.pi*(self.phic<0)
        #elf.phii += 2*np.pi*(self.phii<0)
        if np.min(self.phic) < 0.0:
            self.phic += np.pi
            self.phii += np.pi

    def domain_size(self):
        self.ri = 5.2*np.loadtxt(os.path.join(self.simpath, "domain_y.dat"))[3:-3]
        self.rc = 0.5*(self.ri[1:] + self.ri[:-1])
        return 2*np.max(self.rc)

    def createOpacityFile(self, grain_size, wavelengths):
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

    def write_opacities(self):
        grain_sizes = self.parameters["grain_sizes"]
        l = self.parameters["wavelength"]
        wavelengths = np.array([0.95 * l, l, 1.05*l]) * 1e-4
        self.opacity_files = []
        for a in grain_sizes:
            self.createOpacityFile(a, wavelengths)
            self.opacity_files.append("{:.1e}_micron".format(a*1e4))

    def write_radmc_input_files(self, densities_3d, temperature_3d, data_dtype=float, binary=False):

        lam1     = 0.1e0
        lam2     = 7.0e0
        lam3     = 25.e0
        lam4     = 1.0e4
        n12      = 20
        n23      = 100
        n34      = 30
        lam12    = np.logspace(np.log10(lam1),np.log10(lam2),n12,endpoint=False)
        lam23    = np.logspace(np.log10(lam2),np.log10(lam3),n23,endpoint=False)
        lam34    = np.logspace(np.log10(lam3),np.log10(lam4),n34,endpoint=True)
        lam      = np.concatenate([lam12,lam23,lam34])
        nlam     = lam.size

        with open(os.path.join(self.model_dir, "wavelength_micron.inp"), 'w+') as f:
            f.write('%d\n'%(nlam))
            np.savetxt(f,lam.T,fmt=['%13.6e'])

        nstars = 1
        mstar1 = (self.parameters["stellarmass"] * const.M_sun.cgs).value
        rstar1    = const.R_sun.cgs.value
        tstar1    = self.parameters["stellartemperature"]
        pstar1    = np.array([0, 0., 0.])

        with open(os.path.join(self.model_dir, "stars.inp"), 'w+') as f:
            f.write('2\n')
            f.write('%d %d\n\n'%(nstars, nlam))
            f.write('%13.6e %13.6e %13.6e %13.6e %13.6e\n\n'%(rstar1,mstar1,pstar1[0],pstar1[1],pstar1[2]))
            np.savetxt(f,lam.T,fmt=['%13.6e'])
            f.write('\n%13.6e\n'%(-tstar1))

        with open(os.path.join(self.model_dir, "amr_grid.inp"), 'w+') as f:
            f.write('1\n')                          # iformat
            f.write('0\n')                          # AMR grid style  (0=regular grid, no AMR)
            f.write('100\n')                        # Coordinate system: spherical
            f.write('0\n')                          # gridinfo
            f.write('1 1 1\n')                      # Include r, theta, phi coordinates
            f.write('%d %d %d\n'%(self.nrc, self.nthetac, self.nphic))  # Size of grid
            np.savetxt(f, self.ri.T, fmt=['%13.6e'])      # R coordinates (cell walls)
            np.savetxt(f, self.thetai.T, fmt=['%13.6e'])  # Theta coordinates (cell walls)
            np.savetxt(f, self.phii.T, fmt=['%13.6e'])    # Phi coordinates (cell walls)
            f.write(' ')

        if binary:
            size = 8
            if data_dtype == np.float32:
                size = 4
                for data in densities_3d.values():
                    data = data.astype(np.float32)
            # the binary file expects the following header:
            # Line1: 1 (format, always '1')
            # Line2: sizeof(data) (8 for double, 4 for float. Using double by default.)
            # Line3: Ncells (Nx*Ny*Nz)
            # Line4: Nspecies (1 by default, doubt support for more will ever be needed/added)
            with open(os.path.join(self.model_dir, "dust_density.binp"), 'w') as out:
                header = np.array([1, size, list(densities_3d.values())[0].size, len(densities_3d)], dtype=int) #prep header
                header.tofile(out) #write header 
                for data in densities_3d.values():
                    data = data.swapaxes(0, 2) #dims: (Nr,Nθ,Nφ) -> (Nφ,Nθ,Nr)
                    data.tofile(out)   #write data
        else:
            with open(os.path.join(self.model_dir, "dust_density.inp"), 'w+') as f:
                f.write('1\n')                       # Format number
                f.write('%d\n' % (list(densities_3d.values())[0].size))     # Nr of cells
                f.write('%d\n' % (len(densities_3d)))                       # Nr of dust species
                for data in densities_3d.values():
                    data = data.ravel(order='F')         # Create a 1-D view, fortran-style indexing
                    np.savetxt(f, data.T, fmt=['%13.6e'])  # The data

        if binary:
            data = temperature_3d.swapaxes(0, 2) #dims: (Nr,Nθ,Nφ) -> (Nφ,Nθ,Nr)
            size = 8
            if data_dtype == np.float32:
                size = 4
                data = data.astype(np.float32)
            # see density_3d section for details
            with open(os.path.join(self.model_dir, "dust_temperature.bdat"), 'w') as out:
                header = np.array([1, size, len(data.flatten()), len(densities_3d)], dtype=int) #prep header
                header.tofile(out) #write header 
                for fluid in densities_3d.keys():
                    data.tofile(out)   #write data
        else:
            with open(os.path.join(self.model_dir, "dust_temperature.dat"), 'w+') as f:
                f.write('1\n')
                f.write('%d\n' % (self.nrc * self.nthetac * self.nphic))     # Nr of cells
                f.write('%d\n' % (len(densities_3d)))                       # Nr of dust species
                data = temperature_3d.ravel(order='F')         # Create a 1-D view, fortran-style indexing
                for fluid in densities_3d.keys():
                    np.savetxt(f, data.T, fmt=['%13.6f'])  # The data

        with open(os.path.join(self.model_dir, "dustopac.inp"), 'w+') as f:
            f.write('2               Format number of this file\n')
            f.write('{}               Nr of dust species\n'.format(len(self.opacity_files)))
            f.write('============================================================================\n')

            for opacity_file in self.opacity_files:
                f.write('1               Way in which this dust species is read\n')
                f.write('0               0=Thermal grain\n')
                f.write('{}      Extension of name of dustkappa_***.inp file\n'.format(opacity_file))
                f.write('----------------------------------------------------------------------------\n')

        nphot    = self.parameters["nphot"]
        nphot_scat    = self.parameters["nphot_scat"]
        with open(os.path.join(self.model_dir, "radmc3d.inp"), 'w+') as f:
            f.write('nphot = %d\n'%(int(nphot)))
            f.write('nphot_scat = %d\n'%(int(nphot_scat)))
            f.write('scattering_mode_max = 1\n')
            f.write('iranfreqmode = 1\n')
            f.write('mc_scat_maxtauabs = 5.d0\n')
            f.write('setthreads = {}\n'.format(multiprocessing.cpu_count()-2))
            if binary: f.write('rto_style = 3\n')
