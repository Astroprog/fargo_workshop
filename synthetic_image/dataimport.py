import numpy as np

class DataImport:
    def __init__(self, data, output_number, subsample=1, fluids=None, grain_sizes=None, stride=1, **kwargs):
        """ Reads a simscripts/simdata data object
        input:
            output_number:  number of output snapshot
            subsample:      factor of reduction in resolution
            fluids:         array of strings matching available fluid keys
        """
        self.data = data
        self.dim = "2d"
        self.output_number = output_number

        self.density = dict()
        if fluids is not None:
            fluid_keys = fluids
        else:
            fluid_keys = self.data.fluids.keys()

        if grain_sizes is not None:
            self.grain_sizes = grain_sizes
        else:
            self.grain_sizes = None

        for key in fluid_keys:
            self.density[key] = self.data.fluids[key].get(self.dim, "mass density", self.output_number, stride=stride, **kwargs)
            if key == "gas" and len(self.data.fluids) == 1:
                self.density[key] *= 1e-2
            self.density[key] = subsample_field(self.density[key], subsample)

    def get_grain_sizes(self):
        """ Returns array with grain sizes
        """
        if self.grain_sizes is not None:
            return self.grain_sizes
        else:
            amin = self.data.parameters["amin"]
            amax = self.data.parameters["amax"]
            nfluids = len(self.data.fluids)
            if nfluids == 2:
                grains = np.array([amin])
                selected_sizes = np.array([0])
            else:
                grains = amin * (amax / amin)**(np.arange(nfluids - 1) / (nfluids - 2.0))
                dust_keys = [key for key in self.density.keys() if key != "gas"]
                selected_sizes = [int(key[-1]) - 1 for key in dust_keys]
                sizes = grains[selected_sizes].sort()
            return grains[selected_sizes]

    def get_parameters(self):
        return self.data.parameters

    def get_temperature(self):
        return self.data.fluids[self.fluid].get(self.dim, "temperature", self.output_number)

    def get_densities(self):
        return self.density

    def get_radius(self, interface=False):
        field = list(self.density.values())[0]
        if interface:
            return field.grid.get_interfaces("r")
        else:
            return field.grid.get_centers("r")

    def get_phi(self, interface=False):
        field = list(self.density.values())[0]
        if interface:
            self.phii = field.grid.get_interfaces("phi") + np.pi
        else:
            self.phic = field.grid.get_centers("phi") + np.pi


    def get_coordinates(self, name, interface=False):
        field = list(self.density.values())[0]
        if interface:
            return field.grid.get_interfaces(name)
        else:
            return field.grid.get_centers(name)

    def sample_down(self, factor):
        data_2d = self.sample_2d_density(data_2d, interpolation_factor=8)

    def sample_2d_density(self, data_2d, cached=True, interpolation_factor=1):
        interpolation_factor = int(interpolation_factor)
        points = np.column_stack((self.grid_rc, self.grid_phic))
        self.rc = self.rc[::interpolation_factor]
        self.phic = self.phic[::interpolation_factor]
        self.ri = self.ri[::interpolation_factor]
        self.phii = self.phii[::interpolation_factor]
        self.grid_rc, self.grid_phic = np.meshgrid(self.rc, self.phic)
        self.nrc = len(self.rc)
        self.nphic = len(self.phic)

        try:
            data_2d = np.load("interpolated_grid.npy")
        except FileNotFoundError:
            print(len(points), len(data_2d))
            data_2d = griddata(points, data_2d, (self.grid_rc, self.grid_phic), rescale=True).T
            if cached:
                np.save("interpolated_grid.npy", data_2d)

        return data_2d

def subsample_field(field, step):
    from simdata.field import Field
    # subsample data and grid
    new_grid = subsample_grid(field.grid, step)
    new_data = subsample_data(field.data, step)
    # construct new field
    new_field = Field( new_grid, new_data, field.time, field.name )
    return new_field

def subsample_data(data, step):
    # subsample data with quick and dirty approach
    # to do it better use area weighted averages
    dim = len(data.shape)
    if dim == 2:
        new_data = data[::step, ::step]
    else:
        new_data = data[::step, ::step, ::step]
    return new_data

def subsample_grid(grid, step):
    dim = grid.dim
    # subsample grid
    GridClass = type(grid)
    coord_names = [grid.names[n] for n in range(dim)]
    # sample interfaces
    interfaces = [grid.get_interfaces( coord_names[n])[::step] for n in range(dim)]
    # reconstruct new grid
    grid_args = { coord_names[n]+"_i" : interfaces[n] for n in range(dim) }
    grid_args["active_interfaces"] = grid.active_interfaces
    new_grid = GridClass(**grid_args)
    return new_grid
