import numpy as np
import re

# Constants
UNITLESS = ''
SI_UNIT = 1
C_UNIT_PREFIX = 1e-2
M_UNIT_PREFIX = 1e-3
U_UNIT_PREFIX = 1e-6
N_UNIT_PREFIX = 1e-9
K_UNIT_PREFIX = 1e3
UNIT_PREF_MAP = { 'c' : C_UNIT_PREFIX,
                  'm' : M_UNIT_PREFIX,
                  'u' : U_UNIT_PREFIX,
                  'n' : N_UNIT_PREFIX,
                  'k' : K_UNIT_PREFIX }

CONVERTABLE_UNIT_VEC = ['', 'rad', 'deg', 'm', 's']
CONVERSION_TABLE = np.array([[1, 1, 1, 1, 1],
                             [1, 1, 180.0 / np.pi, -1, -1],
                             [1, np.pi / 180.0, 1, -1, -1],
                             [1, -1, -1, 1, -1],
                             [1, -1, -1, -1, 1]])

def verify_coord_system(csyst : str):
    '''
    Verify a coordinate system name is accepted
    '''
    csyst = csyst.lower()
    if csyst not in ["cylindrical", "cartesian"]:
        raise ValueError(f"{csyst} is not a valid coordinate system."
                         "Only accepts 'cylindrical' or 'cartesian'")
    return

def set_coord_syst_member(obj, csyst : str):
    '''
    Set the coordinate system member of a simulation class
    '''
    csyst = csyst.lower()
    if csyst not in ["cylindrical", "cartesian"]:
        raise ValueError(f"{csyst} is not a valid coordinate system."
                         "Only accepts 'cylindrical' or 'cartesian'")
    obj.coord_syst = csyst
    return

def check_dim_type(dim_type : str, coord_syst = None):
    '''
    Check the dimension type is acceptable for the coordinate system.

    E.g., for cartesian, must be one of 'x', 'y', or 'z'
    '''
    if coord_syst == None:
        if dim_type not in ['t', 'x', 'y', 'z', 'r', 'theta', 'phi']:
            raise ValueError(f"Invalid dimension type: {dim_type}")
    elif coord_syst == "cartesian":
        if dim_type not in ['x', 'y', 'z']:
            raise ValueError(f"Invalid dimension type: {dim_type}")
    else:
        if dim_type not in ['z', 'r', 'theta']:
            raise ValueError(f"Invalid dimension type: {dim_type}")
        
    return

def verify_unit_for_dim_type(dim_type : str, units : str):
    '''
    Verify units make physical sense for dimension type
    '''
    unit_map = { 's' : ['t'],
                 'm' : ['x', 'y', 'z', 'r'],
                 'rad' : ['theta', 'phi'],
                 'deg' : ['theta', 'phi'] }

    if dim_type not in unit_map[units]:
        raise ValueError(f"Invalid unit ({units}) for dimension type {dim_type}")

    return

def convert_cyl_to_cart(r, theta, z):
    '''
    Convert a set of coordinates from cylindrical to cartesian
    '''
    return (r * np.cos(theta), r * np.sin(theta), z)

def convert_cart_to_cyl(x, y, z):
    '''
    Convert a set of coordinates from cartesian to cylindrical
    '''
    theta = np.arctan2(y,x)
    r = np.sqrt(x**2 + y**2)
    return (r, theta, z)


class SimulationVolume:
    """
    Class to made a simulated 3D geometry
    """
    __required_sections = ['geometry', 'x_0', 'x_1', 'x_2', 't']

    def __init__(self, gconfig : dict):
        self._ax_name_map = None
        self.coord_syst = None
        self.dist_vec = None
        self.dist_mag = None

        self._init_space(gconfig)     

    def _verify_config_sections(self, gconfig):
        '''
        Verify the configuration file sections are as expected
        '''
        for sect in self.__required_sections:
            if sect not in gconfig.keys():
                raise RuntimeError(f"Config is missing the {sect} section. Can't parse.")
            if sect == 'geometry':
                if 'coordinate_system' not in gconfig[sect].keys():
                    raise RuntimeError("Geometry section needs to have a coordinate system.")
                for sect2 in self.__required_sections[1:-1]:
                    if sect2 not in gconfig[sect].keys():
                        raise RuntimeError(f"Need the {sect2} element in {sect}")
        return

    def _set_coord_system(self, csyst : str):
        set_coord_syst_member(self, csyst)
        return

    def _set_axes_map(self, x_0_name, x_1_name, x_2_name):
        '''
        Set the map of generic axis names to name based on dimension type
        '''
        uniq_ax_names = sorted(set([x_0_name, x_1_name, x_2_name]))
        if len(uniq_ax_names) != 3:
            raise RuntimeError("Not enough unique axes names provided (need 3 for 3D space).")

        # check the names are expected for the coordinate systems
        if self.coord_syst is None:
            raise RuntimeError("Can't label axes with names without the coordinate system.")
        elif self.coord_syst == 'cylindrical':
            if uniq_ax_names != sorted(set(['r','z','theta'])):
                raise RuntimeError("Axes names need to be 'r', 'z', and 'theta'")
        else:
            if uniq_ax_names != sorted(set(['x','y','z'])):
                raise RuntimeError("Axes names need to be 'r', 'z', and 'theta'")
            
        self._ax_name_map = { 'x_0' : x_0_name,
                              'x_1' : x_1_name,
                              'x_2' : x_2_name }
        self._ax_name_map_inv = { value: key for key, value in self._ax_name_map.items() }

        return            

    def _init_space(self, gconfig : dict):
        # Read the 'geometry' section
        self._verify_config_sections(gconfig)
        g_section = gconfig['geometry']

        # Set up coordinate system and naming of axes
        self._set_coord_system(g_section['coordinate_system'])
        self._set_axes_map(g_section['x_0'], g_section['x_1'],
                           g_section['x_2'])

        self.x_0 = SimulationAxis(g_section['x_0'], **gconfig['x_0'])
        self.x_1 = SimulationAxis(g_section['x_1'], **gconfig['x_1'])
        self.x_2 = SimulationAxis(g_section['x_2'], **gconfig['x_2'])

        self.dim_labels = (self.x_0.dim_type, self.x_1.dim_type,
                           self.x_2.dim_type)
        self.units = (self.x_0.units, self.x_1.units, self.x_2.units)

        # Set up the time axis
        self.x_t = SimulationAxis('t', **gconfig['t'])

        # Get the number of bins and the volume bounds
        self.matrix_shape = (self.x_0.nbins, self.x_1.nbins, self.x_2.nbins)
        self.n_cells = self.x_0.nbins * self.x_1.nbins * self.x_2.nbins

        self.bounds = np.array([[self.x_0.ax_min, self.x_0.ax_max],
                                [self.x_1.ax_min, self.x_1.ax_max],
                                [self.x_2.ax_min, self.x_2.ax_max]])

        # Set arrays defining the cell extents
        self._set_volume_differentials()
        self._centers = None
        self._corners = None
        self._corner_set = ('u','u','u')

    def _set_volume_differentials(self):
        # Get first axis diffs
        if self.coord_syst == 'cartesian':
            self.dV = self.x_0.ax_edges[1:] - self.x_0.ax_edges[:-1]
            self.dV = self.dV[:,np.newaxis] * (self.x_1.ax_edges[1:] - self.x_1.ax_edges[:-1])
            self.dV = self.dV[:,:,np.newaxis] * (self.x_2.ax_edges[1:] - self.x_2.ax_edges[:-1])
        else:
            self.dV = None
            for i, key in enumerate(self._ax_name_map.keys()):
                ax = getattr(self, key)
                if ax.dim_type == 'r':
                    arr = np.pi * (ax.ax_edges[1:]**2 - ax.ax_edges[:-1]**2)
                else:
                    arr = ax.ax_edges[1:] - ax.ax_edges[:-1]
                    if ax.dim_type == 'theta':
                        # Normalize by the fraction of the full 2pi rad space
                        if ax.units == 'rad':
                            arr = arr / (2 * np.pi)
                        else:
                            arr = arr / 360. # normalize by degrees

                if self.dV is None:
                    self.dV = arr
                else:
                    self.dV = np.multiply.outer(self.dV, arr)

        return

    def _set_corner_set(self, c_str : str):
        c_str = c_str.lower()
        if len(c_str) != 3 or re.search('[^ul]', c_str) is not None:
            raise ValueError("Misformatted string for specifying cell corner. "
                             "Must be 3 characters long and only contain u/U (upper) or "
                             "l/L (lower)")
        corner_set = (c_str[0], c_str[1], c_str[2])
        if corner_set != self._corner_set:
            self._corner_set = corner_set
            return True
        else:
            return False

    def get_cell_centers(self):
        '''
        Retrieve the cell corners array
        '''
        if self._centers is None:
            # Get the center of each cell
            self._centers = np.zeros(self.matrix_shape + (3,))
            self._centers[:,:,:,0] = ((self.x_0.ax_edges[:-1] + self.x_0.ax_edges[1:])/2)[:,np.newaxis,np.newaxis]
            self._centers[:,:,:,1] = ((self.x_1.ax_edges[:-1] + self.x_1.ax_edges[1:])/2)[np.newaxis,:,np.newaxis]
            self._centers[:,:,:,2] = ((self.x_2.ax_edges[:-1] + self.x_2.ax_edges[1:])/2)[np.newaxis,np.newaxis,:]

        return self._centers
        
    def get_axes_titles(self):
        '''
        Get titles (dimension types) of the simulated volume
        '''
        return (self.x_0.dim_type, self.x_1.dim_type, self.x_2.dim_type)

    def get_axis(self, dim_label : str):
        ''' Get a single axis by label (dimension type) '''
        if dim_label not in self.dim_labels:
            raise KeyError('Invalid dimension requested for volume.')

        return getattr(key, f'x_{self.dim_labels.index("dim_label")}')
        
    def get_corners(self, c_str : str = ''):
        '''
        Retrieve the array of corner positions of each cell

        c_str is a 'corner string' to specify which corners are used. The format
        is 3 characters each of which can be 'u' or 'l', where 'u' represents upper
        and 'l' lower. So, 'lll' means the corner formed by the lower edge of each cell
        of each axis
        '''
        new_set = False
        if c_str != '':
            new_set = self._set_corner_set(c_str)

        if self._corners is None or new_set == True:
            self._corners = np.zeros(self.matrix_shape + (3,))
            for i, key in enumerate(['x_0', 'x_1', 'x_2']):
                ax = getattr(self, key)
                if self._corner_set[i] == 'u':
                    # Set using the upper edge
                    self._corners[:,:,:,i] = ax.ax_edges[1:]
                else:
                    # Set using the lower edge
                    self._corners[:,:,:,i] = ax.ax_edges[:-1]

        return self._corners

    def get_cell_center(self, i : int, j : int, k : int):
        ''' Get the center position of a cell with indices i,j,k '''
        self.get_cell_centers()
        return self._centers[int(i),int(j),int(k)]

    def get_cell_corner(self, i : int, j : int, k : int, c_str):
        ''' Get the corner position of a cell with indices i,j,k '''
        new_set = False
        if c_str != '':
            new_set = self._set_corner_set(c_str)

        if new_set == False and self._corners is not None:
            return self._corners[i,j,k]
        else:
            out_pos = np.array([0.0,0.0,0.0])
            idxs = (int(i),int(j),int(k))
            for idx, key in enumerate(['x_0', 'x_1', 'x_2']):
                ax = getattr(self, key)
                if self._corner_set[idx] == 'u':
                    # Set using the upper edge
                    out_pos[idx] = ax.ax_edges[idxs[idx]+1]
                else:
                    # Set using the lower edge
                    out_pos[idx] = ax.ax_edges[idxs[idx]]

            return out_pos

    def init_distances(self, pos, cell_position="center"):
        """
        Get distances of all cells from a particular point in space.
        By default, the cell's position is its center point. Otherwise,
        It is a corner.

        Provides distances as cartesian vectors of dimension order x,y,z
        """
        # Get reference to proper array
        if cell_position == "center":
            self.get_cell_centers()
            p_arr = self._centers
        else:
            self.get_corners(cell_position)
            p_arr = self._corners

        # check coordinate system
        if self.coord_syst == "cylindrical":
            # Convert to cartesian for easier distance calculation
            # Use dimension ordering of x,y,z
            r_index = self.dim_labels.index('r')
            t_index = self.dim_labels.index('theta')
            z_index = self.dim_labels.index('z')

            t_arr = p_arr.copy()
            t_arr[:,:,:,0] = p_arr[:,:,:,r_index] * np.cos(p_arr[:,:,:,t_index])
            t_arr[:,:,:,1] = p_arr[:,:,:,r_index] * np.sin(p_arr[:,:,:,t_index])
            t_arr[:,:,:,2] = p_arr[:,:,:,z_index]
            p_arr = t_arr

            # Do the same for the position
            pos.set_coordinate_system("cartesian", dim_order=('x', 'y', 'z'),
                                      units = self.units)
        else:
            pos.reorder_dims(('x','y','z'))

        # Get unitless values
        pos = pos.get_position_unitless_values()

        # Set up the distance data members
        if self.dist_mag is None:
            self.dist_mag = np.zeros(self.matrix_shape, dtype=np.float64)
            self.dist_vec = np.zeros(self.matrix_shape + (3,),
                                     dtype=np.float64)

        # Get the distances
        self.dist_vec[:,:,:,0] -= p_arr[:,:,:,0] - pos[0]
        self.dist_vec[:,:,:,1] -= p_arr[:,:,:,1] - pos[1]
        self.dist_vec[:,:,:,2] -= p_arr[:,:,:,2] - pos[2]
        self.dist_vec *= -1 # Change to point from position to source

        # Get the magnitudes
        self.dist_mag = np.linalg.norm(self.dist_vec, axis=3)

        # Normalize the vector magnitude
        self.dist_vec /= self.dist_mag[:,:,:,np.newaxis]

        return

    def get_info_string(self):
        out_s = f"{' Volume Information ':-^50}\n"
        out_s += f"\tCoorindate System:\t{self.coord_syst}\n"
        out_s += f"\tDimension order:\t{self.x_0.dim_type}"
        out_s += f" {self.x_1.dim_type} {self.x_2.dim_type}\n"
        out_s += f"\tNumber of cells:\t{self.n_cells}\n"
        out_s += f"\tGrid size:\t\t{self.matrix_shape}\n"
        out_s += f"{' Simulated Axis 0 ':-^50}\n"
        out_s += self.x_0.get_info_string()
        out_s += f"{' Simulated Axis 1 ':-^50}\n"
        out_s += self.x_1.get_info_string()
        out_s += f"{' Simulated Axis 2 ':-^50}\n"
        out_s += self.x_2.get_info_string()
        out_s += f"{' Simulated Time Axis ':-^50}\n"
        out_s += self.x_t.get_info_string()

        return out_s

class SimulationAxis:
    def __init__(self, dim_type, ax_min = None,
                 ax_max = None, nbins = None, **kwargs):
        self.dim_type = dim_type.lower()
        self.ax_min = 0
        self.ax_max = 0
        self._unit_scale = 1
        self.bin_width = None
        self.nbins = None
        self.units = None

        check_dim_type(self.dim_type)
        self._set_pos_args(ax_min, ax_max, nbins, **kwargs)
        self._set_spacing(**kwargs)
        self._setup_axis_edges()

    def _set_pos_args(self, ax_min, ax_max, nbins, **kwargs):
        self.units = UNITLESS
        if 'units' in kwargs.keys():
            self.units = kwargs['units']
        
        # Set min
        if ax_min is None:
            if 'min' in kwargs.keys():
                ax_min = float(kwargs['min'])
            else:
                ax_min = 0

        # Get the physical value
        p_val = PhysicalValue(ax_min, self.units)

        # Verify the units makes sense for dimension type
        self.units = p_val.get_units()
        self._unit_scale = p_val.get_unit_scale()
        verify_unit_for_dim_type(self.dim_type, self.units)

        self.ax_min = p_val.get_unitless_value()

        # Set bins
        self.nbins = nbins
        if self.nbins is None:
            if 'num_bins' in kwargs.keys():
                self.nbins = kwargs['num_bins']
            else:
                raise RuntimeError('SimulationAxis required the nbins argument.')
        self.nbins = int(self.nbins)

        # if bin width is provided with nbins, the ax_max is ignored
        # and we have all information we need
        if 'bin_width' in kwargs.keys():
            self.bin_width = float(kwargs['bin_width'])

        # Set max
        self.ax_max = ax_max
        if self.ax_max is None:
            if 'max' in kwargs.keys():
                self.ax_max = float(kwargs['max'])
            elif self.bin_width is not None:
                self.ax_max = self.nbins * self.bin_width

                # Put in right units
                p_val.set_value(self.bin_width)
                self.bin_width = p_val.get_unitless_value()
            else:
                raise RuntimeError('SimulationAxis required the ax_max argument when '
                                   'no bin_width is provided.')

        p_val.set_value(self.ax_max)
        self.ax_max = p_val.get_unitless_value()

        return

    def _set_spacing(self, **kwargs):
        self.spacing = 'linear'
        if 'spacing' in kwargs.keys():
            self.spacing = kwargs['spacing']

        if self.spacing != 'linear' and self.spacing[0:3] != 'log':
            raise ValueError(f"Invalid value for axis spacing: {self.spacing}")

        if self.spacing != 'linear':
            # Get the log base
            self._spacing_base = self.spacing[3:]
            if self._spacing_base not in ['', '2', '10']:
                raise ValueError("Invalid logarithmic base for spacing")
            else:
                if self._spacing_base == '':
                    self._spacing_base = np.exp()
                else:
                    self._spacing_base = float(self._spacing_base)

            # Setup axis first upper edge
            if self.ax_min == 0:
                if 'first_upper_edge' not in kwargs.keys():
                    raise RuntimeError("The 'first_upper_edge' element is needed for axis "
                                       "when spacing is non-linear and the axis min is 0.")
                self._first_uedge = float(kwargs['first_upper_edge']) * self._unit_scale

        

    def _setup_axis_edges(self):
        if self.spacing == 'linear':
            self.ax_edges = np.linspace(self.ax_min, self.ax_max, self.nbins+1)
        else:
            log_func = np.log
            if self._spacing_base == 2.:
                log_func = np.log2
            elif self._spacing_base == 10.:
                log_func = np.log10

            if self._first_uedge is not None:
                self.ax_edges = np.concatenate((np.array([0]),
                                                np.logspace(log_func(self._first_uedge),
                                                            log_func(self.ax_max),
                                                            self.nbins, base=self._spacing_base)))
            else:
                self.ax_edges = np.logspace(log_func(self.ax_min),
                                            log_func(self.ax_max),
                                            self.nbins+1, base=self._spacing_base)
        return

    def get_value_bin(self, val : float):
        arr = np.argwhere(self.ax_edges > val)
        if arr.shape[0] == 0:
            return -1
        elif arr[0,0] == 0:
            return -1

        return arr[0,0] - 1

    def get_info_string(self):
        out_s = f"\tDimension name:\t{self.dim_type}\n"
        out_s += f"\tNumber of bins:\t{self.nbins}\n"
        out_s += f"\tAxis maximum:\t{self.ax_max}\n"
        out_s += f"\tAxis minimum:\t{self.ax_min}\n"
        out_s += f"\tAxis spacing:\t{self.spacing}"
        if self.spacing != 'linear' and self._first_uedge is not None:
            out_s += f"\t(first upper edge: {self._first_uedge})\n"
        else:
            out_s += "\n"
        out_s +="\n"
        return out_s

    def get_bin_centers(self):
        return (self.ax_edges[1:] + self.ax_edges[:-1]) / 2

class SimulationPosition:
    def __init__(self, x_0, x_1, x_2, coord_syst = 'cartesian', units = None,
                 dimension_labels = None):
        self.x_0 = None
        self.x_1 = None
        self.x_2 = None

        # Check the coordinate system
        set_coord_syst_member(self, coord_syst)

        # Set the position physical values
        if units is None:
            if self.coord_syst == 'cartesian':
                units = ('m', 'm', 'm')
            else:
                units = ('m', 'rad', 'm')
        elif isinstance(units, str):
            if self.coord_syst != 'cartesian':
                raise ValueError("Units cannot be provided as a string when not in cartesian coordinates")
            units = (units, units, units)
        elif isinstance(units, tuple):
            if len(units) != 3:
                raise ValueError("Incomplete set of units as a tuple")
        else:
            raise TypeError("Invalid way of providing units.")

        # Now check dimension labels
        if dimension_labels is None:
            dimension_labels = ('r', 'theta', 'z') if self.coord_syst == "cylindrical" else ('x', 'y', 'z')
        elif isinstance(dimension_labels, str):
            if self.coord_syst != 'cartesian':
                raise ValueError("Units cannot be provided as a string when not in cartesian coordinates")
            dimension_labels = (dimension_labels, dimension_labels,
                                dimension_labels)
        elif isinstance(dimension_labels, tuple):
            if len(units) != 3:
                raise ValueError("Incomplete set of units as a tuple")
        else:
            raise TypeError("Invalid way of providing units.")

        # Verify the units are right for the dimension type
        self.dim_labels = dimension_labels
        self.units = units
        for i in range(3):
            check_dim_type(self.dim_labels[i])
            verify_unit_for_dim_type(self.dim_labels[i],
                                     self.units[i])

        # Check that the positions are physical values
        x_in = (x_0, x_1, x_2)
        for i, key in enumerate(['x_0', 'x_1', 'x_2']):
            if isinstance(x_in, PhysicalValue):
                x_i = x_in[i]
                x_i.set_units(units[i])
            else:
                x_i = PhysicalValue(x_in[i], units[i])
            setattr(self, key, x_i)

        return

    def get_string_position(self):
        out_str = ""
        for i, key in enumerate(['x_0', 'x_1', 'x_2']):
            x_i = getattr(self, key)
            x_i_s = x_i.get_string_value()
            out_str += f"\t{self.dim_labels[i]}:\t{x_i_s} \n"
        return out_str

    def get_position_values(self):
        return (self.x_0.get_value(), self.x_1.get_value(),
                self.x_2.get_value())

    def get_position_unitless_values(self):
        return (self.x_0.get_unitless_value(), self.x_1.get_unitless_value(),
                self.x_2.get_unitless_value())

    def get_position(self):
        return (x_0, x_1, x_2)

    def get_distance(self, pos_in):
        convert_coord_syst = True if self.coord_syst == 'cylindrical' else False
        orig_dim_order = self.dim_labels

        if convert_coord_syst:
            # Convert to cartesian
            self.set_coordinate_system("cartesian")

        arr = np.array([self.x_0.get_unitless_value(), self.x_1.get_unitless_value(),
                        self.x_2.get_unitless_value()])

        dist = np.linalg.norm(arr)

        if convert_coord_syst:
            self.set_coordinate_system("cylindrical", dim_order = orig_dim_order)

        return dist


    def get_dim_order(self):
        return self.dim_labels

    def set_coordinate_system(self, coord_syst, units = None,
                              dim_order = None):
        verify_coord_system(coord_syst)

        # Verify the dimension order tuple
        if dim_order is not None:
            if not isinstance(dim_order, tuple):
                raise TypeError("Dimension order parameter must be a tuple")
            if len(dim_order) != 3:
                raise ValueError("Length of dimension order tuple must be 3.")

        # Same for units
        if units is not None:
            if not isinstance(units, tuple):
                raise TypeError("Dimension order parameter must be a tuple")
            if len(units) != 3:
                raise ValueError("Length of dimension order tuple must be 3.")        

        # Perform conversions
        if coord_syst == "cylindrical":
            self._set_cylindrical(units, dim_order)
        else:
            self._set_cartesian(units, dim_order)
        return

    def _set_cylindrical(self, units = None,
                         dim_order = None):
        if self.coord_syst == "cylindrical":
            # nothing to do
            return

        self.coord_syst = "cylindrical"

        # Map the vector elements to coordinates
        pos = (self.x_0, self.x_1, self.x_2)
        x = pos[self.dim_labels.index('x')].get_unitless_value()
        y = pos[self.dim_labels.index('y')].get_unitless_value()
        z = pos[self.dim_labels.index('z')].get_unitless_value()
        
        # Perform conversion
        r, theta, z = convert_cart_to_cyl(x, y, z)

        # Convert units and set a physical values
        if units is None:
            units = ('m', 'rad', 'm')

        self.x_0 = PhysicalValue(r, units[0])
        self.x_1 = PhysicalValue(theta, units[1])
        self.x_2 = PhysicalValue(z, units[2])
        self.dim_labels = ('r', 'theta', 'z')

        # Put coordinates in correct vector order
        dim_order = self.dim_labels if dim_order is None else dim_order
        if dim_order != self.dim_labels:
            self.reorder_dims(dim_order)

        return

    def _set_cartesian(self, units = None, dim_order = None):
        if self.coord_syst == "cartesian":
            # nothing to do
            return

        self.coord_syst = "cartesian"

        # Map the vector elements to coordinates
        pos = (self.x_0, self.x_1, self.x_2)
        r = pos[self.dim_labels.index('r')].get_unitless_value()

        theta = pos[self.dim_labels.index('theta')]
        theta.set_units('rad') # Set units to radians for conversion
        theta = theta.get_unitless_value()

        z = pos[self.dim_labels.index('z')].get_unitless_value()

        # Perform conversion
        x, y, z = convert_cyl_to_cart(r, theta, z)

        # Convert units and set a physical values
        if units is None:
            units = ('m', 'm', 'm')

        self.x_0 = PhysicalValue(x, units[0])
        self.x_1 = PhysicalValue(y, units[1])
        self.x_2 = PhysicalValue(z, units[2])
        self.dim_labels = ('x', 'y', 'z')

        # Put coordinates in correct vector order
        dim_order = self.dim_labels if dim_order is None else dim_order
        if self.dim_labels != dim_order:
            self.reorder_dims(dim_order)

        return

    def reorder_dims(self, dim_order : tuple):
        for i in dim_order:
            check_dim_type(i, self.coord_syst)

        pos = (self.x_0, self.x_1, self.x_2)

        self.x_0 = pos[self.dim_labels.index(dim_order[0])]
        self.x_1 = pos[self.dim_labels.index(dim_order[1])]
        self.x_2 = pos[self.dim_labels.index(dim_order[2])]
        self.dim_labels = dim_order

        return

    def set_unitless_values(self, val_0, val_1, val_2):
        self.x_0.set_unitless_value(val_0)
        self.x_1.set_unitless_value(val_1)
        self.x_2.set_unitless_value(val_2)
        return


class PhysicalValue:
    def __init__(self, value, units):
        self._set_units(units)
        self.set_value(value)
        return

    def _set_units(self, units : str):
        if not isinstance(units, str):
            # Make sure its a string
            raise TypeError("Invalid type for units parameter")

        self.units = units.lower()
        self._unit_scale = 1.0
        self._units = self.units

        # Nothing more to do for unitless
        if self.units == '':
            return

        if self.units not in ['deg', 'rad']:
            unit_pref = None
            if len(self.units) > 1:
                unit_pref = self.units[0]
                self._units = self.units[1:]

            # get the base unit
            if self._units not in ['s', 'm']:
                raise ValueError(f"Unsupported unit: {self._units}")

            # get the unit scale
            if unit_pref is not None:
                if unit_pref not in UNIT_PREF_MAP.keys():
                    raise ValueError(f"Invalid unit scale ({kwargs['units']}) for value.")
                else:
                    self._unit_scale = UNIT_PREF_MAP[unit_pref]

        return

    def set_value(self, value : float):
        try:
            # Make sure the type is right
            value = float(value)
        except:
            raise TypeError("Invalid type for value parameter")

        self._value = value
        self._unitless_value = self._value * self._unit_scale

    def set_unitless_value(self, value : float):
        try:
            # Make sure the type is right
            value = float(value)
        except:
            raise TypeError("Invalid type for value parameter")

        self._unitless_value = value
        self._value = self._unitless_value / self._unit_scale


    def _converted_from(self, old_units):
        i_old = CONVERTABLE_UNIT_VEC.index(old_units)
        i_new = CONVERTABLE_UNIT_VEC.index(self._units)
        conversion_factor = CONVERSION_TABLE[i_old][i_new]
        if conversion_factor == -1:
            raise RuntimeError(f"Invalid conversion from {old_units} to {self._units}")
        elif conversion_factor != 1:
            self._unitless_value *= conversion_factor

        return

    def set_units(self, units):
        if units == self.units:
            # Nothing to do
            return

        # Store the old value
        prev_units = self._units
        prev_unit_scale = self._unit_scale

        # Check and set units
        self._set_units(units)
        if prev_units != self._units:
            self._converted_from(prev_units)
            self._value = self._unitless_value / self._unit_scale
        elif prev_unit_scale != self._unit_scale:
            # Adjust for a scale change
            self._value = self._unitless_value / self._unit_scale

    def get_value(self):
        return self._value

    def get_unitless_value(self):
        return self._unitless_value

    def get_units(self):
        return self._units

    def get_unit_scale(self):
        return self._unit_scale

    def get_string_value(self):
        return f"{self._value} {self.units}"
