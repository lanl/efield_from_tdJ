from read_process_inputs import read_and_analyze_text_file
from simulated_geometry import SimulationVolume, SimulationAxis, SimulationPosition
from write_outputs import write_ascii_output
from efield_math import calculate_e_field
import physical_constants
import calculation_globals as cg

import numpy as np
import sys
import os
import argparse
import configparser
from collections import defaultdict

def config_ini_parse(config_filepath : str):
    """
    Parse an ini config file.

    Parameters
    ----------
    config_filepath : string
        Filepath of the config ini file to parse

    Returns
    -------
    config : dict
        Parsed configuration as a dictionary
    """
    cp = configparser.ConfigParser()
    cp.read(config_filepath)

    try:
        success = cp.read(config_filepath)
        if not success:
            raise configparser.Error
    except Exception:
        raise IOError

    # full dictionary is nested with an entry for each section
    config = defaultdict(dict)
    for sect in cp.sections():
        config[sect] = dict(cp.items(sect))

    return config

def verify_arguments(args):
    """
    Verify that supplied arguments are valid

    Checks that the right options are supplied based on is_index.
    Checks that the geoemetry configuration can be parsed
    Checks that paths are valid
    """

    if args.is_index:
        if args.indices is None:
            raise ValueError("If 'is_index' option is used, 'indices' option "
                             "is required.")
            
        # Make sure indices can be interpretted as floats
        try:
            args.indices = [float(int(i)) for i in args.indices]
        except ValueError:
            raise ValueError(f"An index could not be interpretted as int: {args.indices}")

    else:
        if args.coords is None:
            raise ValueError("The 'coordinates' option is required "
                             "(unless 'is_index' is used)")

        # Make sure coordiantes can be interpretted as floats
        try:
            args.coords = [float(i) for i in args.coords]
        except ValueError:
            raise ValueError(f"A coordinate could not be interpretted as float: {args.coords}")
        

    # Check paths
    args_dict = vars(args)
    for path_key in ['o_dir', 'ifile', 'geometry_conf']:
        fn_path = args_dict[path_key]
        if path_key == 'o_dir':
            fn_path = os.path.split(fn_path)[0]
            if fn_path == '':
                # Nothing to do if current working directory
                continue

        if not os.path.exists(fn_path):
            raise ValueError(f"Invalid path for {path_key}: {fn_path}")

    # Check that the geometry can be parsed
    try:
        args.geometry = config_ini_parse(args.geometry_conf)
    except IOError:
        sys.exit(f"Error parsing the geometry calculation ({args.geometry_conf}).")

    return args

def parse_cmd_line(argv):
    """
    Parse the command-line interface

    Parameters
    ----------
    argv : list
        Command line arguments (probably supplied through sys.argv)

    Returns
    -------
    """
    desc = ''
    parser = argparse.ArgumentParser(description=desc)

    # Add the coordinate system parameter
    parser.add_argument('-s','--probe_coord_syst', dest='coord_syst',
                        metavar='COORD_SYST', default='cartesian',
                        choices=['cartesian', 'cylindrical'],
                        help=('Coordinate system of the provided probe position'
                              'Options: "cartesian" or "cylindrical"'
                              'Default: "cartesian"'))

    # Add the coordinate argument
    parser.add_argument('-c', '--coordinates', dest='coords',
                        nargs=3, metavar=('X_0', 'X_1', 'X_2'),
                        default=None,
                        help=('Coordinates of the probe position in meters from '
                              'the origin. The coordinate system determines how '
                              'these positional arguments are interpetted. If '
                              'COORD_SYST="cartesian", X_0 = X, X_1 = Y, and '
                              'X_2 = Z. If COORD_SYST="cylindrial", X_0 = radius (meters), '
                              'X_1 = theta (degrees), and X_2 = z (meters).'))

    # Flag to indicate the probe position is given as array indices
    parser.add_argument('--is_index', dest='is_index', action='store_true',
                        default=False,
                        help=('If is_index is set, the --indices argument is expected '
                              'for the probe position. Otherwise, the --coordinates '
                              'argument is expected.'))

    # Add arguments for the output time axis
    parser.add_argument('--no-time-offset', dest='out_time_off', action='store_false',
                        default=False,
                        help=('If this flag is set, the origin of the time axis always '
                              'corresponds with the origin of the input time axis. By default, '
                              'the time axis is shifted so the signal starts at the 10th '
                              'time point.'))

    # Add the indices argument
    parser.add_argument('-I', '--indices', dest='indices',
                        nargs=3, metavar=('I_0', 'I_1', 'I_2'),
                        default=None,
                        help=('Indices of the probe position in spatial array. '
                              'When indices are used, the probe is placed in the corner '
                              'of the gridded volume element with index (I_0, I_1, I_2). '
                              'The exact corner is the intersection of the three lower edges '
                              'of the orthogonal 3D grid. The index order corresponds with '
                              'the order of dimensions provided in the geometry config.'))

    # Add the IO options
    parser.add_argument('-o','--output', dest='o_dir',
                        metavar='ODIR', required=True,
                        help='Required full output directory path')
    parser.add_argument('-i','--input', dest='ifile',
                        metavar='IFILE', required=True,
                        help='Required full input filepath')
    parser.add_argument('-n', '--ntime_bins', dest='n_tO_bins', default = None,
                        help=('Number of time bins in the output time axis. The default is 2 times '
                              'the number of bins as the input time axis.'))

    # Add the input and output
    parser.add_argument('-g', '--geometry-config', dest='geometry_conf',
                        help=('Configuration file (INI) providing the simulated '
                              'volume and geometry'),
                        required=True)

    parser.add_argument('--verbose', dest='verbose', action='store_true',
                        default=False,
                        help=('If flag is set, higher verbosity is enabled and '
                              'debugging information is printed.'))

    args = parser.parse_args(argv)
    J_VERBOSE_GLOBAL = args.verbose

    # Check the arguments
    args = verify_arguments(args)

    return args

def verify_and_set_position(in_vec, is_index, coord_syst = None):
    pos = in_vec

    # Build properly ordered position
    if is_index == False:
        # Create a SimulationPosition
        if coord_syst == 'cartesian':
            sim_pos = SimulationPosition(in_vec[0], in_vec[1], in_vec[2])
        else:
            sim_pos = SimulationPosition(in_vec[0], in_vec[1], in_vec[2],
                                         coord_syst = 'cylindrical')

        # Convert the coordinate system if needed
        # This also ensures right ordering of dimensions
        sim_pos.set_coordinate_system(cg.sim_geometry.coord_syst,
                                      dim_order = cg.sim_geometry.get_axes_titles())

        pos_tmp = sim_pos.get_position_unitless_values()

        # Check that the position is not a cell center
        cell_centers = cg.sim_geometry.get_cell_centers()
        is_cell_center = None
        for i in range(cg.sim_geometry.x_0.nbins):
             for j in range(cg.sim_geometry.x_1.nbins):
                  for k in range(cg.sim_geometry.x_2.nbins):
                      if np.array_equal(pos_tmp, cell_centers[i]):
                          is_cell_center = (i,j,k)
                      break

        if is_cell_center is not None:
            # If it is a cell center, we move the location to the
            # corner former by the lower edge of each axis
            new_pos = cg.sim_geometry.get_cell_corner(
                is_cell_center[0], is_cell_center[1], is_cell_center[2],
                'lll')
            print(f"Warning position ({pos_tmp}) is a cell center. "
                  "To avoid singularity issues, position is being shifted to "
                  f"the cell center: {new_pos}")
            sim_pos.set_unitless_values(new_pos[0], new_pos[1], new_pos[2])
            
    else:
        # In the case of an index position, we check that the indices
        # valid
        if pos[0] >= cg.sim_geometry.x_0.nbins or pos[0] < 0 or \
           pos[1] >= cg.sim_geometry.x_1.nbins or pos[1] < 0 or \
           pos[2] >= cg.sim_geometry.x_2.nbins or pos[2] < 0:
            raise ValueError(f"Position by indices ({pos}) is out of range "
                             "for a celled volume of shape ({cg.sim_geometry.matrix_shape}).")

        pos_tmp = cg.sim_geometry.get_cell_corner(pos[0], pos[1], pos[2],
                                               'lll')
        # Convert back to a SimulationPosition
        if coord_syst == 'cartesian':
            sim_pos = SimulationPosition(pos_tmp[0], pos_tmp[1], pos_tmp[2])
        else:
            sim_pos = SimulationPosition(pos_tmp[0], pos_tmp[1], pos_tmp[2],
                                         coord_syst = 'cylindrical')

    return sim_pos
        

def main():

    # Parse the command line
    args = parse_cmd_line(sys.argv[1:])

    # Verify arguments
    verify_arguments(args)

    # Set up the geometry
    cg.init_sim_geometry(args.geometry)

    # Get the probe positions
    position = None
    if args.is_index:
        position = verify_and_set_position(args.indices, True)
    else:
        position = verify_and_set_position(args.coords, False,
                                           args.coord_syst)

    # Log probe positions to console
    print(80 * "*")
    print(f"{' Probe Position ':=^80}")
    p_name = " Probe #0 "
    print(f"{ ' Probe #0 ' :-^80}")
    print(position.get_string_position())
    print(80 * "-")
    print(80 * "=")

    # Log the simulated volume information
    print(f"{' Simulated Volume ':=^80}")
    print(cg.sim_geometry.get_info_string())
    print(80 * "=")

    # Get the output and input files
    o_dir = args.o_dir
    ifile = args.ifile

    # Correct the probe positions
    t_offset = 0 
    if not args.out_time_off:
        # -1 is to offset a meter so signal doesn't start before time origin
        t_offset = (position.get_distance(SimulationPosition(0,0,0)) - 0.5) \
                   / physical_constants.c_m_p_s

    # Set up the output time axis
    n_tO_bins = args.n_tO_bins
    if n_tO_bins is None:
        # By default, its the number of input time bins times 2
        n_tO_bins = cg.sim_geometry.x_t.nbins * 2
        
    tO_ax = SimulationAxis('t', num_bins=n_tO_bins,
                           units ='ns', bin_width=1)

    # Initialize the E-fields
    # Set up the output datasets
    shape_out = (n_tO_bins, 3)

    E_sta = np.zeros(shape_out)
    E_ind = np.zeros(shape_out)
    E_rad = np.zeros(shape_out)

    # Read the input while analyzing and calculating the E-field
    print(f"{' Reading Input ':=^80}")
    read_and_analyze_text_file(ifile, E_sta, E_ind, E_rad,
                               calculate_e_field,
                               position=position, t_offset=t_offset)
    print(f"{' Done ':=^80}")

    # Write to the output file
    print(f"{' Writing Output ':=^80}")
    write_ascii_output(o_dir, 'probe0', tO_ax, E_sta, E_ind, E_rad)
    print(f"{' Done ':=^80}")
    print(80 * "*")

if __name__ == "__main__":
    main()
