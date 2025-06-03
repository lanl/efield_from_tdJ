import struct
import sys
import configparser
from collections import defaultdict
from simulated_geometry import SimulationVolume, PhysicalValue

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


gconfig = config_ini_parse("../config/sample_config.ini")
sim_geometry = SimulationVolume(gconfig)

# Usage python bin_to_csv.py INPUT_BIN_FILE OUTPUT_CSV_FILE
# Converts the binary of format:
#    Byte 0: float (4 bytes) time
#    Byte 4: float (4 bytes) x_0
#    Byte 8: float (4 bytes) x_1
#    Byte 12: float (4 bytes) x_2
#    Byte 16: float (4 bytes) J_0
#    Byte 20: float (4 bytes) J_1
#    Byte 24: float (4 bytes) J_2
#    Length
# to a CSV in that column order

ifile = sys.argv[1]
ofile = sys.argv[2]

with open(ifile, 'rb') as f_in:
    with open(ofile, 'w') as f_out:
        # Write the header
        f_out.write("TIME INDEX,X_0 INDEX,X_1 INDEX,X_2 INDEX,"
                    "J_0 COMPONENT,J_1 COMPONENT, J_2 COMPONENT\n")

        while True:
            bin_dat = f_in.read(28)
            if len(bin_dat) != 28:
                break

            ascii_dat = list(struct.unpack("fffffff", bytearray(bin_dat)))
            ascii_dat[0] = sim_geometry.x_t.get_value_bin(ascii_dat[0])
            ascii_dat[1] = sim_geometry.x_0.get_value_bin(ascii_dat[1])
            theta = PhysicalValue(ascii_dat[2],'rad')
            theta.set_units('deg')
            ascii_dat[2] = sim_geometry.x_1.get_value_bin(theta.get_unitless_value())
            ascii_dat[3] = sim_geometry.x_2.get_value_bin(ascii_dat[3])
            if -1 in ascii_dat:
                print("Problem")

            f_out.write(",".join(str(elem) for elem in ascii_dat))
            f_out.write('\n')            
