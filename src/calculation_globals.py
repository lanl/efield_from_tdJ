from simulated_geometry import SimulationVolume

'''
This file is for storing globally scoped variables to
be shared across python codes.

The most important variable is the sim_geometry SimulationVolume
object. This stores the axes, timing elements, and volume elements
used through the computation.
'''

# Globals
J_VERBOSE_GLOBAL = False

sim_geometry = None
sim_dt = 0

def init_sim_geometry(geo):
    global sim_geometry
    global sim_dt
    sim_geometry = SimulationVolume(geo)
    sim_dt = float(sim_geometry.x_t.bin_width)
