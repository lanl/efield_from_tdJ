import time
import numpy as np
import numpy.ma as ma
import physical_constants
import math
import sys
import calculation_globals as cg

inited_d = False

J_vec = None
dJ_dt_vec = None

def get_distance_differential_array_slice(origin, n_1, n_2, n_3):
    '''
    Get simulation volumes different volume elements and cell distances
    by slicing a number of bins (positive) from an origin

    Parameters
    ----------
    origin : np.array
        A element array giving the cell index of the lower corner of the
        sliced volume
    n_1 : int
        Number of bins along first axis to include in slice (origin[0], origin[0] + n_1)
    n_2 : int
        Number of bins along second axis to include in slice (origin[1], origin[1] + n_2)
    n_3 : int
        Number of bins along first axis to include in slice (origin[2], origin[2] + n_3)

    Returns
    -------
    mags : np.array
        The selected slice of the distance magnitude array
    dists : np.array
        The selected slice of the distance vectors
    dV : np.array
        The  selected slice of the differential volume array
    '''
    if cg.sim_geometry.dist_vec is None:
        raise RuntimeError('Distance array needs to be initialized')

    # Get the span
    span = np.array([[0,n_1],[0,n_2],[0,n_3]])
    span[:,0] += origin
    span[:,1] += origin

    sliced = (slice(span[0,0],span[0,1]),slice(span[1,0],span[1,1]),
              slice(span[2,0],span[2,1]))

    return np.copy(cg.sim_geometry.dist_mag[sliced]), np.copy(cg.sim_geometry.dist_vec[sliced]), \
        np.copy(cg.sim_geometry.dV[sliced])


def get_time_delays(r, t_offset_idx, dt, t, a):
    '''
    Get the number of time steps needed for a signal to travel provided distances

    Index is rounded down using np.floor

    Paramters
    ---------
    r : np.array
        A flat array of distance magnitudes over which the signals with propagate
    t_offset_idx : int
        Number of steps to arbitrarily offset from the result (negative)
    t : int
        The time step at the start of signal propagation
    dt : float
        The time width of a time step in seconds
    a : int
        Number of steps to arbitrarily offset from the result (positive)

    Results
    -------
    t_idx : np.array
        Array of same length of r giving the absolute retarded  time indices 
    '''
    return np.floor((r / physical_constants.c_m_p_s + t) / dt) - t_offset_idx + a

def compute_J_total(Jx, Jy, Jz):
    '''
    Compute the total current density magnitudes from three components

    All component arrays must have the same shape.

    Parameters
    ----------
    Jx : np.array
        First component of the current density
    Jy : np.array
        Second component of the current density
    Jz : np.array
        Third component of the current density

    Return
    ------
    J_tot : ma.masked_array
        Total magnitude of current density. Values of 0 are masked out.
    '''
    J_tot = np.sqrt(Jx*Jx + Jy*Jy + Jz*Jz)

    # If the total J == 0, we mask it
    return ma.masked_array(J_tot, mask = (J_tot == 0))


def set_3d_J_vec(Jx, Jy, Jz, J_vec):
    '''
    Use three component arrays to set a vector array
    '''

    # get vectorized shape
    J_vec[:,:,:,0] = Jx
    J_vec[:,:,:,1] = Jy
    J_vec[:,:,:,2] = Jz

def dot_vector(a_arr, b_arr):
    '''
    Perform doc product on two arrays of vectors
    '''
    return np.sum(a_arr*b_arr, axis=1)

def calculate_e_field(Jx, Jy, Jz, dJx_dt, dJy_dt, dJz_dt, t,
                      E_sta, E_ind, E_rad, **kwargs):
    '''
    Calculate the electric field using current density dependent
    electric field Jefimenko equation. E-field is broken into
    static, induction, and radiation terms. This is for a single
    time step
    
    Parameters
    ----------
    Jx : np.array
        First component of current density
    Jy : np.array
        Second component of current density
    Jz : np.array
        Third component of current density
    dJx_dt : np.array
        First component of current density time derivative
    dJy_dt : np.array
        Second component of current density time derivative
    dJz_dt : np.array
        Third component of current density time derivative
    t : int
        Current time index
    E_sta : np.array
        Output static electric field
    E_ind : np.array
        Output induction electric field
    E_rad : np.array
        Output radiation electric field
    '''
    global J_vec
    global dJ_dt_vec

    # Get keyword args
    off = kwargs['arr_offset'] # Offset of origin of the mesh

    # Offset of time axis origin
    t_offset_idx = int(kwargs['t_offset'] / cg.sim_dt)

    # Get shapes for J + dJ/dt 4 vectors
    n_1 = Jx.shape[0]
    n_2 = Jx.shape[1]
    n_3 = Jx.shape[2]
    shape = (n_1, n_2, n_3)

    #Get the 3 vectors of the J and the partial time derivative.
    shape_3d = (shape[0], shape[1], shape[2], 3)

    if J_vec is None or shape_3d != J_vec.shape:
        J_vec = np.zeros(shape_3d)
        dJ_dt_vec = np.zeros(shape_3d)

    set_3d_J_vec(Jx, Jy, Jz, J_vec)
    set_3d_J_vec(dJx_dt, dJy_dt, dJz_dt, dJ_dt_vec)

    # Get the matrix of distances from the probe position
    r, r_vec, dV = get_distance_differential_array_slice(off, n_1, n_2, n_3)

    # Get the mask (i.e., the cells that have non-zero J at a given time)
    # This is useful in ignoring the cells that don't contribute
    J_tot = compute_J_total(Jx, Jy, Jz)
    dJ_tot = compute_J_total(dJx_dt, dJy_dt, dJz_dt) # Don't forgot cells with non-zero dJ but zero J
    mask = ~J_tot.mask | ~dJ_tot.mask

    # Get the volume of each cell contributing to the E-field
    dv_flat = dV[mask][:]

    # Get the J and J_dot and calculate field
    J_vec_0 = J_vec[mask][:]
    dJ_dt_vec_0 = dJ_dt_vec[mask][:]
    r_vec_0 = r_vec[mask][:]
    r_0 = r[mask][:,np.newaxis]

    # Get the time delays
    t_idx = get_time_delays(r_0, t_offset_idx, cg.sim_dt,
                            (float(t) + 0.5) * cg.sim_dt, 10)
    t_idx = t_idx.astype(int)

    # Determine E-field for position and time
    e_field_at_pos_dt(J_vec_0, dJ_dt_vec_0, r_vec_0, r_0, t_idx, dv_flat,
                      E_sta, E_ind, E_rad)


def e_field_at_pos_dt(J_vec, dJ_dt_vec, r_vec, r, t_idx, dv, E_sta_p,
                      E_ind_p, E_rad_p):
    '''
    Get the eletric field from a single time step at a location

    The J_vec, r_vec, r arrays are flattened arrays where each entry
    is a contribution from a point in space to the E-field at this
    location and time. So is dJ_dt_vec.

    t_idx gives the retarded time index for each entry in the J and r arrays

    dv - the differential volume and has one entry per contributing cell

    E_sta_p is the static electric field at the position
    E_ind_p is the induction electric field at the position
    E_rad_p is the radiation electric field at the position
    '''
    # Get vector products
    JdR = dot_vector(J_vec, r_vec)[:,np.newaxis]
    JxR = np.cross(J_vec, r_vec)
    JxRxR = np.cross(JxR, r_vec)

    # Get time partial derivative products
    dJ_dtxR = np.cross(dJ_dt_vec, r_vec)
    dJ_dtxRxR = np.cross(dJ_dtxR, r_vec) 

    scal = (physical_constants.four_pi_eps / dv)[:,np.newaxis]

    # Calculate each term
    E_sta_i = (2*JdR*r_vec + JxRxR) / (r**3) / scal
    E_ind_i = (2*JdR*r_vec + JxRxR) / (physical_constants.c_m_p_s * (r**2)) \
              / scal
    E_rad_i = dJ_dtxRxR / ((physical_constants.c_m_p_s**2) * r) \
              / scal

    for i in range(3):
        # For each component, add the E-filed contribution of each type to the output
        # e-field time segment
        np.add.at(E_sta_p[:,i], t_idx.flatten(), E_sta_i[:,i])
        np.add.at(E_ind_p[:,i], t_idx.flatten(), E_ind_i[:,i])
        np.add.at(E_rad_p[:,i], t_idx.flatten(), E_rad_i[:,i])

    return
