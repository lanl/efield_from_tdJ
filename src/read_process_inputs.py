import numpy as np
from physical_constants import c_m_p_s as cmps
import calculation_globals as cg
import struct
import sys

def fill_data_for_time(idx, Jx, Jy, Jz, shape, f):
    '''
    Load all data from the unput file for the current time step.

    Parameters
    ----------
    idx : int
        Current time step index
    Jx : np.array
        First component of the current density array for the time step
    Jy : np.array
        Second component of the current density array for the time step
    Jz : np.array
        Third component of the current density array for the time step
    shape : tuple
        Shape of the current density arrays
    f : File
        Opened file handler for the input CSV file
    '''
    
    # Set up the shape and array origin
    idx_l = np.array([idx,0,0,0], dtype=int)

    min_maxs = np.array([[0,shape[0]-1],[0,shape[1]-1],[shape[2],2]])  

    data = True
    while True:
        pos = f.tell()
        line = f.readline()
        if len(line) <= 1:
            idx_l[0] = -1
            data = False
            break

        # get the 4-space indices
        toks = line.split(',')
        idx_l[0] = int(toks[0])
        if idx_l[0] != idx:
            # Rewind a line
            f.seek(pos)
            data = True
            break

        idx_l[1:] = np.array([int(toks[1]), int(toks[2]),
                              int(toks[3])])

        # Report the current density components                                                       
        Jx[idx_l[1], idx_l[2], idx_l[3]] = float(toks[4])
        Jy[idx_l[1], idx_l[2], idx_l[3]] = float(toks[5])
        Jz[idx_l[1], idx_l[2], idx_l[3]] = float(toks[6])

        
        if idx_l[3] < min_maxs[2,0]:
            min_maxs[2,0] = idx_l[3]
        if idx_l[3] > min_maxs[2,1]:
            min_maxs[2,1] = idx_l[3]

    return idx_l[0], min_maxs, data
    

def get_union_span(span_n, span_l):
    """
    Get a slice tuple to represent the union of two array spans
    """
    span = np.copy(span_n)
    # take the union of the spans
    for i in range(3):
        if span_l[i,0] < span[i,0]:
            span[i,0] = span_l[i,0]
        if span_l[i,1] > span[i,1]:
            span[i,1] = span_l[i,1]

    return (slice(span[0,0],span[0,1]+1),slice(span[1,0],span[1,1]+1),
            slice(span[2,0],span[2,1]+1))
    

def verify_csv_file_header(in_file):
    """
    Verify this is a valid file by checking the header and its length
    """
    hdr = in_file.readline()
    toks = hdr.split(',')
    if len(toks) != 7: # 1 time index + 3 space indices + 3 J components
        print(f'--- Invalid number ({len(toks)}) of elements in file header.')
        return False
    
    sort_toks = sorted(set([tok.strip() for tok in toks]))
    exp_hdr = sorted(set(['TIME INDEX', 'X_0 INDEX', 'X_1 INDEX', 'X_2 INDEX',
                          'J_0 COMPONENT', 'J_1 COMPONENT', 'J_2 COMPONENT']))
    if exp_hdr != sort_toks:
        print('--- Header does not reflect the expected input.')
        print(f'---   Expected: {exp_hdr}')
        print(f'---   Got       {sort_toks}')
        return False

    return True

def read_and_analyze_text_file(i_filename, E_sta, E_ind, E_rad,
                               funct, **kwargs):
    '''
    Read the input file and analyze the current densities to get E-field
    using the provided function. kwargs are passed to the provided function

    Parameters
    ----------
    i_filename : str
        Input file name and path
    E_sta : np.array
        Output static electric field result
    E_ind : np.array
        Output induction electric field result
    E_rad : np.array
        Output radiation electric field result    
    funct : function
        Function to use for computing e-field
    kwargs : Key word arguments dict
    '''

    # Open input and verify input file
    f = open(i_filename, 'r')
    if verify_csv_file_header(f) is False:
        print(f'--- Bad input file: {i_filename}')
        return False

    # Initialize the current density (J) arrays and the time derivatives
    max_shape = cg.sim_geometry.matrix_shape
    Jx1 = np.zeros(max_shape, dtype=float)
    Jy1 = np.copy(Jx1)
    Jz1 = np.copy(Jx1)
    Jx2 = np.copy(Jx1)
    Jy2 = np.copy(Jx1)
    Jz2 = np.copy(Jx1)

    cg.sim_geometry.init_distances(kwargs['position'])

    idx_t = 0
    prev_ti = 0

    # open the input
    data = True
    d_idx = 0
    i = 0
    span = np.array([[0,0],[0,0],[0,0]])
    span_l = np.copy(span)
    sl = get_union_span(span, span_l)

    # Get the time differential element
    dt = float(cg.sim_geometry.x_t.bin_width)

    while data:
        # Reset the J arrays
        Jx1[sl] = Jx2[sl]
        Jy1[sl] = Jy2[sl]
        Jz1[sl] = Jz2[sl]
        Jx2[sl] = 0
        Jy2[sl] = 0
        Jz2[sl] = 0

        # Get all data for the current time
        prev_ti = idx_t
        idx_t, span, data = fill_data_for_time(idx_t, Jx2, Jy2, Jz2, max_shape, f)

        # Get the span of this data
        if prev_ti == 0:
            span_l = np.copy(span)

        # The span includes the spatial extention of this timesteps data and that
        # of the prior timestep
        sl = get_union_span(span, span_l)
        span_l = span

        if idx_t == -1:
            idx_t = prev_ti + 1

        # Get the time partial derivatives
        if idx_t - prev_ti == 1:
            dJx_dt = (Jx2[sl] - Jx1[sl]) / dt
            dJy_dt = (Jy2[sl] - Jy1[sl]) / dt
            dJz_dt = (Jz2[sl] - Jz1[sl]) / dt
        else:
            dJx_dt = (Jx2[sl]) / dt
            dJy_dt = (Jy2[sl]) / dt
            dJz_dt = (Jz2[sl]) / dt

        # Evoke the analysis function
        funct(Jx2[sl], Jy2[sl], Jz2[sl], dJx_dt, dJy_dt, dJz_dt,
              idx_t, E_sta, E_ind, E_rad, dt=dt,
              arr_offset=span[:,0],
              **kwargs)

    # Perform the time integration
    for l in range(E_sta.shape[0]-1):
        E_sta[l+1,:] += E_sta[l,:]
    E_sta *= dt

    return
