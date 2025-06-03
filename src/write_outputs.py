from datetime import datetime
import numpy as np
import os

def write_array_to_file(fn, ts, E_field, e_type):
    ''' Get an array for an electric field to a file'''
    column_wid = 30
    full_fn = fn.format(e_type)
    print(f'\tWriting {e_type} field result to: {full_fn}')

    # Create header format string
    fill_fs = f' <{column_wid}'
    header = f"{ 'Time (s)' :{fill_fs}}{ 'E_x (V/m)' :{fill_fs}}"
    header += f"{ 'E_y (V/m)' :{fill_fs}} { 'E_z (V/m)' :{fill_fs}}\n"

    # Create float format string
    sci_fs = f'{fill_fs}.9e'

    # Write output file
    with open(full_fn, 'w') as f:
        # Write the header
        f.write(header)
        for i in range(E_field.shape[0]):
            line = f"{ts[i]:{sci_fs}}{E_field[i,0]:{sci_fs}}"
            line += f"{E_field[i,1]:{sci_fs}}{E_field[i,2]:{sci_fs}}\n"
            f.write(line)

    return

def write_ascii_output(out_path : str, unique_ID : str,
                       time_ax,
                       E_sta : np.array,
                       E_ind : np.array, E_rad : np.array):
    '''
    Write the three electric fields to output files. The output
    is the directory the files are written to. The file names have the format:
       efield_{field_name}_{unique_ID}_{timestamp}.txt, where
       {field_name} - One of "static," "induction," or "radiation."
       {unique_ID} - the value of unique_ID parameter
       {timestamp} - date time of format YYMMSS_hhmmss
    '''
    # Assumes that outpath is valid

    # Get the timestamp
    today_dt = datetime.today()
    timestamp = today_dt.strftime('%Y%m%d_%H%M%S')

    # Get the filename
    fn_stem = f'_{unique_ID}_{timestamp}.txt'
    fn = os.path.join(out_path, 'efield_{0}' + fn_stem)

    # Get times from the axis
    ts = time_ax.get_bin_centers()

    # Write to the file
    write_array_to_file(fn, ts, E_sta, 'static')
    write_array_to_file(fn, ts, E_ind, 'induction')
    write_array_to_file(fn, ts, E_rad, 'radiation')

    return
