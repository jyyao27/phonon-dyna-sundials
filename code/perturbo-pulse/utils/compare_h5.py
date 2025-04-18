# compare_h5.py
"""
A very simple file compare tool specific to perturbo h5 file output

It takes as input two file names, and returns the maximum difference in equivalent values
between the two files.
"""

import numpy as np
import h5py
import sys

def compare_h5(a, b: str, show_progression: bool = False):
    max_diff = 0.
    max_diff_pct = 0.
    max_diff_step = -1
    with h5py.File(a, 'r') as h5_a:
        with h5py.File(b, 'r') as h5_b:
            dr_a = h5_a['dynamics_run_1']
            dr_b = h5_b['dynamics_run_1']
            assert dr_a['num_steps'][()] == dr_b['num_steps'][()]
            for step in range(dr_a['num_steps'][()]):
                sn_a = np.array(dr_a[f'snap_t_{step}'])
                sn_b = np.array(dr_b[f'snap_t_{step}'])
                dif = np.abs(sn_a - sn_b)
                dif_pct = dif/sn_a
                max_idx = np.unravel_index(dif_pct.argmax(), dif_pct.shape)
                m = dif_pct[max_idx]
                if m > max_diff_pct:
                    if show_progression:
                        print(f"New max on step {step + 1} at index {max_idx}: {dif[max_idx]} ({m * 100.}%) ({sn_a[max_idx]} vs. {sn_b[max_idx]}")
                    max_diff_pct = m
                    max_diff = dif[max_idx]
                    max_diff_step = step
    return max_diff, max_diff_pct, max_diff_step

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print("usage: python compare_h5 <h5_file_a> <h5_file_b> [show_progress (bool)]")
        exit(-1)
    show_progress = False
    if len(sys.argv) >= 4:
        show_progress = bool(sys.argv[3])
    max_diff, max_diff_pct, max_diff_step = compare_h5(sys.argv[1], sys.argv[2], show_progress)
    print(f'Max diff: {max_diff} ({max_diff_pct * 100.}%) on step {max_diff_step + 1}')
