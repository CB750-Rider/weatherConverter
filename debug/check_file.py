import argparse
import json
import numpy as np
import pdb
from typing import Tuple

description = """This file reads an output json from the testWeatherCoverter program and helps to identify where siginifcant errors occur. """

parser = argparse.ArgumentParser(description=description)

parser.add_argument('filename')
parser.add_argument('--n_out', default=4)

args = parser.parse_args()

print(f"Checking {args.filename}")

def relative_difference_with_sorting(a: np.ndarray, b: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute element-wise relative difference between two 1D numpy arrays.
    
    Parameters:
    a, b: 1D numpy arrays of floats, same shape
    
    Returns:
    rel_diff: array of relative differences (np.nan where undefined)
    sort_idx: indices that sort from largest "difference" to smallest,
              with NaNs first, then infs, then finite values (largest to smallest)
    """
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    
    if a.shape != b.shape or a.ndim != 1:
        raise ValueError("a and b must be 1D arrays of the same length")
    
    # Avoid division by zero: use np.where to handle zero denominators
    abs_a = np.abs(a)
    abs_b = np.abs(b)
    
    # Relative difference: |a - b| / max(|a|, |b|)
    # This is symmetric and well-behaved when one or both are zero
    rel_diff = np.abs(a - b) / np.maximum(abs_a, abs_b)
    
    # Where both a and b are zero, max(|a|,|b|) = 0 → division by zero
    # In that case, relative difference is 0 (identical values)
    rel_diff = np.where((abs_a == 0) & (abs_b == 0), 0.0, rel_diff)
    
    # Now create sorting keys:
    # 1. NaN anywhere in a or b → highest priority
    # 2. inf anywhere in a or b → second priority  
    # 3. Finite values → sorted by rel_diff descending
    
    has_nan = np.isnan(a) | np.isnan(b)
    has_inf = np.isinf(a) | np.isinf(b)
    finite_idx = np.where(~(has_nan | has_inf))[0]
    nan_idx = np.where(has_nan)[0]
    inf_idx = np.where(has_inf)[0]
    
    # Create a sorting key tuple array
    # Lower tuple values come first in the final sort
    n = len(a)
    key = np.empty(n, dtype=object)
    
    # NaNs first (key 0)
    #key[nan_idx] = (0, 0.0)
    
    # Then infs (key 1)
    #key[inf_idx] = (1, 0.0)
    
    # Then finite values, sorted by rel_diff descending → use -rel_diff
    key[finite_idx] = list(zip(
		np.full(len(finite_idx), 2),
		-rel_diff[finite_idx]))
    
    # Get the sorting indices
    sort_idx = np.argsort(key)
    
    return rel_diff, sort_idx


def print_out(rel_diff, key, idx, set_data, std, clc):
    print(f"  rd {rel_diff[idx]} found for {key} std={std[idx]} vs clc={clc[idx]} starting with", end="")
    for key, vals in set_data.items():
        print(f" {key}={vals[idx]}", end="")
    print("")


def check_vals(key: str, std: np.ndarray, clc: np.ndarray, set_data):
    print(f"Checking {key}")

    rel_diff, sort_idx = relative_difference_with_sorting(std, clc)

    for ri in range(args.n_out):
        print_out(rel_diff, key, sort_idx[ri], set_data, std, clc)		


def check_data(file_data: dict):
    print(file_data.keys())
    standard_data = file_data['standardVariables']
    set_data = file_data['setVariables']
    calculated_data = file_data['calculatedVariables']

    for key, clc_vals in calculated_data.items():
         clc_vals = np.asarray(clc_vals)
         std_vals = np.asarray(standard_data[key])

         check_vals(key, std_vals, clc_vals, set_data)


if __name__ == "__main__":
    with open(args.filename) as infile:
        file_data = json.load(infile)

    check_data(file_data)

    print("Task Complete.")
	
