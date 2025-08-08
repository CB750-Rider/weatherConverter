"""
Make sure and renumber these to look like 
Table8-0.txt
...
Table8-5.txt

Also, remove blank lines.
"""

import os
import json
import numpy as np
from typing import List, Tuple
import matplotlib.pyplot as plt

root_data_path = "/path/to/folders"

table_8 = \
    {'name': 'table_8',
     'folder_name': 'Table8',
     'filename_base': 'Table8-',
     'file_count': 6,
     'title': 'Table VIII, Atmospheric Composition Number Density',
     'columns': ['altitude Z', 'altitude H', 'N2', 'O', 'O2', 'Ar', 'He', 'H'],
     'column_units': ['m', 'm', 'per m^3', 'per m^3', 'per m^3', 'per m^3',
                      'per m^3', 'per m^3'],
     'pages': [225,230],
    }

table_1 = \
    {'name': 'table_1',
     'folder_name': 'Table1',
     'filename_base': 'Table1-',
     'file_count': 6,
     'title': 'Table I, Geometric Altitude, Metric Units',
     'columns': ['altitude Z', 'altitude H', 'Temperature K', 'Temperature C',
                 'Molecular Temperature K', 'Pressure mb', 'Pressure Torr',
                 'Pressure Ratio', 'Mass Density', 'Density Ratio'],
     'column_units': ['m', 'm', 'K', 'C', 'K', 'mb', 'Torr', '', 'kg/m^3', ''],
     'pages': [65,88],
    }

tables = [table_8, table_1]


def _read_page(fname: str) -> np.ndarray:
    return np.genfromtxt(fname, delimiter=" ", skip_header=0,
                         missing_values='I', filling_values=0.0)


def read_pages(table: dict) -> Tuple[List[np.ndarray], List[str]]:
    out = []
    file_names_out = []
    try:
        for fname in table['file_names']:
            out.append(_read_page(os.path.join(root_data_path,
                                               table['folder_name'],
                                               fname)))
            file_names_out.append(fname)
    except KeyError:
        for fi in range(table['file_count']):
            fname = f"{table['filename_base']}{fi}.txt"
            out.append(_read_page(os.path.join(root_data_path,
                                               table['folder_name'],
                                               fname)))
            file_names_out.append(fname)
    return out, file_names_out


def save_data(merged_data, table):

    for idx, column in enumerate(table['columns']):
        table[column] = merged_data[:, idx].tolist()

    with open(table['name'] + '.json', 'w') as outfile:
        json.dump(table, outfile)


if __name__ == "__main__":
    for table in tables:
        print(f"Processing {table['name']}: {table['title']}")
        data_sets, file_names = read_pages(table)
        # for di in range(len(table['columns'])):
        #     merged_data = plot_data(data_sets, table, di, file_names)
        save_data(merged_data, table)
        # plt.show()
