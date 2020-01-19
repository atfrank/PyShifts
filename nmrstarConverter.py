import numpy as np
import pandas as pd
import os.path
import sys


def read_nmrstar(nmrstar_file):
    if not os.path.isfile(nmrstar_file):
        raise FileNotFoundError()
    with open(nmrstar_file, "r") as fi:
        columns = []
        for ln in fi:
            if ln.startswith("_Atom_chem_shift"):
                _, col = ln.rstrip().split('.')
                columns.append(col)
    ncols = len(columns)
    # print(columns, len(columns))
    nmrstar = pd.read_csv(nmrstar_file, skiprows=ncols,
                          header=None, names=columns,
                          delim_whitespace=True)
    # transform column format
    nmrstar[['Atom_ID']] = nmrstar[['Atom_ID']].replace('"', '')
    resname_map = {'G': 'GUA', 'A': 'ADE', 'C': 'CYT', 'U': 'URA'}
    nmrstar[['Comp_ID']] = nmrstar[['Comp_ID']].replace(resname_map)
    nmrstar[['Val_err']] = nmrstar[['Val_err']].replace('?', '.')
    return nmrstar[["Comp_ID", "Comp_index_ID", "Atom_ID", "Val", "Val_err", "Entry_ID"]]


def nmrstar2measure(nmrstar_file, measured_file):
    nmrstar = read_nmrstar(nmrstar_file)
    measured = nmrstar[["Comp_ID", "Comp_index_ID", "Atom_ID", "Val", "Val_err"]]
    measured.to_csv(measured_file, sep=' ', header=None, index=None, float_format="%.2f")
    return 0


def main():
    nmrstar_file = "test/5kh8_cs.str"
    measured_file = "test/measured_shifts_5kh8.dat"
    if len(sys.argv) > 1:
        nmrstar_file = sys.argv[1]
        if len(sys.argv) > 2:
            measured_file = sys.argv[2]
    if not nmrstar2measure(nmrstar_file, measured_file):
        print("Success: Measured shifts written to file", measured_file)
    else:
        print("Failed: Exit on Error.")


if __name__ == '__main__':
    main()
