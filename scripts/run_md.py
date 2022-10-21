import os
from collections import OrderedDict

import numpy as np
import pandas as pd

from epsilon import epsilon
from fireworks import LaunchPad, Workflow
from mispr.lammps.workflows.base import lammps_workflow
from number_molecules import number_molecules
from pymatgen.core.structure import Molecule

lpad = LaunchPad.auto_load()

type_cosolvent = "DMT"
name_cosolvent = "dmt"
df = pd.read_csv("data/data.csv")
number_molecules = number_molecules(type_cosolvent, df, conc_ps=0.25)
eps = epsilon(type_cosolvent, df)
tfsi_charges = [
    -0.66,
    1.02,
    -0.53,
    -0.53,
    0.35,
    -0.16,
    -0.16,
    -0.16,
    1.02,
    -0.53,
    -0.53,
    0.35,
    -0.16,
    -0.16,
    -0.16,
]
ps_charges = [-0.6223, -0.1259, -0.1259, -0.1259, -0.1259, -0.6223, -0.1259, -0.1259]
li_charges = [1]
tfsi_charges = [i / np.sqrt(eps) for i in tfsi_charges]
ps_charges = [i / np.sqrt(eps) for i in ps_charges]
li_charges = [i / np.sqrt(eps) for i in li_charges]

li = Molecule.from_file("li.pdb")
li.set_charge_and_spin(1, 1)
li_param = {
    "Molecule": li,
    "Labels": ["li"],
    "Masses": OrderedDict({"li": 6.941}),
    "Nonbond": [[0.16505739986523993, 1.506000000135153]],
    "Bonds": [],
    "Angles": [],
    "Dihedrals": [],
    "Impropers": [],
    "Improper Topologies": None,
    "Charges": np.asarray(li_charges),
}

ps = Molecule.from_file("ps.pdb")
ps.set_charge_and_spin(-2, 1)
ps_param = {
    "Molecule": ps,
    "Labels": ["sp", "ss", "ss", "ss", "ss", "sp", "ss", "ss"],
    "Masses": OrderedDict({"sp": 32.06, "ss": 32.06}),
    "Nonbond": [
        [0.3440010005854374, 3.5903218336860525],
        [0.3440010005854374, 3.5903218336860525],
    ],
    "Bonds": [
        {"coeffs": [157.87, 2.133], "types": [("sp", "ss")]},
        {"coeffs": [154.36, 2.117], "types": [("ss", "ss")]},
    ],
    "Angles": [
        {"coeffs": [73.56, 112.74], "types": [("sp", "ss", "ss")]},
        {"coeffs": [73.56, 112.74], "types": [("ss", "ss", "ss")]},
    ],
    "Dihedrals": [
        {
            "coeffs": [
                "fourier",
                4,
                5.4317,
                1,
                0,
                -3.59,
                2,
                180,
                1.0556,
                3,
                0,
                0.05832,
                4,
                180,
            ],
            "types": [("sp", "ss", "ss", "ss")],
        },
        {
            "coeffs": [
                "fourier",
                4,
                5.4317,
                1,
                0,
                -3.59,
                2,
                180,
                1.0556,
                3,
                0,
                0.05832,
                4,
                180,
            ],
            "types": [("ss", "ss", "ss", "ss")],
        },
    ],
    "Impropers": [],
    "Improper Topologies": None,
    "Charges": np.asarray(ps_charges),
}

tfsi = Molecule.from_file("tfsi.pdb")
tfsi.set_charge_and_spin(-1, 1)
tfsi_param = {
    "Molecule": tfsi,
    "Labels": [
        "nt",
        "st",
        "ot",
        "ot",
        "ct",
        "ft",
        "ft",
        "ft",
        "st",
        "ot",
        "ot",
        "ct",
        "ft",
        "ft",
        "ft",
    ],
    "Masses": OrderedDict(
        {"nt": 14.01, "st": 32.06, "ot": 16.00, "ct": 12.01, "ft": 19.00}
    ),
    "Nonbond": [
        [0.17, 3.2499985240310356],
        [0.250, 3.549340493071112],
        [0.210, 2.959565541662207],
        [0.066, 3.4994501648552525],
        [0.053, 2.948874757044523],
    ],
    "Bonds": [
        {"coeffs": [441.8, 1.323], "types": [("ct", "ft")]},
        {"coeffs": [235.4, 1.818], "types": [("ct", "st")]},
        {"coeffs": [637.1, 1.442], "types": [("st", "ot")]},
        {"coeffs": [372.0, 1.570], "types": [("st", "nt")]},
    ],
    "Angles": [
        {"coeffs": [93.3, 107.1], "types": [("ft", "ct", "ft")]},
        {"coeffs": [82.90, 111.8], "types": [("ft", "ct", "st")]},
        {"coeffs": [104.0, 102.6], "types": [("ct", "st", "ot")]},
        {"coeffs": [115.8, 118.5], "types": [("ot", "st", "ot")]},
        {"coeffs": [94.30, 113.6], "types": [("ot", "st", "nt")]},
        {"coeffs": [97.50, 100.2], "types": [("ct", "st", "nt")]},
        {"coeffs": [80.20, 125.6], "types": [("st", "nt", "st")]},
    ],
    "Dihedrals": [
        {"coeffs": ["fourier", 1, 0.1734, 3, 0], "types": [("ot", "st", "ct", "ft")]},
        {"coeffs": ["fourier", 1, 0.158, 3, 0], "types": [("nt", "st", "ct", "ft")]},
        {"coeffs": ["fourier", 1, -0.0018, 3, 0], "types": [("ot", "st", "nt", "st")]},
        {
            "coeffs": ["fourier", 3, 3.91646, 1, 0, -1.245, 2, 180, -0.3818, 3, 0],
            "types": [("ct", "st", "nt", "st")],
        },
    ],
    "Impropers": [],
    "Improper Topologies": None,
    "Charges": np.asarray(tfsi_charges),
}

system_species_data = {
    "dol": {
        "molecule": "dol.pdb",
        "molecule_operation_type": "get_from_file",
        "ff_param_method": "get_from_opls",
        "ff_param_data": "dol.pdb",
        "mol_mixture_type": "Solvents",
        "mixture_data": number_molecules["dol"],
        "save_ff_to_file": True,
    },
    name_cosolvent: {
        "molecule": f"{name_cosolvent}.pdb",
        "molecule_operation_type": "get_from_mol",
        "ff_param_method": "get_from_opls",
        "ff_param_data": f"{name_cosolvent}.pdb",
        "mol_mixture_type": "Solvents",
        "mixture_data": number_molecules[type_cosolvent],
        "save_ff_to_file": True,
    },
    "tfsi": {
        "molecule": tfsi,
        "molecule_operation_type": "get_from_mol",
        "ff_param_method": "get_from_dict",
        "ff_param_data": tfsi_param,
        "mol_mixture_type": "Solutes",
        "mixture_data": number_molecules["salt"],
        "save_ff_to_file": True,
    },
    "ps": {
        "molecule": ps,
        "molecule_operation_type": "get_from_mol",
        "ff_param_method": "get_from_dict",
        "ff_param_data": ps_param,
        "mol_mixture_type": "Solutes",
        "mixture_data": number_molecules["ps"],
        "save_ff_to_file": True,
    },
    "li": {
        "molecule": li,
        "molecule_operation_type": "get_from_mol",
        "ff_param_method": "get_from_dict",
        "ff_param_data": li_param,
        "mol_mixture_type": "Solutes",
        "mixture_data": number_molecules["li"],
        "save_ff_to_file": True,
    },
}

wf = lammps_workflow(
    system_species_data=system_species_data,
    system_mixture_type="number of molecules",
    box_data=60.0,
    box_data_type="cubic",
    recipe=[
        ["min", ["template_filename", "min.in"]],
        ["npt", ["template_filename", "npt.in"]],
        ["melt_quench", ["template_filename", "melt_quench.in"]],
        ["nvt", ["template_filename", "nvt.in"]],
    ],
    template_dir="data/templates",
)

lpad.add_wf(wf)
