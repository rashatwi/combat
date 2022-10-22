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

working_dir = os.getcwd()
df = pd.read_csv(f"{working_dir}/data.csv")

# define cation force field parameters
li = Molecule.from_file(f"{working_dir}/pdb/Li.pdb")
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
    "Charges": np.asarray([1]),
}

# define polysulfide force field parameters
ps = Molecule.from_file(f"{working_dir}/pdb/S8.pdb")
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
    "Charges": np.asarray(
        [
            -0.6223,
            -0.1259,
            -0.1259,
            -0.1259,
            -0.1259,
            -0.6223,
            -0.1259,
            -0.1259,
        ]
    ),
}

# define anion force field parameters
tfsi = Molecule.from_file(f"{working_dir}/pdb/TFSI.pdb")
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
        {
            "coeffs": ["fourier", 1, 0.1734, 3, 0],
            "types": [("ot", "st", "ct", "ft")],
        },
        {
            "coeffs": ["fourier", 1, 0.158, 3, 0],
            "types": [("nt", "st", "ct", "ft")],
        },
        {
            "coeffs": ["fourier", 1, -0.0018, 3, 0],
            "types": [("ot", "st", "nt", "st")],
        },
        {
            "coeffs": ["fourier", 3, 3.91646, 1, 0, -1.245, 2, 180, -0.3818, 3, 0],
            "types": [("ct", "st", "nt", "st")],
        },
    ],
    "Impropers": [],
    "Improper Topologies": None,
    "Charges": np.asarray(
        [
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
    ),
}

# define DOL force field parameters
dol = Molecule.from_file(f"{working_dir}/pdb/DOL.pdb")
dol.set_charge_and_spin(0, 1)
dol_param = {
    "Molecule": dol,
    "Labels": ["os", "os", "ct", "ct", "co", "hc", "hc", "hc", "hc", "hc", "hc"],
    "Masses": OrderedDict({"os": 16.00, "ct": 12.01, "co": 12.01, "hc": 1.008}),
    "Nonbond": [[0.14, 2.9], [0.066, 3.5], [0.066, 3.5], [0.03, 2.5]],
    "Bonds": [
        {"coeffs": [320.0, 1.41], "types": [("ct", "os")]},
        {"coeffs": [320.0, 1.38], "types": [("co", "os")]},
        {"coeffs": [268.0, 1.529], "types": [("ct", "ct")]},
        {"coeffs": [340.0, 1.09], "types": [("ct", "hc")]},
        {"coeffs": [340.0, 1.09], "types": [("co", "hc")]},
    ],
    "Angles": [
        {"coeffs": [60.0, 109.5], "types": [("co", "os", "ct")]},
        {"coeffs": [50.0, 109.5], "types": [("ct", "ct", "os")]},
        {"coeffs": [35.0, 109.5], "types": [("hc", "ct", "os")]},
        {"coeffs": [37.5, 110.7], "types": [("ct", "ct", "hc")]},
        {"coeffs": [33.0, 107.8], "types": [("hc", "ct", "hc")]},
        {"coeffs": [92.6, 111.55], "types": [("os", "co", "os")]},
        {"coeffs": [35.0, 109.5], "types": [("hc", "co", "os")]},
        {"coeffs": [33.0, 109.5], "types": [("hc", "co", "hc")]},
    ],
    "Dihedrals": [
        {
            "coeffs": ["fourier", 3, 0.0625, 1, 0, -0.011, 2, 180, 0.2805, 3, 0],
            "types": [("ct", "os", "ct", "ct")],
        },
        {
            "coeffs": ["fourier", 3, 0.0625, 1, 0, -0.011, 2, 180, 0.2805, 3, 0],
            "types": [("co", "os", "ct", "ct")],
        },
        {
            "coeffs": ["fourier", 1, 0.3705, 3, 0],
            "types": [("ct", "os", "ct", "hc")],
        },
        {
            "coeffs": ["fourier", 1, 0.3705, 3, 0],
            "types": [("co", "os", "ct", "hc")],
        },
        {
            "coeffs": ["fourier", 3, -0.5335, 1, 0, -0.7835, 2, 180, 0.3375, 3, 0],
            "types": [("ct", "os", "co", "os")],
        },
        {
            "coeffs": ["fourier", 1, 0.3705, 3, 0],
            "types": [("ct", "os", "co", "hc")],
        },
        {
            "coeffs": ["fourier", 3, 1.119, 1, 0, -1.1635, 2, 180, -0.3415, 3, 0],
            "types": [("os", "ct", "ct", "os")],
        },
        {
            "coeffs": ["fourier", 1, 0.234, 3, 0],
            "types": [("hc", "ct", "ct", "os")],
        },
        {"coeffs": ["fourier", 1, 0.15, 3, 0], "types": [("hc", "ct", "ct", "hc")]},
    ],
    "Impropers": [],
    "Improper Topologies": None,
    "Charges": np.asarray(
        [
            -0.4000,
            -0.4000,
            0.1400,
            0.1400,
            0.2000,
            0.0300,
            0.0300,
            0.0300,
            0.0300,
            0.1000,
            0.1000,
        ]
    ),
}

for solvent in df["Abbreviation"]:
    if solvent != "DOL":
        type_cosolvent = solvent
        name_cosolvent = solvent.lower()
        num_mols = number_molecules(type_cosolvent, df, conc_ps=0.25)
        charge_scaling_factor = 1 / np.sqrt(epsilon(type_cosolvent, df))
        system_species_data = {
            "dol": {
                "molecule": dol,
                "molecule_operation_type": "get_from_mol",
                "ff_param_method": "get_from_dict",
                "ff_param_data": dol_param,
                "mol_mixture_type": "Solvents",
                "mixture_data": num_mols["dol"],
                "save_ff_to_file": True,
            },
            name_cosolvent: {
                "molecule": f"{working_dir}/pdb/{name_cosolvent}.pdb",
                "molecule_operation_type": "get_from_file",
                "ff_param_method": "get_from_opls",
                "ff_param_data": f"{working_dir}/pdb/{name_cosolvent}.pdb",
                "mol_mixture_type": "Solvents",
                "mixture_data": num_mols[type_cosolvent],
                "save_ff_to_file": True,
            },
            "tfsi": {
                "molecule": tfsi,
                "molecule_operation_type": "get_from_mol",
                "ff_param_method": "get_from_dict",
                "ff_param_data": tfsi_param,
                "mol_mixture_type": "Solutes",
                "mixture_data": num_mols["salt"],
                "save_ff_to_file": True,
            },
            "ps": {
                "molecule": ps,
                "molecule_operation_type": "get_from_mol",
                "ff_param_method": "get_from_dict",
                "ff_param_data": ps_param,
                "mol_mixture_type": "Solutes",
                "mixture_data": num_mols["ps"],
                "save_ff_to_file": True,
            },
            "li": {
                "molecule": li,
                "molecule_operation_type": "get_from_mol",
                "ff_param_method": "get_from_dict",
                "ff_param_data": li_param,
                "mol_mixture_type": "Solutes",
                "mixture_data": num_mols["li"],
                "save_ff_to_file": True,
            },
        }

        wf = lammps_workflow(
            system_species_data=system_species_data,
            system_mixture_type="number of molecules",
            box_data=60.0,
            box_data_type="cubic",
            scale_charges=True,
            charge_scaling_factor=charge_scaling_factor,
            working_dir=f"{working_dir}/{name_cosolvent}",
            recipe=[
                ["min", ["template_filename", "min.in"]],
                ["npt", ["template_filename", "npt.in"]],
                ["melt_quench", ["template_filename", "melt_quench.in"]],
                ["nvt", ["template_filename", "nvt.in"]],
            ],
            template_dir=f"{working_dir}/md_templates",
        )

        lpad.add_wf(wf)
