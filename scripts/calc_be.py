import os
import pandas as pd

from fireworks import LaunchPad
from pymatgen.core.structure import Molecule

from mispr.gaussian.workflows.base.binding_energy import get_binding_energies


def select_site(sites_dict):
    if sites_dict["O"]:
        selected_site = sites_dict["O"][0]
    elif sites_dict["S"]:
        selected_site = sites_dict["S"][0]
    elif sites_dict["Si"]:
        selected_site = sites_dict["Si"][0]
    elif sites_dict["N"]:
        selected_site = sites_dict["N"][0]
    elif sites_dict["C"]:
        selected_site = sites_dict["C"][0]
    else:
        raise ValueError(
            "Could not find a correct site to bind the ligand to. Exiting!"
        )
    return selected_site


lpad = LaunchPad.auto_load()

opt_gaussian_inputs = {
    "functional": "B3LYP",
    "basis_set": "6-31+G*",
    "route_parameters": {
        "Opt": "(calcfc, tight)",
        "SCF": "(tight, xqc)",
        "int": "ultrafine",
        "NoSymmetry": None,
        "test": None,
        "EmpiricalDispersion": "GD3",
    },
    "link0_parameters": {"%chk": "opt.chk", "%mem": "45GB", "%NProcShared": "28"},
}

freq_gaussian_inputs = {
    "functional": "B3LYP",
    "basis_set": "6-31+G*",
    "route_parameters": {
        "Freq": None,
        "iop(7/33=1)": None,
        "NoSymmetry": None,
        "test": None,
        "EmpiricalDispersion": "GD3",
    },
    "link0_parameters": {"%chk": "freq.chk", "%mem": "45GB", "%NProcShared": "28"},
}

working_dir = os.getcwd()
df = pd.read_csv(f"{working_dir}/data.csv")

for solvent in df["Abbreviation"]:
    name = solvent.lower()
    os.makedirs(f"{working_dir}/{name}", exist_ok=True)
    mol = Molecule.from_file(f"{working_dir}/pdb/{name}.pdb")
    sites = {"O": [], "S": [], "Si": [], "N": [], "C": []}
    for ind, site in enumerate(mol.sites):
        if str(site.specie) == "O":
            sites["O"].append(ind)
        elif str(site.specie) == "S":
            sites["S"].append(ind)
        elif str(site.specie) == "Si":
            sites["Si"].append(ind)
        elif str(site.specie) == "N":
            sites["N"].append(ind)
        elif str(site.specie) == "C":
            sites["C"].append(ind)
    site = select_site(sites)
    for lg_site, ligand in zip([9, 15], ["li2s8", "litfsi"]):
        wf = get_binding_energies(
            mol_operation_type=["get_from_file", "get_from_mol"],
            mol=[f"{working_dir}/pdb/{ligand}.pdb", mol],
            index=[lg_site, site],
            process_mol_func=False,
            mol_name=[ligand, name],
            working_dir=f"{working_dir}/{name}/{name}_{ligand}",
            opt_gaussian_inputs=opt_gaussian_inputs,
            freq_gaussian_inputs=freq_gaussian_inputs,
            save_to_db=True,
            save_to_file=True,
            additional_prop_doc_fields={"solvent": name, "ligand": ligand},
            wall_time=172800,
            tag="lis",
        )
        lpad.add_wf(wf)
