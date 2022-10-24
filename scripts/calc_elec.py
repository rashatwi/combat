import os

import pandas as pd

from fireworks import LaunchPad

from mispr.gaussian.workflows.base.esp import get_esp_charges

lpad = LaunchPad.auto_load()

opt_gaussian_inputs = {
    "functional": "B3LYP",
    "basis_set": "6-31+G*",
    "route_parameters": {
        "Opt": "(calcfc, tight)",
        "SCF": "Tight",
        "int": "ultrafine",
        "NoSymmetry": None,
        "test": None,
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
        "Polar": None,
        "test": None,
    },
    "link0_parameters": {"%chk": "freq.chk", "%mem": "45GB", "%NProcShared": "28"},
}

esp_gaussian_inputs = {
    "functional": "B3LYP",
    "basis_set": "6-31+G*",
    "route_parameters": {
        "pop": "MK",
        "iop(6/50=1)": None,
        "NoSymmetry": None,
        "Density": None,
        "test": None,
    },
    "link0_parameters": {"%chk": "esp.chk", "%mem": "45GB", "%NProcShared": "28"},
}

working_dir = os.getcwd()
df = pd.read_csv(f"{working_dir}/data.csv")

for solvent in df["Abbreviation"]:
    name = solvent.lower()
    wf = get_esp_charges(
        mol_operation_type="get_from_file",
        mol=f"{working_dir}/pdb/{name}.pdb",
        process_mol_func=False,
        mol_name=name,
        opt_gaussian_inputs=opt_gaussian_inputs,
        freq_gaussian_inputs=freq_gaussian_inputs,
        esp_gaussian_inputs=esp_gaussian_inputs,
        format_chk=True,
        save_to_db=True,
        save_to_file=True,
        additional_prop_doc_fields={"solvent": name},
        tag="lis",
    )
    lpad.add_wf(wf)
