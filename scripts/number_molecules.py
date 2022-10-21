import numpy as np

AVOGADRO = 6.023 * (10 ** 23)
CONVERSION_CONSTANT = 10 ** (-24)

rho = {
    "LiTFSI": 1330,
    "DOL": 1060,
    "Li2S8": 1660,
}

mw = {
    "LiTFSI": 287.085518,
    "DOL": 74.08,
    "Li2S8": 270.402,
}


def number_molecules(cosolvent, prop, conc_ps=0.25, cosolvent_ratio=0.5, length=6):
    """
    Calculates the number of species of each type for a solution of 1 M LiTFSI,
    0.25 M Li2S8, in 1:1 DOL:cosolvent mixture.
    Parameters
    ----------
    cosolvent: str
        name of the cosolvent.
    prop: pd.DataFrame
        Dataframe of solvent properties including name, MW, and density.
    conc_ps: float
        Concentration of PS to use. Defaults to 0.25 M.
    cosolvent_ratio: float
        volumetric ratio of FLS co-solvent.
    length: int
        Simulation box in nm.

    Returns
    -------
    num_mols: dictionary
        Number of molecules of each type
    """
    cosolvent_prop = prop[prop["Abbreviation"] == cosolvent]
    volume = (length ** 3) * CONVERSION_CONSTANT
    number_ps = np.round(conc_ps * AVOGADRO * volume, 0)
    number_salt = np.round(1 * AVOGADRO * volume)  # conc_salt = 1
    volume_salt = number_salt * mw["LiTFSI"] / (AVOGADRO * rho["LiTFSI"])
    volume_ps = number_ps * mw["Li2S8"] / (AVOGADRO * rho["Li2S8"])

    remaining_volume = volume - volume_salt - volume_ps

    volume_cosolvent = remaining_volume * cosolvent_ratio
    volume_dol = remaining_volume * (1 - cosolvent_ratio)
    number_dol = np.round(volume_dol * AVOGADRO * rho["DOL"] / mw["DOL"])
    number_cosolvent = np.round(
        volume_cosolvent
        * AVOGADRO
        * cosolvent_prop["Density"].values[0] * 1000
        / cosolvent_prop["MolecularWeight"].values[0]
    )
    number_li = number_ps * 2 + number_salt

    num_mols = {
        "li": number_li,
        "salt": number_salt,
        "ps": number_ps,
        cosolvent: number_cosolvent,
        "dol": number_dol,
    }
    return num_mols