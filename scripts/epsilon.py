def epsilon(cosolvent, data):
    """
    Calculates the number of species of each type for a solution of 1 M LiTFSI,
    0.25 M Li2S8, in 1:1 DOL:cosolvent mixture.
    Parameters
    ----------
    cosolvent: str
        name of the cosolvent.
    data: pd.DataFrame
        Dataframe of solvent properties including name and refractive index.

    Returns
    -------
    e: float

    """
    e1 = data[data["Abbreviation"] == "DOL"]["Refraction"].values[0]
    e2 = data[data["Abbreviation"] == cosolvent]["Refraction"].values[0]
    e = ((e1 + e2) / 2) ** 2
    return e
