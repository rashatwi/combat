def epsilon(cosolvent, data):
    """
    Calculates the electronic high-frequency dielectric permittivity of
    the solvent mixture composed of DOL and a cosolvent using the
    refractive index of each solvent

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
