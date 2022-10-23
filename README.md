# ComBat Database 

<img src=logo.png>

## Contents
* [Overview](#overview-)
* [Contents](#contents-)
* [Updates](#updates-)
* [Citation](#citation-)
* [Licensing](#licensing-)
* [Contact](#contact-)

## Overview [↑](#overview)
The Computational Database for Li-S Batteries (ComBat) is a publicly 
available dataset of quantum-chemical and molecular dynamics properties 
for Li-S electrolytes containing solvents spanning 16 different chemical 
classes. 

The data includes DFT-optimized geometries (in PDB format) and several 
properties, such as binding energies, partial atomic charges, 
polarizabilities, dipole moments, as well as MD-derived properties, such
as radial distribution functions, coordination numbers, diffusion coefficients,
and more. More details about each property can be found in the ReadME file 
associated with each dataset. 

## Contents [↑](#contents)
* data.csv: A CSV file containing the solvent metadata and the following properties:

| Column           | Description                                                                                                                        |
|------------------|------------------------------------------------------------------------------------------------------------------------------------|
| Abbreviation     | Common abbreviation of the solvent                                                                                                 |
| Type             | Chemical class of the solvent                                                                                                      |
| InChI Key        | InChI Key of the solvent                                                                                                           |
| InChI            | InChI of the solvent                                                                                                               |
| CanonicalSMILES  | SMILES of the solvent                                                                                                              |   
| IUPACName        | IUPAC name of the solvent                                                                                                          |
| CAS              | CAS number of the solvent                                                                                                          |
| CID              | PubChem CID of the solvent                                                                                                         |
| MolecularFormula | Molecular formula of the solvent                                                                                                   |
| MolecularWeight  | Molecular weight of the solvent                                                                                                    |
| Density          | Density of the solvent (from public databases)                                                                                     |
| DN               | Donor number of the solvent (from public databases & literature)                                                                   |
| De               | Dielectric constant of the solvent (from public databases & literature)                                                            |
| Viscosity        | Viscosity of the solvent (from public databases & literature)                                                                      |
| Refraction       | Refraction index of the solvent (from ChemSpider)                                                                                  |
| FlashPoint       | Flash point of the solvent in °C (from public databases & literature)                                                              |             | 
| MeltingPoint     | Melting point of the solvent in °C (from public databases & literature)                                                            |
| BoilingPoint     | Boiling point of the solvent in °C (from public databases & literature)                                                            |
| DensityT         | Temperature in °C at which the density was measured (from public databases)                                                        |
| ViscosityT       | Temperature in °C at which the viscosity was measured (from public databases & literature)                                         |
| DeT              | Temperature in °C at which the dielectric constant was measured (from public databases & literature)                               |
| Diffusion        | Diffusion coefficient of the pure solvent (from public databases)                                                                  |
| DiffusionT       | Temperature in °C at which the diffusion coefficient was measured (from public databases)                                          |
| BE_Salt          | Binding energy of the solvent with LiTFSI in kcal/mol (from DFT)                                                                   |  
| BE_PS            | Binding energy of the solvent with Li2S8 in kcal/mol (from DFT)                                                                    |
| Li - O (DOL) | Coordination number of Li+ with O atoms of DOL (from MD simulations of 1 M LiTFSI, 0.25 M Li2S8 in DOL/solvent (1/1, v/v))         |
| Li - X (Solvent) | Coordination number of Li+ with X sites of the solvent (from MD simulations of 1 M LiTFSI, 0.25 M Li2S8 in DOL/solvent (1/1, v/v)) |
| Li - O (TFSI) | Coordination number of Li+ with O atoms of TFSI (from MD simulations of 1 M LiTFSI, 0.25 M Li2S8 in DOL/solvent (1/1, v/v))        |
| Li - S (PS) | Coordination number of Li+ with terminal S atoms of PS (from MD simulations of 1 M LiTFSI, 0.25 M Li2S8 in DOL/solvent (1/1, v/v)) |
| Li - Li | Coordination number of Li+ with Li+ (from MD simulations of 1 M LiTFSI, 0.25 M Li2S8 in DOL/solvent (1/1, v/v))                    |
| Dissociation | Fraction of Li+ dissociated from the PS (from MD simulations of 1 M LiTFSI, 0.25 M Li2S8 in DOL/solvent (1/1, v/v))                 |
| Bridging | Fraction of Li+ with more than one PS in their solvation shell (from MD simulations of 1 M LiTFSI, 0.25 M Li2S8 in DOL/solvent (1/1, v/v)) |

## Updates [↑](#updates)
For a description of the updates corresponding to each version, 
see [updates.md](https://github.com/rashatwi/combat/blob/main/updates.md).

## Citation [↑](#citation)
MISPR software, which was used to produce the ComBat database, was
introduced [here](https://doi.org/10.1038/s41598-022-20009-w). 
Please include this citation where relevant. 

## Licensing [↑](#licensing)
The ComBat Database is made publicly available under a CC BY 4.0 license. 
You are free to copy, share, and adapt the material for any purpose, 
provided that you give appropriate credit and indicate any changes made.

## Contact [↑](#contact)
If you have any questions, you can reach the corresponding author at the 
e-mail [here](https://www.rashatwi.com).