from schrodinger import structure
import numpy, pandas, os
from pathlib import Path

# Global settings
# Definition of the coordinate system. The numbering is as found in the molecular viewer.
coordinate_system = {
    "center": ("A", 208),
    "x": ("A", 205),
    "y": ("A", 210),
    "z": ("E", 115)
}


def compute(folder_path):
    '''
    Function calculating the coordinates of all the ligand atoms for all the given structures.
    The function calculates the angles to the three coordinate axes and the distance of the atom to the center.
    :param folder_path: String or path-like object delineating the folder with PDB files to be analyzed
    :return: pandas.DataFrame with coordinates for all the atoms for all the analyzed files.  
    '''
    
    folder_path = Path(folder_path)
    files = [x for x in os.listdir(folder_path) if x.endswith(".pdb")]
    s = structure.StructureReader.read(files[0])

    # Get the coordinates of the alpha carbons defining the coordinate system
    center = s.findResidue(":".join([str(x) for x in coordinate_system['center']])).getAlphaCarbon()
    axes = {
        "alpha": s.findResidue(":".join([str(x) for x in coordinate_system['x']])).getAlphaCarbon(),
        "beta":  s.findResidue(":".join([str(x) for x in coordinate_system['y']])).getAlphaCarbon(),
        "gamma": s.findResidue(":".join([str(x) for x in coordinate_system['z']])).getAlphaCarbon()
    }

    # Constructing the results table
    ligand_chain = [x for x in s.chain if len(x.residue) == 1][0]

    columns = ["Complex"]
    for ligand_atom in ligand_chain.atom:
        columns.extend(["_".join((ligand_atom.property["s_m_pdb_atom_name"].strip(), x)) for x in list(axes) + ["distance"]])
    results = pandas.DataFrame({c: v for c, v in zip(columns, [[x.split(".")[0] for x in files]] + [[0 for _ in files] for _ in columns[1:]])}, columns=columns)
    results.index = results["Complex"]

    for fname in files:
        complex_name = fname.split(".")[0]
        s = structure.StructureReader.read(fname)
        ligand_chain = [x for x in s.chain if len(x.residue) == 1][0]

        for ligand_atom in ligand_chain.atom:
            ligand_atom_name = ligand_atom.property["s_m_pdb_atom_name"].strip()
            # Calculating the angles to the axes
            for axis_name, axis_atom in axes.items():
                results.loc[complex_name, "_".join((ligand_atom_name, axis_name))] = s.measure(ligand_atom, center, axis_atom)
            # Adding distance from the center atom
            results.loc[complex_name, f"{ligand_atom_name}_distance"] = s.measure(ligand_atom, center)

    return results