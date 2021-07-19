from pathlib import Path
import os
import inspect
from subprocess import run
import re
from collections import defaultdict

_location = str(Path(inspect.getfile(inspect.currentframe())).parents[0].absolute())
_script = Path(_location, r"svl_scripts\Molecule_builder.svl")

def construct_list(tag, l):
    result = [tag] * (len(l)*2)
    result[1::2] = l
    return result

def generate(folder, total_substitutions, substitutions, sidechains, prohibited, script=_script):
    """
    Constructs the molecule files.
    param folder: Path to the target folder.
    param total_substitutions: The number of positions on the molecule that can be substituted at the same time.
    param substitutions: The (hidrogen) atom codes that should be substituted.
    param sidechains: Sidechains that should be built. Note these are dependent on the script that gets called!!
    param prohibited: Defaultdict with sets as default values. It indicates which sidechains should not be attached at certain positions (i.e. "H25": {"F",}).
    param script: The SVL script used for the generation of the molecules. Default is the one included in the package. 
    """
    for i in range(total_substitutions):
        all = set([x for x in os.listdir(folder) if x.endswith(".mol2")])
        if i == 0:
            organized = [all]
        else:
            organized.append(all - set([y for x in organized for y in x]))
        for ligand in organized[-1]:
            substitutions_current = {s for s in substitutions if not re.findall(f"[A-z]+{s[1:]}[_|\.]", ligand)}
            for position in substitutions_current:
                sidechains_current = sidechains - prohibited[position]
                run(["moebatch",
                        "-run", str(script.absolute()),
                        "-i", str(folder.absolute()),
                        "-t", ligand] +
                        construct_list("-s", substitutions_current) +
                        construct_list("-c", sidechains_current), shell=True)

if __name__ == "__main__":
    folder = Path(r"C:\Users\Jure\Desktop\Test")
    total_substitutions = 2
    substitutions = ["H25", "H26", "H27", "H30", "H19", "H20", "H1"]
    sidechains = {'F', 'Cl', 'Br', 'Me', 'NH', 'Et', 'Pr', 'iPr', 'tBu', 'CCH', 'CF', 'COMe', 'OMe', 'CN', 'COOH', 'COOMe', 'COOEt', 'NO'}
    prohibited = defaultdict(set, {"H25": {"F",}}) # To be filled

    generate(folder, total_substitutions, substitutions, sidechains, prohibited)
