from pathlib import Path
import os
from subprocess import run

script = Path(r"C:\Users\Jure\Documents\GitHub\Docking\Docking\Molecule_builder.svl")

def adjust_script(substitutions, template, directory):
    with open(script, "r") as script_file_handle:
        script_file_text = script_file_handle.read().split("\n")
    if substitutions:
        line_parts = script_file_text[151].split("=")
        line_parts[1] = " [\'" + "\', \'".join(substitutions).replace("\"", "\'") + "\'];"
        script_file_text[151] = "=".join(line_parts)
    if template:
        line_parts = script_file_text[153].split("\'")
        line_parts[1] = template
        script_file_text[153] = "\'".join(line_parts)
    if directory:
        line_parts = script_file_text[152].split("\'")
        line_parts[1] = str(directory.absolute()).replace(os.sep, "/")
        script_file_text[152] = "\'".join(line_parts)
    with open(script, "w") as script_file_handle:
        script_file_handle.write("\n".join(script_file_text))

if __name__ == "__main__":
    folder = Path(r"C:\Users\Jure\Desktop\Test")
    total_substitutions = 3
    all = set([x for x in os.listdir(folder) if x.endswith(".mol2")])
    organized = [all]

    substitutions = ["H25", "H26", "H27"]

    # First line
    adjust_script(substitutions, list(all)[0], folder)
    run(["moebatch", "-run", str(script.absolute())], shell=True)

    for i in range(total_substitutions-1):
        all = set([x for x in os.listdir(folder) if x.endswith(".mol2")])
        organized.append(all - set([y for x in organized for y in x]))

        for x in organized[-1]:
            adjust_script([s for s in substitutions if s[1:] not in x], x, folder)
            run(["moebatch", "-run", str(script.absolute())], shell=True)