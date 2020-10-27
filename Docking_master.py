import os, time
from subprocess import Popen, PIPE
from pathlib import Path
from io import TextIOWrapper

settings = [{"root": r"Z:\Ongoing work\Projects\a6 non-review review\Docking\Runs\a1b3\6HUP",
             "ligand_path": r"Z:\Ongoing work\Projects\a6 non-review review\Docking\Templates\Ligands",
             "configuration": "gold_template.conf",
             "template_protein": "Template_protein.mol2",
             "template_ligand": "Template_ligand.mol2"},
             
            {"root": r"Z:\Ongoing work\Projects\a6 non-review review\Docking\Runs\a6b3",
             "ligand_path": r"Z:\Ongoing work\Projects\a6 non-review review\Docking\Templates\Ligands",
             "configuration": "gold_template.conf",
             "template_protein": "Template_protein.mol2",
             "template_ligand": "Template_ligand.mol2"}]
"""
settings = [{"root": r"Z:\Ongoing work\Projects\a6 non-review review\Docking\Runs\Testing",
             "ligand_path": r"Z:\Ongoing work\Projects\a6 non-review review\Docking\Templates\Ligands",
             "configuration": "gold_template.conf",
             "template_protein": "Template_protein.mol2",
             "template_ligand": "Template_ligand.mol2"}]
"""
# Global settings
process_limit = 5
n_poses = 10
processes = []
processes_ligands = []

ligands = ["XHEII087C.mol2"]

for current_settings in settings:
    #ligands = [x for x in os.listdir(current_settings["ligand_path"]) if x.endswith("mol2")]
    root = Path(current_settings["root"])
    ligand_path = Path(current_settings["ligand_path"])

    for ligand in ligands:
        if not os.path.exists(Path(root, ligand.split(".")[0])):
            processes_ligands.append(ligand)
            processes.append(Popen(["python", "Docking_script_1.py",
                                    str(root.absolute()),
                                    current_settings["configuration"],
                                    current_settings["template_protein"],
                                    current_settings["template_ligand"],
                                    str(ligand_path.absolute()),
                                    ligand], shell=True, stdout=PIPE))
            time.sleep(10)
            if len(processes) >= 1:
                print(f"Current number of processes: {len(processes_ligands)}")
                if len(processes) >= process_limit:
                    while True:
                        stdout = TextIOWrapper(processes[0].stdout, encoding="utf-8").read()
                        if "Done" in stdout:
                            processes.pop(0)
                            break
                else:
                    while len([x for x in os.listdir(Path(root, processes_ligands[0].split(".")[0])) if x.startswith("gold_soln")]) < max((n_poses/(len(processes)*process_limit), 5)):
                        time.sleep(60)
 