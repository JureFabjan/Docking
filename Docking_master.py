import os, time
from subprocess import Popen, PIPE
from pathlib import Path
from io import TextIOWrapper

settings = [{"root": r"Z:\Ongoing work\Projects\a6 non-review review\Docking\Runs\Replication experiment",
             "ligand_path": r"Z:\Ongoing work\Projects\a6 non-review review\Docking\Templates\Ligands",
             "configuration": "gold_template.conf",
             "template_protein": "Template_protein.mol2",
             "template_ligand": "Template_ligand.mol2"}
             ]

# Global settings
process_limit = 10
n_poses = 200
processes = []
processes_ligands = []
repetitions = 100

ligands = ["PZII029.mol2"]

for current_settings in settings:
    # ligands = [x for x in os.listdir(current_settings["ligand_path"]) if x.endswith("mol2")]
    root = Path(current_settings["root"])
    ligand_path = Path(current_settings["ligand_path"])
    output_dirs = []

    for ligand in ligands:
        for repetition in range(repetitions):
            out_dir = "_".join((ligand.split(".")[0], str(repetition)))
            if not os.path.exists(Path(root, out_dir)):
                processes_ligands.append(ligand)
                output_dirs.append(out_dir)
                processes.append(Popen(["python", "Docking_script_1.py",
                                        str(root.absolute()),
                                        current_settings["configuration"],
                                        current_settings["template_protein"],
                                        current_settings["template_ligand"],
                                        str(ligand_path.absolute()),
                                        ligand,
                                        out_dir], shell=True, stdout=PIPE))
                time.sleep(10)
                if len(processes) >= 1:
                    print(f"Current number of processes: {len(processes_ligands)}\n" + "Current ligands: " + "\n".join(processes_ligands))
                    if len(processes) >= process_limit:
                        while True:
                            stdout = TextIOWrapper(processes[0].stdout, encoding="utf-8").read()
                            if "Done" in stdout:
                                processes.pop(0)
                                print(f"{processes_ligands.pop(0)} run finished")
                                output_dirs.pop(0)
                                break
                    else:
                        while len([x for x in os.listdir(Path(root, output_dirs[0])) if x.startswith("gold_soln")]) < max((n_poses/process_limit*len(processes), 5)):
                            time.sleep(60)
    