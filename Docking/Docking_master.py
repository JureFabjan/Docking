import os, time, inspect
from subprocess import Popen, PIPE
from pathlib import Path
from io import TextIOWrapper


_location = str(Path(inspect.getfile(inspect.currentframe())).parents[0].absolute())
class Master:
    """
    The master manager of all docking threads.
    """
    def __init__(self, settings, ligands=None, process_limit=10, n_poses=200, repetitions=1):
        """
        Class used for constructing a paralelized docking run. Takes the settings and uses Docking_script to spawn individual processes. Manages the concurrent number of processes.
        :param settings: A list of dictionaries containing all the settings. The dicts have to contain the following entries:
            - root: the root directory in which docking is run
            - ligand_path: path to the dorectory containing all the ligand input files
            - configuration: path to the configuration file used for the docking
            - template_protein: path to the protein file used for docking
            - template_ligand: path to the ligand file used for defining the binding site
            Optional entries:
            - center, x_axis, y_axis, z_axis: Lists defining the atoms for coordinate system definition. Lists should have the following structure: [chain_number, aa_code, residue_number, atom_name]. 
        :param ligands:
            - List of the ligand files that should be run. If left out, all the ligands in the ligand folder
        :specified in settings will be used.
        :param process_limit: The limit on the number of processes that can be spawned.
        :param n_poses: Number of generated poses in a single run.
        :param repetitions: Number of times every ligand needs to be run.
        """
        self.settings = settings
        self.ligands = ligands
        self.process_limit = process_limit
        self.n_poses = n_poses
        self.repetitions = repetitions
    
    def run(self):
        """
        Runs the dockings.
        """
        processes_ligands = []
        processes = []
        output_dirs = []

        if self.ligands is None:
            ligand_check = True
        else:
            ligand_check = False
        for current_settings in self.settings:
            root = Path(current_settings["root"])

            # Transformation of the root and ligand_path to be sure it is system-compatible
            current_settings["root"] = str(root.absolute())
            current_settings["ligand_path"] = str(Path(current_settings["ligand_path"]).absolute())
            if ligand_check:
                self.ligands = self.ligands = [x for x in os.listdir(current_settings["ligand_path"]) if x.endswith(".mol2")]
           
            for ligand in self.ligands:
                # Adding the ligand to the current settings
                current_settings["ligand"] = ligand

                for repetition in range(self.repetitions):
                    out_dir = "_".join((ligand.split(".")[0], str(repetition)))
                    # Adding out_dir to the current settings
                    current_settings["out_dir"] = out_dir
                    if not os.path.exists(Path(root, out_dir)):
                        processes_ligands.append(ligand)
                        output_dirs.append(Path(root, out_dir))                        
                        processes.append(Popen(["python", str(Path(_location, "Docking_script.py").absolute())] + package_args(current_settings),
                                               shell=True, stdout=PIPE))
                        time.sleep(10)

                        if len(processes) >= 1:
                            print(f"Current number of processes: {len(processes_ligands)}\n" + "Current ligands: " + "\n".join(processes_ligands))
                        if len(processes) >= self.process_limit:
                            while True:
                                stdout = TextIOWrapper(processes[0].stdout, encoding="utf-8").read()
                                if "Done" in stdout:
                                    processes.pop(0)
                                    print(f"{processes_ligands.pop(0)} run finished")
                                    output_dirs.pop(0)
                                    break
                        else:
                            while len([x for x in os.listdir(output_dirs[0]) if x.startswith("gold_soln")]) < max((self.n_poses/self.process_limit*len(processes), 5)):
                                time.sleep(60)

def package_args(args_dict):
    """
    Takes a dictionary and creates a list where keys are transformed into command line flags, followed by the values.
    If the value is a list, the function unpacks it and separates the elements with appropriate flag.
    :param args_dict: Dictionary that should be transformed.
    :return: List ready to be used for construction of a shell command.
    """
    args_list = []
    for key, value in args_dict.items():
        if isinstance(value, list):
            for subvalue in value:
                args_list.append("-" + key)
                args_list.append(str(subvalue))
        else:
            args_list.append("-" + key)
            args_list.append(str(value))
    return args_list
