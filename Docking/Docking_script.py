from . import Docking
import os, sys, time
from pathlib import Path

def dock(root, configuration, template_protein, template_ligand, ligand_path, ligand, out_dir):
    """
    Initializing the docking procedure.
    :param root: The root directory in which the docking should be done,
    :param configuration: Path to the configuration file.
    :param template_protein: Path to the template protein file.
    :param template_ligand: Path to the template ligand file.
    :param ligand_path: Path to the ligand.
    :param out_dir: Path to the output directory.
    """
    os.chdir(root)

    run = Docking.Dock(template_protein,
                        str(Path(ligand_path, ligand).absolute()),
                        template_ligand=template_ligand,
                        fitness_fun="chemscore", rescore_fun="asp",
                        site_radius=7,
                        output_dir=str(Path(root, out_dir).absolute()),
                        ndocks=200, autoscale=200, early_termination=False,
                        configuration=configuration, overwrite_protein=False)
    results = Docking.Results(run.settings.conf_file)
    results.save(save_complex=True, clean_complex=True, extract_all_positions=True)

    print("Done")

if __name__ == "__main__":
    # Arguments are:
    #   0. Path to execution
    #   1. Configuration file name (presumably in the execution path)
    #   2. Protein file
    #   3. Template ligand file
    #   4. Path to ligands
    #   5. Name of the ligand, which should be run
    #   6. Output directory

    arguments = sys.argv[1:]

    root, configuration, template_protein, template_ligand, ligand_path, ligand, out_dir, *_ = arguments
    ligand_path = Path(ligand_path)
    dock(root, configuration, template_protein, template_ligand, ligand_path, ligand, out_dir)