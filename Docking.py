import os
from pathlib import Path

from ccdc.docking import Docker
from ccdc.io import MoleculeReader, EntryWriter
from ccdc.protein import Protein


class Dock:
    """
    The class wrapper for the data, used in the docking experiment.
    """
    def __init__(self, protein_file, target_ligand,
                 siteing_method="ligand", template_ligand="",
                 fitness_fun="goldscore", rescore_fun="chemscore",
                 output_dir="", ndocks=10, site_radius=12.0,
                 autoscale=10.0, early_termination=True, configuration="",
                 split_output=True):
        """

        :param protein_file: The file with the template protein. If it already contains ligands, they will be removed.
        The preferred file type is .mol2.
        :param target_ligand: The file with the ligand to be docked. The preferred file type is .mol2.
        :param siteing_method: Which method should be used for defining the binding site. The options are:
        'ligand': the ligand from the template_ligand or extracted from template protein file is used.
        :param template_ligand: The file with the ligand, which is used for defining the binding site. Preferred
        file type is .mol2.
        :param fitness_fun: Scoring function
        :param rescore_fun: Rescoring function
        :param output_dir: The directory to which the results will be saved. If left empty, 'Results_XX' gets created
        in the working directory.
        :param ndocks: Number of docking runs to be done.
        :param site_radius: The radius tu be used when defining the binding site.
        :param autoscale: The speed of the docking - 100 is the slowest, while numbers close to 0 will result in faster
        runs.
        :param early_termination: True is the computation can be terminated early.
        :param configuration: gold.conf file, which can be imported. Mainly intended when applying soft potentials and
        flexible side chains.
        :param split_output: If true, output docking poses will be saved in separate files.
        """
        # Initiating the docker and the settings
        if configuration:
            self.settings = Docker.Settings.from_file(configuration)
            self.docker = Docker(settings=self.settings)
        else:
            self.docker = Docker()
            self.settings = self.docker.settings

        # Setting directory and file names
        if not output_dir:
            current_dir = Path(".")
            existing = [str(x)[-2:] for x in current_dir.glob("Results*")]
            n = 0
            while f"{n:02}" in existing:
                n += 1
            output_dir = Path(current_dir, f"Results_{n:02}")
            os.mkdir(output_dir)
            self.settings.output_directory = str(output_dir.absolute())
        else:
            if not os.path.isdir(output_dir):
                os.mkdir(output_dir)
            self.settings.output_directory = output_dir

        # Adding the settings file to the output directory
        self.settings.output_format = str(Path(output_dir, "gold.conf").absolute)

        # Setting the file names
        self.settings.output_format = "mol2"
        if split_output:
            self.settings.output_file = ""
        else:
            self.settings.output_file = "Docking_results.mol2"

        # Checks if both scoring functions are the same (should not be).
        if fitness_fun == rescore_fun:
            if fitness_fun == "chemscore":
                rescore_fun = "goldcore"
            else:
                rescore_fun = "chemscore"
            print(f"Exception: Rescore function equals fitness function. Changing rescore function to {rescore_fun}.")
        self.f_fitness = fitness_fun
        self.f_rescore = rescore_fun

        # Other settings
        # Speed of the computation (lower => faster)
        self.settings.autoscale = autoscale
        self.settings.early_termination = early_termination

        # Prepares the protein file and extracts the ligands present in it
        self.settings.clear_protein_files()
        self.protein_file = protein_file
        self.ligands = self.prepare_protein()

        # Creating the coordinates of the binding site
        self.template_ligand = template_ligand
        self.site_radius = site_radius
        {
            "ligand": self.reference_ligand
         }[siteing_method]()

        # Setting the target ligand
        self.settings.clear_ligand_files()
        self.settings.add_ligand_file(target_ligand, ndocks)

        self.results = self.docker.dock()

    def prepare_protein(self):
        """
        Prepares the protein structure for docking.
        :return:
        """
        protein = Protein.from_file(self.protein_file)
        protein.remove_all_waters()
        protein.remove_unknown_atoms()
        protein.add_hydrogens()

        ligands = protein.ligands
        for p_ligand in ligands:
            protein.remove_ligand(p_ligand.identifier)

        clean_protein_file = os.path.join(
            self.settings.output_directory, f"{protein.identifier}_clean.mol2"
        )

        with EntryWriter(clean_protein_file) as writer:
            writer.write(protein)
        self.settings.add_protein_file(clean_protein_file)
        
        return ligands
    
    def reference_ligand(self):
        """
        Takes the ligand, which defines the binding site, into the settings.
        :return:
        """
        if not self.template_ligand:
            if not self.ligands:
                raise RuntimeError("No ligand file given and no ligands in the protein.")
            self.template_ligand = os.path.join(self.settings.output_directory, "protein_ligands.mol2")
            with EntryWriter(self.template_ligand) as writer:
                for p_ligand in self.ligands:
                    writer.write(p_ligand)
        self.settings.binding_site = self.settings.BindingSiteFromLigand(self.settings.proteins[0],
                                                                         MoleculeReader(self.template_ligand)[0],
                                                                         self.site_radius)


if __name__ == "__main__":
    _protein_file = "Protein.mol2"
    _protein_ligand_file = "Ligand.mol2"
    _ligand_file = "Molecule.mol2"
    _dock = Dock(_protein_file,
                 _ligand_file,
                 template_ligand=_protein_ligand_file,
                 ndocks=100,
                 autoscale=100,
                 configuration="api_gold_UI.conf",
                 split_output=False)
