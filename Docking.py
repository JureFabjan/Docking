import os
from pathlib import Path

from ccdc.docking import Docker
from ccdc.io import MoleculeReader, MoleculeWriter, EntryWriter
from ccdc.protein import Protein
from pandas import DataFrame


class Dock:
    """
    The class wrapper for the data, used in the docking experiment.
    """
    def __init__(self, protein_file, target_ligand,
                 siteing_method="ligand", template_ligand="",
                 fitness_fun="goldscore", rescore_fun="chemscore",
                 output_dir="", ndocks=10, site_radius=12.0,
                 autoscale=10.0, early_termination=True, configuration="",
                 split_output=True, overwrite_protein=True):
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
        :param overwrite_protein: If configuration file is used and this parameter is true, theused  protein in the
        configuration file will be ignored.
        """
        # Initiating the docker and the settings
        if configuration:
            self.settings = Docker.Settings.from_file(configuration)
            self.docker = Docker(settings=self.settings)
        else:
            self.docker = Docker()
            self.settings = self.docker.settings

        # Setting directory and file names
        n = 0
        if not output_dir:
            current_dir = Path(".")
            existing = [str(x)[-2:] for x in current_dir.glob("Results*")]
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
        self.settings.fitness_function = fitness_fun
        self.settings.rescore_function = rescore_fun

        # Other settings
        # Speed of the computation (lower => faster)
        self.settings.autoscale = autoscale
        self.settings.early_termination = early_termination

        # Prepares the protein file and extracts the ligands present in it
        if configuration and not overwrite_protein:
            # If the configuration is present and the protein in it should be used
            self.protein_file = self.settings.protein_files[0]
            self.ligands = []
        else:
            # If there is no configuration or there is configuration but the protein in it should not be used
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

        # Changing the configuration file name; n was created when naming the output folder, so n should equal
        # to the final number of the fulder with the results
        self.results = self.docker.dock(file_name=f"api_gold_{n}.conf")

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


class Results:
    def __init__(self, settings_file):
        """
        Results operation class.
        :param settings_file: The .conf settings file, used to perform docking.
        """
        self.settings = Docker.Settings.from_file(settings_file)
        self.results = Docker.Results(self.settings)
        self.ligands = [x for x in self.results.ligands]

    def save(self, end_notation=True):

        # Collecting the scoring and clusters
        scores = self.ligand_score_extraction()
        clusters = self.clusters_extraction()

        # Saving the scores
        scores.to_csv(Path(self.settings.output_directory, "Ligand scores.csv"), index=False)

        # Reading the output file with the docked molecules and extracting the ligands themselves
        # The ligands are already ordered by score
        with MoleculeReader(self.settings.output_file) as docked_ligands:
            ligands = [ligand for ligand in docked_ligands]

        # Renaming of the ligands to accommodate the clusters
        for i, cluster in enumerate(clusters):
            for ligand in cluster:
                ligand = int(ligand)
                if end_notation:
                    # Appends the cluster number to the end of the identifier
                    ligands[ligand-1].identifier = ligands[ligand-1].identifier + f"|{i+1}"
                else:
                    # Original name: Template_name|Molecule|mol2|1|dockN
                    # Constructed name: Templ_name|Molecule|mol2|cluster|dockN
                    name = ligands[ligand-1].identifier.split("|")
                    ligands[ligand-1].identifier = "|".join(name[:-2] + [str(i+1)] + name[-1:])

        with MoleculeWriter(self.settings.output_file) as docked_ligands:
            for ligand in ligands:
                docked_ligands.write(ligand)

    def ligand_score_extraction(self):
        """
        Extracts the scores of fitness and rescore functions from the results.
        :return: Pandas DataFrame with ligand ID, fitness function result and rescore function results as columns.
        """
        return DataFrame({"Identifier": [ligand.identifier for ligand in self.ligands],
                          self.settings.fitness_function: [ligand.fitness(self.settings.fitness_function) for ligand in self.ligands],
                          self.settings.rescore_function: [ligand.fitness(self.settings.rescore_function) for ligand in self.ligands]})

    def clusters_extraction(self, treshold=3.0):
        """
        Extracts the clusters of docking poses from the ligand log file. Note! The numbers returned are not the
        numbers of docking poses, but rather the ranking of the score.
        :param treshold: The distance for the clusters to be used. The first distance lower than provided treshold
        will be taken from the results.
        :return: A list of lists of poses in the same clusters.
        """
        ligand_log = self.results.ligand_log(0)
        # Extracting the start of the table
        ligand_log = ligand_log[ligand_log.index("Distance | Clusters"):].split("\n")
        cluster_line = ""
        for line in ligand_log[::-1]:
            distance = 0
            try:
                distance = float(line.split("|")[0].strip())
            except ValueError:
                pass
            if distance and distance < treshold:
                cluster_line = line
                break
        if not cluster_line:
            raise ValueError("No cluster found with a distance lower than treshold.")

        distance, *clusters = cluster_line.split("|")
        return [cluster.strip().split() for cluster in clusters]

    def save_ligand_complexes(self):
        """
        Makes a complex of protein with each pose of the ligand in the results with the side chains adjustments and
        saves them as individual .pdb files in a subfolder Complexes inside the output directory. The files are named
        Pose_XXX.pdb, where XXX denotes the number of the pose.
        :return:
        """
        # Create base output directory
        output = Path(self.settings.output_directory, "Complexes")
        os.mkdir(output)
        # Loop though all the ligands
        for ligand in self.ligands:
            # Extract the ligand docking number from its identifier
            n = int(ligand.identifier.split("dock")[-1].split("|")[0])
            # Save
            with EntryWriter(str(Path(output, f"Pose_{n:03}.pdb"))) as protein_writer:
                protein_writer.write(self.results.make_complex(ligand))


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
                 split_output=False,
                 overwrite_protein=False)

    _results = Results(_dock.settings.conf_file)
    _results.save()
