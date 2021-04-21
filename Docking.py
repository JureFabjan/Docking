import inspect
import os
import re
from pathlib import Path
from subprocess import run

from Bio.PDB import PDBIO
from Bio.PDB import PDBParser
from ccdc.docking import Docker
from ccdc.io import MoleculeReader, MoleculeWriter, EntryWriter
from ccdc.protein import Protein
from pandas import DataFrame, read_csv, merge
from shutil import rmtree

# Just getting the script location
_location = str(Path(inspect.getfile(inspect.currentframe())).parents[0].absolute())


class Dock:
    """
    The class wrapper for the data, used in the docking experiment.
    """
    def __init__(self, protein_file, target_ligand,
                 siteing_method="ligand", template_ligand="",
                 fitness_fun="goldscore", rescore_fun="chemscore",
                 output_dir="", ndocks=10, site_radius=12.0,
                 autoscale=10.0, early_termination=True, configuration="",
                 split_output=True, overwrite_protein=True, output_format="mol2"):
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
        :param output_format: The desired format of the output files. Mol2 advised.
        """
        # Initiating the docker and the settings
        if configuration:
            self.docker = Docker(settings=Docker.Settings.from_file(configuration))
        else:
            self.docker = Docker()
        self.settings = self.docker.settings

        # Setting directory and file names
        n = 0
        if not output_dir:
            current_dir = Path(".")
            n = max((int(str(x)[-2:]) for x in current_dir.glob("Results*"))) + 1
            output_dir = Path(current_dir, f"Results_{n:02}")
            os.mkdir(output_dir)
            self.settings.output_directory = str(output_dir.absolute())
        else:
            if not os.path.isdir(output_dir):
                os.mkdir(output_dir)
            self.settings.output_directory = output_dir
            # Setting n:
            #   if the name of the dir ends with n, then that number is used
            numbers = re.findall(r"\d+\s*$", output_dir)
            if numbers:
                n = int(numbers[0])

        # Setting the file names
        self.settings.output_format = output_format
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
        self.protein_file = protein_file
        if configuration and not overwrite_protein:
            # If the configuration is present and the protein in it should be used
            self.protein_file = self.settings.protein_files[0].file_name
            # If protein file already exists, it is assumed it is already clean. Else the file gets created.
            if not Path(self.protein_file).exists():
                self.ligands = self.prepare_protein(output=self.protein_file)
            else:
                self.ligands = []
        else:
            # If there is no configuration or there is configuration but the protein in it should not be used
            if configuration:
                # If the protein should not be overwritten, then the correspondence of the path is checked and adjusted.
                if Path(self.protein_file).absolute() != Path(protein_file).absolute():
                    self.protein_file = str(Path(protein_file).absolute())
            else:
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
        # to the final number of the folder with the results
        self.results = self.docker.dock(file_name="gold_{}_{}.conf".format(target_ligand.split("\\")[-1].split("/")[-1].split(".")[0],
                                                                           n))

    def prepare_protein(self, output=""):
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

        if output:
            clean_protein_file = output
        else:
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

    def save(self, end_notation=True, cluster_threshold=3.0, save_complex=False, clean_complex=False, extract_distances=False,
             extract_positions=False, extract_all_positions=False):
        """
        Saves the scores of the docking poses, checks the clustering results and renames the ligands in the
        results .mol2 file to include the clusters. If specified, the ligand-protein complexes of all ligand poses
        will be saved in individual .pdb files.
        :param extract_distances: Boolean specifying if the protein-ligand distances should be extracted. Note this
        will be performed only in case save_complex is True too.
        :param end_notation: If the ligands should be renamed so that the cluster number is included at the end. Else
        a number before the docking number is exchanged for the cluster number.
        :param cluster_threshold: The threshold for cluster extraction in angstroms.
        :param save_complex: True if the ligand-protein complexes should also be generated.
        :param clean_complex: Clean the complexes from doubles of chains.
        :param extract_positions: Extract the positions of ligands in the relative coordinate system.
        :param extract_all_positions: Extract the positions of all ligand atoms in the relative coordinate system.
        :return:
        """
        # Collecting the scoring and clusters
        try:
            # Trying to extrtact the scores; this works with the API only the first time
            scores = self.ligand_score_extraction()
            scores.to_csv(Path(self.settings.output_directory, "Ligand scores.csv"), index=False)
        except RuntimeError as _:
            # If there is no score available, the scores try to be extracted from the possibly saved csv
            scores = read_csv(Path(self.settings.output_directory, "Ligand scores.csv"))
            scores = scores[scores.columns[0:3]]
    
        clusters = self.clusters_extraction(threshold=cluster_threshold)

        # Adding the clusters into the scores DataFrame
        scores["Cluster"] = 1
        # Making sure the scores are ordered by the second column (fitness function score)
        scores.sort_values(by=self.settings.fitness_function, axis=0, ascending=False, inplace=True, ignore_index=True)
        for cluster_index, indices in enumerate(clusters):
            scores.loc[[int(i)-1 for i in indices], "Cluster"] = cluster_index+1

        # Reading the output file with the docked molecules and extracting the ligands themselves
        # The ligands are already ordered by score
        if self.settings.output_file:
            # In case the output file was not split
            with MoleculeReader(self.settings.output_file) as docked_ligands:
                ligands = [ligand for ligand in docked_ligands]
        else:
            # In case there is no single output file in the docking
            ligand_list = [fname for fname in os.listdir(self.settings.output_directory) if fname.startswith("gold_soln_")]
            ligands = [MoleculeReader(str(Path(self.settings.output_directory,
                                               ligand).absolute()))[0] for ligand in ligand_list]

        # Renaming of the ligands to accommodate the clusters
        for i, cluster in enumerate(clusters):
            for ligand in cluster:
                ligand = int(ligand)
                if end_notation:
                    # Appends the cluster number to the end of the identifier
                    name = ligands[ligand-1].identifier + f"|{i+1}"
                else:
                    # Original name: Template_name|Molecule|mol2|1|dockN
                    # Constructed name: Templ_name|Molecule|mol2|cluster|dockN
                    name = ligands[ligand-1].identifier.split("|")
                    name = "|".join(name[:-2] + [str(i+1)] + name[-1:])
                ligands[ligand-1].identifier = name

        if self.settings.output_file:
            # Saving ligands into a single file, if the output_file is defined
            with MoleculeWriter(self.settings.output_file) as docked_ligands:
                for ligand in ligands:
                    docked_ligands.write(ligand)
        else:
            # Ligands were initially saved into separate files, so we save them into appropriate files again.
            # A merged file is also created and filled with the ligands
            with MoleculeWriter(str(Path(self.settings.output_directory, "Merged.mol2").absolute())) as docked_ligands:
                for fname, ligand in zip(ligand_list, ligands):
                    docked_ligands.write(ligand)
                    with MoleculeWriter(str(Path(self.settings.output_directory, fname).absolute())) as current_ligand:
                        current_ligand.write(ligand)

        if save_complex:
            self.save_ligand_complexes(clean_complex)

            # Adding complexes file names to the scoring data frame
            scores["Complex"] = scores["Identifier"].str.split("dock").str[-1].apply(lambda x: f"Pose_{int(x):03}")

            if extract_distances or extract_positions or extract_all_positions:
                # Making the MOE database from the extracted complexes
                self.moe_complex_import()

            if extract_distances:
                # Extracting the distances between the ligands and protein from the database
                self.moe_distance_extract()

                # Opening the distances file
                distances = read_csv(Path(self.settings.output_directory, "Complexes", "Results.txt"), sep="\t")
                distances.rename({"File": "Complex"}, axis=1, inplace=True)
                scores = merge(scores, distances, on="Complex")

            if extract_positions:
                # Extracting the ligand positions
                self.moe_position_extract()

                # Opening the positions file
                positions = read_csv(Path(self.settings.output_directory, "Complexes", "Results.txt"), sep="\t")
                positions.rename({"File": "Complex"}, axis=1, inplace=True)
                scores = merge(scores, positions, on="Complex")

            if extract_all_positions:
                # Extracting the positions of all ligand atoms
                self.moe_all_positions_extract()

                # Opening the positions file
                positions = read_csv(Path(self.settings.output_directory, "Complexes", "Results.txt"), sep="\t")
                positions.rename({"File": "Complex"}, axis=1, inplace=True)
                scores = merge(scores, positions, on="Complex") 

        # Saving the scores. This is done last because of potential additions to the scores files (i.e. file names,
        # distances etc.)
        scores.to_csv(Path(self.settings.output_directory, "Ligand scores.csv"), index=False)

    def ligand_score_extraction(self):
        """
        Extracts the scores of fitness and rescore functions from the results.
        :return: Pandas DataFrame with ligand ID, fitness function result and rescore function results as columns.
        """
        return DataFrame({"Identifier": [ligand.identifier for ligand in self.ligands],
                          self.settings.fitness_function: [ligand.fitness(self.settings.fitness_function) for ligand in self.ligands],
                          self.settings.rescore_function: [ligand.fitness(self.settings.rescore_function) for ligand in self.ligands]})

    def clusters_extraction(self, threshold=3.0):
        """
        Extracts the clusters of docking poses from the ligand log file. Note! The numbers returned are not the
        numbers of docking poses, but rather the ranking of the score.
        :param threshold: The distance for the clusters to be used. The first distance lower than provided threshold
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
            if distance and distance < threshold:
                cluster_line = line
                break
        if not cluster_line:
            raise ValueError("No cluster found with a distance lower than threshold.")

        distance, *clusters = cluster_line.split("|")
        return [cluster.strip().split() for cluster in clusters]

    def save_ligand_complexes(self, clean_complex):
        """
        Makes a complex of protein with each pose of the ligand in the results with the side chains adjustments and
        saves them as individual .pdb files in a subfolder Complexes inside the output directory. The files are named
        Pose_XXX.pdb, where XXX denotes the number of the pose.
        :param clean_complex: Boolean specifying if the complexes should be cleaned from chain doubles.
        :return:
        """
        # Create base output directory
        output = Path(self.settings.output_directory, "Complexes")
        if os.path.exists(output):
            rmtree(output)
        os.mkdir(output)
        # Loop though all the ligands
        for ligand in self.ligands:
            # Extract the ligand docking number from its identifier
            n = int(ligand.identifier.split("dock")[-1].split("|")[0])
            # Save
            with EntryWriter(str(Path(output, f"Pose_{n:03}.pdb"))) as protein_writer:
                complex = self.results.make_complex(ligand)
                complex.remove_unknown_atoms()
                protein_writer.write(complex)
        
        if clean_complex:
            # The PDBParser cannot parse residues with the same numbering, so it just retains the correct ones.
            # This means that to clean the complexes we just read the files and save them.
            parser = PDBParser(PERMISSIVE=1)
            saver = PDBIO()
            for n in range(1, len(self.ligands)+1):
                structure_path = str(Path(output, f"Pose_{n:03}.pdb"))
                structure = parser.get_structure(f"Pose_{n:03}", structure_path)
                saver.set_structure(structure)
                saver.save(structure_path)

    def moe_complex_import(self):
        """
        In essence the function activates an MOE script, which imports all the .pdb files in a folder into a
        database called 'database.mdb', which is created in the same folder. In the context of the results,
        this function creates a database with all the generated complex files.
        :return:
        """
        script = Path(_location, "db_Import.svl")

        # Preparing the script by including the output path
        self.adjust_script(script, 10)

        run(["moebatch", "-run", str(script.absolute())], shell=True)

    def moe_distance_extract(self):
        """
        Runs a script, which extracts the specific protein-ligand distances from a MOE database.
        Currently the distances are specified in the script, but this will be more flexible in the future.
        :return:
        """
        script = Path(_location, "db_Distance.svl")

        # Preparing the script by including the output path
        self.adjust_script(script, 11)

        run(["moebatch", "-run", str(script.absolute())], shell=True)

    def moe_position_extract(self):
        """
        Runs a script, which extracts the protein-ligand distances and angles in the defined coordinate
        system. Works on an MOE database.
        Currently the distances are specified in the script, but this will be made to be more flexible in the future.
        :return:
        """
        script = Path(_location, "db_Position.svl")

        # Preparing the script by including the output path
        self.adjust_script(script, 11)

        run(["moebatch", "-run", str(script.absolute())], shell=True)

    def moe_all_positions_extract(self):
        """
        Runs a script, which extracts the protein-ligand distances and angles in the defined coordinate
        system for all ligand atoms. Works on an MOE database.
        :return:
        """

        script = Path(_location, "db_AllAtomPosition.svl")

        self.adjust_script(script, 11)

        run(["moebatch", "-run", str(script.absolute())], shell=True)

    def adjust_script(self, script_file, line):
        """
        Abstraction of the script adjustment - it adjusts the path in the specified script, on a specified line.
        :param script_file: Path object, pointing to the script.
        :param line: Integer, pointing to the line, which needs to be changed (numbering starts from 0).
        :return:
        """
        output = Path(self.settings.output_directory, "Complexes")
        with open(script_file, "r") as script_file_handle:
            script_file_text = script_file_handle.read().split("\n")
        line_parts = script_file_text[line].split("\'")
        script_file_text[line] = "\'".join([line_parts[0], str(output.absolute()).replace(os.sep, "/"),
                                            line_parts[2]])
        with open(script_file, "w") as script_file_handle:
            script_file_handle.write("\n".join(script_file_text))

    def adjust_coordinate_acids(self, script, center=[], x=[], y=[], z=[]):
        """
        Changes the amino acids used for the coordinate system to the defined values. Just the specified positions
        will be changed, the rest will be kept as is. The coordinates should be given in a list of format
        [Chain, Name, N, Atom], where Chain is the intiger of chain number, Name is a 3-letter string name of the
        amino acid (all in capital), N is the intiger delineating the amino acid's position in the chain, and Atom
        is a string with the name of the amino acid's atom used.
        :param script: Name of the script (together with the '.svl' ending)
        :param center: Definition of the central atom
        :param x: Definition of the atom delineating together with the center the x-axis
        :param y: Definition of the atom delineating together with the center the y-axis
        :param z: Definition of the atom delineating together with the center the z-axis
        """
        with open(Path(_location, script), "r") as file:
            script_text = file.read()
        if len(center) == 4:
            script_text = re.sub(r"local aCenter = \[\d+,.+,\s\d+,.+\]",
                                 f"local aCenter = {str(center)}", script_text)
        if len(x) == 4:
            script_text = re.sub(r"local aX_axis = \[\d+,.+,\s\d+,.+\]",
                                 f"local aX_axis = {str(center)}", script_text)
        if len(y) == 4:
            script_text = re.sub(r"local aY_axis = \[\d+,.+,\s\d+,.+\]",
                                 f"local aY_axis = {str(center)}", script_text)
        if len(z) == 4:
            script_text = re.sub(r"local aZ_axis = \[\d+,.+,\s\d+,.+\]",
                                 f"local aZ_axis = {str(center)}", script_text)
        with open(Path(_location, script), "w") as file:
            file.write(script_text)


if __name__ == "__main__":
    os.chdir(Path(".", "a1g2"))
    _protein_file = "Source_protein.mol2"
    _protein_ligand_file = "Source_molecule.mol2"
    _ligand_file = "DCBSPU19.mol2"
    _dock = Dock(_protein_file,
                 _ligand_file,
                 template_ligand=_protein_ligand_file,
                 ndocks=10,
                 autoscale=10,
                 # configuration="gold.conf",
                 early_termination=False,
                 split_output=False)
                 # overwrite_protein=True)

    _results = Results(_dock.settings.conf_file)
    _results.save(save_complex=True,
                  clean_complex=True,
                  extract_distances=True,
                  extract_positions=True)
