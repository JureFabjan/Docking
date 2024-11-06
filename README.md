# Docking

## Description
Collection of tools for docking and some extras.
The tool requires CCDC python API. Parts work with Schrödinger python API or are written in MOE SVL, though not having those will not make the package unusable.

## Project setup
First create the python package by running the command:
```
python -m build
```

Then move to `dist` subfolder and run:
```
pip install docking-0.0.2-py3-none-any.whl
```

Package requirements:
- bipython
- numpy
- pandas
- ccdc
- MOE for using the included SVL scripts
- optionally schrodinger

## Usage

### Basic docking

The main docking object is `Docking.Docking_main.Dock` class, which works on CCDC Gold interfacing. It either takes either the parameters for generating a new docking configuration, or an already existing configuration file. If new parameters are passed to it, it will overrite the value for these from the configuration. It also manages the preparation steps of the protocol. This means it prepares the protein for the docking and extracts any ligans present in it. Next it creates the coordinates of the binding site based on the specified siteing method. Finally, it will create a new configuration file for the current experiment to guarantee the reproducibility and run the experiment. Note, this all runs when invoking the class.

The second class in the same module is `Docking.Docking_main.Results`, which is meant for reading and processing the results of the docking run. We initiate it by passing to it the path to the configuration file of an already run experiment. Through it's functions the class can process the results by performing the following tasks:
- extracting the docking scores and assigned clusters of the poses
- saving the poses individually in protein-ligand complexes as individual pdb files with adjusted flexible side chains
- saving all the ligand poses into a single file
- extracting distances between prespecified ligand atoms and defined protein residues (through `Scripts/db_Distance.svl`, where the atoms should be specified as of now)
- extracting the positions of all the atoms for every pose in a spherical coordinate system, defined through 4 points (defining the center and 3 axes of the system through 4 atoms of the protein - by setting `extract_all_positions` to `True` and passing the coordinate system dictionary with `center`, `x_axis`, `y_axis` and `z_axis` entries; each of these should be a list [chain number, three-letter code of the AA, sequential number of the AA, atom code])

### Automatization and parallelization

The docking procedure can be automated and paralellized by using the docking master in `Docking.Docking_master.Master`. This class takes a list of settings dictionaries, with one dictionary per run we want to perform. The dictionaries contain:
- `root`: root folder in which to run an individual run
- `ligand_path`: path to the folder containing all the ligands that should be used
- `configuration`: path to the configuration file that should be used
- `template_protein`: path to the protein file that will be used in the docking
- `template_ligand`: path to the ligand file that is used to define the binding site
- optional entries are the lists for defining the coordinate system from above

Additionally the class accepts:
- a list of ligands if we don't wish to use all the ligands in the folders defined in the settings dictionaries,
- the maximum number of processes that should be spawned,
- number of poses that should be retained,
- number of times every ligand needs to be run (mainly intended to test the reproducibility of the experiment)

After defining the master, we run `master.run()` function, which runs the experiment through spawning new processes. These first run `Dock` class which performs the docking procedure, followed by the `Results` class which processes the output.
