from Docking import Docking_main
import os, sys
from pathlib import Path
from collections import defaultdict

def dock(root, configuration, template_protein, template_ligand, ligand_path, ligand, out_dir, coordinate_system=None):
    """
    Initializing the docking procedure.
    :param root: The root directory in which the docking should be done,
    :param configuration: Path to the configuration file.
    :param template_protein: Path to the template protein file.
    :param template_ligand: Path to the template ligand file.
    :param ligand_path: Path to the ligand.
    :param out_dir: Path to the output directory.
    :param coordinate_system: A dictionary defining the coordinate system axesa. Keys should be ('center', 'x_axis', 'y_axis', 'z_axis') and as values lists [chain_number, aa_code, residue_number, atom_name].
    """
    os.chdir(root)

    run = Docking_main.Dock(template_protein,
                            str(Path(ligand_path, ligand).absolute()),
                            template_ligand=template_ligand,
                            fitness_fun="chemscore", rescore_fun="asp",
                            site_radius=7,
                            output_dir=str(Path(root, out_dir).absolute()),
                            ndocks=200, autoscale=200, early_termination=False,
                            configuration=configuration, overwrite_protein=False)
    results = Docking_main.Results(run.settings.conf_file)
    results.save(save_complex=True, clean_complex=True, extract_all_positions=True, coordinate_system=coordinate_system)

    print("Done")

def arg_parser(arguments):
    """
    Parses and organizes the arguments from a sys.argv-like list.
    :param arguments: A sys.argv-like list of arguments.
    :return: A defaultdict that can be used for passing to dock() function.
    """
    arg_dict = defaultdict(list)
    skip = False    # A bool telling if the next iteration of the loop should be skipped over because the previous
                    # iteration contained a flag for a1
    for a1, a2 in zip(arguments, arguments[1:]):
        if not skip:
            if a1.startswith("-"):
                value = a2
                key = a1[1:]
                skip = True
            else:
                # The non-flagged arguments get stacked into "else" key in order of the appearance.
                value = a1
                key = "else"

            try:
                # If the second argument can be converted to intiger it should be.
                # In case the conversion returns ValueError we write the string.
                arg_dict[key].append(int(value))
            except ValueError:
                arg_dict[key].append(value)
                skip = True

        else:
            skip = False
    # Unpacking lists for keys where lists have only one item
    for key, value in arg_dict.items():
        if len(value) == 1:
            arg_dict[key] = value[0]
    return arg_dict

if __name__ == "__main__":
    arg_dict = arg_parser(sys.argv[1:])
    try:
        ligand_path = Path(arg_dict["ligand_path"])
    except TypeError as e:
        raise Exception(str(arg_dict))
    if "else" in arg_dict.keys():
        unflagged = arg_dict.pop("else")
    if "center" in arg_dict.keys():
        arg_dict["coordinate_system"] = {"center": arg_dict.pop("center"),
                                         "x_axis": arg_dict.pop("x_axis"),
                                         "y_axis": arg_dict.pop("y_axis"),
                                         "z_axis": arg_dict.pop("z_axis")}
    dock(**arg_dict)
