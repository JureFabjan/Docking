"""
Module for streamlining the download of data from SRA.

Usage:
The tool works through STAtoolkit. If there is no binary downloaded on
the used machine first run sratoolkit_install(). This will install the
toolkit into a subfolder of the current working directory.

Once the SRAtoolkit is available, just run download function supplying
a list of SRR codes to be downloaded.

Possible command line arguments:
-c: single run code to be downloaded
-f: path to the folder in which the data should be put
-t: path to temporary storage of the data before it is transfered
"""
import sys
import subprocess
import os
from pathlib import Path
import tempfile
from time import sleep
from datetime import datetime
from collections import defaultdict

class Command:
    """
    Wrapper for subprocess.Popen.
    Mainly written to enable easier switching between parallel and
    non-parallel commands.
    """
    def __init__(self, command_list, parallel=True):
        """
        Initalization of the class instance. It already opens a command instance.
        
        Args:
            command_list (list): A list of strings for composing the total command.
            See the Popen reference for further detail.
            parallel (bool, optional): Should the command run in parallel to further code. If True,
            class initialization will not wait for the command to complete. Defaults to True.
        """
        self._subprocess = subprocess.Popen(command_list,
                                            stdout=subprocess.PIPE,
                                            stderr=subprocess.PIPE)
        
        self.stderr = ""
        self.stdout = ""
        if not parallel:
            self.stdout, self.stderr = self._subprocess.communicate()
        self.parallel = parallel
        
    def state(self):
        """
        Returns the state of the job. If None, the job is still running, else a return code
        is returned.

        Returns:
            int or None: Return code or None if the code is not available.
        """
        return self._subprocess.poll()
    
    def get_returns(self):
        """
        Fetches the stdout and stderror of the job.

        Returns:
            tuple: Tuple with stdout and stderr values as strings.
        """
        if not (self.stderr or self.stdout):
            self.stdout, self.stderr = self._subprocess.stdout.readlines(), self._subprocess.stderr.readlines()
        return self.stdout, self.stderr


def sratoolkit_install():
    """
    Function for setting up SRAtoolkit in the local folder.
    It downloads the latest version of the toolkit, untars it and starts the interactive configuration.
    Does not return anything. 
    """
    Command(["wget", "--output-document", "sratoolkit.tar.gz", "http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz"],
            parallel=False)
    Command(["tar", "-zxvf", "sratoolkit.tar.gz"],
            parallel=False)
    os.remove("sratoolkit.tar.gz")
    folder_name = [x for x in os.listdir() if os.path.isdir(x) and x.startswith("sratoolkit")][0]
    os.system(f"{folder_name}/bin/vdb-config -i")

def download(codes, final_path="", temp_path="", toolkit_path=""):
    """
    Function for downloading the data from SRA.

    Args:
        codes (list): List of strings containing SRR codes of runs that should be downloaded.
        final_path (str, optional): String containing path pointing to where the data should be downloaded.
        If it is not provided or an empty string, the files will be downloaded to the current folder.
        Defaults to "".
        temp_path (str, optional): Path to the temporary file. If left empty defaults to the current working directory.
        Defaults to "".
        toolkit_path (str, optional):String containing path to the SRAtoolkit. If not provided it searches 
        for the toolkit in the current working folder. 
        Defaults to "".

    Returns:
        list: A list of completion codes for individual provided codes.
    """
    if temp_path:
        temp_path = Path(temp_path)
    else:
        temp_path = Path(".").absolute()
    if final_path:
        final_path = Path(final_path)
    else:
        final_path = Path(".").absolute()
    if not toolkit_path:
        toolkit_path = [x for x in os.listdir() if os.path.isdir(x) and x.startswith("sratoolkit")][0]
    running = []
    for code in codes:
        running.append(Command([toolkit_path + "/bin/fasterq-dump", "--split-files", code, "-O", str(temp_path)],
                                parallel=True))
    completion_codes = [x.state() for x in running]
    t_beginning = datetime.now()
    while None in completion_codes:
        sleep(10)
        print(f"Running time: {str(datetime.now() - t_beginning)}")
        completion_codes = [x.state() for x in running]
        print("Completion codes:\n", "\n".join([f"{x}: {y}" for x, y in zip(codes, completion_codes)]))
    if final_path != temp_path:
        for file in os.listdir(temp_path):
            print("Transferring the data")
            transfer = Command(["pipe", "storage", "mv", temp_path, "-rf"],
                               parallel=False)
            print(transfer.get_returns()[0])
    return completion_codes

def arg_parser(arguments):
    """
    Takes a list of arguments passed from the command line and
    collects them into a dictionary.

    Args:
        arguments (list): List of arguments ordered as they were specified in the command line.

    Returns:
        dict: Dictionary where the keys are the command line arguments and the values are lists of values. 
    """
    arguments = " ".join(arguments).split("-")
    arg_dict = defaultdict(list)
    for arg in arguments:
        code, value = arg.split(" ")
        arg_dict[code].append(value)
    return arg_dict

def cleanup():
    """
    Removal of the SRAtoolkit.
    """
    target_folder = [x for x in os.listdir() if os.path.isdir(x) and x.startswith("sratoolkit")]
    if target_folder:
        for folder_path in target_folder:
            os.rmdir(folder_path)
    
        
if __name__ == "__main__":
    codes = arg_parser(sys.argv[1:])
    sratoolkit_install()
    for x in ("p", "t"):
        if not codes[x]:
            codes[x] = ""
    download(codes["c"],
             final_path=codes["p"],
             temp_path=codes["t"])
    cleanup()
