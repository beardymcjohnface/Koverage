import os
import subprocess


def profile_self(output_file):
    """Run py-spy on self

    Args:
        output_file (str): output SVG file for saving
    """
    pid = str(os.getpid())
    command = ["py-spy", "record", "-s", "-o", output_file, "--pid", pid]
    subprocess.Popen(command, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
