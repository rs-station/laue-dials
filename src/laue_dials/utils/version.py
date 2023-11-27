# Get package versioning information
from importlib.metadata import version
import os

def laue_version():
    try:
        width = os.get_terminal_size().columns 
    except:
        width = 65
    dials_version = version('dials')
    laue_dials_version = version('laue-dials')

    print('-' * width)
    print(f"DIALS version " + dials_version + ".")
    print(f"laue-dials version " + laue_dials_version + ".")
    print('-' * width)
