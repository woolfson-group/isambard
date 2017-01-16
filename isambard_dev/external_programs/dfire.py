import os
import subprocess
import tempfile

from settings import global_settings
from tools.isambard_warnings import check_availability


def test_dfire():
    """Test to see if dfire is available using the path to the binary in global settings."""
    is_dfire_available = False
    if os.path.isfile(global_settings['dfire']['path']):
        try:
            subprocess.check_output([global_settings['dfire']['path']], stderr=subprocess.DEVNULL)
        except subprocess.CalledProcessError:
            is_dfire_available = True
    return is_dfire_available


@check_availability('dfire', test_dfire, global_settings)
def run_dfire(pdb, path=True):
    """Run dfire with a path to a pdb file or the contents of a PDB file as a string.

    Parameters
    ----------
    pdb : str
        PDB file or the contents of a PDB file as a string
    path : bool
        Indicates whether the `pdb` is a path or PDB format string.

    Returns
    -------
    dfire_std_our : str
        Output from dfire as a unicode string.
    """
    if path:
        with open(pdb, 'r') as inf:
            pdb = inf.read()
    pdb = pdb.encode()
    try:
        with tempfile.NamedTemporaryFile(delete=False) as dfire_tmp:
            dfire_tmp.write(pdb)
            dfire_tmp.seek(0)  # Resets the buffer back to the first line
            dfire_std_out = subprocess.check_output(
                [global_settings['dfire']['path'], global_settings['dfire']['lib'], dfire_tmp.name])
    finally:
        os.remove(dfire_tmp.name)
    return dfire_std_out.decode()


def parse_dfire_out(dfire_std_out):
    """Parses the an output string produced by dfire.

    Parameters
    ----------
    dfire_std_out : str
        Output from dfire as a unicode string.

    Returns
    -------
    dfire_score : float
        The dfire score as a float.
    """
    return float(dfire_std_out.split()[1])


def calculate_dfire_score(pdb, path=True):
    """Calculates the dfire score of the pdb file.

    Parameters
    ----------
    pdb : str
        PDB file or the contents of a PDB file as a string
    path : bool
        Indicates whether the `pdb` is a path or PDB format string.

    Returns
    -------
    dfire_score : float
        The dfire score as a float.
    """
    dfire_out = run_dfire(pdb, path=path)
    return parse_dfire_out(dfire_out)
