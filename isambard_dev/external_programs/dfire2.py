import os
import subprocess
import tempfile

from optimisation.optimizer import BaseScore, OptCMAES, OptDE, OptGA,  OptPSO
from settings import global_settings
from tools.isambard_warnings import check_availability


def test_dfire2():
    """Test to see if dfire2 is available using the path to the binary in global settings."""
    is_dfire2_available = False
    if os.path.isfile(global_settings['dfire2']['path']):
        try:
            subprocess.check_output([global_settings['dfire2']['path']], stderr=subprocess.DEVNULL)
        except subprocess.CalledProcessError:
            is_dfire2_available = True
    return is_dfire2_available


@check_availability('dfire2', test_dfire2, global_settings)
def run_dfire2(pdb, path=True):
    """Run dfire2 with a path to a pdb file or the contents of a PDB file as a string.

    Parameters
    ----------
    pdb : str
        PDB file or the contents of a PDB file as a string
    path : bool
        Indicates whether the `pdb` is a path or PDB format string.

    Returns
    -------
    dfire2_std_our : str
        Output from dfire2 as a unicode string.
    """
    if path:
        with open(pdb, 'r') as inf:
            pdb = inf.read()
    pdb = pdb.encode()
    try:
        with tempfile.NamedTemporaryFile(delete=False) as dfire2_tmp:
            dfire2_tmp.write(pdb)
            dfire2_tmp.seek(0)  # Resets the buffer back to the first line
            dfire2_std_out = subprocess.check_output(
                [global_settings['dfire2']['path'], global_settings['dfire2']['lib'], dfire2_tmp.name])
    finally:
        os.remove(dfire2_tmp.name)
    return dfire2_std_out.decode()


def parse_dfire2_out(dfire2_std_out):
    """Parses the an output string produced by dfire2.

    Parameters
    ----------
    dfire2_std_out : str
        Output from dfire2 as a unicode string.

    Returns
    -------
    dfire2_score : float
        The dfire2 score as a float.
    """
    return float(dfire2_std_out.split()[1])


def calculate_dfire2_score(pdb, path=True):
    """Calculates the dfire2 score of the pdb file.

    Parameters
    ----------
    pdb : str
        PDB file or the contents of a PDB file as a string
    path : bool
        Indicates whether the `pdb` is a path or PDB format string.

    Returns
    -------
    dfire2_score : float
        The dfire2 score as a float.
    """
    dfire2_out = run_dfire2(pdb, path=path)
    return parse_dfire2_out(dfire2_out)


# Dfire2 Optimiser


def dfire2_eval(params):
    """Builds and evaluates Dfire2 energy of model

    Parameters
    ----------
    params: list
        Tuple containing the specification to be built, the sequence, and the parameters for model building.

    Returns
    -------
    model.bude_score: float
        BUFF score for model to be assigned to particle fitness value.
    """
    specification, sequence, parsed_ind = params
    model = specification(*parsed_ind)
    model.build()
    model.pack_new_sequences(sequence)
    return calculate_dfire2_score(model.pdb, path=False)


class Dfire2Score(BaseScore):
    """
    Assigns Dfire2 score as fitness to individuals in optimization
    """
    evaluation_function = staticmethod(dfire2_eval)


class CMAES_Dfire2_Opt(OptDE, Dfire2Score):
    """
    Class for DE algorithm optimizing Dfire2 fitness
    """
    def __init__(self, specification, **kwargs):
        super().__init__(**kwargs)
        self._params['specification'] = specification


class DE_Dfire2_Opt(OptCMAES, Dfire2Score):
    """
    Class for DE algorithm optimizing Dfire2 fitness
    """
    def __init__(self, specification, **kwargs):
        super().__init__(**kwargs)
        self._params['specification'] = specification


class GA_Dfire2_Opt(OptGA, Dfire2Score):
    """
    Class for DE algorithm optimizing Dfire2 fitness
    """
    def __init__(self, specification, **kwargs):
        super().__init__(**kwargs)
        self._params['specification'] = specification


class PSO_Dfire2_Opt(OptPSO, Dfire2Score):
    """
    Class for DE algorithm optimizing Dfire2 fitness
    """
    def __init__(self, specification, **kwargs):
        super().__init__(**kwargs)
        self._params['specification'] = specification



