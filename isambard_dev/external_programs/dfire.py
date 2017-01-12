import os
import subprocess
import tempfile
import warnings

from settings import global_settings
from tools.isambard_warnings import DependencyNotFoundWarning


def check_dfire_avail():
    is_dfire_available = False
    if os.path.isfile(global_settings['dfire']['path']):
        try:
            subprocess.check_output([global_settings['dfire']['path']], stderr=subprocess.DEVNULL)
        except subprocess.CalledProcessError:
            is_dfire_available = True
    else:
        warning_string = ('\n\nDfire not found and so cannot be used. Check that the path to the Dfire binary'
                          ' in `.isambard_settings` is correct.\n'
                          'Suggestion:\n'
                          'You might want to try running isambard.settings.configure() after importing ISAMBARD in a\n'
                          'Python interpreter or running `configure.py` in the module folder.')
        warnings.warn(warning_string, DependencyNotFoundWarning)
    return is_dfire_available


global_settings['dfire']['available'] = check_dfire_avail()


def run_dfire(pdb, path=True):
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
    return float(dfire_std_out.split()[1])


def calculate_score(pdb, path=True):
    if global_settings['dfire']['available'] is None:
        global_settings['dfire']['available'] = check_scwrl_avail()
    if not global_settings['dfire']['available']:
        warning_string = ('Dfire not found, score has not been calculated.\n'
                          'Check that the path to the Scwrl binary in `.isambard_settings` is correct.\n'
                          'You might want to try rerunning `configure.py`\n')
        warnings.warn(warning_string, DependencyNotFoundWarning)
        return
    dfire_out = run_dfire(pdb, path=path)
    return parse_dfire_out(dfire_out)
