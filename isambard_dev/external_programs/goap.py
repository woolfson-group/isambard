import subprocess
import tempfile
import os
from pathlib import Path
from shutil import copyfile

from settings import global_settings
from tools.isambard_warnings import check_availability


def test_goap():
    is_goap_available = False
    if os.path.isfile(global_settings['goap']['goap_exe']):
        is_goap_available = True
    return is_goap_available


@check_availability('goap', test_goap, global_settings)
def run_goap(input_file, path=True):
    if path:
        input_path = Path(input_file)
        if not input_path.exists():
            print('No file found at', path)
            return None

        pathf = tempfile.NamedTemporaryFile(dir=os.getcwd())
        copyfile(input_file, pathf.name)
        file_path = pathf.name
        input_path = Path(file_path)
    else:
        pathf = tempfile.NamedTemporaryFile(dir=os.getcwd())
        encoded_input = input_file.encode()
        pathf.write(encoded_input)
        pathf.seek(0)
        file_path = pathf.name
        input_path = Path(file_path)

    goap_dir = Path(global_settings['goap']['folder'])
    goap_exe = Path(global_settings['goap']['goap_exe'])

    goap_fh = tempfile.NamedTemporaryFile(dir=os.getcwd())
    goap_input = "{}\n{}\n".format(str(goap_dir), str(input_path.name))

    encoded_goap_input = goap_input.encode()

    goap_fh.write(encoded_goap_input)
    goap_fh.seek(0)
    try:
        goap_output = subprocess.run(str(goap_exe), stdin=goap_fh, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    except FileNotFoundError as e:
        print (e, '\nSomething went wrong')
        return None

    try:
        goap_results = goap_output.stdout.decode()

    except ValueError:
        print ("No result")
        return None

    goap_results.rstrip()
    scores = GoapScore(goap_results)
    return scores


class GoapScore(object):
    def __init__(self, scores):
        """Object containing all the different scores calculated by GOAP.

        Parameters
        ----------

        scores: [(str, str, str)]
            List of GOAP scores
        """
        scores.rstrip()
        score_words = scores.split()
        goap_score = float(score_words[2])
        dfire_score = float(score_words[3])
        goap_ag_score = float(score_words[4])

        self.goap = goap_score
        self.dfire = dfire_score
        self.goap_ag = goap_ag_score

    def __repr__(self):
        return "<GOAP Score {:.2f}: | DFIRE Score {:.2f} | GOAP_AG Score {:.2f}>".format(
            self.goap, self.dfire, self.goap_ag)
