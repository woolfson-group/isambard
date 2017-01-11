import subprocess
import tempfile
from pathlib import Path

from settings import global_settings

def run_goap(input_file,path=True):

    if path:
        input_path=Path(input_file)
        if not input_path.exists():
            print('No file found at', path)
            return None,None
    else:
        pathf = tempfile.NamedTemporaryFile()
        encoded_input = input_file.encode()
        pathf.write(encoded_input)
        pathf.seek(0)
        file_path = pathf.name
        input_path=Path(file_path)

    goap_dir = "/Users/chgjb/goap-alone"
    goap_exe = "/Users/chgjb/goap-alone/goap"
    goap_fh = tempfile.NamedTemporaryFile()
    goap_input = "{}\n{}\n".format(goap_dir,input_path.name)
    encoded_goap_input = goap_input.encode()
    goap_fh.write(encoded_goap_input)
    goap_fh.seek(0)

    try:
        goap_output = subprocess.run([str(goap_exe)],stdin=goap_fh,stdout=subprocess.PIPE,stderr=subprocess.PIPE)

    except FileNotFoundError as e:
        print (e, '\nSomething went wrong')
        return None

    try:
        goap_results = goap_output.decode()

    except ValueError:
        print ("No result")
        return None

    return goap_results