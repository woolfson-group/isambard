import re
import os
import subprocess
import tempfile

from settings import global_settings


# TODO: Separate this into a section that generates the PDB files for BUDE and running BUDE
def run_bude(pdb_str, omp_threads=1, experimental=False):
    """Runs BUDE on a given PDB.

    Notes
    -----
    This requires an input PDB str with different chain identifiers for running BUDE
    and also a bude_directory specified in settings.json for output files.

    Parameters
    ----------

    pdb_str : str
        Input PDB string for BUDE
    omp_threads : int
        Set number of omp threads

    Returns
    -------
    sum(bude_scores) : float
        Sum of all bude scores

    len(active_bude_chains) : int

        Number of secondary structure elements

    Raises
    ------

    IOError
        If BUDE fails to run

    """

    bude_chains = {chr(x): '' for x in range(33, 127)}

    os.environ['OMP_NUM_THREADS'] = str(omp_threads)

    bude_scores = []
    pdb_str_lines = pdb_str.split('\n')
    bude_directory = global_settings['bude']['cmd_files']
    active_bude_chains = []
    for line in pdb_str_lines:
        line = line.rstrip()
        find_chain = re.search(r'\w+\s+\d+\s+\w+\s+\w\w\w\s+(\w)\s+\d+', line)
        if find_chain:
            active_chain = find_chain.group(1)
            active_bude_chains.append(active_chain)
            bude_chains[active_chain] += line + ' ' * 80 + '\n'
        else:
            pass
    active_bude_chains = set(active_bude_chains)
    for x in active_bude_chains:
        bude_chains[x] += 'TER' + ' ' * 80 + '\n'
    for lig in active_bude_chains:
        with tempfile.NamedTemporaryFile(mode='w+') as lig_tmp, tempfile.NamedTemporaryFile(mode='w+') as rec_tmp, \
                tempfile.NamedTemporaryFile() as out_tmp:
            lig_tmp.write(bude_chains[lig])
            lig_tmp.write(('END' + ' ' * 80 + '\n'))
            rec_chains = [i for i in set(active_bude_chains) if i != lig]
            for rec in rec_chains:
                rec_tmp.write(bude_chains[rec])
            rec_tmp.write(('END' + ' ' * 80 + '\n'))
            lig_tmp.seek(0)
            rec_tmp.seek(0)
            bude_command = ['bude', '-f', os.path.join(bude_directory, 'boostlessbude.bctl'),
                            '--transformations-filename', os.path.join(bude_directory, 'ligand_trans_rot.bltr'),
                            '--forcefield-filename', (
                                os.path.join(bude_directory, 'force_field_experimental.bhff')) if experimental else (
                                os.path.join(bude_directory, 'force_field.bhff')),
                            '--emc-filename', os.path.join(bude_directory, 'protein_docking.bemc'),
                            '--zero-generation-filename', os.path.join(bude_directory, 'zeroGenSingleEnergy.bzgn'),
                            '--receptor-coordinates-filename', '{0}'.format(rec_tmp.name),
                            '--ligand-coordinates-filename', '{0}'.format(lig_tmp.name),
                            '--output-file-base', '{0}'.format(out_tmp.name),
                            '--output-energy-filename', '{0}'.format(out_tmp.name)]
            try:
                bude_std_out = subprocess.check_output(bude_command).decode()
            except subprocess.CalledProcessError:
                raise RuntimeError('BUDE failed to run, if BUDE usually runs then check the input pdb, otherwise make '
                                   'sure that the isambard/settings.json bude_directory path is correct.')
            try:
                bude_out_path = out_tmp.name + '_R01_S0001_C0001_L00001_S0001_C0001_G0001.btle'
                with open(bude_out_path, 'r') as outf:
                    bude_result = outf.read()
                os.remove(bude_out_path)
            except IOError:
                raise IOError('\nERROR: BUDE failed to run on the given input files. BUDE error:\n{}'.format(
                    bude_std_out))

        pull_score = re.search(r'(\s+\d){6}\s+([-0-9\.]+)', bude_result).group(2)
        bude_scores.append(float(pull_score))
    return bude_scores


def run_bude_internal_energy(pdb_str, omp_threads=1, experimental=False):
    """Runs BUDE in internal energy mode, rather than the default interaction energy mode.

    Parameters
    ----------
    pdb_str : str
        String of the pdb file for a protein.
    omp_threads : int
        Number of OpenMP threads that BUDE is allowed to use.

    Returns
    -------
    bude_ie_score : float
        The internal energy of the protein.
    """
    os.environ['OMP_NUM_THREADS'] = str(omp_threads)
    bude_directory = global_settings['bude']['cmd_files']
    with tempfile.NamedTemporaryFile(mode='w+') as s_tmp, tempfile.NamedTemporaryFile(mode='w+') as out_tmp:
        s_tmp.write(pdb_str)
        s_tmp.write(('END' + ' ' * 80 + '\n'))
        s_tmp.seek(0)
        bude_command = [global_settings['bude']['internal_energy_binary'],
                        '-f', os.path.join(bude_directory, 'boostlessbude.bctl'),
                        '--transformations-filename', os.path.join(bude_directory, 'ligand_trans_rot.bltr'),
                        '--forcefield-filename', (
                            os.path.join(bude_directory, 'force_field_experimental.bhff')) if experimental else (
                            os.path.join(bude_directory, 'force_field.bhff')),
                        '--emc-filename', os.path.join(bude_directory, 'protein_docking.bemc'),
                        '--zero-generation-filename', os.path.join(bude_directory, 'zeroGenSingleEnergy.bzgn'),
                        '--receptor-coordinates-filename', '{0}'.format(s_tmp.name),
                        '--ligand-coordinates-filename', '{0}'.format(s_tmp.name),
                        '--output-file-base', '{0}'.format(out_tmp.name),
                        '--output-energy-filename', '{0}'.format(out_tmp.name)]
        try:
            bude_std_out = subprocess.check_output(bude_command).decode()
        except subprocess.CalledProcessError:
            raise RuntimeError('BUDE failed to run, if BUDE usually runs then check the input pdb, otherwise make sure'
                               ' that the isambard/settings.json bude_directory path is correct.')
        try:
            match = re.search(r'Conf: \d+ has an energy value of: ([0-9-.]+)', bude_std_out)
            bude_ie_score = float(match.groups()[0])
        except AttributeError:
            raise IOError('\nERROR: BUDE failed to run on the given input files. BUDE error:\n{}'.format(
                bude_std_out))
    return bude_ie_score


def parse_bude_out(pdb_str, mode='additive', experimental=False):
    bude_scores = run_bude(pdb_str, experimental=experimental)
    if mode == 'additive':
        return sum(bude_scores)
    elif mode == 'average':
        return sum(bude_scores)/len(bude_scores)
    elif mode == 'raw':
        return bude_scores
    else:
        raise ValueError('Parse mode "{}" not found. Please choose from "additive", "average" or "raw".'.format(mode))


def run_bude_additive(pdb_str):
    bude_scores = run_bude(pdb_str)
    return sum(bude_scores)


def run_bude_average(pdb_str):
    bude_scores = run_bude(pdb_str)
    return sum(bude_scores)/len(bude_scores)

__author__ = 'Christopher W. Wood'
