import os
import re
import subprocess
import sys

from add_ons.filesystem import FileSystem

from settings import global_settings


# TODO Think about how to deal with modified/obsolete pdb codes.
# TODO Are we going to want to run stats on the clusterings? Maybe for logging?
# TODO Clever ways to create false fasta files for subsets of the whole sequence.

def get_fasta_filepath(code):
    """Return filepath to fasta file if it exists.


    Parameters
    ----------
    code : str
        A pdb code, e.g. '2ebo'.


    Returns
    -------
    filepath : str or None
        Filepath for the fasta file associated with the pdb code if it exists.
    """
    data_dir = global_settings["structural_database"]["path"]
    filepath = ("{0}{1}/{2}/fasta/{2}.fasta".format(data_dir, code[1:3], code))
    if os.path.exists(filepath):
        return filepath
    else:
        return None


def get_fastas(pdb_list):
    """ Get list of fasta filepaths for codes in pdb_list.

    Parameters
    ----------
    pdb_list: list of str
        List of PDB codes.

    Returns
    -------
    List of str
        Fasta filepaths associated with codes in pdb_list.
    """
    fastas = []
    for code in pdb_list:
        fasta = get_fasta_filepath(code)
        if fasta is None:
            fs = FileSystem(code)
            fasta = fs.fasta
        fastas.append(fasta)
    return fastas


def get_all_fastas():
    """ Get list of fasta filepaths for all fasta files in the database.

    Returns
    -------
    List of str
        Fasta filepaths for all PDB codes in the database.
    """
    data_dir = global_settings["structural_database"]["path"]
    fastas = []
    for root, dirs, files in os.walk(data_dir):
        for file in files:
            if file.endswith(".fasta"):
                 fastas.append(os.path.join(root, file))
    return fastas


def concatenate_fastas(file_list, outputfile):
    """ Takes list of files and concatenates their contents into outputfile.

    Parameters
    ----------
    file_list : list of str
        Filepath strings for fasta files.
    outputfile : str
        Filepath for concatenated file containing contents of all files in file_list.

    Returns
    -------
    None
    """
    with open(outputfile, 'w') as outfile:
        for fname in file_list:
            with open(fname) as infile:
                for line in infile:
                    outfile.write(line)
    return


def cdhit_choose_word_size(cval):
    """ Get appropriate -n flag, given -c flag.

    CD-Hit specific function used to ensure the input parameters are appropriately chosen.
    The CD-Hit user guide gives appropriate values for the word size to use (-n flag), depending on the comparison
    threshold (-c flag). This function returns the appropriate word size to use when using a particular comparison
    threshold.

    Parameters
    ----------
    cval : float
        Comparison threshold (-c flag) to use in CD-Hit.

    Returns
    -------
    nval : float
        Word size (-n flag) to use in CD-Hit, or None if cval is invalid.

    Raises
    ------
    ValueError
        If cval is less than 0.4 or greater than 1.0.

    """
    if 0.7 <= cval <= 1.0:
        word_size = 5
    elif 0.6 <= cval < 0.7:
        word_size = 4
    elif 0.5 <= cval < 0.6:
        word_size = 3
    elif 0.4 <= cval < 0.5:
        word_size = 2
    else:
        raise ValueError("CD-HIT needs -c flag to be in [0.4, 1.0]")
    return word_size


def run_cascading_cdhit(fasta_file, run_svals=True):
    """ Runs CD-Hit hierarchically on fasta_file for different values of '-c', '-s'.

    See CD-Hit manual for details of -c and -s flags.
    If run_svals, then run at different values of '-s', otherwise just vary '-c' flag.

    Parameters
    ----------
    fasta_file : str
        Input file of sequences to be clustered, in fasta format.
    run_svals : bool, optional
        If True then run CD-Hit at different values of '-s' (default is True).

    Returns
    -------
    None
    """
    cdhit = global_settings["cd_hit"]["path"]
    # Set range of values for '-c' and '-s' flag
    cvals = [0.9, 0.8, 0.7, 0.6]
    svals = [0.9, 0.8, 0.7, 0.6]
    cval_names = [int(cval * 100) for cval in cvals]
    # Values for '-n' flag determined from '-c' flag (see CD-Hit manual for more detail).
    nvals = [cdhit_choose_word_size(cval) for cval in cvals]

    for i in range(len(cvals)):
        if i == 0:
            cdhit_infile = fasta_file
        else:
            cdhit_infile = "{0}{1}".format(fasta_file, cval_names[i - 1])
        cdhit_outfile = "{0}{1}".format(fasta_file, cval_names[i])
        cmd = [cdhit, '-i', cdhit_infile, '-o', cdhit_outfile, '-c', str(cvals[i]), '-n', str(nvals[i])]
        print(subprocess.check_output(cmd))

        # For each value of '-c', run at each value of '-s'.
        if run_svals:
            for sval in svals:
                if i == 0:
                    cdhit_infile = fasta_file
                else:
                    cdhit_infile = "{0}{1}s{2}".format(fasta_file, cval_names[i - 1], int(sval * 100))
                cdhit_outfile = "{0}{1}s{2}".format(fasta_file, cval_names[i], int(sval * 100))
                cmd = [cdhit, '-i', cdhit_infile, '-o', cdhit_outfile, '-c', str(cvals[i]),
                       '-s', str(sval), '-n', str(nvals[i])]
                print(subprocess.check_output(cmd))
    return


def get_representative_codes(cluster_file):
    """ Return PDB codes of the representative structure from each cluster.

    Parameters
    ----------
    cluster_file : str
        Filestring that contains output of CD-Hit clustering.

    Returns
    -------
    None
    """
    representative_pdbs = []
    with open(cluster_file, "r") as foo:
        for line in foo.readlines():
            line = line.rstrip()
            if re.search("\*$", line):
                line = line.split(">")
                code = line[1][:4].lower()
                if code not in representative_pdbs:
                    # validate it's a pdb code
                    original_fasta = get_fasta_filepath(code)
                    if original_fasta is not None:
                        representative_pdbs.append(code)
    return


def merge_clustering(old_clustering, new_fasta_file, sval=0):
    """ Merges the sequences new sequences into an existing clustering.

    See CD-Hit manual for more information. Specifically, the section entitled "Incremental clustering".
    Link to CD-Hit manual: http://weizhongli-lab.org/cd-hit/wiki/doku.php?id=cd-hit_user_guide.

    This function runs the CD-Hit commands outlined in the manual, merging the new sequences into the old clustering
    system. The clusters are then renumbered.

    Finally, some perl scripts from CD-Hit are run to get some stats on the clustering, and to generate clustering
    table files that may be useful for putting redundancy data into sql tables.

    old_clustering is a fasta file that has already been run through CD-Hit hierarchically using each value in cvals.

    During the running of this script, many files are created and then destroyed.
    If successful, the new clustering will simply overwrite the old files, and new_fasta_file will be deleted.

    Notes
    -----
        Sequences in new that are longer than any in old are considered to be non-redundant.
        This has potential to be inaccurate, and is more likely to be so after long periods of time.
        For this reason, it is recommended to do a 'big update' every so often, where the clustering is done on the
        whole database from scratch. I will set this 'big update' to run in a cron script, every 6(?) months?

    Parameters
    ----------
    old_clustering : str
        The fasta file that has already been clustered. (Not the cluster file that resulted from the clustering).
        old_clustering should already have been run through CD-Hit hierarchically using run_cascading_cdhit.
    new_fasta_file : str
        The fasta file of new sequences to be merged into the old clustering.
    sval : float, optional
        Runs CD-Hit with the '-s' flag set to sval. Defaults to 0.

    Returns
    -------
    None
    """
    # Parameters for running CD-Hit.
    cvals = [0.9, 0.8, 0.7, 0.6]
    nvals = [cdhit_choose_word_size(x) for x in cvals]
    # Input and output files for running CD-Hit. Oldfiles should exist already.
    if sval == 0:
        oldfiles = ["{0}{1}".format(old_clustering, int(x * 100)) for x in cvals]
        newfiles = ["{0}{1}".format(new_fasta_file, int(x * 100)) for x in cvals]
    else:
        oldfiles = ["{0}{1}s{2}".format(old_clustering, int(x * 100), int(sval * 100)) for x in cvals]
        newfiles = ["{0}{1}".format(new_fasta_file, int(x * 100), int(sval * 100)) for x in cvals]
    for fname in oldfiles:
        if not os.path.exists(fname):
            print("{0} does not exist. Run run_cascading_cdhit on {1} first.".format(fname, old_clustering))
            sys.exit()

    # These are temporary filepaths that will be written to and later removed.
    old_clustering_dir = os.path.dirname(old_clustering)
    temp_new = "{0}/temp_new_file".format(old_clustering_dir)
    temp_clstr = "{0}/temporary.clstr".format(old_clustering_dir)
    updated_file = "{0}/updated_db".format(old_clustering_dir)
    if sval == 0:
        updated_files = ["{0}{1}.clstr".format(updated_file, int(x * 100)) for x in cvals]
    else:
        updated_files = ["{0}{1}.clstr".format(updated_file, int(x * 100), int(sval * 100)) for x in cvals]

    # Handles to the CD-Hit programs that will be run.
    cdhit = global_settings["cd_hit"]["path"]
    cdhit2d = '{0}-2d'.format(cdhit)
    cdhit_dir = os.path.dirname(cdhit)
    cluster_merge = os.path.join(cdhit_dir, 'clstr_merge.pl')
    cluster_renumber = os.path.join(cdhit_dir, 'clstr_renumber.pl')

    for i in range(len(cvals)):
        # Run CD-Hit commands as outlined in the CD-Hit manual.
        if i == 0:
            i2 = new_fasta_file
        else:
            i2 = newfiles[i - 1]
        cmd = [cdhit2d, '-i', oldfiles[i], '-i2', i2, '-o', temp_new, '-c', str(cvals[i]), '-s', str(sval), '-n', str(nvals[i])]
        print(subprocess.check_output(cmd))
        cmd = ([cdhit, '-i', temp_new, '-o', newfiles[i], '-c', str(cvals[i]), '-s', str(sval), '-n', str(nvals[i])])
        print(subprocess.check_output(cmd))
        cmd = ['cat', newfiles[i]]
        contents = subprocess.check_output(cmd)
        with open(oldfiles[i], 'a') as foo:
            foo.write(contents)
        cmd = ([cluster_merge, "{0}.clstr".format(oldfiles[i]), "{0}.clstr".format(temp_new)])
        cmd_out = subprocess.check_output(cmd)
        with open(temp_clstr, 'w') as foo:
            foo.write(cmd_out)
        cmd_one = ['cat', temp_clstr]
        cmd_two = ['cat', "{0}.clstr".format(newfiles[i])]
        cmd_one_out = subprocess.check_output(cmd_one)
        cmd_two_out = subprocess.check_output(cmd_two)
        with open(updated_files[i], 'a') as foo:
            foo.write(cmd_one_out)
            foo.write(cmd_two_out)
        cmd = [cluster_renumber, updated_files[i]]
        cmd_out = subprocess.check_output(cmd)
        with open(updated_files[i], 'w') as foo:
            foo.write(cmd_out)

    # Cleaning up: delete all temporary files. Replace old_clustering files with the updated clustering files.
    delete_files = []
    for i in range(len(cvals)):
        replaced_file = "{0}.clstr".format(oldfiles[i])
        cmd = ['mv', updated_files[i], replaced_file]
        subprocess.check_output(cmd)
        delete_files.append(newfiles[i])
        delete_files.append("{0}.clstr".format(newfiles[i]))
    other_removals = [new_fasta_file, temp_clstr, temp_new, "{0}.clstr".format(temp_new)]
    for f in other_removals:
        delete_files.append(f)
    for f in delete_files:
        try:
            os.remove(f)
        except:
            pass

    # Generate sql_table, histogram, stats and lenplot files from .clstr files. May be useful later!
    for fname in oldfiles:
        fname += ".clstr"
        run_cdhit_sql_table(clstr_file=fname)
        run_cdhit_perl_scripts(clstr_file=fname)

    return


def run_cdhit_sql_table(clstr_file):
    """Create a .table file with cluster information from .clstr file.

    See CD-Hit manual for more information on the perl script clstr_sql_tbl.pl that this function runs.

    Parameters
    ----------
    clstr_file : str
        Filepath for a .clstr file output from CD-Hit.

    Returns
    -------
    None
    """
    cdhit = global_settings['cd_hit']['path']
    cdhit_dir = os.path.dirname(cdhit)
    cdhit_sql_table = os.path.join(cdhit_dir, 'clstr_sql_tbl.pl')

    if os.path.exists(cdhit_sql_table):
        table_name = re.sub('.clstr$', '.table', clstr_file)
        cmd = [cdhit_sql_table, clstr_file, table_name]
        subprocess.check_output(cmd)
    return


def run_cdhit_perl_scripts(clstr_file):
    """ Run three perl scripts from within CD-Hit on a .clstr file.

    See CD-Hit manual for more information.
    The perl scripts from CD-Hit implemented here are:
        1) clstr_size_histogram.pl
        2) clstr_size_stat.pl
        3) plot_len1.pl
    The perl scripts are run and output saved to files.
    The output files have the same main name a clstr_file,
        with 1) .hist, 2) .stat and 3) .lenplot
        as their respective filename extensions.

    Parameters
    ----------
    clstr_file : str
        Filepath for a .clstr file output from CD-Hit.

    Returns
    -------
    None
    """
    cdhit_dir = os.path.dirname(global_settings["cd_hit"]["path"])
    cdhit_hist = os.path.join(cdhit_dir, 'clstr_size_histogram.pl')
    cdhit_stat = os.path.join(cdhit_dir, 'clstr_size_stat.pl')
    cdhit_plotlen = os.path.join(cdhit_dir, 'plot_len1.pl')
    perl_scripts = {cdhit_hist: ".hist", cdhit_stat: ".stat", cdhit_plotlen: ".lenplot"}
    for perl_script, file_ext in perl_scripts.items():
        if os.path.exists(perl_script):
            savefile_name = re.sub('.clstr$', file_ext, clstr_file)
            cmd = [perl_script, clstr_file]
            if perl_script == cdhit_plotlen:
                clstr_sizes = "1,2-4,5-9,10-19,20-49,50-99,100-499,500-99999"
                seq_lengths = "10-59,60-149,150-499,500-1999,2000-999999"
                cmd.append(clstr_sizes)
                cmd.append(seq_lengths)
            cmd_out = subprocess.check_output(cmd)
            with open(savefile_name, 'w') as foo:
                foo.write(cmd_out)
    return


def run_cdhit_full_database(subdir_name="full_db", file_name="all_fasta"):
    """ Runs CD-Hit on all sequences in the database.

    Create a fasta file from all the fasta files in the database and runs CD-Hit hierarchically.
    Run subsequent CD-Hit analysis using CD-Hit perl scripts.
    Genenerates a series of .clstr, .table, .hist and .stat files.

    Parameters
    ----------
    subdir_name : str, optional
        Subdirectory within redundancy_dir in which to store files.
    file_name : str, optional
        Main name used in all of the files generated.
        Each file has a different file extension according to file type.

    Returns
    -------
    None
    """
    data_dir = global_settings["structural_database"]["path"]
    full_db = os.path.join(data_dir, 'redundancy', subdir_name)
    allf = get_all_fastas()
    all_fasta = os.path.join(full_db, file_name)
    concatenate_fastas(file_list=allf, outputfile=all_fasta)
    run_cascading_cdhit(fasta_file=all_fasta)
    for fname in os.listdir(full_db):
        if re.search('.clstr$', fname):
            fname = os.path.join(full_db, fname)
            run_cdhit_sql_table(clstr_file=fname)
            run_cdhit_perl_scripts(clstr_file=fname)
    return


def get_cdhit_dict(fasta_name):
    """ Get dictionary of redundancy information from .table files.

    This dictionary will be used to populate an sql table with redundancy information.
    Dictionary keys are the sequence ids from the fasta file.
    Format for each key is PDB_CODE:CHAIN_ID e.g. 1CWA:A.
    Each value has pdb_code as first entry, and a list of Booleans as the second entry.
    There is one Boolean per cval (-c flag used in CD-Hit).
    Boolean is 1 if that sequence is the representative for a cluster.
    Boolean is 0 if that sequence is not the representative for a cluster.

    Parameters
    ----------
    fasta_name : str
        Filepath for the full fasta file used to generate the .clstr and .table files.

    Returns
    -------
    cdhit_dict : dict
        Keys are sequence ids from the fasta file (str).
            Each key has format PDB_CODE:CHAIN_ID, e.g. 1CWA:A.
        Values are list of lists.
        The first list has two elements:
            1) pdb_code (str)
            2) chain_id (str)
        The second list is a list of bool.
        These are 0 if representative, 1 otherwise. 
    """
    cvals = [0.9, 0.8, 0.7, 0.6]
    cval_names = [str(int(x * 100)) for x in cvals]
    table_files = []
    for cval_name in cval_names:
        fname = "{0}{1}.table".format(fasta_name, cval_name)
        table_files.append(fname)

    cdhit_dict = {}
    with open(table_files[0], 'r') as foo:
        for line in foo.readlines():
            line = line.split()
            seqid = line[0].split("|")[0]
            if seqid not in list(cdhit_dict.keys()):
                cdhit_dict[seqid] = []
                pdb_code = seqid.split(":")[0].lower()
                chain = seqid.split(":")[1]
                cdhit_dict[seqid].append([pdb_code, chain])
                # Set each Boolean in the list as 0. If representative, change to one.
                cdhit_dict[seqid].append([0] * len(table_files))
                rep = int(line[-1])
                if rep:
                    cdhit_dict[seqid][1][0] = 1

    # Remaining .table files only contain entries for representative sequences of .table files at higher -c.
    for i in range(1, len(table_files)):
        with open(table_files[i], 'r') as foo:
            for line in foo.readlines():
                line = line.split()
                seqid = line[0].split("|")[0]
                rep = int(line[-1])
                if rep:
                    cdhit_dict[seqid][1][i] = 1

    return cdhit_dict


def get_pdb_cdhit_dict(fasta_name):
    """ Get dictionary of redundancy information from .table files.

    The dictionary associates PDB codes with information on their sequence redundancy.
    This will be used to populate an sql table with redundancy information.

    Parameters
    ----------
    fasta_name : str
        Full fasta file used to generate the .clstr and .table files are based on.

    Returns
    -------
    pdb_cdhit_dict : dict
        CD-Hit redundancy dictionary.
        Keys are PDB codes.
        Values are list of bool.
            bool is 1 if at least one of the chains in that sequence is the representative for a cluster.
            bool is 0 if none of the chains in the sequence is representative for a cluster.
    """
    cvals = [0.9, 0.8, 0.7, 0.6]
    cval_names = [str(int(x * 100)) for x in cvals]
    table_files = []
    for cval_name in cval_names:
        fname = "{0}{1}.table".format(fasta_name, cval_name)
        table_files.append(fname)

    pdb_cdhit_dict = {}
    with open(table_files[0], 'r') as foo:
        for line in foo.readlines():
            line = line.split()
            code = line[0].split("|")[0][:4].lower()
            if code not in list(pdb_cdhit_dict.keys()):
                pdb_cdhit_dict[code] = [0, 0, 0, 0]
            if pdb_cdhit_dict[code][0] == 0:
                rep = int(line[-1])
                if rep:
                    pdb_cdhit_dict[code][0] += 1

    for i in range(1, len(table_files)):
        with open(table_files[i], 'r') as foo:
            for line in foo.readlines():
                line = line.split()
                code = line[0].split("|")[0][:4].lower()
                if pdb_cdhit_dict[code][i] == 0:
                    rep = int(line[-1])
                    if rep:
                        pdb_cdhit_dict[code][i] = 1

    return pdb_cdhit_dict


__author__ = 'Jack W. Heal'
