from bs4 import BeautifulSoup
import parmed
from pathlib import Path
from urllib.request import Request, urlopen
from urllib.error import URLError


def dict_from_mmcif(mmcif, path=True):
    """Parse mmcif file into a dictionary.

    Notes
    -----
    Full list of keys/value types, and further information on them can be viewed here:
            http://mmcif.wwpdb.org/docs/pdb_to_pdbx_correspondences.html
    All values in the returned dict are str or list(str).
    This means that some of the data values are string representations of integers
    - parse these outside of this function if desired.
    An alternative approach to this can be found in Biopython (via the function Bio.PDB.MMCIF2Dict.MMCIF2Dict).
    mmcif files are subject to the usual "here be dragons" problems of the PDB and difficult file formats.
    As such, this function is likely to be in a permanent state of flux as more dragons are found.

    Parameters
    ----------
    mmcif : str
        mmcif string or a path to an mmcif file.
    path : bool
        True if mmcif is a path.

    Returns
    -------
    cif_data : dict
        Keys are cif data names, e.g. '_struct_keywords.text'.
        Values are str or list(str).

    """
    if path:
        with open(mmcif, 'r') as foo:
            lines = foo.readlines()
    else:
        lines = mmcif.splitlines()
    lines = [' '.join(x.strip().split()) for x in lines]

    # Some of the data in a .cif files are stored between 'loop_' to initiate a loop, and '#' to terminate it.
    # The variable 'loop' is a flag to keep track of this behaviour.
    loop = False

    # Set up the dictionary to populate as the lines of the .cif file are iterated over.
    cif_data = {}

    for i, line in enumerate(lines):

        if not line:
            continue

        # hash signifies end of a loop. Ensure loop flag is set to False.
        if line == '#':
            loop = False
            continue

        if not loop:
            # This line initiates a loop section, in which keys are listed first,
            # followed by lines of data in which the values are listed in the same order as the above keys.
            # The values in the loop section will be stored as lists - there are multiple values for one key.
            # An example of this type of data is the 'REVDAT' section, which stores details on the (potentially
            # numerous) various revisions made to the PDB file during its history.

            if line[:5] == 'loop_':
                loop = True
                key_list = []
                continue

            # Lines beginning with '_' start with data names, i.e. keys in the cif_data dictionary.
            elif line[0] == '_':
                # If line consists only of a key, then subsequent lines may contain the associated value.
                if len(line.split()) == 1:
                    current_key = line
                    count = 1
                    while True:
                        # Look forward until a key is found, keeping count of the number of lines in between.
                        try:
                            if lines[i + count][0] != '_':
                                count += 1

                            # prevent infinite loop.
                            elif i + count > len(lines):
                                break

                            else:
                                if count > 1:
                                    try:
                                        cif_data[current_key] = ' '.join(lines[i + 1: i + count])
                                    except IndexError:
                                        cif_data[current_key] = None
                                else:
                                    cif_data[current_key] = None
                                break
                        except IndexError:
                            break
                    continue

                # Simplest case. Line is a key-value pair, with the key identified by its first character, '_'.
                elif len(line.split()) > 1:
                    line = line.split()
                    try:
                        cif_data[line[0]] = ' '.join(line[1:])
                    except IndexError:
                        cif_data[line[0]] = None
                    continue

            # Line is one of multiple lines that are combined into a value in the while True: loop above.
            else:
                continue

        else:
            # Within a loop section, keys are identified by their first character '_'.
            # Add them to the list of keys in the loop.
            if line[0] == '_':
                if len(line.split()) == 1:
                    key_list.append(line)
                    if line not in cif_data.keys():
                        cif_data[line] = []

            # Within a loop section, the values are listed within a single space-separated line in the same order
            # that the keys were listed at the start of the loop.
            else:
                # Cannot do a simple split if any of the values themselves are strings containing at least one space.
                if '\"' in line and line.count('\"') % 2 == 0:
                    line_parts = [x.strip() for x in line.split('\"') if x]
                    line = []
                    for part in line_parts:
                        if line_parts.index(part) % 2 == 0:
                            for x in part.split():
                                line.append(x)
                        else:
                            line.append(part)

                elif '\'' in line and line.count('\'') % 2 == 0:
                    line = [x.strip() for x in line.split('\'') if x]

                elif len(key_list) == len(line.split()):
                    line = line.split()

                if len(key_list) == len(line):
                    for j, v in enumerate(line):
                        cif_data[key_list[j]].append(line[j])
                else:
                    # CURRENTLY THERE IS A PROBLEM WITH REALLY LONG LOOPS eg _pdbx_refine_tls*, _pdbx_struct_oper_list*
                    # The values span multiple lines, and therefore do not satisfy
                    # the condition of the above 'if' statement.
                    # A correction for this needs to keep track of the value count on subsequent lines,
                    # until the 'if' condition is met.
                    continue

    return cif_data


def get_protein_dict_from_parmed(cif_file):
    s = parmed.formats.CIFFile.parse(cif_file, skip_bonds=True)
    protein_dict = {
        'keywords': s.keywords,
        'space group': s.space_group,
        'experimental_method': s.experimental,
        'resolution': s.resolution,
        'atoms_total': len(s.atoms),
        'title': s.title,
        'journal': s.journal,
        'authors': s.authors,
        'journal_year': s.year,
        'journal_volume': s.volume,
        'journal_volume_page': s.volume_page,
        'journal_page': s.page,
        'doi': s.doi,
    }
    return protein_dict


def parse_PISCES_output(pisces_output, path=False):
    """ Takes the output list of a PISCES cull and returns in a usable dictionary.

    Notes
    -----
    Designed for outputs of protein sequence redundancy culls conducted using the PISCES server.
    http://dunbrack.fccc.edu/PISCES.php
    G. Wang and R. L. Dunbrack, Jr. PISCES: a protein sequence culling server. Bioinformatics, 19:1589-1591, 2003.

    Parameters
    ----------
    pisces_output : str or path
        Output list of non-redundant protein chains from PISCES, or path to text file.
    path : bool
        True if path given rather than string.

    Returns
    -------
    pisces_dict : dict
        Data output by PISCES in dictionary form.
    """
    pisces_dict = {}
    if path:
        pisces_path = Path(pisces_output)
        pisces_content = pisces_path.read_text().splitlines()[1:]
    else:
        pisces_content = pisces_output.splitlines()[1:]
    for line in pisces_content:
        pdb = line.split()[0][:4].lower()
        chain = line.split()[0][-1]
        pdb_dict = {'length': line.split()[1],
                    'method': line.split()[2],
                    'resolution': line.split()[3],
                    'R-factor': line.split()[4],
                    'R-free': line.split()[5]}
        if pdb in pisces_dict:
            pisces_dict[pdb]['chains'].append(chain)
        else:
            pdb_dict['chains'] = [chain]
            pisces_dict[pdb] = pdb_dict
    return pisces_dict


def download_decode(URL, encoding='utf-8', verbose=True):
    """ Downloads data from URL and returns decoded contents."""
    if verbose:
        print("Downloading data from " + URL)
    req = Request(URL)
    try:
        with urlopen(req) as u:
            decoded_file = u.read().decode(encoding)
    except URLError as e:
        if hasattr(e, 'reason'):
            print('Server could not be reached.')
            print('Reason: ', e.reason)
        elif hasattr(e, 'code'):
            print('The server couldn\'t fulfill the request.')
            print('Error code: ', e.code)
        return None
    return decoded_file


def olderado_best_model(pdb_id):
    """ Checks the Olderado web server and returns the most representative conformation for PDB NMR structures.

    Notes
    -----
    Uses OLDERADO from the EBI.
    See http://www.ebi.ac.uk/pdbe/nmr/olderado/ and citations therein.

    Parameters
    ----------
    pdb_id : str
        The 4-character PDB code for the NMR structure of interest.

    Returns
    -------
    model_no : int
        The conformation number of the most-representative conformation.

    Raises
    ------
    ValueError
        If the model number it finds is not an integer. This might indicate that the website format has changed.
    """
    pdb_code = pdb_id[:4].lower()
    olderado_url = 'http://www.ebi.ac.uk/pdbe/nmr/olderado/searchEntry?pdbCode=' + pdb_code
    olderado_page = download_decode(olderado_url, verbose=False)
    if olderado_page:
        parsed_page = BeautifulSoup(olderado_page, 'html.parser')
    else:
        return None
    try:
        best_model = parsed_page.find_all('td')[1]
    except IndexError:
        print("No model info could be found for {0} - ensure that it's an NMR structure.".format(pdb_id))
        return None
    try:
        model_no = int(best_model.string)
    except ValueError as v:
        print("Did not find a number for best model.")
        raise v
    return model_no


__author__ = 'Jack W. Heal, Kieran L. Hudson'
