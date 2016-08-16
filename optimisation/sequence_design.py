import json
import os
import random
from collections import Counter

from settings import global_settings

_palette_path = os.path.join(global_settings['package_path'], 'optimisation', 'residue_palettes')
available_palettes = {x[:-5]: os.path.join(_palette_path, x) for x in os.listdir(_palette_path)}


def natural_similarity(sequence, threshold=0.1):
    """ Tests if the amino acid usage in a sequence is similar to that observed in proteins

    Parameters
    ----------
    sequence : str
        Input sequence to be tested.
    threshold : float
        A percentage variation that's tolerated away from the ideal, i.e. 0.02
        is a 2% variation from the listed value in aa_freq_50.

    Returns
    -------
    bool
        True if it passes the naturality test.
    """
    aa_freq_50 = {
        'A': 0.07,
        'C': 0.02,
        'D': 0.04,
        'E': 0.05,
        'F': 0.05,
        'G': 0.06,
        'H': 0.03,
        'I': 0.07,
        'K': 0.06,
        'L': 0.10,
        'M': 0.04,
        'N': 0.04,
        'P': 0.05,
        'Q': 0.04,
        'R': 0.06,
        'S': 0.07,
        'T': 0.05,
        'V': 0.06,
        'W': 0.01,
        'Y': 0.03
    }

    res_passed = {}
    res_bd = Counter(sequence)
    for k, v in list(res_bd.items()):
        aa_fq = aa_freq_50[k]
        res_percentage = v/len(sequence)
        if res_percentage < (aa_fq+threshold):
            res_passed[k] = True
        else:
            res_passed[k] = False
    if all(res_passed.values()):
        return True
    else:
        return False


def measure_charge(seq):
    residues = Counter(seq)
    charge = 0
    try:
        charge += residues['K']
    except KeyError:
        pass
    try:
        charge += residues['R']
    except KeyError:
        pass
    try:
        charge -= residues['D']
    except KeyError:
        pass
    try:
        charge -= residues['E']
    except KeyError:
        pass
    return charge


def weighted_choice(weights):
    totals = []
    running_total = 0

    for w in weights:
        running_total += w
        totals.append(running_total)

    rnd = random.random() * running_total
    for i, total in enumerate(totals):
        if rnd < total:
            return i


class SequenceDesigner:
    def __init__(self, residue_pattern, palette='unweighted'):

        """Converts a sequence of residue types to a semi-random amino acid based of amino acid weightings.

        Lowercase characters represent different groups of amino acids while upper case residues dictate an exact
        residue type. The exact group and relative weightings of amino acids that a lower case letter represents depends
        on which 'palette' method is use to generate the sequence."""

        self.palette = self.load_palette(palette)
        self.residue_p = residue_pattern

    def __repr__(self):
        return "<Sequence Designer: Residue Pattern={}>".format(str(self.residue_p))

    @staticmethod
    def load_palette(palette):
        with open(available_palettes[palette], 'r') as inf:
            aa_palette = json.loads(inf.read())
        return aa_palette

    def make_sequence(self):
        cannonical = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W',
                      'Y']
        working_seq = []
        for aa in self.residue_p:
            if aa.istitle():
                if aa in cannonical:
                    working_seq.append(aa)
                else:
                    raise ValueError(
                        'Residue label {0} in non-cannonical, please use standard amino acid labels.'.format(aa))
            else:
                choice_i = weighted_choice([x[1] for x in self.palette[aa]])  # 1st element is the residue weight
                working_seq.append(self.palette[aa][choice_i][0])  # 0th element is the res string
        return str(''.join(working_seq))


__author__ = 'Christopher W. Wood'
