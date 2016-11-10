import os

import numpy

from tools.geometry import angle_between_vectors, dihedral, distance
from tools.tools_pdb import extract_all_ss_dssp, gen_parsed_pdb_dict, find_ss_regions


# TODO: Rewrite using APM structure
def ss_pattern(fragments):
    """Creates a string that describes the pattern of secondary structure in a set of fragments."""
    ss_types = {
        'H': 'a',
        'E': 'b',
        ' ': 'l',
        'T': 'l',
        'S': 'l',
        'B': 'l',
        'G': '3',
        'I': 'p'
    }
    pattern = []
    for frag in fragments:
        pattern.append(ss_types[frag[0][1]])
    return ''.join(pattern)


def gather_loop_information(pdb_path, dssp_path):
    """Generates a Pandas dataframe with information about all the loops in a pdb file and its DSSP file.

    Parameters
    ----------
    pdb_path : str
        Path to PDB file.
    dssp_path : str
        Path to DSSP file.

    Returns
    -------
    loop_data : pandas.dataframe
        Pandas dataframe containing information on every loop region. Columns:
            "loop_type" : str 3 letters showing preceding and following ss elements e.g. 'ala' = alpha-loop-alpha
            "res_start" : int starting residue number
            "res_end" : int ending residue number
            "chain" : str Chain identifier
            "sequence" : str Sequence
            "distance" : float Distance between C of preceding sse and N of following sse
            "angle" : float Angle between vectors CA-C of preceding sse and N-CA of following sse
            "dihedral" : float Dihedral from CA1-C1-N2-CA2
            "backbone_coordinates" : str([triple]) Backbone coordinates of loop
            "entering_residue" : str([triple]) Preceding residue backbone coordinates
            "exiting_residue" : str([triple]) Following residue backbone coordinates
    """
    pdb_code = os.path.basename(pdb_path).split('.')[0]
    ss_eles = extract_all_ss_dssp(dssp_path)
    frags = find_ss_regions(ss_eles)
    sspat = ss_pattern(frags)
    frames = [sspat[i - 1:i + 2] for i in range(1, len(sspat) - 1)]
    frag_groups = [frags[i - 1:i + 2] for i in range(1, len(frags) - 1)]
    frames_frags = list(zip(frames, frag_groups))
    loops = []
    for frame, frag in frames_frags:
        if frame[1] == 'l':
            loops.append((frame, frag[1][0][0], frag[1][-1][0], frag[1][0][2], ''.join([x[3] for x in frag[1]])))
    parsed_pdb = gen_parsed_pdb_dict(pdb_path)
    loop_d = []
    for loop in loops:
        try:
            if loop[1] > loop[2]:
                continue  # This is to stop the end of one chain being paired with the start of another
            coords = []
            for res in range(loop[1] - 1, loop[2] + 2):
                coords.extend(parsed_pdb[loop[3]][(res, '')][1][:4])
            entering = [numpy.array(x) for x in coords[0:4]]
            exiting = [numpy.array(x) for x in coords[-4:]]
            ca_c_v = entering[2] - entering[1]
            n_ca_v = exiting[1] - exiting[0]
            ca_ca_distance = distance(entering[1], exiting[1])
            c_n_distance = distance(entering[2], exiting[0])
            end_angle = angle_between_vectors(ca_c_v, n_ca_v)
            end_dihedral = dihedral(entering[1], entering[2], exiting[0], exiting[1])
            loop_d.append([pdb_code] + list(loop) + [len(loop[-1]), ca_ca_distance, c_n_distance,
                                                     end_angle, end_dihedral,
                                                     str(coords[4:-4]), str(coords[0:4]), str(coords[-4:])])
        except KeyError:
            loop_d.append(None)
    return loop_d


__author__ = 'Christopher W. Wood'
