"""loop_closure.py contains code for performing kinematic loop closure."""

import copy
import math
import random
import sys

from ampal.secondary_structure.ta_polypeptide import TAPolypeptide
from tools.geometry import distance


k_boltz = 1.987206504191549E-003  # Boltzmann constant


def mutate_ta(in_tas, d_ta=10):
    """Randomly alters a single torsion angle in a list.

    Parameters
    ----------
    in_tas: [[float, float, float]]
        List of input torsion angles.
    d_ta: float
        Maximum allowed variance in torsion angle. Uniformly distributed +/-.

    Returns
    -------
    tas: [[float, float, float]]
        List of modified torsion torsion angles.
    """
    tas = in_tas[:]
    ta_i = random.choice(range(len(tas)))
    mut_t = random.choice(range(1, 3))
    angle = random.uniform(-d_ta, d_ta)
    ta = tas[ta_i][:]
    ta[mut_t] += angle
    tas[ta_i] = ta
    return tas


def calc_align_max_dist(frag1, frag2):
    """Returns the maximum distance between a pair of atoms in equivalent sets.

    Parameters
    ----------
    frag1: Polypeptide or Residue
        Reference section of protein for alignment.
    frag2: Polypeptide or Residue
        Target section of protein for alignment.

    Returns
    -------
    max_distance: float
        Maximum distance between a pair of like atoms.
    """
    frag1_atoms = frag1.backbone.get_atoms()
    frag2_atoms = frag2.backbone.get_atoms()
    return max([distance(x, y) for x, y in zip(frag1_atoms, frag2_atoms)])


def check_move(new, old, t=298.15):
    """Determines if a torsion angle move will be accepted.

    Uses Boltzmann distribution scaled by temperature."""
    return math.exp(-(new - old)/(k_boltz*t))


def fit_loop_between(polypeptide, target_monomers, loop_length, rounds=10000, ts=[297.0], print_fit=True):
    """Attempts to fit a residue of a set length between two regions of protein.

    Parameters
    ----------
    polypeptide: Polypeptide
        Polypeptide that the loop will join C terminally.
    target_monomers: Polypeptide or Residue
        Target region of protein that the loop will attempt to fit to.
    loop_length: int
        Length of loop to be fitted.
    rounds: int
        Number of rounds of moves per temperature.
    ts: [float]
        Range of temperatures to be used during fitting. The best fit from
        a particular temperature will be used as the starting point for the
        next temperature.
    print_fit: bool
        If True will print fit progress to standard out.

    Returns
    -------
    loop: Polypeptide
        Fitted loop.
    """
    polypeptide.tag_torsion_angles()
    target_monomers[0].ampal_parent.tag_torsion_angles()
    working_polypeptide = copy.deepcopy(polypeptide)
    best_tas = None
    best_rmsd = None
    for t in ts:
        best_tas, best_rmsd = loop_move(working_polypeptide, target_monomers, loop_length,
                                        rounds, t=t, starting_tas=best_tas, starting_rmsd=best_rmsd,
                                        print_fit=print_fit)
    loop = TAPolypeptide(best_tas)
    return loop


def loop_move(polypeptide, target_monomers, loop_length, rounds, t=297.0,
              starting_tas=None, starting_rmsd=None, move_func=mutate_ta, print_fit=True):
    """

    Parameters
    ----------
    polypeptide: Polypeptide
        Polypeptide that the loop will join C terminally.
    target_monomers: Polypeptide or Residue
        Target region of protein that the loop will attempt to fit to.
    loop_length: int
        Length of loop to be fitted.
    rounds: int
        Number of rounds of moves per temperature.
    t: float
        Temperature used during moves. This will affect the probability that a
        high-energy move will be accepted.
    starting_tas: [[float, float, float]] or None
        Used if the loop_move has to be started from a particular conformation.
    starting_rmsd: float or None
        Used if the loop_move has to be started from a particular conformation.
    move_func: function
        A function that perfroms a MC move on a set of torsion angles.
    print_fit: bool
        If True will print fit progress to standard out.

    Returns
    -------
    best_tas: [[float, float, float]] or None
        Set of torsion angles selected during minimisation.
    best_rmsd: float
        Best RMSD between loop reference and target section.
    """
    if not starting_tas:
        current_tas = [[178, 100, 100] for _ in range(loop_length)]
        best_tas = None
    else:
        current_tas = starting_tas[:]
        best_tas = current_tas
    target_tas = [list(x.tags['tas']) for x in target_monomers]
    if None in target_tas[0]:
        raise ValueError("All torsion angles must be present in first residue.\n")
    if starting_rmsd:
        current_rmsd = starting_rmsd
        best_rmsd = starting_rmsd
    else:
        current_rmsd = 10000
        best_rmsd = 10000
    search = True
    max_rounds = rounds
    current_round = 0
    while search:
        tas = move_func(current_tas)
        loop = TAPolypeptide(tas + target_tas)
        polypeptide.c_join(loop)
        rmsd = calc_align_max_dist(polypeptide[-len(target_monomers):-1], target_monomers[:-1])
        if rmsd < current_rmsd:
            current_tas = tas
            current_rmsd = rmsd
            if current_rmsd < best_rmsd:
                best_rmsd = current_rmsd
                best_tas = current_tas
        else:
            if check_move(rmsd, current_rmsd, t=t) > random.uniform(0, 1):
                current_rmsd = rmsd
                current_tas = tas
        current_round += 1
        if not (current_round % 10) and print_fit:
            sys.stdout.write("\rRMSD ({}): {} (best {}), t={}".format(
                move_func, *[float_f(x) for x in (current_rmsd, best_rmsd)], t))
            sys.stdout.flush()
        if current_round == max_rounds:
            search = False
        del (polypeptide._monomers[-(loop_length + len(target_monomers)):])
    return best_tas, best_rmsd


def float_f(f):
    """Formats a float for printing to std out."""
    return '{:3.3f}'.format(f).rjust(7)


__author__ = "Christopher W. Wood"
