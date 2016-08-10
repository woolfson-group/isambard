from pathlib import Path

from ampal import Atom, Ligand
from ampal.base_ampal import write_pdb
from tools.components import get_alignment_dict


# TODO replicate alignment functionality in Isambard, ie redefining coordinates around common origin and axes
def output_with_alignment_tags(monomer_list, tagged_monomer=None, output_path='', alt_states=False, strip_states=False):
    """ Outputs a PDB file with a monomer side chain atoms tagged with primes for alignment.

    Notes
    -----
    Outputs minimal PDB files with atoms tagged for alignment using align_all_to_all script in PyMOL.

    Parameters
    ----------
    monomer_list : list
        List of AMPAL monomers to be included in output file
    tagged_monomer : AMPAL monomer
        Monomer to have side chain atoms tagged for alignment
    output_path : str or pathlib Path
        Destination for output file
    alt_states : bool
        True to include only one state if multiple exist for any monomers
    strip_states : bool
        True to remove state labels to facilitate alignment of monomers

    Returns
    -------
    output_file : pathlib Path
        Location of output file
    """
    if output_path:
        output_file = Path(output_path)
    elif tagged_monomer:
        output_file = Path(tagged_monomer.mol_code + tagged_monomer.id + '_tagged.pdb')
        print('No output path specified, outputting to', output_file)
    else:
        print('Please specify an output path and/or a monomer to tag, or just use write_pdb')
        return None
    if output_file.exists():
        return output_file
    if tagged_monomer:
        temp_atoms = get_alignment_dict(tagged_monomer.mol_code)
        if tagged_monomer not in monomer_list:
            monomer_list.append(tagged_monomer)
    else:
        print('No monomer specified for tagging - outputting without tags')
        temp_atoms = {}
    for atom in temp_atoms:
        if atom not in tagged_monomer.atoms:
            print(atom, "not found in", tagged_monomer.mol_code, tagged_monomer.id, tagged_monomer.ampal_parent.id,
                  tagged_monomer.atoms.keys())
            continue
        tagged_monomer.atoms[temp_atoms[atom]] = tagged_monomer.atoms.pop(atom)
    output_pdb = write_pdb(monomer_list, alt_states=alt_states, strip_states=strip_states)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    output_file.write_text(output_pdb)
    for atom in temp_atoms:
        if temp_atoms[atom] in tagged_monomer.atoms:
            tagged_monomer.atoms[atom] = tagged_monomer.atoms.pop(temp_atoms[atom])
    return output_file


def output_environment(monomer, output_path='', tag_monomer=False, strip_states=False, include_neighbours=True,
                       include_ligands=True, include_solvent=True, inter_chain=True, cutoff=4, side_chain=False):
    """ Output a monomer with its nearby residues as a PDB file.

    Parameters
    ----------
    monomer : AMPAL Monomer object
        Reference Monomer to export the environment of
    output_path : str or pathlib Path
        Destination for output file
    tag_monomer : bool
        If true, tags select atoms of the monomer with primes to facilitate identification/alignment
    strip_states : bool
        True to remove state labels to facilitate alignment of monomers
    include_neighbours : bool
        If false, does not include monomers at i-1, i+1 positions in same chain as Monomer
    inter_chain : bool
        If false, only includes nearby monomers in the same chain as the Monomer
    include_ligands : bool
        If true, Monomers classed as ligands but not identified as solvent will be included in the environment
    include_solvent : bool
        If true, Monomers classed as categorised as solvent will be included in the environment
    cutoff : float
        Maximum inter-atom distance for nearby Monomer to be included

    Returns
    -------
    output_file : pathlib Path
        Location of output file
    """
    if output_path:
        output_file = Path(output_path)
        if output_file.exists():
            return output_file
    if side_chain:
        if monomer.mol_code == 'GLY':
            monomer_list = [monomer]
        else:
            monomer_list = monomer.side_chain_environment(cutoff=cutoff, include_neighbours=include_neighbours,
                                                          inter_chain=inter_chain, include_ligands=include_ligands,
                                                          include_solvent=include_solvent)
    else:
        monomer_list = monomer.environment(cutoff=cutoff, include_self=True, include_neighbours=include_neighbours,
                                           inter_chain=inter_chain, include_ligands=include_ligands,
                                           include_solvent=include_solvent)
    if tag_monomer:
        tagged_monomer = monomer
    else:
        tagged_monomer = None
    output_file = output_with_alignment_tags(monomer_list=monomer_list, tagged_monomer=tagged_monomer,
                                             output_path=output_path, strip_states=strip_states)
    return output_file


def output_environment_spheres(monomer, output_path='', cutoff=4, include_neighbours=True, inter_chain=True,
                               strip_states=False, tag_chpi_acceptors=False, tag_chpi_donors=False, acceptor_codes=None,
                               donor_codes=None, side_chain=False):
    """ Output a monomer with the centres of nearby residue side-chains represented as non-bonded atoms.

    Parameters
    ----------
    monomer : [Monomer]
        Reference Monomer to export with nearby residues.
    output_path : str or [Path]
        Destination for output file.
    cutoff : float
        Maximum inter-atom distance for nearby Monomer to be included.
    include_neighbours : bool
        If False, does not include monomers at i-1, i+1 positions in same chain as Monomer.
    inter_chain : bool
        If False, only includes nearby monomers in the same chain as the Monomer.
    strip_states : bool
        True to remove state labels to facilitate alignment of monomers.
    tag_chpi_acceptors : bool
        If True, dummy atoms for residues acting as CH-pi acceptors to the monomer will have the atom name 'ACC'.
    tag_chpi_donors : bool
        If True, dummy atoms for residues acting as CH-pi donors to the monomer will have the atom name 'DON'.
    acceptor_codes : list
    donor_codes : list
        Lists of CH-pi donors and acceptors to include if they are tagged.


    Returns
    -------
    output_file : [Path]
        Location of output file.
    """
    if output_path:
        output_file = Path(output_path)
    else:
        output_file = Path(monomer.mol_code + monomer.id + '_environment.pdb')
        print("No output path specified, outputting to {0}".format(output_file))
    if output_file.exists():
        return output_file
    if side_chain:
        environment_residues = monomer.side_chain_environment(cutoff=cutoff, include_neighbours=include_neighbours,
                                                              inter_chain=inter_chain)
        if monomer in environment_residues:
            del(environment_residues[environment_residues.index(monomer)])
    else:
        environment_residues = monomer.environment(cutoff=cutoff, include_self=True,
                                                   include_neighbours=include_neighbours, inter_chain=inter_chain)
    acceptors = []
    donors = []
    if tag_chpi_acceptors:
        if not donor_codes or monomer.mol_code in donor_codes:
            for interaction in monomer.CH_pi_interactions(as_acceptor=False, acceptor_codes=acceptor_codes):
                if interaction.acceptor_monomer not in acceptors:
                    acceptors.append(interaction.acceptor_monomer)
    if tag_chpi_donors:
        if not acceptor_codes or monomer.mol_code in acceptor_codes:
            for interaction in monomer.CH_pi_interactions(as_donor=False, donor_codes=donor_codes):
                if interaction.donor_monomer not in donors:
                    donors.append(interaction.donor_monomer)
    dummy_environment = [monomer]
    for residue in environment_residues:
        side_chain_centre = residue._visual_scc
        if side_chain_centre is None:
            continue
        dummy_atom_parameters = {'coordinates': side_chain_centre,
                                 'element': 'X',
                                 'atom_id': '0'}
        if residue in acceptors:
            dummy_atom_parameters['res_label'] = 'ACC'
        elif residue in donors:
            dummy_atom_parameters['res_label'] = 'DON'
        else:
            dummy_atom_parameters['res_label'] = 'UNK'
        dummy_atom = Atom(**dummy_atom_parameters)
        dummy_ligand = Ligand(atoms={'A': {dummy_atom.res_label: dummy_atom}}, monomer_id=residue.id,
                              mol_code=residue.mol_code, ampal_parent=residue.ampal_parent)
        dummy_environment.append(dummy_ligand)
    output_pdb = write_pdb(dummy_environment, strip_states=strip_states)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    output_file.write_text(output_pdb)
    return output_file


def output_chpis_donor(monomer, output_path='', strip_states=False, acceptor_codes=None, inter_chain=True,
                       include_neighbours=True):
    """ Output a monomer with residues that act as CH-pi acceptors to it as a PDB file.

    Parameters
    ----------
    monomer : [Monomer]
        Monomer to export.
    output_path : str or [Path]
        Destination for output file.
    strip_states : bool
        True to remove state labels to facilitate alignment of monomers.
    acceptor_codes : list
        List of mol codes to be queried as CH-pi acceptors.
    inter_chain : bool
        If false, only includes monomers in the same chain as the Monomer as potential CH-pi acceptors.
    include_neighbours : bool
        If false, does not include monomers at i-1, i+1 positions in same chain as Monomer as potential CH-pi acceptors.

    Returns
    -------
    output_file : pathlib Path
        Location of output file
    """
    monomer_list = []
    for interaction in monomer.CH_pi_interactions(as_acceptor=False, acceptor_codes=acceptor_codes,
                                                  inter_chain=inter_chain):
        if interaction.acceptor in monomer_list:
            continue
        if interaction.acceptor.ampal_parent == monomer.ampal_parent and not include_neighbours:
            if int(interaction.acceptor.id) == int(monomer.id) + 1 or int(interaction.acceptor.id) == int(
                    monomer.id) - 1:
                continue
        monomer_list.append(interaction.acceptor)
    if len(monomer_list) == 0:
        return None
    output_path = output_with_alignment_tags(monomer_list=monomer_list, tagged_monomer=monomer, output_path=output_path,
                                             strip_states=strip_states)
    return output_path


def output_chpis_acceptor(monomer, output_path='', strip_states=False, donor_codes=None, inter_chain=True,
                          include_neighbours=True):
    """ Output a monomer with residues that act as CH-pi donors to it as a PDB file.

    Parameters
    ----------
    monomer : [Monomer]
        Monomer to export.
    output_path : str or [Path]
        Destination for output file.
    strip_states : bool
        True to remove state labels to facilitate alignment of monomers.
    donor_codes : list
        List of mol codes to be queried as CH-pi donors.
    inter_chain : bool
        If false, only includes monomers in the same chain as the Monomer as potential CH-pi acceptors.
    include_neighbours : bool
        If false, does not include monomers at i-1, i+1 positions in same chain as Monomer as potential CH-pi acceptors.

    Returns
    -------
    output_file : pathlib Path
        Location of output file
    """
    monomer_list = []
    for interaction in monomer.CH_pi_interactions(as_donor=False, donor_codes=donor_codes, inter_chain=inter_chain):
        if interaction.donor in monomer_list:
            continue
        if interaction.donor.ampal_parent == monomer.ampal_parent and not include_neighbours:
            if int(interaction.donor.id) == int(monomer.id) + 1 or int(interaction.donor.id) == int(monomer.id) - 1:
                continue
        monomer_list.append(interaction.donor)
    if len(monomer_list) == 0:
        return None
    output_path = output_with_alignment_tags(monomer_list=monomer_list, tagged_monomer=monomer, output_path=output_path,
                                             strip_states=strip_states)
    return output_path


__author__ = 'Kieran L. Hudson'
