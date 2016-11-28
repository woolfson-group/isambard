from pathlib import Path

from databases.general_tools import get_or_create
from databases.interactions.interaction_tables import *
from tools.file_parsing import parse_PISCES_output
from tools.components import component_category


def populate_from_pisces(pisces_output, path=False, session=interaction_session):
    """ Adds the PDB and chain data from a PISCES output to the interactions database.

    Parameters
    ----------
    pisces_output : str or Path
        Text output from PISCES as a string, or path to output file.
    path : bool
        True if giving a path to a file.
    session : SQLAlchemy session
        Session for interacting with the interactions database.

    Returns
    -------
    outputs : dict
        Dict of chain entries for each PDB entry in the Interactions database from the PISCES output.
    """
    outputs = {}
    pisces_dict = parse_PISCES_output(pisces_output, path=path)
    for pdb in pisces_dict:
        pdb_args = {'pdb': pdb,
                    'resolution': pisces_dict[pdb]['resolution']}
        pdb_entry, new_pdb = get_or_create(PdbTable, session=session, **pdb_args)
        for chain in pisces_dict[pdb]['chains']:
            chain_args = {'chain': chain,
                          'pdb_id': pdb_entry.id}
            chain_entry, new_chain = get_or_create(ChainTable, session=session, **chain_args)
            if pdb_entry in outputs:
                outputs[pdb_entry].append(chain_entry)
            else:
                outputs[pdb_entry] = [chain_entry]
    return outputs


def populate_privateer(privateer_output, path=True, session=interaction_session, only_ligands=True):
    """Populates the Privateer table in the ligand database from a Privateer output

    Notes
    -----
    Also back-populates Ligand, Monomer, Component, Chain, and Pdb tables if required.

    Parameters
    ----------
    privateer_output: path to Privateer output file or str
        The tab-delimited list of Privateer outputs.
    path: bool
        True if path to file, false if opened file as str.
    session: Session
        An SQLAlchemy session for the database.
    only_ligands : bool
        True if all residues in Privateer output are ligands.

    Returns
    -------
    new_privateer: list
        A list of the new entries for Privateer table.
    """
    if path:
        privateer_path = Path(privateer_output)
        privateer_list = privateer_path.read_text().splitlines()
    else:
        privateer_list = privateer_output.splitlines()
    outputs = []
    for line in privateer_list:
        items = line.split()
        pdb_dict = {'pdb': items[0].lower(),
                    'resolution': float(items[2])}
        pdb_entry, new_pdb = get_or_create(PdbTable, session=session, **pdb_dict)
        sugar_id = items[1].split('-')
        chain_dict = {'chain': sugar_id[1][0],
                      'pdb_id': pdb_entry.id}
        chain_entry, new_chain = get_or_create(ChainTable, session=session, **chain_dict)
        component_dict = {'code': sugar_id[0],
                          'category': component_category(sugar_id[0])}
        component_entry, new_component = get_or_create(ComponentTable, session=session, **component_dict)
        insertion_code = ' '
        if sugar_id[-2] == '':
            res_no = '-' + sugar_id[-1]
        elif len(sugar_id[-1].split(':')) == 2:
            res_no = sugar_id[-1].split(':')[0]
            insertion_code = sugar_id[-1].split(':')[-1]
        else:
            res_no = sugar_id[-1]
        monomer_dict = {'monomer': res_no,
                        'chain_id': chain_entry.id,
                        'component_id': component_entry.id,
                        'insertion_code': insertion_code}
        monomer_entry, new_monomer = get_or_create(MonomerTable, session=session, **monomer_dict)
        if only_ligands:
            ligand_dict = {'monomer_id': monomer_entry.id}
            get_or_create(LigandTable, session=session, **ligand_dict)
        type_detected = items[7].split('-')
        privateer_dict = {'monomer_id': monomer_entry.id,
                          'q': float(items[3]),
                          'phi': float(items[4]),
                          'theta': False if items[5] == '--' else float(items[5]),
                          'rscc': float(items[6]),
                          'anomer': type_detected[0],
                          'configuration': type_detected[1],
                          'cho_type': type_detected[2][:4] + 'se',
                          'ring': type_detected[2][4:],
                          'conformation': items[8],
                          'mean_density': float(items[9]),
                          'mean_b_factor': float(items[10]),
                          'bond_length_deviation': float(items[11]),
                          'bond_angle_deviation': float(items[12]),
                          'environment': items[13][1],
                          'diagnostic': items[14]}
        privateer_entry, new_privateer = get_or_create(PrivateerTable, session=session, **privateer_dict)
        outputs.append(privateer_entry)
    return outputs


__author__ = 'Kieran L. Hudson'
