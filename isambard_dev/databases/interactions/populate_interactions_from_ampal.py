from subprocess import CalledProcessError

from databases.general_tools import get_or_create
from databases.interactions.interaction_tables import *
from tools.components import component_category


def get_or_create_pdb(protein, session=interaction_session):
    """ Adds a protein imported as an AMPAL Polymer to the PDB table of the interactions database.

        Parameters
        ----------
        protein : AMPAL Polymer
        session : SQLAlchemy session
            A session to interface with the interactions database.

        Returns
        -------
        pdb_entry : PdbTable object
            The new or existing entry in the PDB table.
        new_pdb : bool
            True if entry had been newly added.
        """
    pdb_id = protein.id[:4].lower()
    pdb_entry = session.query(PdbTable).filter(PdbTable.pdb == pdb_id).one_or_none()
    if not pdb_entry:
        resolution = protein.structural_information['resolution']
        pdb_dict = {'pdb': pdb_id,
                    'resolution': resolution}
        pdb_entry, new_pdb = get_or_create(PdbTable, session=session, **pdb_dict)
        return pdb_entry, new_pdb
    if not pdb_entry.resolution:
        try:
            resolution = protein.structural_information['resolution']
            pdb_entry.resolution = resolution
        except AttributeError:
            pass
    return pdb_entry, False


def get_or_create_chain(chain, chain_id=None, session=interaction_session):
    """ Adds an AMPAL Chain to the Chain table of the interactions database.

        Notes
        -----
        Back-populates Pdb tables if required.

        Parameters
        ----------
        chain : AMPAL Chain
        chain_id : str
            Specify the chain ID if it can't be found from AMPAL, e.g., for a ligand.
        session : SQLAlchemy session
            A session to interface with the interactions database.

        Returns
        -------
        chain_entry : ChainTable object
            The new or existing entry in the Chain table.
        new_chain : bool
            True if entry had been newly added.
        """
    pdb_entry, new_pdb = get_or_create_pdb(chain.ampal_parent, session=session)
    chain_dict = {'chain': chain_id if chain_id else chain.id,
                  'pdb_id': pdb_entry.id}
    chain_entry, new_chain = get_or_create(ChainTable, session=session, **chain_dict)
    return chain_entry, new_chain


def get_or_create_component(monomer, session=interaction_session, ask_unknown=False):
    """ Adds component details of an AMPAL Monomer to the Component table of the interactions database.

        Parameters
        ----------
        monomer : AMPAL Monomer
        session : SQLAlchemy session
            A session to interface with the interactions database.
        ask_unknown : Bool
            Ask category of component if it is not known. Set as true if you later want to query by category.

        Returns
        -------
        comp_entry : ComponentTable object
            The new or existing entry in the Component table.
        new_comp : bool
            True if entry had been newly added.
        """
    comp_dict = {'category': component_category(monomer.mol_code, ask_unknown=ask_unknown),
                 'code': monomer.mol_code}
    comp_entry, new_comp = get_or_create(ComponentTable, session=session, **comp_dict)
    return comp_entry, new_comp


def get_or_create_monomer(monomer, session=interaction_session):
    """Adds an AMPAL Monomer to the Monomer table of the interactions database.

    Notes
    -----
    Back-populates Component, Chain, and Pdb tables if required.

    Parameters
    ----------
    monomer : AMPAL Monomer
    session : SQLAlchemy session
        A session to interface with the interactions database.

    Returns
    -------
    monomer_entry : MonomerTable object
        The new or existing entry in the Monomer table.
    new_monomer : bool
        True if entry had been newly added.
    """
    chain_id = monomer.ampal_parent.id
    chain_entry, new_chain = get_or_create_chain(monomer.ampal_parent, chain_id=chain_id, session=session)
    comp_entry, new_comp = get_or_create_component(monomer, session=session)
    monomer_dict = {'monomer': monomer.id,
                    'chain_id': chain_entry.id,
                    'component_id': comp_entry.id,
                    'insertion_code': monomer.insertion_code}
    monomer_entry, new_monomer = get_or_create(MonomerTable, session=session, **monomer_dict)
    return monomer_entry, new_monomer


def add_monomer_environment(monomer, cutoff=4, include_neighbours=True, inter_chain=True, include_ligands=False,
                            include_solvent=False, session=interaction_session):
    """ Adds details of an AMPAL Monomer's environment to the interactions database."""
    monomer_entry, new_monomer = get_or_create_monomer(monomer, session=session)
    for res in monomer.environment(cutoff=cutoff, include_neighbours=include_neighbours, inter_chain=inter_chain,
                                   include_ligands=include_ligands, include_solvent=include_solvent):
        res_entry, new_res = get_or_create_monomer(res, session=session)
        if res_entry not in monomer_entry.environment:
            monomer_entry.environment.append(res_entry)
    return monomer_entry.environment


def add_residue_sc_environment(residue, cutoff=4, include_neighbours=True, inter_chain=True, include_ligands=False,
                               include_solvent=False, session=interaction_session, sc_only=False):
    """ Adds details of an AMPAL Residue's side-chain environment to the interactions database."""
    residue_entry, new_residue = get_or_create_monomer(residue, session=session)
    for res in residue.side_chain_environment(cutoff=cutoff, include_neighbours=include_neighbours,
                                              inter_chain=inter_chain, include_ligands=include_ligands,
                                              include_solvent=include_solvent, sc_only=sc_only):
        if res is residue:
            continue
        res_entry, new_res = get_or_create_monomer(res, session=session)
        if res_entry not in residue_entry.environment:
            residue_entry.environment.append(res_entry)
    return residue_entry.environment


def get_or_create_amino_acid(residue, session=interaction_session):
    """ Adds an AMPAL Residue to the AminoAcid table of the interactions database.

        Notes
        -----
        Requires secondary structure, main-chain, and side-chain torsion angles to have been tagged.
        Back-populates Monomer, Component, Chain, and Pdb tables if required.

        Parameters
        ----------
        residue : [Residue]
            AMPAL Residue object
        session : SQLAlchemy session
            A session to interface with the interactions database.

        Returns
        -------
        amino_acid_entry : AminoAcidTable object
            The new or existing entry in the Ligand table.
        new_aa : bool
            True if entry had been newly added.

        Raises
        ------
        KeyError
            If the main-chain torsion angles have not been tagged prior to adding.
        """
    monomer_entry, new_monomer = get_or_create_monomer(residue)
    try:
        residue_dict = {'monomer_id': monomer_entry.id,
                        'phi': residue.tags['phi'],
                        'psi': residue.tags['psi'],
                        'omega': residue.tags['omega'],
                        'secondary_structure': residue.tags['secondary_structure'] if 'secondary_structure' in residue.tags
                                               else ' '}
        for x in range(0, len(residue.tags['chi_angles'])):
            residue_dict['chi_' + str(x+1)] = residue.tags['chi_angles'][x]
    except KeyError as k:
        print("One or more tags are missing. Ensure you have run tag_torsion_angles and tag_secondary_structure.")
        raise k
    amino_acid_entry, new_aa = get_or_create(AminoAcidTable, session=session, **residue_dict)
    return amino_acid_entry, new_aa


def add_chain_secondary_structure(chain, session=interaction_session):
    """ Add all amino acids in an AMPAL chain to the AminoAcid table, with requisite tags."""
    chain_aas = []
    chain.tag_torsion_angles()
    try:
        chain.tag_secondary_structure()
    # Catch for error message raised in Windows thwarting check_output when DSSP is run
    except CalledProcessError as c:
        print("DSSP could not be run. Secondary structure tags will be missing.\n", c)
    chain.tag_sidechain_dihedrals()
    for residue in chain:
        if residue.is_hetero:
            print(residue, 'is not a natural amino acid - skipping')
            continue
        amino_acid, new_aa = get_or_create_amino_acid(residue, session=session)
        chain_aas.append(amino_acid)
    return chain_aas


def get_or_create_ligand(ligand, session=interaction_session):
    """ Adds an AMPAL Ligand to the Ligand table of the interactions database.

        Notes
        -----
        Back-populates Monomer, Component, Chain, and Pdb tables if required.

        Parameters
        ----------
        ligand : AMPAL Ligand
        session : SQLAlchemy session
            A session to interface with the interactions database.

        Returns
        -------
        ligand_entry : LigandTable object
            The new or existing entry in the Ligand table.
        new_ligand : bool
            True if entry had been newly added.
        """
    monomer_entry, new_monomer = get_or_create_monomer(ligand, session=session)
    ligand_dict = {'monomer_id': monomer_entry.id}
    ligand_entry, new_ligand = get_or_create(LigandTable, session=session, **ligand_dict)
    return ligand_entry, new_ligand


def get_or_create_interaction(interaction, session=interaction_session):
    """Adds an AMPAL Interaction to the Interaction table of the interactions database.

    Notes
    -----
    Back-populates other tables as required through add_monomer.

    Parameters
    ----------
    interaction : AMPAL Interaction
    session : SQLAlchemy session
        A session to interface with the interactions database.

    Returns
    -------
    interaction_entry : InteractionTable object
        The new or existing entry in the Interaction table.
    new_interaction : bool
        True if entry had been newly added.
    """
    donor_entry, new_donor = get_or_create_monomer(interaction.donor_monomer, session=session)
    acceptor_entry, new_acceptor = get_or_create_monomer(interaction.acceptor_monomer, session=session)
    interaction_dict = {'donor_id': donor_entry.id,
                        'acceptor_id': acceptor_entry.id}
    interaction_entry, new_interaction = get_or_create(InteractionTable, session=session, **interaction_dict)
    return interaction_entry, new_interaction


def get_or_create_CH_pi(CH_pi, session=interaction_session):
    """ Adds an AMAPAL CH_pi object to the Interactions database, with full back-population.

    Parameters
    ----------
    CH_pi : AMPAL CH_pi Interaction
    session : SQLAlchemy session
        A session to interface with the interactions database.

    Returns
    -------
    CH_pi_entry : CHpiTable object
        The new or existing entry in the CH-pi table.
    new_CHpi : bool
        True if entry had been newly added.
    """
    interaction_entry, new_interaction = get_or_create_interaction(CH_pi, session=session)
    CH_pi_dict = {'donor_C': CH_pi.c,
                  'donor_H': CH_pi.h,
                  'pi_system': CH_pi.pi_system,
                  'distance': CH_pi.distance,
                  'proj_dist': CH_pi.proj_dist,
                  'angle': CH_pi.angle,
                  'interaction_id': interaction_entry.id}
    CH_pi_entry, new_CHpi = get_or_create(CHpiTable, session=session, **CH_pi_dict)
    return CH_pi_entry, new_CHpi


def get_or_create_H_bond(H_bond, session=interaction_session):
    """ Adds an AMAPAL H_bond object to the Interactions database, with full back-population.

        Parameters
        ----------
        H_bond : AMPAL CH_pi Interaction
        session : SQLAlchemy session
        A session to interface with the interactions database.

        Returns
        -------
        H_bond_entry : HbondTable object
        The new or existing entry in the H-bond table.
        new_Hbond : bool
        True if entry had been newly added.
        """
    interaction_entry, new_interaction = get_or_create_interaction(H_bond, session=session)
    H_bond_dict = {'acceptor_atom': H_bond.acceptor.res_label,
                   'donor_atom': H_bond.donor.res_label,
                   'distance': H_bond.dist,
                   'angle_a': H_bond.ang_a,
                   'angle_d': H_bond.ang_d,
                   'interaction_id': interaction_entry.id}
    H_bond_entry, new_Hbond = get_or_create(HbondTable, session=session, **H_bond_dict)
    return H_bond_entry, new_Hbond


__author__ = 'Kieran L. Hudson, Ali Scott'
