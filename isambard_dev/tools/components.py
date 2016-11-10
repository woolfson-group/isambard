from pathlib import Path
from sqlalchemy import create_engine, Column, Integer, String, Float, ForeignKey, Boolean, Table
from sqlalchemy.orm import sessionmaker, relationship
from sqlalchemy.ext.declarative import declarative_base

from settings import global_settings
from tools.file_parsing import dict_from_mmcif, download_decode
from databases.general_tools import get_or_create


try:
    db_path = Path(global_settings['package_path'], 'tools', 'component_properties.db')
    engine_name = 'sqlite:///' + str(db_path)
    component_engine = create_engine(engine_name, echo=False)
    ComponentBase = declarative_base()
    ComponentSession = sessionmaker(bind=component_engine)
    component_session = ComponentSession()
finally:
    component_session.close()


# Used for H-bond identification
hbond_heteroatoms = ('N', 'O', 'S', 'F')


modified_table = Table('modified_component', ComponentBase.metadata,
                       Column('modified_id', Integer, ForeignKey('component.id')),
                       Column('parent_id', Integer, ForeignKey('component.id')))


class ComponentTable(ComponentBase):
    __tablename__ = 'component'

    id = Column(Integer, primary_key=True)
    code = Column(String(3), nullable=False, unique=True)
    name = Column(String(255))
    charge = Column(Integer)
    formula = Column(String(255))
    formula_weight = Column(Float)
    cif_category = Column(String(255))
    category = Column(String(255))
    hbond_donors_added = Column(Boolean)
    hbond_acceptors_added = Column(Boolean)
    chbond_donors_added=Column(Boolean)

    parent = relationship('ComponentTable', secondary=modified_table, uselist=False,
                          primaryjoin=modified_table.c.modified_id == id,
                          secondaryjoin=modified_table.c.parent_id == id)

    def __repr__(self):
        return '<Component(code={0}, name={1}, category={2})>'.format(self.code, self.name, self.category)


class PiSystemTable(ComponentBase):
    __tablename__ = 'pi_system'

    id = Column(Integer, primary_key=True, unique=True)
    component_id = Column(Integer, ForeignKey('component.id'), nullable=False)
    pi_system = Column(String(10), nullable=False)

    component = relationship('ComponentTable', backref='pi_systems')

    def __repr__(self):
        return '<Pi system(label={0}, component={1})>'.format(self.pi_system, self.component.code)


class ElementTable(ComponentBase):
    __tablename__ = 'element'

    id = Column(Integer, primary_key=True, unique=True)
    symbol = Column(String(3), unique=True, nullable=False)
    electronegativity = Column(Float)
    CPK = Column(String(10))
    group = Column(Integer, nullable=False)
    name = Column(String(255), nullable=False)
    atomic_mass = Column(Float)
    density = Column(Float)
    EA = Column(Integer)
    ion_radius = Column(Integer)
    ion_charge = Column(Integer)
    bonding_type = Column(String(255))
    electronic_config = Column(String(255), nullable=False)
    metal = Column(String(255))
    period = Column(Integer, nullable=False)
    IE = Column(Integer)
    atomic_number = Column(Integer, nullable=False)
    standard_state = Column(String(255))
    vdW_radius = Column(Integer)
    boiling_point = Column(Integer)
    year_discovered = Column(String(255))
    atomic_radius = Column(Integer)
    melting_point = Column(Integer)

    def __repr__(self):
        return '<Element(symbol={0}, atomic number={1})>'.format(self.symbol, self.atomic_number)


class AtomTable(ComponentBase):
    __tablename__ = 'atom'

    id = Column(Integer, primary_key=True, unique=True)
    component_id = Column(Integer, ForeignKey('component.id'), nullable=False)
    element_id = Column(Integer, ForeignKey('element.id'), nullable=False)
    atom_label = Column(String(10), nullable=False)
    alt_atom_label = Column(String(4))
    charge = Column(Integer)
    aromatic = Column(Boolean)
    leaving_atom = Column(Boolean)
    stereo_config = Column(String(1))

    component = relationship('ComponentTable', backref='atoms')
    element = relationship('ElementTable', backref='atoms')

    def __repr__(self):
        return '<Atom(atom_id={0}, component={1})>'.format(self.atom_label, self.component.code)


class BondTable(ComponentBase):
    __tablename__ = 'bond'

    id = Column(Integer, primary_key=True, unique=True)
    component_id = Column(Integer, ForeignKey('component.id'), nullable=False)
    atom_1_id = Column(Integer, ForeignKey('atom.id'), nullable=False)
    atom_2_id = Column(Integer, ForeignKey('atom.id'), nullable=False)
    order = Column(Integer, nullable=False)
    aromatic = Column(Boolean)
    stereo_config = Column(String(1))

    component = relationship('ComponentTable', backref='covalent_bonds')
    atom_1 = relationship('AtomTable', foreign_keys=[atom_1_id], backref='covalent_bonds_1')
    atom_2 = relationship('AtomTable', foreign_keys=[atom_2_id], backref='covalent_bonds_2')

    def __repr__(self):
        return '<Covalent bond(atom 1={0}, atom 2={1}, component={2})>'.format(
            self.atom_1.atom_label, self.atom_2.atom_label, self.component.code)


class AlignmentTable(ComponentBase):
    __tablename__ = 'alignment'

    id = Column(Integer, primary_key=True, unique=True)
    component_id = Column(Integer, ForeignKey('component.id'), nullable=False)
    atom_id = Column(Integer, ForeignKey('atom.id'), nullable=False)

    component = relationship('ComponentTable', backref='alignment_atoms')
    atom = relationship('AtomTable', backref='alignment')

    def __repr__(self):
        return '<Alignment atom(atom_id={0}, component={1})>'.format(self.atom.atom_label, self.component.code)


class SideChainCentreTable(ComponentBase):
    __tablename__ = 'side_chain_centre'

    id = Column(Integer, primary_key=True, unique=True)
    component_id = Column(Integer, ForeignKey('component.id'), nullable=False)
    atom_id = Column(Integer, ForeignKey('atom.id'), nullable=False)

    component = relationship('ComponentTable', backref='side_chain_centre_atoms')
    atom = relationship('AtomTable', backref='side_chain_centre')

    def __repr__(self):
        return '<Side chain centre atom(atom_id={0}, component={1})>'.format(self.atom.atom_label, self.component.code)


class PiSystemAtomTable(ComponentBase):
    __tablename__ = 'pi_system_atoms'

    id = Column(Integer, primary_key=True, unique=True)
    pi_system_id = Column(Integer, ForeignKey('pi_system.id'), nullable=False)
    atom_id = Column(Integer, ForeignKey('atom.id'), nullable=False)

    atom = relationship('AtomTable', backref='pi_systems')
    pi_system = relationship('PiSystemTable', backref='atoms')

    def __repr__(self):
        return '<Pi system atom(atom={0}, pi system={1}, component={2})>'.\
            format(self.atom.atom_label, self.pi_system.pi_system, self.pi_system.component.code)


class HBondDonorTable(ComponentBase):
    __tablename__ = 'hbond_donor'

    id = Column(Integer, primary_key=True, unique=True)
    component_id = Column(Integer, ForeignKey('component.id'), nullable=False)
    proton_id = Column(Integer, ForeignKey('atom.id'), nullable=False)
    heteroatom_id = Column(Integer, ForeignKey('atom.id'), nullable=False)

    component = relationship('ComponentTable', backref='hbond_donors')
    proton = relationship('AtomTable', foreign_keys=[proton_id], backref='hbond_donor')
    heteroatom = relationship('AtomTable', foreign_keys=[heteroatom_id], backref='donor_protons')

    def __repr__(self):
        return '<H-bond donor(proton={0}, heteroatom={1}, component={2})>'.format(
            self.proton.atom_label, self.heteroatom.atom_label, self.component.code)

class CHBondDonorTable(ComponentBase):

    __tablename__ = 'chbond_donor'

    id = Column(Integer,primary_key=True,unique=True)
    component_id = Column(Integer, ForeignKey('component.id'), nullable=False)
    proton_id = Column(Integer,ForeignKey('atom.id'), nullable=False)
    heteroatom_id = Column(Integer,ForeignKey('atom.id'),nullable=False)

    component = relationship('ComponentTable',backref='chbond_donors')
    proton = relationship('AtomTable',foreign_keys=[proton_id], backref='chbond_donor')
    heteroatom = relationship('AtomTable', foreign_keys=[heteroatom_id],backref='chbond_donor_protons')

    def __repr__(self):
        return '<CH-bond donor(proton={0}, heteroatom={1}, component={2})>'.format(
            self.proton.atom_label, self.heteroatom.atom_label, self.component.code)

class HBondAcceptorTable(ComponentBase):
    __tablename__ = 'hbond_acceptor'

    id = Column(Integer, primary_key=True, unique=True)
    component_id = Column(Integer, ForeignKey('component.id'), nullable=False)
    heteroatom_id = Column(Integer, ForeignKey('atom.id'), nullable=False)
    bound_atom_id = Column(Integer, ForeignKey('atom.id'))

    component = relationship('ComponentTable', backref='hbond_acceptors')
    heteroatom = relationship('AtomTable', foreign_keys=[heteroatom_id], backref='hbond_acceptor')
    bound_atom = relationship('AtomTable', foreign_keys=[bound_atom_id], backref='bound_hbond_acceptor')

    def __repr__(self):
        return '<H-bond acceptor(heteroatom={0}, bound atom={1}, component={2})>'.format(
            self.heteroatom.atom_label, self.bound_atom.atom_label if self.bound_atom else None, self.component.code)


def download_cif(ligand_code):
    """ Downloads ligand data from PDB website as mmcif, if it exists."""
    url = 'http://ligand-expo.rcsb.org/reports/' + ligand_code[0] + '/' + ligand_code + '/' + ligand_code + '.cif'
    cif_file = download_decode(url)
    return cif_file


def get_cif(mol_code, cif_path=''):
    """ Retrieves the mmcif for a component, either downloaded or from the chemical component dictionary.

    Notes
    -----
    Makes use of the chemical component dictionary as a mmCIF file if it's in the default location or at the specified
    path. This can be downloaded from http://ligand-expo.rcsb.org/dictionaries/Components-pub.cif.
    Otherwise tries to download from the PDB website.

    Parameters
    ----------
    mol_code : str
        Standard (up to) three-letter code for components (amino acid, ligand, solvent, etc.).
    cif_path : str or pathlib.Path
        Location of component mmcif file, if not default.

    Returns
    -------
    cif_file : str
        The mmcif file with the data for a component from the PDB.
    """
    if cif_path:
        all_component_cif = Path(cif_path)
    else:
        all_component_cif = Path(global_settings['structural_database']['path'], 'Components-pub.cif')
    if all_component_cif.exists():
        all_components = all_component_cif.read_text().split('\n\n#')
        all_components_dict = {x.split()[0][5:]: x for x in all_components}
        try:
            cif_file = ('#' + all_components_dict[mol_code.upper()])
        except KeyError:
            print("{0} is not in the database. Try updating the components cif file at {1}."
                  .format(mol_code, all_component_cif))
            cif_file = download_cif(mol_code.upper())
    else:
        print("The cif database was not found at {0}".format(all_component_cif))
        cif_file = download_cif(mol_code.upper())
    if not cif_file:
        print("The mmcif file corresponding to component {0} could not be found.".format(mol_code))
        return None
    return cif_file


def get_component_dict(mol_code, cif_path=''):
    """ Retrieves relevant data for a component from the PDB ligand explorer.

    Notes
    -----
    Only needs to be run for new labels not in component_properties.db.
    Data for known components should be obtained by querying that table.

    Parameters
    ----------
    mol_code : str
        Standard (up to) three-letter code for components (amino acid, ligand, solvent, etc.).
    cif_path : str or pathlib.Path
        Location of component mmcif file, if not default.

    Returns
    -------
    component_dict : dict
        A subset of the data in the component cif file.
    """
    cif_file = get_cif(mol_code, cif_path=cif_path)
    if not cif_file:
        return {}
    full_dict = dict_from_mmcif(cif_file, path=False)
    component_dict_names = {'name': '_chem_comp.name',
                            'charge': '_chem_comp.pdbx_formal_charge',
                            'formula_weight': '_chem_comp.formula_weight',
                            'cif_category': '_chem_comp.type',
                            'formula': '_chem_comp.formula',
                            'modified': '_chem_comp.mon_nstd_parent_comp_id'}
    component_dict = {x: full_dict[component_dict_names[x]].strip('\"') for x in component_dict_names if
                      component_dict_names[x] in full_dict}
    component_dict['code'] = mol_code
    if '_chem_comp_atom.comp_id' not in full_dict:
        component_dict['atoms'] = []
    elif isinstance(full_dict['_chem_comp_atom.comp_id'], str):
        component_dict['atoms'] = [{'atom_label': full_dict['_chem_comp_atom.atom_id'].strip('\"'),
                                    'alt_atom_label': full_dict['_chem_comp_atom.alt_atom_id'].strip('\"'),
                                    'element': full_dict['_chem_comp_atom.type_symbol'],
                                    'charge': full_dict['_chem_comp_atom.charge'],
                                    'aromatic': True if full_dict['_chem_comp_atom.pdbx_aromatic_flag'] == 'Y'
                                                else False,
                                    'leaving_atom': True if full_dict['_chem_comp_atom.pdbx_leaving_atom_flag'] == 'Y'
                                                    else False,
                                    'stereo_config': None if full_dict['_chem_comp_atom.pdbx_stereo_config'] == 'N'
                                                     else full_dict['_chem_comp_atom.pdbx_stereo_config']}]
    else:
        component_dict['atoms'] = [{'atom_label': full_dict['_chem_comp_atom.atom_id'][x].strip('\"'),
                                    'alt_atom_label': full_dict['_chem_comp_atom.alt_atom_id'][x].strip('\"'),
                                    'element': full_dict['_chem_comp_atom.type_symbol'][x],
                                    'charge': full_dict['_chem_comp_atom.charge'][x],
                                    'aromatic': True if full_dict['_chem_comp_atom.pdbx_aromatic_flag'][x] == 'Y'
                                                else False,
                                    'leaving_atom': True if full_dict['_chem_comp_atom.pdbx_leaving_atom_flag'][x] == 'Y'
                                                    else False,
                                    'stereo_config': None if full_dict['_chem_comp_atom.pdbx_stereo_config'][x] == 'N'
                                                    else full_dict['_chem_comp_atom.pdbx_stereo_config'][x]}
                                   for x in range(0, len(full_dict['_chem_comp_atom.comp_id']))]
    if '_chem_comp_bond.comp_id' not in full_dict:
        component_dict['bonds'] = []
    elif isinstance(full_dict['_chem_comp_bond.comp_id'], str):
        component_dict['bonds'] = [{'atom_1': full_dict['_chem_comp_bond.atom_id_1'].strip('\"'),
                                    'atom_2': full_dict['_chem_comp_bond.atom_id_2'].strip('\"'),
                                    'order': 1 if full_dict['_chem_comp_bond.value_order'] == 'SING'
                                             else 2 if full_dict['_chem_comp_bond.value_order'] == 'DOUB'
                                             else 3 if full_dict['_chem_comp_bond.value_order'] == 'TRIP'
                                             else None,
                                    'aromatic': True if full_dict['_chem_comp_bond.pdbx_aromatic_flag'] == 'Y'
                                                else False,
                                    'stereo_config': None if full_dict['_chem_comp_bond.pdbx_stereo_config'] == 'N'
                                                     else full_dict['_chem_comp_bond.pdbx_stereo_config']}]
    else:
        component_dict['bonds'] = [{'atom_1': full_dict['_chem_comp_bond.atom_id_1'][x].strip('\"'),
                                    'atom_2': full_dict['_chem_comp_bond.atom_id_2'][x].strip('\"'),
                                    'order': 1 if full_dict['_chem_comp_bond.value_order'][x] == 'SING'
                                             else 2 if full_dict['_chem_comp_bond.value_order'][x] == 'DOUB'
                                             else 3 if full_dict['_chem_comp_bond.value_order'][x] == 'TRIP'
                                             else None,
                                    'aromatic': True if full_dict['_chem_comp_bond.pdbx_aromatic_flag'][x] == 'Y'
                                                else False,
                                    'stereo_config': None if full_dict['_chem_comp_bond.pdbx_stereo_config'][x] == 'N'
                                                     else full_dict['_chem_comp_bond.pdbx_stereo_config'][x]}
                                    for x in range(0, len(full_dict['_chem_comp_bond.comp_id']))]
    return component_dict


def add_component(mol_code, session=component_session):
    """ Adds a component to the component properties database if it's new.

    Parameters
    ----------
    mol_code : str
    session : SQLAlchemy session
        Session linking to database storing component properties.

    Returns
    -------
    entry : ComponentTable entry
        Entry added or retrieved to the component properties table.
    bond_entries : list or None
        List of covalent bonds added for the component, if any.
    """
    component_dict = get_component_dict(mol_code)
    if not component_dict:
        return (None, False), [], []
    short_comp_dict = dict(component_dict)
    if 'modified' in short_comp_dict:
        if short_comp_dict['modified'] != '?':
            parent_comp_entry = session.query(ComponentTable).filter(
                ComponentTable.code == short_comp_dict['modified']).one_or_none()
            if not parent_comp_entry:
                (parent_comp_entry, new_parent_comp), all_parent_comp_atoms, all_parent_comp_bonds =\
                    add_component(short_comp_dict['modified'], session=session)
            if parent_comp_entry:
                short_comp_dict['parent'] = parent_comp_entry
                short_comp_dict['category'] = parent_comp_entry.category
        del(short_comp_dict['modified'])
    del (short_comp_dict['atoms'], short_comp_dict['bonds'])
    component_entry, new_component = get_or_create(ComponentTable, session=session, **short_comp_dict)
    all_atoms = []
    for atom_dict in component_dict['atoms']:
        atom_dict['component_id'] = component_entry.id
        element_entry = session.query(ElementTable).filter(ElementTable.symbol == atom_dict['element'].title()).one()
        atom_dict['element_id'] = element_entry.id
        del (atom_dict['element'])
        atom_entry, new_atom = get_or_create(AtomTable, session=session, **atom_dict)
        all_atoms.append((atom_entry, new_atom))
    all_bonds = []
    for bond_dict in component_dict['bonds']:
        bond_dict['component_id'] = component_entry.id
        atom_1_entry = session.query(AtomTable).filter(AtomTable.atom_label == bond_dict['atom_1']) \
            .join(ComponentTable).filter(ComponentTable.code == mol_code).one()
        bond_dict['atom_1_id'] = atom_1_entry.id
        atom_2_entry = session.query(AtomTable).filter(AtomTable.atom_label == bond_dict['atom_2']) \
            .join(ComponentTable).filter(ComponentTable.code == mol_code).one()
        bond_dict['atom_2_id'] = atom_2_entry.id
        del (bond_dict['atom_1'], bond_dict['atom_2'])
        bond_entry, new_bond = get_or_create(BondTable, session=session, **bond_dict)
        all_bonds.append((bond_entry, new_bond))
    return (component_entry, new_component), all_atoms, all_bonds


def component_category(mol_code, session=component_session, ask_unknown=True):
    """ Queries the component database for manually entered category, or asks user to assign one.

    Parameters
    ----------
    mol_code : str
    session : SQLAlchemy session
        Session for interacting with component properties database.
    ask_unknown : Bool
        If true will ask user to assign category to compounds for which one isn't known.

    Returns
    -------
    component_category : str or None
        Found category, entered category, or None if none was found or entered.
    """
    component = session.query(ComponentTable).filter(ComponentTable.code == mol_code).one_or_none()
    if not component:
        (component, new_component), all_atoms, all_bonds = add_component(mol_code, session=session)
        session.commit()
    if not component.category and ask_unknown:
        known_categories = [x for y in session.query(ComponentTable.category).distinct() for x in y if x is not None]
        category = input("Category unknown. What category is " + component.name + "? Known categories: " +
                         ', '.join(known_categories) + ". The cif file says " + component.cif_category +
                         ". (Leave blank to exit.): ").lower()
        if not category:
            return None
        component.category = category
        session.commit()
    return component.category


def component_covalent_bonds(mol_code, atom_a='', atom_b='', session=component_session):
    """ Returns a dictionary of covalent bonds in a component, between particular atoms if specified.

    Notes
    -----
    Specify atom_a and atom_b to find specific covalent bonds, e.g. atom_a='C' and atom_b='H' gives all C-H bonds, or
    leave blank for all covalent bonds.
    Partial queries also work - just atom_a='C' will return all bonds involving carbon atoms, just atom_b='H' will give
    all atoms bound to a proton, etc.
    Will try and populate component properties database if the component is missing.

    Parameters
    ----------
    mol_code : str
    atom_a : optional str
        Start of atoms to be returned as keys, e.g. 'C' for carbons or 'CA'.
    atom_b : optional str
        Start of atoms to be returned in value lists.
    session : SQLAlchemy session
        Session linking to database storing component properties.

    Returns
    -------
    bond_dict : dict
        Dictionary of covalent bonds queried.
    """
    bond_dict = {}
    component = session.query(ComponentTable).filter(ComponentTable.code == mol_code).one_or_none()
    if not component:
        (component, new_component), all_atoms, all_bonds = add_component(mol_code, session=session)
        session.commit()
        if not component:
            print("No component corresponding to {0} is known.".format(mol_code))
            return bond_dict
    if not component.covalent_bonds:
        return bond_dict
    for bond in component.covalent_bonds:
        if bond.atom_1.atom_label.startswith(atom_a) and bond.atom_2.atom_label.startswith(atom_b):
            if bond.atom_1.atom_label in bond_dict:
                bond_dict[bond.atom_1.atom_label].append(bond.atom_2.atom_label)
            else:
                bond_dict[bond.atom_1.atom_label] = [bond.atom_2.atom_label]
        elif bond.atom_2.atom_label.startswith(atom_a) and bond.atom_1.atom_label.startswith(atom_b):
            if bond.atom_2.atom_label in bond_dict:
                bond_dict[bond.atom_2.atom_label].append(bond.atom_1.atom_label)
            else:
                bond_dict[bond.atom_2.atom_label] = [bond.atom_1.atom_label]
    return bond_dict


def ch_bond_dict(mol_codes=None, session=component_session):
    if not mol_codes:
        mol_codes = [x for y in session.query(ComponentTable.code).all() for x in y]
    ch_bonds = {x: component_covalent_bonds(x, atom_a='C', atom_b='H') for x in mol_codes}
    return(ch_bonds)


def get_alignment_dict(mol_code, session=component_session):
    """ Returns the atoms to be labelled for alignment in output files."""
    mol_code = mol_code.upper()
    component = session.query(ComponentTable).filter(ComponentTable.code == mol_code).one_or_none()
    if not component:
        print("{0} is not a known component.".format(mol_code))
        return {}
    if not component.alignment_atoms:
        print("No known alignment atoms for {0}.".format(mol_code))
        return {}
    alignment_dict = {x.atom.atom_label: x.atom.atom_label + "\'" for x in component.alignment_atoms}
    return alignment_dict


def side_chain_centre_atoms(mol_code, session=component_session):
    """ Returns the atoms deemed functional in amino acid side chain."""
    mol_code = mol_code.upper()
    component = session.query(ComponentTable).filter(ComponentTable.code == mol_code).one_or_none()
    if not component:
        print("{0} is not a known component.".format(mol_code))
        return []
    if not component.side_chain_centre_atoms:
        print("No known side chain centre atoms for {0}.".format(mol_code))
        return []
    scc_atoms = {x.atom.atom_label for x in component.side_chain_centre_atoms}
    return scc_atoms


def component_pi_systems(mol_code, session=component_session):
    """ Returns pi system labels and list of atoms for components."""
    pi_systems = session.query(PiSystemTable).join(ComponentTable).\
                 filter(ComponentTable.code == mol_code).all()
    return {x.pi_system: [y.atom.atom_label for y in x.atoms] for x in pi_systems}


def known_pi_systems(session=component_session):
    """ Gives all of the known pi systems in the component database, with the constituent atoms."""
    comp_with_pi = session.query(ComponentTable).join(PiSystemTable).all()
    return {x.code: {y.pi_system: [z.atom. atom_label for z in y.atoms] for y in x.pi_systems} for x in comp_with_pi}


def add_hbond_donors(mol_code, session=component_session):
    """ Populate the hbond_donor table for a given component.

    Note
    ----
    This must be run for new components to allow H-bonds to be found for that component.

    Parameters
    ----------
    mol_code : str
        Three-letter code for a component.
    session : sqlalchemy.orm.session.Session
        SQLAlchemy session connecting to components database.

    Returns
    -------
    hbond_entries : list
        List of H-bond entries for component.
    """
    hbond_entries = []
    comp_entry = session.query(ComponentTable).filter(ComponentTable.code == mol_code).one_or_none()
    if not comp_entry:
        (comp_entry, new_component), all_atoms, all_bonds = add_component(mol_code, session=session)
        session.commit()
        if not comp_entry:
            print("{0} not a known component.".format(mol_code))
            return hbond_entries
    donor_bonds = component_covalent_bonds(mol_code, 'H', hbond_heteroatoms)
    for bond in donor_bonds:
        proton_entry = session.query(AtomTable).filter(AtomTable.component_id == comp_entry.id).\
                       filter(AtomTable.atom_label == bond).one()
        heteroatom_entry = session.query(AtomTable).filter(AtomTable.component_id == comp_entry.id).\
                           filter(AtomTable.atom_label == donor_bonds[bond][0]).one()
        bond_dict = {'component_id': comp_entry.id,
                     'proton_id': proton_entry.id,
                     'heteroatom_id': heteroatom_entry.id}
        hbond_donor_entry, new_hbond_donor = get_or_create(HBondDonorTable, session=session, **bond_dict)
        hbond_entries.append(hbond_donor_entry)
    comp_entry.hbond_donors_added = True
    return hbond_entries

def add_chbond_donors(mol_code,session=component_session):
    """Populate the chbond donor table for a given component
    Notes
    -----
    This must be run for new components to allow CH-bonds to be found for that component.

    Parameters
    ----------
    mol_code: str
        Three-letter code for a component.
    session : sqlalchemy.orm.session.Session
        SQLAlchemy session connecting to components database

    Returns
    -------
    chbond_entries : list
        List of CH-bond entries for that component
    """

    chbond_entries = []
    chbond_heteroatoms = ('C')
    comp_entry = session.query(ComponentTable).filter(ComponentTable.code==mol_code).one_or_none()
    if not comp_entry:
        (comp_entry, new_component), all_atoms, ll_bonds = add_component(mol_code,session=session)
        session.commit()
        if not comp_entry:
            print ("{0} not a known component.".format(mol_code))
            return chbond_entries
    donor_chbonds = component_covalent_bonds(mol_code,'H',chbond_heteroatoms)
    for bond in donor_chbonds:
        proton_entry = session.query(AtomTable).filter(AtomTable.component_id == comp_entry.id).\
            filter(AtomTable.atom_label==bond).one()
        heteroatom_entry = session.query(AtomTable).filter(AtomTable.component_id == comp_entry.id).\
            filter(AtomTable.atom_label==donor_chbonds[bond][0]).one()
        chbond_dict={'component_id' : comp_entry.id,
                     'proton_id' : proton_entry.id,
                     'heteroatom_id' : heteroatom_entry.id}

        chbond_donor_entry, new_chbond_donor = get_or_create(CHBondDonorTable, session=session, **chbond_dict)
        chbond_entries.append(chbond_donor_entry)
    comp_entry.chbond_donors_added = True
    return chbond_entries

def get_hbond_donors(mol_code, session=component_session):
    """ Retrieve or populate and return the H-bond donors for a component."""
    comp_entry = session.query(ComponentTable).filter(ComponentTable.code == mol_code).one_or_none()
    if not comp_entry:
        (comp_entry, new_component), all_atoms, all_bonds = add_component(mol_code, session=session)
        session.commit()
        if not comp_entry:
            print("{0} not a known component.".format(mol_code))
            return {}
    if not bool(comp_entry.hbond_donors_added):
        print("Finding H-bond donors for new component {0}.".format(mol_code))
        donors = add_hbond_donors(mol_code, session=session)
    else:
        donors = comp_entry.hbond_donors
    return {x.proton.atom_label: x.heteroatom.atom_label for x in donors}

def get_chbond_donors(mol_code, session=component_session):
    """Retrieve or populate and return CH-bond donors for a component."""
    comp_entry = session.query(ComponentTable).filter(ComponentTable.code == mol_code).one_or_none()
    if not comp_entry:
        (comp_entry, new_component), all_atoms, all_bonds = add_component(mol_code, session=session)
        session.commit()
        if not comp_entry:
            print("{0} not a known component.".format(mol_code))
            return{}
    if not bool(comp_entry.chbond_donors_added):
        print ("Finding CH-bond donors for new component {0}.".format(mol_code))
        donors = add_chbond_donors(mol_code,session=session)
    else:
        donors = comp_entry.chbond_donors

    return{x.proton.atom_label: x.heteroatom.atom_label for x in donors}


def hbond_donor_dict(mol_codes=[], session=component_session):
    """ Return the H-bond donors for all components in the database."""
    if not mol_codes:
        mol_codes = [x for y in session.query(ComponentTable.code).all() for x in y]
    donor_dict = {x: get_hbond_donors(x) for x in mol_codes}
    session.commit()
    return donor_dict

def chbond_donor_dict(mol_codes=[],session=component_session):
    """ Return CH-bond donors for all components in the database """
    if not mol_codes:
        mol_codes = [x for y in session.query(ComponentTable.code).all() for x in y]
    donor_dict = {x: get_chbond_donors(x) for x in mol_codes}
    session.commit()
    return donor_dict

def add_hbond_acceptors(mol_code, session=component_session):
    """ Populate the hbond_acceptor table for a given component.

        Note
        ----
        This must be run for new components to allow H-bonds to be found for that component.

        Parameters
        ----------
        mol_code : str
            Three-letter code for a component.
        session : sqlalchemy.orm.session.Session
            SQLAlchemy session connecting to components database.

        Returns
        -------
        hbond_acceptor_entries : list
            List of H-bond acceptor entries for component.
        """
    acceptor_entries = []
    comp_entry = session.query(ComponentTable).filter(ComponentTable.code == mol_code).one_or_none()
    if not comp_entry:
        (comp_entry, new_component), all_atoms, all_bonds = add_component(mol_code, session=session)
        session.commit()
        if not comp_entry:
            print("{0} not a known component.".format(mol_code))
            return acceptor_entries
    heteroatoms = [x for x in comp_entry.atoms if x.element.symbol in hbond_heteroatoms]
    acceptor_dict = {}
    for atom in heteroatoms:
        bound_atoms = [x.atom_2.atom_label for x in atom.covalent_bonds_1] +\
                      [x.atom_1.atom_label for x in atom.covalent_bonds_2]
        acceptor_dict[atom.atom_label] = bound_atoms
    for aa in acceptor_dict:
        aa_entry = session.query(AtomTable).filter(AtomTable.atom_label == aa).\
                   filter(AtomTable.component_id == comp_entry.id).one()
        if acceptor_dict[aa]:
            for ba in acceptor_dict[aa]:
                ba_entry = session.query(AtomTable).filter(AtomTable.atom_label == ba).\
                           filter(AtomTable.component_id == comp_entry.id).one()
                hbond_acceptor_dict = {'component_id': comp_entry.id,
                                       'heteroatom_id': aa_entry.id,
                                       'bound_atom_id': ba_entry.id}
                hbond_acceptor_entry, new_hbond_acceptor = get_or_create(HBondAcceptorTable,
                                                                         session=session, **hbond_acceptor_dict)
                acceptor_entries.append(hbond_acceptor_entry)
        else:
            hbond_acceptor_dict = {'component_id': comp_entry.id,
                                   'heteroatom_id': aa_entry.id,
                                   'bound_atom_id': None}
            hbond_acceptor_entry, new_hbond_acceptor = get_or_create(HBondAcceptorTable, session=session,
                                                                     **hbond_acceptor_dict)
            acceptor_entries.append(hbond_acceptor_entry)
    comp_entry.hbond_acceptors_added = True
    return acceptor_entries


def get_hbond_acceptors(mol_code, session=component_session):
    """ Retrieve or populate and return the H-bond acceptors for a component."""
    comp_entry = session.query(ComponentTable).filter(ComponentTable.code == mol_code).one_or_none()
    if not comp_entry:
        (comp_entry, new_component), all_atoms, all_bonds = add_component(mol_code, session=session)
        session.commit()
        if not comp_entry:
            print("{0} not a known component.".format(mol_code))
            return {}
    if not bool(comp_entry.hbond_acceptors_added):
        print("Finding H-bond acceptors for new component {0}.".format(mol_code))
        acceptors = add_hbond_acceptors(mol_code, session=session)
    else:
        acceptors = comp_entry.hbond_acceptors
    acceptor_dict = {}
    for acceptor in acceptors:
        if not acceptor.bound_atom:
            acceptor_dict[acceptor.heteroatom.atom_label] = []
        elif acceptor.heteroatom.atom_label in acceptor_dict:
            acceptor_dict[acceptor.heteroatom.atom_label].append(acceptor.bound_atom.atom_label)
        else:
            acceptor_dict[acceptor.heteroatom.atom_label] = [acceptor.bound_atom.atom_label]
    return acceptor_dict


def hbond_acceptor_dict(mol_codes=None, session=component_session):
    """ Return the H-bond acceptors for all components in the database."""
    if not mol_codes:
        mol_codes = [x for y in session.query(ComponentTable.code).all() for x in y]
    acceptor_dict = {x: get_hbond_acceptors(x) for x in mol_codes}
    session.commit()
    return acceptor_dict


def get_hbond_dicts(mol_codes=None, session=component_session):
    """ Return the dictionaries of H-bond donor and acceptor atoms for a list of mol codes."""
    donor_dict = hbond_donor_dict(mol_codes=mol_codes, session=session)
    acceptor_dict = hbond_acceptor_dict(mol_codes=mol_codes, session=session)
    return donor_dict, acceptor_dict

def get_chbond_dict(mol_codes=None, session=component_session):
    donor_dict = chbond_donor_dict(mol_codes=mol_codes,session=session)
    return donor_dict


def atomic_mass(symbol, session=component_session):
    """ Returns the atomic mass of an element."""
    element_entry = session.query(ElementTable.atomic_mass).filter(ElementTable.symbol == symbol).one_or_none()
    if not element_entry:
        print("{0} is not a known element.".format(symbol))
        return None
    return element_entry.atomic_mass


__author__ = 'Kieran L. Hudson'
