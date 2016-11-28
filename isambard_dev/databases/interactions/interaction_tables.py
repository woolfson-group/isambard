from pathlib import Path
from sqlalchemy import create_engine, Column, Integer, String, Float, ForeignKey, Table
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import backref, sessionmaker, relationship

from settings import global_settings


try:
    db_dir = Path(global_settings['interactions_database']['folder'])
    print("Interactions database directory (this can be changed in the settings.json under 'interactions_db_dir'):",
          str(db_dir), "\nExisting databases:", ', '.join([x.stem for x in db_dir.glob('*.db')]))
    db_name = input("Name of database to use or create: ") + '.db'
    db_path = Path(db_dir, db_name)
    engine_name = 'sqlite:///' + str(db_path)
    interaction_engine = create_engine(engine_name, echo=False)
    InteractionBase = declarative_base()
    InteractionSession = sessionmaker(bind=interaction_engine)
    interaction_session = InteractionSession()
finally:
    interaction_session.close()


class PdbTable(InteractionBase):
    __tablename__ = 'pdb'

    id = Column(Integer, primary_key=True, unique=True)
    pdb = Column(String(4), nullable=False, unique=True)
    # additional column covered by other tables in Coeus
    resolution = Column(Float)

    def __repr__(self):
        return '<Pdb(pdb={0})>'.format(self.pdb)


class CAZyTable(InteractionBase):
    __tablename__ = 'cazy'

    pdb_id = Column(Integer, ForeignKey('pdb.id'), primary_key=True, nullable=False)
    family = Column(String(10), nullable=False)
    kingdom = Column(String(255))
    organism = Column(String(255))
    trivial_name = Column(String(255))
    gene_name = Column(String(255))
    locus_tag_name = Column(String(255))
    standard_name = Column(String(255))

    pdb = relationship('PdbTable', backref=backref("cazy", uselist=False))

    def __repr__(self):
        return '<CAZy(pdb={0}, protein={1}, family={2}>'.format(self.pdb.pdb,  self.trivial_name, self.family)


# similar to ChainDB in coeus, but skips PDBe table
class ChainTable(InteractionBase):
    __tablename__ = 'chain'

    id = Column(Integer, primary_key=True, unique=True)
    pdb_id = Column(Integer, ForeignKey('pdb.id'), nullable=False)
    chain = Column(String(4), nullable=False)

    pdb = relationship('PdbTable', backref='chains')

    def __repr__(self):
        return '<Chain(pdb={0}, chain={1})>'.format(self.pdb.pdb, self.chain)


class ComponentTable(InteractionBase):
    __tablename__ = 'component'

    id = Column(Integer, primary_key=True, unique=True)
    code = Column(String(3), nullable=False, unique=True)
    category = Column(String(255))

    def __repr__(self):
        return '<Component(code={0})>'.format(self.code)


environment_table = Table('environment', InteractionBase.metadata,
                          Column('left_id', Integer, ForeignKey('monomer.id')),
                          Column('right_id', Integer, ForeignKey('monomer.id')))


class MonomerTable(InteractionBase):
    __tablename__ = 'monomer'

    id = Column(Integer, primary_key=True, unique=True)
    chain_id = Column(Integer, ForeignKey('chain.id'), nullable=False)
    monomer = Column(Integer, nullable=False)
    component_id = Column(Integer, ForeignKey('component.id'), nullable=False)
    insertion_code = Column(String(1))

    chain = relationship('ChainTable', backref='monomers')
    component = relationship('ComponentTable', backref='monomers')
    environment = relationship('MonomerTable', secondary=environment_table,
                               primaryjoin=environment_table.c.left_id == id,
                               secondaryjoin=environment_table.c.right_id == id)

    def __repr__(self):
        return '<Monomer(code={0}, id={1}, chain={2}, pdb={3})>'.format(self.component.code, self.monomer,
                                                                        self.chain.chain, self.chain.pdb.pdb)


class AminoAcidTable(InteractionBase):
    __tablename__ = 'amino_acid'

    monomer_id = Column(Integer, ForeignKey('monomer.id'), primary_key=True, unique=True)
    phi = Column(Float)
    psi = Column(Float)
    omega = Column(Float)
    secondary_structure = Column(String(1))
    chi_1 = Column(Float)
    chi_2 = Column(Float)
    chi_3 = Column(Float)
    chi_4 = Column(Float)

    monomer = relationship('MonomerTable', backref=backref("amino_acid", uselist=False))

    def __repr__(self):
        return '<Amino acid(pdb={0}, chain={1}, residue={2}>)'.format(self.monomer.chain.pdb.pdb,
            self.monomer.chain.chain, self.monomer.monomer)


class LigandTable(InteractionBase):
    __tablename__ = 'ligand'

    monomer_id = Column(Integer, ForeignKey('monomer.id'), primary_key=True, unique=True)

    monomer = relationship('MonomerTable', backref=backref("ligand", uselist=False))

    def __repr__(self):
        return '<Ligand(pdb={0}, chain={1}, residue={2}>)'.format(self.monomer.chain.pdb.pdb, self.monomer.chain.chain,
                                                                  self.monomer.monomer)


class PrivateerTable(InteractionBase):
    __tablename__ = 'privateer'

    monomer_id = Column(Integer, ForeignKey('monomer.id'), primary_key=True, unique=True)
    q = Column(Float)
    phi = Column(Float)
    theta = Column(Float)
    rscc = Column(Float)
    anomer = Column(String(5))
    configuration = Column(String(1))
    cho_type = Column(String(6))
    ring = Column(String(10))
    conformation = Column(String(3))
    mean_density = Column(Float)
    mean_b_factor = Column(Float)
    bond_length_deviation = Column(Float)
    bond_angle_deviation = Column(Float)
    environment = Column(String(1))
    diagnostic = Column(String(10))

    monomer = relationship('MonomerTable',  backref=backref("privateer", uselist=False))

    def __repr__(self):
        return"<Privateer data(pdb={0}, residue={1}-{2}-{3})>".\
            format(self.monomer.chain.pdb.pdb, self.monomer.component.code,
                   self.monomer.monomer, self.monomer.chain.chain)


class InteractionTable(InteractionBase):
    __tablename__ = 'interaction'

    id = Column(Integer, primary_key=True, unique=True)
    donor_id = Column(Integer, ForeignKey('monomer.id'), nullable=False)
    acceptor_id = Column(Integer, ForeignKey('monomer.id'), nullable=False)

    donor = relationship('MonomerTable', foreign_keys=[donor_id], backref='donor_interactions')
    acceptor = relationship('MonomerTable', foreign_keys=[acceptor_id], backref='acceptor_interactions')

    def __repr__(self):
        return '<Interaction(donor={0}, acceptor={1})>'.format(self.donor, self.acceptor)


class CHpiTable(InteractionBase):
    __tablename__ = 'CH-pi'

    id = Column(Integer, primary_key=True, unique=True)
    interaction_id = Column(Integer, ForeignKey('interaction.id'), nullable=False)
    donor_C = Column(String(4), nullable=False)
    donor_H = Column(String(4), nullable=False)
    pi_system = Column(String(10))
    distance = Column(Float)
    proj_dist = Column(Float)
    angle = Column(Float)

    interaction = relationship('InteractionTable', backref='CH_pis')

    def __repr__(self):
        return '<CH-pi interaction(donor={0} and {1} of {2}, acceptor={3} of {4})>'.\
               format(self.donor_C, self.donor_H, self.interaction.donor, self.pi_system, self.interaction.acceptor)

class HbondTable(InteractionBase):
    __tablename__ = 'H-bond'

    id = Column(Integer, primary_key=True, unique=True)
    interaction_id = Column(Integer, ForeignKey('interaction.id'), nullable=False)
    donor_atom = Column(String(5), nullable=False)
    acceptor_atom = Column(String(5), nullable=False)
    distance = Column(Float)
    angle_a = Column(Float)
    angle_d = Column(Float)

    interaction = relationship('InteractionTable', backref='H_bonds')

    def __repr__(self):
        return '<H-bond interaction(donor={0} of {1}, acceptor={2} of {3)>'.format(
            self.donor_atom, self.interaction.donor, self.acceptor_atom, self.interaction.acceptor)

__author__ = 'Kieran L. Hudson, Ali Scott'
