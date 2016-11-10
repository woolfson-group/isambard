# coding: utf-8
from sqlalchemy import Column, Date, Float, ForeignKey, UniqueConstraint,\
    Integer, Numeric, SmallInteger, Boolean, String, create_engine
from sqlalchemy.dialects.mysql import TINYINT
from sqlalchemy.orm import relationship, sessionmaker
from sqlalchemy.ext.declarative import declarative_base

from settings import global_settings

try:
    db_name = 'coeus'
    _user = global_settings['coeus']['user']
    _host = global_settings['coeus']['host']
    _pwd = global_settings['coeus']['pwd']
    engine_name = 'mysql+pymysql://{0}:{1}@{2}/{3}'.format(_user, _pwd, _host, db_name)
    coeus_engine = create_engine(engine_name)
    CoeusBase = declarative_base()
    coeus_metadata = CoeusBase.metadata
    CoeusSession = sessionmaker(bind=coeus_engine)
    coeus_session = CoeusSession()
finally:
    coeus_session.close()


class AminoAcidDB(CoeusBase):
    __tablename__ = 'amino_acid'
    __table_args__ = (UniqueConstraint('code', 'letter'), {'mysql_engine': 'InnoDB'})

    id = Column(Integer, primary_key=True)
    letter = Column(String(1), nullable=False)
    code = Column(String(5), nullable=False)
    modified = Column(String(255))
    description = Column(String(255))

    residues = relationship('ResidueDB', back_populates='amino_acid')

    def __repr__(self):
        return '<AminoAcidDB(code={0}, letter={1})>'.format(self.code, self.letter)


class AtlasDB(CoeusBase):
    __tablename__ = 'atlas'
    __table_args__ = {'mysql_engine': 'InnoDB'}

    id = Column(Integer, primary_key=True)
    nodes = Column(SmallInteger, nullable=False)
    edges = Column(SmallInteger, nullable=False)
    name = Column(String(30), nullable=False, unique=True)
    two_core_id = Column(Integer, ForeignKey('atlas.id', ondelete='CASCADE'))

    graphs = relationship('GraphDB', back_populates='atlas')
    two_core = relationship('AtlasDB', cascade='all, delete-orphan', passive_deletes=True)

    def __repr__(self):
        return '<AtlasDB(name={0}, nodes={1}, edges={2})>'.format(self.name, self.nodes, self.edges)


class CdhitFullDB(CoeusBase):
    __tablename__ = 'cdhit_full'
    __table_args__ = {'mysql_engine': 'InnoDB'}

    id = Column(Integer, primary_key=True)
    c90 = Column(Boolean, default=False)
    c80 = Column(Boolean, default=False)
    c70 = Column(Boolean, default=False)
    c60 = Column(Boolean, default=False)
    pdb_id = Column(ForeignKey('pdb.id', ondelete='CASCADE'), nullable=False, index=True)

    pdb = relationship('PdbDB', back_populates='cdhit_full')

    def __repr__(self):
        return '<CdhitFullDB(pdb={0}, c90={1})>'.format(self.pdb.pdb, self.c90)


class ChainDB(CoeusBase):
    __tablename__ = 'chain'
    __table_args__ = (UniqueConstraint('pdbe_id', 'chain'), {'mysql_engine': 'InnoDB'})

    id = Column(Integer, primary_key=True)
    chain = Column(String(4), nullable=False)
    pdbe_id = Column(ForeignKey('pdbe.id', ondelete='CASCADE'), nullable=False, index=True)

    pdbe = relationship('PdbeDB', back_populates='chains')
    helices = relationship('HelixDB', back_populates='chain', cascade='all, delete-orphan', passive_deletes=True)
    residues = relationship('ResidueDB', back_populates='chain', cascade='all, delete-orphan', passive_deletes=True)

    def __repr__(self):
        return '<ChainDB(pdb={0}, chain={1})>'.format(self.pdbe.pdb.pdb, self.chain)


class CutoffDB(CoeusBase):
    __tablename__ = 'cutoff'
    __table_args__ = (UniqueConstraint('scut', 'kcut'), {'mysql_engine': 'InnoDB'})

    id = Column(Integer, primary_key=True)
    scut = Column(Numeric(4, 2), nullable=False)
    kcut = Column(Integer, nullable=False)

    graphs = relationship('GraphDB', back_populates='cutoff')

    def __repr__(self):
        return '<CutoffDB(scut={0}, kcut={1})>'.format(self.scut, self.kcut)


class GraphDB(CoeusBase):
    __tablename__ = 'graph'
    __table_args__ = {'mysql_engine': 'InnoDB'}

    id = Column(Integer, primary_key=True)
    connected_component = Column(SmallInteger, nullable=False, index=True)
    atlas_id = Column(ForeignKey('atlas.id'), index=True, nullable=False)
    cutoff_id = Column(ForeignKey('cutoff.id'), nullable=False, index=True)
    pdbe_id = Column(ForeignKey('pdbe.id', ondelete='CASCADE'), nullable=False, index=True)

    atlas = relationship('AtlasDB', back_populates='graphs')
    cutoff = relationship('CutoffDB', back_populates='graphs')
    pdbe = relationship('PdbeDB', back_populates='graphs')
    helices = relationship('GraphHelixDB', back_populates='graph', cascade='all, delete-orphan', passive_deletes=True)

    def __repr__(self):
        return '<GraphDB(pdb={0}, name={1})>'.format(self.pdbe.pdb.pdb, self.atlas.name)


# TODO add __repr__ and check deletion behaviour.
class GraphHelixDB(CoeusBase):
    __tablename__ = 'graph_helix'
    __table_args__ = {'mysql_engine': 'InnoDB'}

    helix_id = Column(ForeignKey('helix.id', ondelete='CASCADE'), primary_key=True)
    graph_id = Column(ForeignKey('graph.id', ondelete='CASCADE'), primary_key=True)
    degree = Column(Integer, nullable=False)

    graph = relationship('GraphDB', back_populates='helices')
    helix = relationship('HelixDB', back_populates='graphs')

    def __repr__(self):
        return '<GraphHelixDB(pdb={0}, name={1})>'.format(self.graph.pdbe.pdb.pdb, self.graph.atlas.name)


class HelixDB(CoeusBase):
    __tablename__ = 'helix'
    __table_args__ = (UniqueConstraint('chain_id', 'number'), {'mysql_engine': 'InnoDB'})

    id = Column(Integer, primary_key=True)
    number = Column(SmallInteger, nullable=False)
    length = Column(SmallInteger, nullable=False)
    chain_id = Column(ForeignKey('chain.id', ondelete='CASCADE'), nullable=False, index=True)

    chain = relationship('ChainDB', back_populates='helices')
    residues = relationship('ResidueDB', back_populates='helix', cascade='all, delete-orphan', passive_deletes=True)
    graphs = relationship('GraphHelixDB', back_populates='helix', cascade='all, delete-orphan', passive_deletes=True)

    def __repr__(self):
        return '<HelixDB(pdb={0}, number={1})>'.format(self.chain.pdbe.pdb.pdb, self.number)


class KnobDB(CoeusBase):
    __tablename__ = 'knob'
    __table_args__ = {'mysql_engine': 'InnoDB'}

    id = Column(Integer, primary_key=True)
    residue_id = Column(ForeignKey('residue.id', ondelete='CASCADE'), nullable=False)
    packing_angle = Column(Numeric(6, 3))
    knob_type = Column(SmallInteger)
    max_cv_dist = Column(Numeric(5, 2), nullable=False)

    holes = relationship('HoleDB', back_populates='knob')
    residue = relationship('ResidueDB', back_populates='knob')

    def __repr__(self):
        return '<KnobDB(max_cv_dist={0})>'.format(self.max_cv_dist)


# TODO add __repr__ and check deletion behaviour.
class HoleDB(CoeusBase):
    # Association table between KnobDB and ResidueDB.
    __tablename__ = 'hole'
    __table_args__ = {'mysql_engine': 'InnoDB'}

    knob_id = Column(ForeignKey('knob.id', ondelete='CASCADE'), primary_key=True)
    residue_id = Column(ForeignKey('residue.id', ondelete='CASCADE'), primary_key=True)
    # Hole residues are numbered from 0-3. (see socket reference).
    hole_res_type = Column(TINYINT, nullable=False, index=True)

    knob = relationship('KnobDB', back_populates='holes')
    residue = relationship('ResidueDB', back_populates='holes')


class PdbDB(CoeusBase):
    __tablename__ = 'pdb'
    __table_args__ = {'mysql_engine': 'InnoDB'}

    id = Column(Integer, primary_key=True)
    pdb = Column(String(4), nullable=False, unique=True)

    protein = relationship('ProteinDB', back_populates='pdb', cascade='all, delete-orphan', passive_deletes=True)
    cdhit_full = relationship('CdhitFullDB', back_populates='pdb', cascade='all, delete-orphan', passive_deletes=True)
    pdbes = relationship('PdbeDB', back_populates='pdb', cascade='all, delete-orphan', passive_deletes=True)

    def __repr__(self):
        return '<PdbDB(pdb={0})>'.format(self.pdb)


class PdbeDB(CoeusBase):
    __tablename__ = 'pdbe'
    __table_args__ = (UniqueConstraint('pdb_id', 'mmol'), {'mysql_engine': 'InnoDB'})

    id = Column(Integer, primary_key=True)
    # specific mysql data type, maximum value = 255.
    mmol = Column(TINYINT(unsigned=True), nullable=False)
    preferred = Column(Boolean, nullable=False, index=True)
    pdb_id = Column(ForeignKey('pdb.id', ondelete='CASCADE'), nullable=False, index=True)

    pdb = relationship('PdbDB', back_populates='pdbes')
    chains = relationship('ChainDB', back_populates='pdbe', cascade='all, delete-orphan', passive_deletes=True)
    graphs = relationship('GraphDB', back_populates='pdbe', cascade='all, delete-orphan', passive_deletes=True)

    def __repr__(self):
        return '<PdbeDB(pdb={0}, mmol={1}, preferred={2})>'.format(self.pdb.pdb, self.mmol, bool(self.preferred))


# TODO Should this be linked to PdbeDB instead of PdbDB? Depends what data we're storing.
class ProteinDB(CoeusBase):
    __tablename__ = 'protein'
    __table_args__ = {'mysql_engine': 'InnoDB'}

    id = Column(Integer, primary_key=True)
    keywords = Column(String(255))
    header = Column(String(100))
    space_group = Column(String(15))
    experimental_method = Column(String(100))
    crystal_growth = Column(String(255))
    resolution = Column(Numeric(4, 2))
    r_value_obs = Column(Numeric(6, 4))
    atoms_protein = Column(Integer)
    atoms_solvent = Column(Integer)
    atoms_ligand = Column(Integer)
    atoms_nucleic_acid = Column(Integer)
    atoms_total = Column(Integer)
    title = Column(String(255))
    pdb_descriptor = Column(String(255))
    model_details = Column(String(255))
    model_type_details = Column(String(255))
    casp_flag = Column(String(255))
    deposition_date = Column(Date)
    release_date = Column(Date)
    last_modified_date = Column(Date)
    ncbi_taxonomy = Column(Integer)
    ncbi_taxonomy_host_org = Column(Integer)
    pdb_id = Column(ForeignKey('pdb.id', ondelete='CASCADE'), nullable=False, index=True, unique=True)

    pdb = relationship('PdbDB', back_populates='protein')

    def __repr__(self):
        return '<ProteinDB(pdb={0}, title={1})>'.format(self.pdb.pdb, self.title)


# TODO Calculate approximate data size of each row, and therefore the PDB. Is it worth making floats smaller numerics?
# http://dev.mysql.com/doc/refman/5.7/en/storage-requirements.html
class ResidueDB(CoeusBase):
    __tablename__ = 'residue'
    __table_args__ = (UniqueConstraint('chain_id', 'resno', 'icode_int', 'amino_acid_id',
                                       name='uix_chain_id_resno_icode_int_amino_acid_id'),
                      {'mysql_engine': 'InnoDB'})

    id = Column(Integer, primary_key=True)
    acc = Column(Numeric(4, 1), nullable=True)
    phi = Column(Numeric(4, 1), nullable=True)
    psi = Column(Numeric(4, 1), nullable=True)
    omega = Column(Numeric(4, 1), nullable=True)
    resno = Column(Integer)
    icode_int = Column(TINYINT(unsigned=True), nullable=False, default=0)
    dssp_ss = Column(String(1))
    rpt = Column(Float(asdecimal=True))
    rpr = Column(Float(asdecimal=True))
    roc = Column(Float(asdecimal=True))
    amino_acid_id = Column(ForeignKey('amino_acid.id'), nullable=False, index=True)
    chain_id = Column(ForeignKey('chain.id', ondelete='CASCADE'), nullable=False, index=True)
    helix_id = Column(ForeignKey('helix.id', ondelete='CASCADE'), nullable=True, index=True)

    amino_acid = relationship('AminoAcidDB', back_populates='residues')
    chain = relationship('ChainDB', back_populates='residues')
    holes = relationship('HoleDB', back_populates='residue',
                        cascade='all, delete-orphan', passive_deletes=True)
    knob = relationship('KnobDB', back_populates='residue',
                        cascade='all, delete-orphan', passive_deletes=True)
    helix = relationship('HelixDB', back_populates='residues')

    def __repr__(self):
        return '<ResidueDB(pdb={0}, chain={1}, resno={2})>'.format(self.chain.pdbe.pdb.pdb, self.chain.chain, self.resno)


__author__ = 'Jack W. Heal'
__status__ = 'Development'
