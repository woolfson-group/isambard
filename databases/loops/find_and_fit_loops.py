import copy
import os
import ast

import numpy
from numpy import array
from sqlalchemy import create_engine, Column, Float, Integer, String
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker

from settings import global_settings
from tools.geometry import unit_vector, rotation_matrix, distance, angle_between_vectors, dihedral

#TODO change eval methods to ast, alter structure to handle new data structures

# Start loop_db
try:
    engine = create_engine('sqlite:///' + os.path.join(
            global_settings['package_path'], 'databases', 'loops', 'loops.db'), echo=False)
    Base = declarative_base()
    Session = sessionmaker(bind=engine)
    loop_db_session = Session()
finally:
    loop_db_session.close()


# Loop database entry class
class Loops(Base):
    __tablename__ = 'loops'

    id = Column(Integer, primary_key=True)
    pdb = Column(String)
    loop_type = Column(String)
    start_res = Column(Integer)
    end_res = Column(Integer)
    chain = Column(String)
    sequence = Column(String)
    length = Column(Integer)
    n_n_distance = Column(Float)
    n_ca_distance = Column(Float)
    n_c_distance = Column(Float)
    ca_n_distance = Column(Float)
    ca_ca_distance = Column(Float)
    ca_c_distance = Column(Float)
    c_c_distance = Column(Float)
    c_ca_distance = Column(Float)
    c_n_distance = Column(Float)
    dihedral = Column(Float)
    bb_coords = Column(String)
    entering_res = Column(String)
    exiting_res = Column(String)
    resolution = Column(Float)

    def __repr__(self):
        return "<Loops(pdb='{0}', loop_type='{1}', sequence='{2})>".format(
            self.pdb, self.loop_type, self.sequence)


# Generate loop database query
def query_loop_db(n_n_distance, n_ca_distance, n_c_distance,
                  ca_n_distance, ca_ca_distance, ca_c_distance,
                  c_n_distance, c_ca_distance, c_c_distance,
                  min_resolution, distance_threshold, ss_before, ss_after, min_length, max_length):
    """Queries loop database for loops that match input parameters.

    Parameters
    ----------
    ca_ca_distance : float
        Distance between CA atoms of the entering and exiting residues of secondary structure before and after the loop.
    c_n_distance : float
        Distance between the C and the N atoms of the entering and exiting residues of secondary structure respectively,
        before and after the loop.
    min_resolution : float
        Minimum resolution of crystal structures included in query, measured in angstroms.
    distance_threshold : float
        Tolerance of deviation away from stated ca_ca_distance and c_n_distance. Default value = 0.2
    ss_before : str
        String containing secondary structure type leading into loop, 'a', '3', 'p' and 'b' are permitted values.
    ss_after : str
        String containing secondary structure type immediately after the loop, 'a', '3', 'p' and 'b' are permitted
        values.

    Returns
    -------
    loop_matches : [Loops]
        Returns a list of Loops database objects.
    """
    loop_matches = loop_db_session.query(Loops).filter(
        Loops.n_n_distance > n_n_distance - distance_threshold,
        Loops.n_n_distance < n_n_distance + distance_threshold,
        # Loops.n_ca_distance > n_ca_distance - distance_threshold,
        # Loops.n_ca_distance < n_ca_distance + distance_threshold,
        # Loops.n_c_distance > n_c_distance - distance_threshold,
        # Loops.n_c_distance < n_c_distance + distance_threshold,

        Loops.ca_n_distance > ca_n_distance - distance_threshold,
        Loops.ca_n_distance < ca_n_distance + distance_threshold,
        Loops.ca_ca_distance > ca_ca_distance - distance_threshold,
        Loops.ca_ca_distance < ca_ca_distance + distance_threshold,
        Loops.ca_c_distance > ca_c_distance - distance_threshold,
        Loops.ca_c_distance < ca_c_distance + distance_threshold,

        # Loops.c_n_distance > c_n_distance - distance_threshold,
        # Loops.c_n_distance < c_n_distance + distance_threshold,
        # Loops.c_ca_distance > c_ca_distance - distance_threshold,
        # Loops.c_ca_distance < c_ca_distance + distance_threshold,
        Loops.c_c_distance > c_c_distance - distance_threshold,
        Loops.c_c_distance < c_c_distance + distance_threshold,

        Loops.loop_type.like('{0}l{1}'.format(ss_before, ss_after)),
        Loops.resolution < min_resolution,
        Loops.length >= min_length,
        Loops.length <= max_length).all()
    return sorted(loop_matches, key=lambda x: x.length)


def find_loops(ent_res, ext_res, min_resolution=3, distance_threshold=0.2, ss_before='%', ss_after='%',
               min_length=0, max_length=50):
    """Finds loops that will fit between the two residues provided.

    Parameters
    ----------
    ent_res : base.ampal.Residue
        Last residue of secondary structure entering the loop.
    ext_res : base.ampal.Residue
        First residue of secondary structure exiting the loop.
    min_resolution : float (optional)
        Lowest resolution allowed in loop filtering i.e. 3.0 will give all loops with resolution of 3.0 or below,
        measured in angstroms.
    distance_threshold : float (optional)
        Tolerance of deviation away from e_x_dihedral. Default value = 10.0
    angle_threshold : float (optional)
        Tolerance of deviation away from e_x_angle. Default value = 10.0
    dihedral_threshold : float (optional)
        Tolerance of deviation away from e_x_dihedral. Default value = 10.0
    ss_before : str (optional)
        String containing secondary structure type leading into loop, 'a', '3', 'p' and 'b' are permitted values.
    ss_after : str (optional)
        String containing secondary structure type immediately after the loop, 'a', '3', 'p' and 'b' are permitted
        values.

    Returns
    -------
    potential_loops : [Loops]
        A list of database entries for loops that match the specified criteria.
    fit_coords : (numpy.array, numpy.array, numpy.array, numpy.array)
        Tuple containing the coordinates for the atomic positions used for the fitting, these are:
            [0] Entering CA
            [1] Entering C
            [2] Exiting N
            [3] Exiting CA
    """

    match_n_n_distance = distance(ent_res.atoms['N'], ext_res.atoms['N'])
    match_n_ca_distance = distance(ent_res.atoms['N'], ext_res.atoms['CA'])
    match_n_c_distance = distance(ent_res.atoms['N'], ext_res.atoms['C'])

    match_ca_n_distance = distance(ent_res.atoms['CA'], ext_res.atoms['N'])
    match_ca_ca_distance = distance(ent_res.atoms['CA'], ext_res.atoms['CA'])
    match_ca_c_distance = distance(ent_res.atoms['CA'], ext_res.atoms['C'])

    match_c_n_distance = distance(ent_res.atoms['C'], ext_res.atoms['N'])
    match_c_ca_distance = distance(ent_res.atoms['C'], ext_res.atoms['CA'])
    match_c_c_distance = distance(ent_res.atoms['C'], ext_res.atoms['C'])

    potential_loops = query_loop_db(match_n_n_distance, match_n_ca_distance, match_n_c_distance,
                                    match_ca_n_distance, match_ca_ca_distance, match_ca_c_distance,
                                    match_c_n_distance, match_c_ca_distance, match_c_c_distance,
                                    min_resolution,
                                    distance_threshold=distance_threshold,
                                    ss_before=ss_before, ss_after=ss_after,
                                    min_length=min_length, max_length=max_length)
    return potential_loops


def loop_db_entry_to_ampal(loop_db_entry):
    """Converts a database object to a ampal Chain object.

    Parameters
    ----------
    loop_db_entry : Loops
        Raw database entry object.

    Returns
    -------
    loop_chain : base_ampal.Chain
        Chain object generated from a database object, extra info on the loop source is in tags.['loop_data'].
    """
    from ampal.protein import flat_list_to_polymer

    l_ent = eval(loop_db_entry.entering_res)
    l_ext = eval(loop_db_entry.exiting_res)
    l_bb = eval(loop_db_entry.bb_coords)
    full_loop_coords = l_ent + l_bb + l_ext
    loop_chain = flat_list_to_polymer(full_loop_coords)

    loop_chain.tags['loop_data'] = {}
    loop_chain.tags['loop_data']['pdb'] = loop_db_entry.pdb
    loop_chain.tags['loop_data']['loop_type'] = loop_db_entry.loop_type
    loop_chain.tags['loop_data']['start_res'] = loop_db_entry.start_res
    loop_chain.tags['loop_data']['end_res'] = loop_db_entry.end_res
    loop_chain.tags['loop_data']['chain'] = loop_db_entry.chain
    loop_chain.tags['loop_data']['sequence'] = loop_db_entry.sequence
    loop_chain.tags['loop_data']['resolution'] = loop_db_entry.resolution
    return loop_chain


def rmsd_of_loop_fit(ent_res, ext_res, loop):
    """Returns the RMSD of the input residues and the first and last residue of the loop."""
    deltas = []
    deltas.append(distance(ent_res.atoms['N'], loop[0].atoms['N']))
    deltas.append(distance(ent_res.atoms['CA'], loop[0].atoms['CA']))
    deltas.append(distance(ent_res.atoms['C'], loop[0].atoms['C']))
    deltas.append(distance(ext_res.atoms['N'], loop[-1].atoms['N']))
    deltas.append(distance(ext_res.atoms['CA'], loop[-1].atoms['CA']))
    deltas.append(distance(ext_res.atoms['C'], loop[-1].atoms['C']))
    loop_rms = sum(deltas)/float(len(deltas))
    loop.tags['fit_rmsd'] = loop_rms
    return loop_rms


def fit_loop(ent_res, ext_res, in_loop):
    """Takes a chain object containing information of a loop and fits

    Parameters
    ----------
    ent_res : ampal.base_ampal.Residue
        Residue that will be entering the loop.
    ext_res : ampal.base_ampal.Residue
        Residue that will be exiting the loop.
    in_loop : ampal.base_ampal.Chain
        Loop to fit between the entering and exiting residues.

    Returns
    -------
    loop : ampal.base_ampal.Chain
        Post fit loop.
    """
    loop = copy.deepcopy(in_loop)
    ss_e_centre = (ent_res.atoms['N'].array + ent_res.atoms['CA'].array)/2
    ss_x_centre = (ext_res.atoms['CA'].array + ext_res.atoms['C'].array)/2
    ss_centre = (ss_e_centre + ss_x_centre)/2
    l_e_centre = (loop[0]['N'].array + loop[0]['CA'].array)/2
    l_x_centre = (loop[-1]['CA'].array + loop[-1]['C'].array)/2
    l_centre = (l_e_centre + l_x_centre)/2
    # Translate loop to centre of secondary structure elements
    trans_vector = l_centre - ss_centre
    for atom in loop.get_atoms():
        atom._vector -= trans_vector
    # Align C/N vectors
    l_e_centre = (loop[0]['N'].array + loop[0]['CA'].array)/2
    rot_angle = numpy.deg2rad(angle_between_vectors(ss_e_centre - ss_centre, l_e_centre - ss_centre))
    rot_vector = unit_vector(numpy.cross(l_e_centre - ss_centre, ss_e_centre - ss_centre))
    align_rm = rotation_matrix(rot_angle, rot_vector, ss_centre)
    for atom in loop.get_atoms():
        atom._vector = numpy.dot(align_rm, numpy.transpose(numpy.append(atom.array, 1)))[:3]
    # Rotate around ss_e and ss_x vector, first aligning the e_ns and then half the x_ns
    l_e_centre = (loop[0]['N'].array + loop[0]['CA'].array)/2
    if distance(ss_centre, ss_e_centre) > distance(ss_centre, l_e_centre):
        dihedral_angle1 = -dihedral(ent_res.atoms['N'].array, ss_e_centre,
                                    l_e_centre, loop[0].atoms['N'].array, radians=True)
    else:
        dihedral_angle1 = dihedral(ent_res.atoms['N'].array, ss_e_centre,
                                   l_e_centre, loop[0].atoms['N'].array, radians=True)
    dihedral_vector = unit_vector(ss_x_centre - ss_e_centre)
    dihedral_rm = rotation_matrix(dihedral_angle1, dihedral_vector, ss_centre)
    for atom in loop.get_atoms():
        atom._vector = numpy.dot(dihedral_rm, numpy.transpose(numpy.append(atom.array, 1)))[:3]
    l_x_centre = (loop[-1]['CA'].array + loop[-1]['C'].array)/2
    if distance(ss_centre, ss_x_centre) > distance(ss_centre, l_x_centre):
        dihedral_angle2 = -dihedral(ext_res.atoms['C'].array, ss_x_centre,
                                    l_x_centre, loop[-1].atoms['C'].array, radians=True)
    else:
        dihedral_angle2 = dihedral(ext_res.atoms['C'].array, ss_x_centre,
                                   l_x_centre, loop[-1].atoms['C'].array, radians=True)
    dihedral_vector = unit_vector(ss_x_centre - ss_e_centre)
    dihedral_rm = rotation_matrix(-dihedral_angle2/2, dihedral_vector, ss_centre)
    for atom in loop.get_atoms():
        atom._vector = numpy.dot(dihedral_rm, numpy.transpose(numpy.append(atom.array, 1)))[:3]
    rmsd_of_loop_fit(ent_res, ext_res, loop)
    return loop


__author__ = 'Christopher W. Wood'
