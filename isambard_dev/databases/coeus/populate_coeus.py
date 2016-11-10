import itertools
from sqlalchemy.exc import IntegrityError
import numpy
import decimal
import networkx

from databases.coeus.coeus_tables import *
from databases.general_tools import populate_model
from tools.graph_theory import list_of_graphs, two_core_names, sorted_connected_components, store_graph,\
    get_graph_name, get_unknown_graph_list, add_two_core_name_to_json
from tools.amino_acids import amino_acids_dict, get_aa_info, add_amino_acid_to_json
from ampal.assembly import Assembly, AmpalContainer
from ampal.protein import Polypeptide
from add_ons.parmed_to_ampal import convert_cif_to_ampal
from add_ons.filesystem import FileSystem
from add_ons.knobs_into_holes import KnobGroup


class CoeusAmpalContainer(AmpalContainer):
    def __init__(self, code):
        super(CoeusAmpalContainer, self).__init__()
        self.code = code
        self.filesystem = FileSystem(self.code)
        for cif_file in self.filesystem.cifs.values():
            a = convert_cif_to_ampal(cif=cif_file, path=True, assembly_id=self.code)
            if isinstance(a, AmpalContainer):
                a = a[0]
            a = Assembly([x for x in a if isinstance(x, Polypeptide)])
            if len(a) > 0:
                self.append(a)
        self.tag_for_coeus()

    def tag_mmols(self):
        for i, a in enumerate(self):
            mmol_number = i + 1
            if mmol_number == self.filesystem.preferred_mmol:
                is_preferred = True
            else:
                is_preferred = False
            a.tags['mmol_number'] = mmol_number
            a.tags['preferred_mmol'] = is_preferred
        return

    def tag_knob_group(self, cutoff=10.0, force=False, tag_name='knob_group'):
        for a in self:
            tagged = [tag_name in x.tags.keys() for x in a.get_monomers(ligands=False)]
            if (not any(tagged)) or force:
                knob_group = KnobGroup.from_helices(assembly=a, cutoff=cutoff)
                if knob_group is not None:
                    a.tags[tag_name] = knob_group
        return

    def tag_for_coeus(self):
        self.tag_mmols()
        self.tag_knob_group()
        for a in self:
            a.tag_secondary_structure()
            a.tag_torsion_angles()
            # tag helix membership information
            a.tag_dssp_solvent_accessibility()
            a.tag_ca_geometry()
            for helix_number, h in enumerate(a.helices):
                for r in h:
                    r.tags['helix_number'] = helix_number
        return


def populate_amino_acid():
    """Populate AminoAcidDB model from the json file of amino acids amino_acids.json.

    Notes
    -----
    json file is loaded into memory on import of isambard.tools.amino_acids as 'aas'.
    This function does not run coeus_session.commit().

    Returns
    -------
    created_objs : list, or None
        List of Django model objects created, if any.
    """
    # Start with a clean session.
    coeus_session.rollback()
    # Get amino acid codes already in the database.
    amino_acid_codes = [x[0] for x in coeus_session.query(AminoAcidDB.code).all()]
    # Get AminoAcidDB instances for all amino acid codes not yet in the database.
    amino_acid_dbs = [AminoAcidDB(code=k, **v) for k, v in amino_acids_dict.items() if k not in amino_acid_codes]
    if not amino_acid_dbs:
        return 0
    # add and commit new instances to the session.
    coeus_session.add_all(amino_acid_dbs)
    try:
        coeus_session.commit()
    except IntegrityError:
        coeus_session.rollback()
        return 0
    return 1


def populate_cutoff():
    """ Populate CutoffDB using internally-defined range of kcuts and scuts.

    Returns
    -------
    created_objs : list, or None
        List of Django model objects created, if any.
    """
    # clear session before starting - good practice.
    coeus_session.rollback()
    # Create cutoff objects from lists of kcut and scut values.
    kcuts = list(range(4))
    scuts = [7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0]
    cutoffs = [CutoffDB(scut=scut, kcut=kcut) for kcut, scut in itertools.product(kcuts, scuts)]
    # Add cutoffs to session and commit.
    coeus_session.add_all(cutoffs)
    try:
        coeus_session.commit()
    except IntegrityError:
        coeus_session.rollback()
        return 0
    return 1


def populate_atlas():
    """ Populate AtlasDB with graphs from the extended list_of_graphs.

    Notes
    -----
    Depends on any new graphs being in the unknown_graphs shelf and in the two_core_names dictionary.
    Running isambard.tools.graph_theory.store_graph(g) for each new graph g will accomplish this.
    Runs session.commit() at the end, so any new graphs will have been added to the atlas table.
    First adds new graphs and then subsequently fills the two_core_id.
    Skips over graphs already in AtlasDB.
    Run this function after new unknown graphs are added to the unknown_graph_shelf.

    Returns
    -------
    None

    """
    # clear session before starting - good practice.
    coeus_session.rollback()
    # Full list of all graphs in Atlas and/or encountered to date.
    graphs = list_of_graphs(unknown_graphs=True)
    # Get AtlasDB instances for graphs not currently in the coeus.atlas table. Add and commit them.
    atlas_names = [x[0] for x in coeus_session.query(AtlasDB.name).all()]
    atlas_dbs = [AtlasDB(nodes=g.number_of_nodes(), name=g.name, edges=g.number_of_edges())
                 for g in graphs if g.name not in atlas_names]
    coeus_session.add_all(atlas_dbs)
    coeus_session.commit()
    # Fill in the two_core_id column for the recently added graphs.
    for atlas_db in atlas_dbs:
        try:
            two_core_name = two_core_names[atlas_db.name]
        except KeyError:
            g = next(filter(lambda x: x.name == atlas_db.name, get_unknown_graph_list()))
            two_core_name = get_graph_name(networkx.k_core(g, 2))
            add_two_core_name_to_json(atlas_db.name, two_core_name, force_add=False)
        two_core_db = coeus_session.query(AtlasDB).filter_by(name=two_core_name).one()
        atlas_db.two_core_id = two_core_db.id
    coeus_session.add_all(atlas_dbs)
    try:
        coeus_session.commit()
    except IntegrityError:
        coeus_session.rollback()
        return 0
    return 1


def simplified_add_object(instance):
    if instance not in coeus_session:
        coeus_session.add(instance)
    try:
        coeus_session.commit()
    except IntegrityError:
        print('Integrity error for instance: {0}. Rolling back...'.format(instance))
    coeus_session.rollback()
    return


def amino_acid_dbs_from_protein(protein):
    """ Dictionary of all AminoAcidDB instances from the database that are encountered in Protein.

    Parameters
    ----------
    protein : An ampal Protein object.

    Returns
    -------
    amino_acid_dbs : dict
        Keys are tuples (amino_acid_code, amino_acid_letter)
        Values are the corresponding AminoAcidDB instances.
    """
    coeus_session.rollback()
    amino_acid_dbs = {}
    # Get distinct amino acids from the Protein, and their corresponding AminoAcidDB instances.
    for code, letter in set([(residue.mol_code, residue.mol_letter)
                             for residue in protein.get_monomers(ligands=False)]):
        amino_acid_db = coeus_session.query(AminoAcidDB).filter_by(code=code, letter=letter).one_or_none()
        # If amino acid is not in the database:
        if not amino_acid_db:
            # Get information for new amino acid.
            aa_dict = get_aa_info(code=code)
            # Use this information to create AminoAcidDB instance.
            amino_acid_db = AminoAcidDB(**aa_dict)
            # Add new amino acid to json file that populates aas.
            add_amino_acid_to_json(**aa_dict)
            # Add this new amino acid to the amino_acid_table.
            populate_amino_acid()
        amino_acid_dbs[(code, letter)] = amino_acid_db
    coeus_session.rollback()
    return amino_acid_dbs


def populate_structural_tables(coeus_ampal_container):
    if not coeus_ampal_container:
        return
    if not len(set([a.id for a in coeus_ampal_container])) == 1:
        raise IOError('all Assembly objects in db_ac should have the same id - the pdb code')
    code = coeus_ampal_container.code
    coeus_session.rollback()
    pdb_db = coeus_session.query(PdbDB).filter_by(pdb=code).one_or_none()
    if not pdb_db:
        pdb_db = PdbDB(pdb=code)
    else:
        return 0
    pdbe_dbs = [PdbeDB(pdb=pdb_db, mmol=a.tags['mmol_number'], preferred=a.tags['preferred_mmol'])
                for a in coeus_ampal_container]
    cutoff_dbs = coeus_session.query(CutoffDB).all()
    atlas_dbs = coeus_session.query(AtlasDB).all()
    # Start a series of loops over aspects of the structure in order to populate the structural tables in the database.
    for i, a in enumerate(coeus_ampal_container):
        amino_acid_dbs = amino_acid_dbs_from_protein(a)
        helix_count = 0
        # Helix dict with helix numbers as keys and the corresponding HelixDB instances as values.
        helix_dbs = {}
        r_dbs = {}
        g = a.tags['knob_group'].graph
        pdbe_db = pdbe_dbs[i]
        for c in a:
            chain_db = ChainDB(pdbe=pdbe_db, chain=c.id)
            for h in c.helices:
                helix_dbs[helix_count] = (HelixDB(chain=chain_db, number=helix_count, length=len(h)))
                helix_count += 1
            icode = 0
            for r in c:
                # grab the AminoAcidDB instance corresponding to the residue. (None if no such amino acid exists).
                r_amino_acid = amino_acid_dbs[(r.mol_code, r.mol_letter)]
                # Compile residue information for the database from the residue tags.
                torsion_angles = dict(phi=r.tags['phi'], psi=r.tags['psi'], omega=r.tags['omega'])
                for k, v in torsion_angles.items():
                    if type(v) == str:
                        torsion_angles[k] = None
                    elif (type(v) == float) or (type(v) == numpy.float64):
                        torsion_angles[k] = float('{0:.1f}'.format(v))
                if (r.insertion_code == ' ') or (r.insertion_code == ''):
                    icode = 0
                else:
                    icode += 1
                if 'helix_number' in r.tags.keys():
                    r_helix_db = helix_dbs[r.tags['helix_number']]
                else:
                    r_helix_db = None
                # Create a ResidueDB instance from the properties compiled above.
                try:
                    dssp_acc = r.tags['dssp_acc']
                except KeyError:
                    dssp_acc = None
                try:
                    dssp_ss = r.tags['secondary_structure']
                except KeyError:
                    dssp_ss = None
                r_db = ResidueDB(acc=dssp_acc,
                                 phi=torsion_angles['phi'],
                                 psi=torsion_angles['psi'],
                                 omega=torsion_angles['omega'],
                                 resno=int(r.id),
                                 icode_int=icode,
                                 dssp_ss=dssp_ss,
                                 rpr=r.tags['rise_per_residue'],
                                 rpt=r.tags['residues_per_turn'],
                                 roc=r.tags['radius_of_curvature'],
                                 amino_acid=r_amino_acid,
                                 chain=chain_db,
                                 helix=r_helix_db
                                 )

                # r_dbs is a dict with unique_ids as keys (an ampal Residue property),
                # and ResidueDB instances as values.
                r_dbs[r.unique_id] = r_db
        if 'knob_group' in a.tags:
            kg = a.tags['knob_group']
        for kih in kg:
            packing_angle = kih.packing_angle
            if packing_angle is not None:
                packing_angle = float('{0:.3f}'.format(packing_angle))
            knob_db = KnobDB(packing_angle=packing_angle,
                             knob_type=kih.knob_type(cutoff=7.0),
                             max_cv_dist=float('{0:.2f}'.format(kih.max_kh_distance)),
                             residue=r_dbs[kih.knob_residue.unique_id])
            for hole_index in range(len(kih.hole)):
                hole_uid = kih.hole_residues[hole_index].unique_id
                HoleDB(residue=r_dbs[hole_uid], knob=knob_db, hole_res_type=hole_index)
        for cutoff_db in cutoff_dbs:
            connected_components = sorted_connected_components(kg.filter_graph(g=g,
                                                                               cutoff=cutoff_db.scut,
                                                                               min_kihs=cutoff_db.kcut))
            for cc_num, cc in enumerate(connected_components):
                storage_changed = store_graph(cc)
                cc_name = get_graph_name(cc)
                atlas_db = None
                if not storage_changed:
                    atlas_db = next(filter(lambda x: x.name == cc_name, atlas_dbs), None)
                if atlas_db is None:
                    populate_atlas()
                    atlas_db = coeus_session.query(AtlasDB).filter_by(name=cc_name).one()
                graph_db = GraphDB(connected_component=cc_num, cutoff=cutoff_db, atlas=atlas_db, pdbe=pdbe_db)
                for node, degree in cc.degree().items():
                    graph_helix_db = GraphHelixDB(graph=graph_db, helix=helix_dbs[node], degree=degree)
    coeus_session.add(pdb_db)
    try:
        coeus_session.commit()
    except IntegrityError:
        coeus_session.rollback()
        return 0
    coeus_session.rollback()
    return 1


def populate_pdb(a):
    """Populate coeus table 'pdb' from an Assembly.

    Parameters
    ----------
    a : Assembly

    Returns
    -------
    obj : coeus PdbDB object, or None.
    """
    code = a.id
    pdb_db = PdbDB(pdb=code)
    simplified_add_object(pdb_db)
    return


def populate_pdbe(a):
    """Populate PdbeDB from an Assembly.

    Parameters
    ----------
    a : Assembly

    Returns
    -------
    obj : coeus PdbeDB object, or None.

    Raises
    ------
    ValueError
        If tags 'mmol' and 'preferred' are not both present.
    sqlalchemy.orm.exc.NoResultFound
        If no pdb associated with pdbe is found.
    sqlalchemy.orm.exc.MultipleResultsFound
        If more than one pdb associated with pdbe is found.
    """
    code = a.id
    try:
        mmol = a.tags['mmol_number']
        preferred = int(a.tags['preferred_mmol'])
    except KeyError:
        raise ValueError('Assembly must have tags \'mmol\' and \'preferred\' in order to proceed')
    pdb = PdbDB(pdb=code)
    pdbe_db = PdbeDB(pdb_id=pdb.id, mmol=mmol, preferred=preferred)
    simplified_add_object(pdbe_db)
    return


def populate_chain(a):
    """Populate ChainDB from an Assembly.

    Parameters
    ----------
    a : Assembly

    Returns
    -------
    created_objs : list, or None
        List of Django model ChainDB objects created, if any.
    """
    code = a.id
    try:
        mmol = a.tags['mmol_number']
    except KeyError:
        raise ValueError('Assembly must have tag \'mmol\' in order to proceed')
    pdbe = coeus_session.query(PdbeDB).filter(PdbDB.pdb == code, PdbeDB.mmol == mmol).one()
    # List of dictionaries containing kwargs for feeding to populate_model.
    data_dicts = [dict(pdbe=pdbe, chain=c.id) for c in a._polymers]
    return populate_model(model=ChainDB, data_dicts=data_dicts)


def populate_protein(code):
    """ Populate ProteinDB model from a pdb code.

    Notes
    -----
    Should only be run after the PdbDB instance corresponding to code has been added to the database.
    Requires functions from isambard.add_ons.filesystem.

    Parameters
    ----------
    code : str
        pdb code

    Returns
    -------
    bool :
        1 if object is created.
        0 if object is not created due to an IntegrityError (it is already in the database).
    """
    # Start with a clean session.
    coeus_session.rollback()
    fs = FileSystem(code)
    mmcif_dict = fs.protein_info
    pdb_db = coeus_session.query(PdbDB).filter_by(pdb=code).one_or_none()
    if not pdb_db:
        coeus_session.rollback()
        return
    table_columns = ProteinDB.__table__.columns
    kwargs = {}
    for col_name, col in table_columns.items():
        # mmcif_dict keys must be named correctly - i.e. have the same names as in ProteinDB.
        if col_name in mmcif_dict.keys():
            v = mmcif_dict[col_name]
            if not v:
                continue
            elif col.type.python_type == str:
                # Avoid annoying sql problems caused by quotes within strings.
                v = v.replace("'", "")
                # Ensure that string is not too long for the table.
                max_len = col.type.length
                v = v[:max_len]
            elif col.type.python_type == decimal.Decimal:
                # Format floats to the correct precision level for the table.
                scale = col.type.scale
                v = float('{0:.{1}f}'.format(v, scale))
            kwargs[col_name] = v
    # Create ProteinDB instance. Add to session and commit if possible. If not, rollback.
    protein_db = ProteinDB(pdb=pdb_db, **kwargs)
    coeus_session.add(protein_db)
    try:
        coeus_session.commit()
    except IntegrityError:
        coeus_session.rollback()
        return 0
    return 1


def populate_helix(a):
    """ Populate ProteinDB model from an Assembly.

    Parameters
    ----------
    a : Assembly

    Returns
    -------
    created_objs : list, or None
        List of Django model HelixDB objects created, if any.
    """
    code = a.id
    try:
        mmol = a.tags['mmol_number']
    except KeyError:
        raise ValueError('Assembly must have tag \'mmol\' in order to proceed')
    data_dicts = []
    pdbe_db = coeus_session.query(PdbeDB).filter(PdbDB.pdb == code, PdbeDB.mmol == mmol).one()
    helix_count = 0
    for c in a._polymers:
        helices = c.helices
        # if there are any helices in the chain.
        if len(helices) > 0:
            chain_db = coeus_session.query(ChainDB).filter_by(pdbe_id=pdbe_db.id, chain=c.id).one()
            for helix in helices:
                add_dict = {'chain': chain_db, 'number': helix_count, 'length': len(helix)}
                data_dicts.append(add_dict)
                helix_count += 1
    return populate_model(HelixDB, data_dicts=data_dicts)


def populate_residue(a):
    """Populate ResidueDB from an Assembly.

    Parameters
    ----------
    a : Assembly

    Returns
    -------
    created_objs : list, or None
        List of Django model ResidueDB objects created, if any.
    """
    code = a.id
    try:
        mmol = a.tags['mmol_number']
    except KeyError:
        raise ValueError('Assembly must have tag \'mmol\' in order to proceed')

    pdbe = coeus_session.query(PdbeDB).filter(PdbDB.pdb == code, PdbeDB.mmol == mmol).one()
    chain_dict = {c.id: coeus_session.query(ChainDB).filter(ChainDB.pdbe==pdbe,
                                                            ChainDB.chain==c.id).one()
                  for c in a._polymers}
    chain_ids = [x.id for x in chain_dict.values()]
    helix_dbs = coeus_session.query(HelixDB).filter(HelixDB.chain_id.in_(chain_ids)).all()
    residue_dbs = []
    for residue in a.get_monomers(ligands=False):
        amino_acid = coeus_session.query(AminoAcidDB).filter_by(code=residue.mol_code, letter=residue.mol_letter).one()
        add_dict = dict(chain=chain_dict[residue.ampal_parent.id],
                        amino_acid=amino_acid,
                        acc=residue.tags['dssp_acc'],
                        dssp_ss=residue.tags['secondary_structure'],
                        rpr=residue.tags['rise_per_residue'],
                        rpt=residue.tags['residues_per_turn'],
                        roc=residue.tags['radius_of_curvature'],
                        resno=int(residue.id),
                        )
        torsion_angles = dict(phi=residue.tags['phi'], psi=residue.tags['psi'], omega=residue.tags['omega'])
        for k, v in torsion_angles.items():
            if type(v) == str:
                torsion_angles[k] = None
            elif (type(v) == float) or (type(v) == numpy.float64):
                torsion_angles[k] = float('{0:.1f}'.format(v))
        add_dict.update(torsion_angles)
        if residue.insertion_code == ' ':
            icode = 0
        else:
            icode += 1
        add_dict['icode_int'] = icode
        if 'helix_number' in residue.tags.keys():
            helix_db = next((h for h in helix_dbs if h.number==residue.tags['helix_number']), None)
            add_dict['helix'] = helix_db
        residue_db = ResidueDB(**add_dict)
        if residue_db not in coeus_session:
            coeus_session.add(residue_db)
        else:
            continue
    coeus_session.commit()
    coeus_session.rollback()
    return


def populate_graph_helix(dbi):
    """Populate GraphHelixDB from a DatabaseItem.

    Parameters
    ----------
    dbi : isambard.DatabaseItem
        An instance of the DatabaseItem class.

    Returns
    -------
    created_objs : list, or None
        List of Django model GraphHelixDB objects created, if any.
    """
    if not dbi.mmols:
        return
    data_dicts = []
    for mmol in dbi.mmols:
        if not mmol.graphs:
            continue
        pdbe = PdbeDB.objects.filter(pdb__pdb=dbi.code, mmol=mmol.number).get()
        chain_dict = {}
        for c in mmol.chains:
            chain = ChainDB.objects.filter(pdbe=pdbe, chain=c).get()
            chain_dict[c] = chain
        helix_dict = {}
        for h in mmol.helices:
            helix = HelixDB.objects.filter(chain=chain_dict[h.chain], number=h.number).get()
            helix_dict[h.number] = helix
        for g in mmol.graphs:
            cutoff = CutoffDB.objects.filter(kcut=g.kcut, scut=g.scut).get()
            graph = GraphDB.objects.filter(cutoff=cutoff, pdbe=pdbe, connected_component=g.connected_component).get()
            for helix_number in g.helix_numbers:
                helix = helix_dict[helix_number]
                degree = g.degree[helix_number]
                add_dict = {'graph': graph, 'helix': helix, 'degree': degree}
                data_dicts.append(add_dict)
    return populate_model(model=GraphHelixDB, data_dicts=data_dicts)


def populate_cdhit_full(pdb_cdhit_dict):
    """Populate CdhitFullDB from the dictionary output from get_pdb_cdhit_dict().

    Parameters
    ----------
    pdb_cdhit_dict : dict
        Output of get_pdb_cdhit_dict(). (See that function for more information).
        Keys are PDB codes.
        Values are booleans indicating whether the sequences are representative.

    Returns
    -------
    created_objs : list, or None
        List of Django model CdhitFullDB objects created, if any.
    """
    data_dicts = []
    for key, val in pdb_cdhit_dict.items():
        try:
            pdb = PdbDB.objects.filter(pdb=key).get()
            c90, c80, c70, c60 = val
            add_dict = {'pdb': pdb, 'c90': c90, 'c80': c80, 'c70': c70, 'c60': c60}
            data_dicts.append(add_dict)
        except:
            continue
    return populate_model(model=CdhitFullDB, data_dicts=data_dicts)


def check_new_aa(dbi):
    """Check for new amino acids. If found, add them to AminoAcidDB.

    Cycle through each aa_code, aa_letter pairing for all residues in the dbi and check for new amino acids.
    If a new amino acid is found, add it to the json file by running get_aa_info()
     and then add_amino_acid_to_json().
    Populate the model AminoAcidDB with any new amino acids.
    Print new amino acid codes.

    Parameters
    ----------
    dbi : isambard.DatabaseItem
        An instance of the DatabaseItem class.

    Raises
    ------
    ValueError
        If a 3-letter amino acid code is paired with a different 1-letter code in the database, compared to in the dbi.

    Returns
    -------
    None
    """
    for mmol in dbi.mmols:
        rset = set([(r.aa_code, r.aa_letter) for r in mmol.residues])
        for code, letter in rset:
            # Check for entry in the amino_acid table.
            count = AminoAcidDB.objects.filter(code=code, letter=letter).count()
            if count == 0:
                # Check that it's not anything weird happening with the one-letter code
                test = AminoAcidDB.objects.filter(code=code).count()
                if test != 0:
                    raise ValueError("Existing code {0} paired with letter {1}. Check for problems.".format(code,
                                                                                                            letter))
                else:
                    # Get information associated with amino acid that is not in amino_acid_table.
                    aa_dict = get_aa_info(new_aa_code=code)
                    # Add this information to the amino_acids.json file.
                    add_amino_acid_to_json(**aa_dict)
                    # Add the information to the amino_acid table.
                    data_dicts = [{'code': code, 'letter': letter,
                                   'modified': aa_dict['modified'], 'description': aa_dict['description']}]
                    populate_model(model=AminoAcidDB, data_dicts=data_dicts)
                    print("Added amino acid {0}".format(code))
    return


__author__ = 'Jack W. Heal'
__status__ = 'Development'
