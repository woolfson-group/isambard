import numpy

# import the relevant secondary structure
from ampal.assembly import Assembly
from ampal.specifications.polymer_specs.helix import Helix
from tools.geometry import dihedral


# TODO Remove rosetta mode, move average_2_points to geometry

def gen_octahedron(el):
    """Generates vertices of an octahedron with a defined edge length.

    Parameters
    ----------
    el : float
        The edge length of the polyhedron in arbitrary units.

    Returns
    -------
    vertices : [triple]
        List containing coordinates of the vertices of the polyhedron.
    """
    xy = numpy.sin(numpy.pi / 4) * el
    vertices = [
        (0, 0, xy),  # a
        (0, xy, 0),  # b
        (xy, 0, 0),  # c
        (0, -xy, 0),  # d
        (-xy, 0, 0),  # e
        (0, 0, -xy),  # f
    ]
    return vertices


def gen_snub_disphenoid(el):
    """Generates vertices of a snub-nosed disphenoid with a defined edge length.

    Parameters
    ----------
    el : float
        The edge length of the polyhedron in arbitrary units.

    Returns
    -------
    vertices : [triple]
        List containing coordinates of the vertices of the polyhedron.
    """
    # -z1 to move gh vector off x axis
    x2 = 0.644584 * el
    z1 = 0.578369 * el
    z2 = 0.989492 * el
    z3 = 1.56786 * el

    vertices = [
        (0, el / 2, z3 - z1),  # a
        (0, -el / 2, z3 - z1),  # b
        (0, -x2, z1 - z1),  # c
        (-x2, 0, z2 - z1),  # d
        (0, x2, z1 - z1),  # e
        (x2, 0, z2 - z1),  # f
        (el / 2, 0, 0 - z1),  # g
        (-el / 2, 0, 0 - z1),  # h
    ]
    return vertices


def gen_gyro_square_bipyramid(el):
    """Generates vertices of a gyroelongated bipyramid with a defined edge length.

    Parameters
    ----------
    el : float
        The edge length of the polyhedron in arbitrary units.

    Returns
    -------
    vertices : [(float, float, float)]
        List containing coordinates of the vertices of the polyhedron.
    """
    rl = (0.5 * el) / numpy.sin(numpy.pi / 4)
    zs = (numpy.sin(numpy.pi / 3) * el)
    pl = rl - (el / 2)
    z1 = numpy.sqrt((zs ** 2) - (pl ** 2)) / 2
    rxy = el / 2
    theta = numpy.arccos(rl / el)
    z2 = numpy.sin(theta) * el
    z3 = z1 + z2
    vertices = [
        (0, 0, z3),  # a
        (0, rl, z1),  # b
        (rl, 0, z1),  # c
        (0, -rl, z1),  # d
        (-rl, 0, z1),  # e
        (-rxy, rxy, -z1),  # f
        (rxy, rxy, -z1),  # g
        (rxy, -rxy, -z1),  # h
        (-rxy, -rxy, -z1),  # k
        (0, 0, -z3),  # l
    ]
    return vertices


def gen_icosahedron(el):
    """Generates vertices of an icosahedron with a defined edge length.

    Parameters
    ----------
    el : float
        The edge length of the polyhedron in arbitrary units.

    Returns
    -------
    vertices : [(float, float, float)]
        List containing coordinates of the vertices of the polyhedron.
    """
    rl = (el / 2) / numpy.sin(numpy.pi / 5)
    x2 = numpy.cos(numpy.pi / 2 - (2 * numpy.pi / 5)) * rl
    y2 = numpy.sin(numpy.pi / 2 - (2 * numpy.pi / 5)) * rl
    x3 = numpy.sin(numpy.pi - 2 * (2 * numpy.pi / 5)) * rl
    y3 = numpy.cos(numpy.pi - 2 * (2 * numpy.pi / 5)) * rl
    x4 = numpy.sin(numpy.pi / 5) * rl
    y4 = numpy.cos(numpy.pi / 5) * rl
    x5 = numpy.cos((3 * numpy.pi / 5) - (numpy.pi / 2)) * rl
    y5 = numpy.sin((3 * numpy.pi / 5) - (numpy.pi / 2)) * rl
    zs = numpy.sqrt(el ** 2 - (el / 2) ** 2)
    z1 = numpy.sqrt(zs ** 2 - (rl - y4) ** 2) / 2
    z2 = numpy.sqrt(el ** 2 - rl ** 2)
    vertices = [
        (0, 0, z1 + z2),  # a
        (x2, y2, z1),  # b
        (x3, -y3, z1),  # c
        (-x3, -y3, z1),  # d
        (-x2, y2, z1),  # e
        (0, rl, z1),  # f
        (x4, y4, -z1),  # g
        (x5, -y5, -z1),  # h
        (0, -rl, -z1),  # k
        (-x5, -y5, -z1),  # l
        (-x4, y4, -z1),  # m
        (0, 0, -z1 + -z2),  # n
    ]
    return vertices


class DeltaProt(Assembly):
    """Generates a deltahedral protein with specified rib orientation and length.

    Notes
    -----
    If centre_helices is True then the helices will be translated along the ribs so the centre of the helices is
    equal to the centre of the rib. They will also be rotated so that the middle most residue (rounded-up) will be
    rotated to face the centre of the assembly.

    All distances in angstroms.

    Parameters
    ----------
    conformation : string
        String describing the deltaprot form to be modelled. See the keys of the rib_orient dictionary for details of
        the various forms available.
    aa : int
        Number of amino acids in the individual helices of the assembly. Keyword argument, default value = 8.
    centre_helices : bool
        If centre_helices is True then the helices will be translated along the ribs so the centre of the helices is
        equal to the centre of the rib. They will also be rotated so that the middle most residue (rounded-up) will be
        rotated to face the centre of the assembly. Keyword argument, default value = True.
    centred_ca : int
        Residue on which the helices are centred. Keyword argument, default value = 1.

    Attributes
    ----------
    ax_trans_v : [float]

    assembly : [[(float, float, float)]]
        List containing a list for each chain in the assembly, each of which contain the coordinates of the mainchain/cb
        atoms.
    """
    rib_orientations = {
        # 3
        'b3iii': [[(0, 1), (2, 3), (4, 5)], [0.0, 270.0, 180.0]],  # ab, cd, ef
        'b3nnn': [[(0, 2), (1, 4), (3, 5)], [0.0, 90.0, 180.0]],  # ac, be, cf
        # 4
        'b4iiiix': [[(0, 5), (1, 3), (2, 7), (4, 6)], [319.0, 319.0, 221.0, 221.0]],  # af, bd, ch, eg
        'b4iiiiy': [[(0, 1), (2, 5), (3, 4), (6, 7)], [0.0, 90.0, 270.0, 180.0]],  # ab, cf, de, gh
        'b4nnnnx': [[(0, 3), (1, 5), (2, 6), (4, 7)], [41.0, 41.0, 139.0, 139.0]],  # ad, bf, cg, eh
        'b4nnnny': [[(0, 1), (2, 3), (4, 5), (6, 7)], [0.0, 270.0, 270.0, 180.0]],  # ab, cd, ef, gh
        'b4iiin': [[(0, 4), (1, 5), (2, 6), (3, 7)], [0.0, 32.0, 127.0, 180.0]],  # ae, bf, cg, dh
        'b4innn': [[(0, 5), (1, 2), (3, 7), (4, 6)], [319.0, 0.0, 180.0, 221.0]],  # af, bc, dh, eg
        'b4inin': [[(0, 4), (1, 2), (3, 7), (5, 6)], [0.0, 0.0, 180.0, 180.0]],  # ae, bc, dh, fg
        'l4iin': [[(0, 1), (2, 6), (3, 7), (4, 5)], [0.0, 139.0, 180.0, 270.0]],  # ab, cg, dh, ef
        'l4inn': [[(0, 1), (2, 5), (3, 7), (4, 6)], [0.0, 90.0, 180.0, 221.0]],  # ab, cf, dh, eg
        'h4i.n': [[(0, 5), (1, 3), (2, 6), (4, 7)], [319.0, 319.0, 139.0, 139.0]],  # af, bd, cg, eh
        # 5
        'b5iiiin': [[(0, 4), (1, 2), (3, 8), (5, 6), (7, 9)], [0.0, 310.0, 270.0, 230.0, 180.0]],
        # ae, bc, dk, fg, hl
        'b5innnn': [[(0, 1), (2, 6), (3, 4), (5, 8), (7, 9)], [0.0, 90.0, 310.0, 130.0, 180.0]],
        # ab, cg, de, fk, hl
        'b5iinin': [[(0, 4), (1, 5), (2, 6), (3, 8), (7, 9)], [0.0, 90.0, 90.0, 270.0, 180.0]],
        # ae, bf, cg, dk, hl
        'b5ininn': [[(0, 1), (2, 6), (3, 8), (4, 5), (7, 9)], [0.0, 90.0, 270.0, 270.0, 180.0]],
        # ab, cg, dk, ef, hl
        'l5iiin': [[(0, 2), (1, 6), (3, 4), (5, 8), (7, 9)], [0.0, 270.0, 310.0, 130.0, 180.0]],
        # ac, bg, de, fk, hl
        'l5innn': [[(0, 3), (1, 4), (2, 6), (5, 8), (7, 9)], [0.0, 50.0, 90.0, 130.0, 180.0]],  # ad, be, cg, fk, hl
        'l5inni': [[(0, 2), (1, 6), (3, 8), (4, 5), (7, 9)], [0.0, 270.0, 270.0, 270.0, 180.0]],
        # ac, bg, dk, ef, hl
        'l5niin': [[(0, 3), (1, 5), (2, 6), (4, 8), (7, 9)], [0.0, 90.0, 90.0, 90.0, 180.0]],  # ad, bf, cg, ek, hl
        'h5i.i': [[(0, 1), (2, 3), (4, 8), (5, 6), (7, 9)], [0.0, 310.0, 90.0, 230.0, 180.0]],  # ab, cd, ek, fg, hl
        'h5n.n': [[(0, 4), (1, 6), (2, 3), (5, 8), (7, 9)], [0.0, 270.0, 310.0, 130.0, 180.0]],
        # ae, bg, cd, fk, hl
        # 6
        'b6ininin': [[(0, 3), (1, 7), (2, 8), (4, 9), (5, 10), (6, 11)], [0.0, 270.0, 270.0, 90.0, 90.0, 180.0]],
        # ad, bh, ck, el, fm, gn
        'b6iiniin': [[(0, 2), (1, 7), (3, 8), (4, 9), (5, 10), (6, 11)], [0.0, 270.0, 90.0, 90.0, 90.0, 180.0]],
        # ac, bh, dk, el, fm, gn
        'b6inninn': [[(0, 3), (1, 2), (4, 9), (5, 10), (6, 11), (7, 8)], [0.0, 302.0, 90.0, 90.0, 180.0, 238.0]],
        # ad, bc, el, fm, gn, hk
        'l6innni': [[(0, 1), (2, 3), (4, 9), (5, 10), (6, 11), (7, 8)], [0.0, 302.0, 90.0, 90.0, 180.0, 238.0]],
        # ab, cd, el, fm, gn, hk
        'l6niiin': [[(0, 1), (2, 7), (3, 8), (4, 9), (5, 10), (6, 11)], [0.0, 90.0, 90.0, 90.0, 90.0, 180.0]],
        # ab, ch, dk, el, fm, gn
        'h6i.i.i': [[(0, 4), (1, 2), (3, 9), (5, 6), (7, 8), (10, 11)], [0.0, 302.0, 270.0, 270.0, 238.0, 180.0]],
        # ae, bc, dl, fg, hk, mn
        'h6n.n.n': [[(0, 1), (2, 7), (3, 4), (5, 10), (6, 11), (8, 9)], [0.0, 90.0, 302.0, 90.0, 180.0, 238.0]],
        # ab, ch, de, fm, gn, kl
        's6': [[(0, 5), (1, 7), (2, 3), (4, 9), (6, 10), (8, 11)], [0.0, 270.0, 302.0, 90.0, 122.0, 180.0]],
        # af, bh, cd, el, gm, kn
    }

    choose_delta = {
        3: gen_octahedron,
        4: gen_snub_disphenoid,
        5: gen_gyro_square_bipyramid,
        6: gen_icosahedron,
    }

    # TODO: Separate out into start_and_ends and build().
    def __init__(self, conformation, aa=8, centre_helices=True, centred_ca=1):
        super(DeltaProt, self).__init__()

        conformation = conformation.lower()
        if conformation not in self.rib_orientations.keys():
            raise ValueError('Invalid deltaprot conformation {}'.format(conformation))
        else:
            self.conformation = conformation

        # Move choose_delta and rib_orient outside class??
        self.rib_num = int(self.conformation[1])
        self.rib_len = 11.0
        self.aa = aa

        # antiparallel flag
        self.ap = [0] * self.rib_num

        self.centred_ca = centred_ca - 1
        # Alter ax_trans to position helix at the centre of the primitives if centre_helices is selected
        if centre_helices:
            self.ax_trans_adjust = ((self.rib_len / 2.0) - (((self.aa - 1) * 1.52) / 2.0))
        else:
            self.ax_trans_adjust = 0

        self.build()

    def conformation_edges(self):
        dv = self.deltahedron_vertices()
        # vertices to choose specific to the conformation
        rib_vertices = self.rib_orientations[self.conformation][0]
        edges = []
        for i, vertices in enumerate(rib_vertices):
            v1, v2 = vertices
            # flip if antiparallel
            if self.ap[i]:
                edges.append((dv[v2], dv[v1]))
            else:
                edges.append((dv[v1], dv[v2]))
        return edges

    def deltahedron_vertices(self):
        # Number and length of ribs determine the vertices of the deltaprot shape
        return self.choose_delta[self.rib_num](self.rib_len)

    @property
    def centre(self):
        dv = self.deltahedron_vertices()
        centre = sum([numpy.array(x) for x in dv]) / len(dv)
        return centre

    def build(self):
        """Uses input parameters to build model in selected OldDeltaProt topology.

        Rerunning build overwrites previous assembly."""
        polymers = []
        for start, end in self.conformation_edges():
            start = numpy.array(start)
            end = numpy.array(end)
            helix = Helix(aa=self.aa)
            helix.move_to(start=start, end=end)
            helix.translate(self.ax_trans_adjust * helix.axis.unit_tangent)
            ax_rot = dihedral(self.centre, end, start, helix[self.centred_ca]['CA']._vector)
            helix.rotate(angle=ax_rot, axis=helix.axis.unit_tangent, point=helix.axis.midpoint)
            polymers.append(helix)

        # TODO make Assembly method reset_ampal_parents. (ampal_children?) May need quivalent at polymer level.
        self._molecules = polymers[:]
        for polymer in self._molecules:
            polymer.ampal_parent = self
            for monomer in polymer._monomers:
                monomer.ampal_parent = polymer
        self.relabel_polymers()  # relabel to give each a chain label
        self.relabel_atoms()
        return


__author__ = 'Christopher W. Wood'
__status__ = 'Development'
