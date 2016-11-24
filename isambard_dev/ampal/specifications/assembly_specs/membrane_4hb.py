"""Specifications for building models of helical bundles in membranes."""

from ampal.assembly import Assembly
from ampal.specifications.polymer_specs import Helix

import numpy


def find_to_mem_plane_vector(helix, membrane_hp_thickness):
    """Finds the vector required to move the helix centre to the membrane plane.

    Parameters
    ----------
    helix: ampal.secondary_structure.Helix
        The reference helix.
    membrane_hp_thickness: float
        The width of the hydrophobic region of the membrane.

    Returns
    -------
    adjust_v: (float, float, float)
        The vector require to move the helix centre to the membrane plane.
    """
    adjust_v = (0, 0, (membrane_hp_thickness/2) - helix.centre_of_mass[2])
    return adjust_v


def find_tilt_to_mem_angle(helix, membrane_hp_thickness):
    """Find the angle required to keep the helix in the hydrophobic membrane region.

    Parameters
    ----------
    helix: ampal.secondary_structure.Helix
        The reference helix.
    membrane_hp_thickness: float
        The width of the hydrophobic region of the membrane.

    Returns
    -------
    tilt_angle: float
        The tilt angle require
    """
    helix_length = helix.rise_per_residue * len(helix)  # In Ångströms
    if helix_length > membrane_hp_thickness:
        tilt_angle = numpy.arccos(membrane_hp_thickness/helix_length)
        return tilt_angle
    return 0


class Membrane4HelixBundle(Assembly):
    default_sh_rot = (0, 90, 180, 270)

    def __init__(self, aas=(20, 20, 20, 20), radii=(5, 5, 5, 5), ah_rotations=(0, 0, 0, 0),
                 sh_rotations=(0, 0, 0, 0), zshifts=(0, 0, 0, 0), tilt_angles=(0, 0, 0, 0),
                 lmp_rotations=(0, 0, 0, 0),
                 auto_build=True):
        """Generates a model of a membrane 4-helix bundle.

        Parameters
        ----------
        aas: [int, int, int, int]
            The number of amino acids to be modelled for each helix.
        radii: [float, float, float, float]
            The distance between the midpoint of the helices and the
            super-helical axis.
        ah_rotations: [float, float, float, float]
            The rotation, in degrees of the helices around their local
            a-helical axis.
        sh_rotations: [float, float, float, float]
            The rotation, in degrees of the helices about the super-helical
            axis.
        auto_build: bool
            Will automatically run the build method if true.
        """
        super().__init__()
        self.aas = list(aas)
        self.radii = list(radii)
        self.ah_rotations = list(ah_rotations)
        self.sh_rotations = [x + y for x, y in zip(sh_rotations, self.default_sh_rot)]
        self.zshifts = list(zshifts)
        self.tilt_angles = list(tilt_angles)
        self.lmp_rotations = list(lmp_rotations)
        if auto_build:
            self.build()

    def build(self):
        """Builds a model using the parameters defined during the __init__.

        Raises
        ------
        ValueError
            Raises a value error if all of the lists of parameters are not the
            same length.

        """
        self._molecules = []
        parameters = (self.aas, self.radii, self.ah_rotations,
                      self.sh_rotations, self.zshifts, self.tilt_angles, self.lmp_rotations)
        if not all([len(x) == 4 for x in parameters]):
            raise ValueError('All parameter lists must have 4 values.')
        for aa, radius, ah_rot, sh_rot, zshift, tilt_angle, lmp_rot in zip(*parameters):
            helix = Helix(aa)
            helix.rotate(ah_rot, (0, 0, 1))
            helix.translate((radius, 0, zshift))
            helix.rotate(sh_rot, (0, 0, 1), (0, 0, 0))
            tilt_rotate_vector = helix.centre_of_mass - (0, 0, helix.centre_of_mass[2])
            helix.rotate(tilt_angle, tilt_rotate_vector, helix.centre_of_mass)
            # Apply the local membrane perpendicular rotation
            helix.rotate(lmp_rot, (0, 0, 1), helix.centre_of_mass)
            self._molecules.append(helix)
        self.relabel_all()
        return


__author__ = "Christopher W. Wood"
