from .solenoid import HelixPair, Solenoid
from ..polymer_specs.nucleic_acid_strand import _helix_parameters


_repeating_unit_parameters = {
    'axis_dist': 5.45,
    'z_shift': 3.99,
    'ia_1': 32.06,
    'ia_2': -169.61,
    'splay': 12.65,
    'off_plane': 160.60
}


_tale_parameters = {
    'radius': 17.08,
    'rise': 3.46,
    'angular_offset': 33.81
}


def z_rot_adjust(z):
    """Calculate angular adjustment required for shift up the z-axis."""
    bdna_parameters = _helix_parameters['b_dna']
    pitch = bdna_parameters[0] * bdna_parameters[1]
    return 360.0*(z/pitch)


class TaleHelixPair(HelixPair):
    """Helix pair class with reduced number of parameters."""
    def __init__(self):
        rad = _repeating_unit_parameters['axis_dist']
        zshift = _repeating_unit_parameters['z_shift']
        phi1 = _repeating_unit_parameters['ia_1']
        phi2 = _repeating_unit_parameters['ia_2']
        splay = _repeating_unit_parameters['splay']
        off_plane_rot = _repeating_unit_parameters['off_plane']
        super().__init__(aas=(9, 8), axis_distances=(-rad, rad), z_shifts=(0, zshift),
                         phis=(phi1, phi2), splays=(0, splay), off_plane=(0, off_plane_rot))
        self.relabel_all()


class Tale(Solenoid):
    def __init__(self, repeats, dna_shift, up=True):
        """

        Parameters
        ----------
        repeats
        dna_shift
        up
        """
        ru = TaleHelixPair()
        # These operations rotate the repeating unit into the correct orientation
        ru.rotate(125.974, (1, 0, 0))
        ru.rotate(-65.135, (0, 1, 0))
        ru.rotate(156.975, (0, 0, 1))
        # Solenoid para
        rad = _tale_parameters['radius']
        rise = _tale_parameters['rise']
        rot_ang = _tale_parameters['angular_offset']
        # Rotate TALE into phase with DNA
        self.dna_shift = dna_shift
        if up:
            super().__init__(ru, repeats, rad, rise, rot_ang, 'right')
            self.rotate(-138 + z_rot_adjust(self.dna_shift), (0, 0, 1))
            self.translate((0, 0, self.dna_shift))
        else:
            super().__init__(ru, repeats, rad, -rise, -rot_ang, 'right')
            self.rotate(55 - z_rot_adjust(self.dna_shift), (0, 0, 1))
            self.rotate(180, (1, 0, 0))
            self.translate((0, 0, self.dna_shift))


__author__ = 'Christopher W. Wood'
__status__ = 'Development'
