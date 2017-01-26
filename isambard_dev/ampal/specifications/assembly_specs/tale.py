"""Specifications for generating models of TALE proteins."""

from ampal import Assembly
from .solenoid import HelixPair, Solenoid
from .nucleic_acid_duplex import DNADuplex
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
    def __init__(self, repeats):
        """Builds a model of a TALE protein.

        Super-helical parameters are based on 3v6t.

        Parameters
        ----------
        repeats: int
            Number of repeating units to model.
        """
        ru = TaleHelixPair()
        # These operations rotate the repeating unit into the correct orientation
        ru.rotate(125.974, (1, 0, 0))
        ru.rotate(-65.135, (0, 1, 0))
        ru.rotate(156.975, (0, 0, 1))
        # Solenoid initialisation
        rad = _tale_parameters['radius']
        rise = _tale_parameters['rise']
        rot_ang = _tale_parameters['angular_offset']
        super().__init__(ru, repeats, rad, rise, rot_ang, 'right')


class TaleDNA(Assembly):
    def __init__(self, up_repeats, down_repeats, up_shift, down_shift, dna_sequence, build=True):
        """Build a model of DNA with two bound TALE proteins.

        Parameters
        ----------
        up_repeats: int
            Number of repeats in the up TALE.
        down_repeats: int
            Number of repeats in the down TALE.
        up_shift: float
            Distance, in angstoms, that the up TALE is shifted along the DNA.
        down_shift: float
            Distance, in angstoms, that the down TALE is shifted along the DNA.
        dna_sequence: str
            DNA sequence of the forward strand of the DNA duplex. The reverse
            compliment will be generated automatically.
        """
        super().__init__()
        self.up_repeats = up_repeats
        self.up_shift = up_shift
        self.down_repeats = down_repeats
        self.down_shift = down_shift
        self.dna_sequence = dna_sequence
        if build:
            self.build()

    def build(self):
        """Builds a model using the parameters provided during initialisation."""
        self._molecules = []
        # Up TALE
        up_tale = Tale(self.up_repeats)
        up_tale.rotate(-138 + z_rot_adjust(self.up_shift), (0, 0, 1))
        up_tale.rotate(55 - z_rot_adjust(self.up_shift), (0, 0, 1))
        up_tale.translate((0, 0, self.up_shift))
        # Down TALE
        down_tale = Tale(self.down_repeats)
        down_tale.rotate(55 - z_rot_adjust(self.down_shift), (0, 0, 1))
        down_tale.rotate(180, (1, 0, 0))
        down_tale.rotate(180, (1, 0, 0))
        down_tale.translate((0, 0, self.down_shift))
        # DNA
        dna = DNADuplex.from_sequence(self.dna_sequence)
        # Relabel for convenience
        self._default_tale_labelling(up_tale, 'A')
        self._default_tale_labelling(down_tale, 'B')
        dna.relabel_polymers(['C', 'D'])
        self.extend(up_tale)
        self.extend(down_tale)
        self.extend(dna)

    @staticmethod
    def _default_tale_labelling(assembly, chain_label):
        """Labels a tale as a single chain."""
        for mol in assembly:
            mol.id = chain_label
        monomers = list(assembly.get_monomers())
        for res, label in zip(monomers, range(1, len(monomers) + 1)):
            res.id = str(label)
        return


__author__ = 'Christopher W. Wood'
__status__ = 'Development'
