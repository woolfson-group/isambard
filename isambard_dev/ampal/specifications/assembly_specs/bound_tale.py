from ampal.specifications.assembly_specs.solenoid import HelixPair


_repeating_unit_parameters = {
    'axis_dist': 5.45,
    'z_shift': 3.99,
    'ia_1': 32.06,
    'ia_2': -169.61,
    'splay': 12.65,
    'off_plane': 160.60
}


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
