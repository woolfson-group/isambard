import unittest

from hypothesis import given, settings
from hypothesis.strategies import integers, floats, text

import isambard_dev as isambard


class TestBoundTaleModelling(unittest.TestCase):
    def test_tale_helix_pair(self):
        thp = isambard.specifications.assembly_specs.tale.TaleHelixPair()
        parameters = isambard.specifications.assembly_specs.tale._repeating_unit_parameters
        self.assertAlmostEqual(thp.axis_distances[0], -parameters['axis_dist'], places=3)
        self.assertAlmostEqual(thp.axis_distances[1], parameters['axis_dist'], places=3)
        self.assertEqual(thp.z_shifts[0], 0)
        self.assertAlmostEqual(thp.z_shifts[1], parameters['z_shift'])
        self.assertAlmostEqual(thp.phis[0], parameters['ia_1'], places=3)
        self.assertAlmostEqual(thp.phis[1], parameters['ia_2'], places=3)
        self.assertEqual(thp.splays[0], 0)
        self.assertAlmostEqual(thp.splays[1], parameters['splay'], places=3)
        self.assertEqual(thp.off_plane[0], 0)
        self.assertAlmostEqual(thp.off_plane[1], parameters['off_plane'], places=3)

    @given(integers(max_value=50))
    @settings(max_examples=20)
    def test_tale_model(self, repeats):
        if repeats < 1:
            with self.assertRaises(ValueError):
                tale = isambard.specifications.assembly_specs.tale.Tale(repeats)
        else:
            tale = isambard.specifications.assembly_specs.tale.Tale(repeats)
            self.assertEqual(len(tale), repeats*2)

    @given(integers(min_value=1, max_value=50), integers(min_value=1, max_value=50), floats(), floats(), text('ATGC'))
    @settings(max_examples=10)
    def test_tale_dna_model(self, up_repeats, down_repeats, up_shift, down_shift, dna_sequence):
        if len(dna_sequence) < 1:
            with self.assertRaises(ValueError):
                tale_dna = isambard.specifications.assembly_specs.tale.TaleDNA(
                    up_repeats, down_repeats, up_shift, down_shift, dna_sequence)
        else:
            tale_dna = isambard.specifications.assembly_specs.tale.TaleDNA(
                up_repeats, down_repeats, up_shift, down_shift, dna_sequence)
            if up_repeats > 0:
                self.assertEqual(len(tale_dna['A']), up_repeats * 17)  # Number of residues
            else:
                self.assertEqual(len(tale_dna['A']), 0)
            if down_repeats > 0:
                self.assertEqual(len(tale_dna['B']), down_repeats * 17)
            else:
                self.assertEqual(len(tale_dna['B']), 0)
            self.assertEqual(tale_dna['C'].sequence, 'D' + ' D'.join(dna_sequence))
            reverse_compliment = (
                isambard.specifications.assembly_specs.nucleic_acid_duplex.generate_antisense_sequence(dna_sequence))
            self.assertEqual(
                tale_dna['D'].sequence, 'D' + ' D'.join(reverse_compliment))
