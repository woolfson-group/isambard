import unittest
import os

import isambard_dev as isambard


class KnobGroupTestCase(unittest.TestCase):
    """Tests for class KnobGroup"""
    def setUp(self):
        test_file = os.path.join('unit_tests', 'testing_files', '2ebo_1.mmol')
        self.test_assembly = isambard.ampal.convert_pdb_to_ampal(test_file)

    def test_number_of_kihs(self):
        """ Test there are 38 kihs found at 7.0 cutoff for 2ebo. """
        a = self.test_assembly
        kg = isambard.add_ons.knobs_into_holes.KnobGroup.from_helices(a, min_helix_length=0)
        self.assertTrue(len(kg) == 38)
        kihs_locations = ''
        for i in range(len(a.helices)):
            kihs_locations += str(len([x for x in kg if x.knob_helix.number == i]))
        self.assertTrue(int(kihs_locations) == 914804714)
