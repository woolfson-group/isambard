import os
import unittest

import isambard_dev as isambard


class CovalentBondFinderTestCase(unittest.TestCase):
    """Test covalent bond assignment."""

    def count_bb_bond(self, ampal):
        for chain in ampal:
            ch_len = len(chain)
            bonds = isambard.ampal.interactions.find_covalent_bonds(chain.backbone)
            self.assertEqual(len(bonds), (ch_len*4)-1)
        return

    def test_3qy1(self):
        """Find covalent bonds in backbone of 3qy1"""
        test_path = os.path.join(os.path.dirname(isambard.__file__), 'unit_tests', 'testing_files', '3qy1.pdb')
        test_pdb = isambard.external_programs.reduce.assembly_plus_protons(test_path)
        self.count_bb_bond(test_pdb)

    '''def test_2ht0(self):
        test_path = os.path.join(os.path.dirname(isambard.__file__), 'unit_tests', 'testing_files', '2ht0.pdb')
        test_pdb = assembly_plus_protons(test_path)
        self.count_bb_bond(test_pdb)'''

    def test_1ek9(self):
        """Find covalent bonds in backbone of 1ek9"""
        test_path = os.path.join(os.path.dirname(isambard.__file__), 'unit_tests', 'testing_files', '1ek9.pdb')
        test_pdb = isambard.external_programs.reduce.assembly_plus_protons(test_path)
        self.count_bb_bond(test_pdb)


class HBondFinderTestCase(unittest.TestCase):
    """Test hydrogen bond finding mechanisms."""

    def test_3qy1(self):
        """Find number of hydrogen bonds in 3qy1"""
        test_path = os.path.join(os.path.dirname(isambard.__file__), 'unit_tests', 'testing_files', '3qy1.pdb')
        test_pdb = isambard.external_programs.reduce.assembly_plus_protons(test_path)
        self.assertEqual(len(isambard.ampal.interactions.find_hydrogen_bonds(
            test_pdb, dist_range=(1.5, 2.7), angular_cutoff=90.0)), 961)


if __name__ == '__main__':
    unittest.main()