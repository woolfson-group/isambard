import os
import unittest

import isambard_dev as isambard

class RadiusOfGyrationTestCase(unittest.TestCase):

    def test_3qy1(self):

        file_path=os.path.join(os.path.dirname(isambard.__file__), 'unit_tests', 'testing_files', '3qy1.pdb')
        m = isambard.ampal.convert_pdb_to_ampal(file_path)
        r_o_g = m.radius_of_gyration
        self.assertAlmostEqual(r_o_g,20.948377)


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

class NPiStarInteractionFinderTestCase(unittest.TestCase):
    """Test n-to-pi-star interaction finding mechanism."""

    def test_3qy1(self):
        test_path=os.path.join(os.path.dirname(isambard.__file__),'unit_tests','testing_files','3qy1.pdb')
        test_pdb = isambard.ampal.convert_pdb_to_ampal(test_path)
        self.assertEqual(len(isambard.ampal.interactions.find_N_pis(test_pdb)), 239)

class SaltBridgeFinderTestCase(unittest.TestCase):
    """Test salt bridge finding mechanism"""

    def test_3qy1(self):
        """Find number of salt bridges in 3qy1"""
        test_path = os.path.join(os.path.dirname(isambard.__file__), 'unit_tests', 'testing_files', '3qy1.pdb')
        test_pdb = isambard.ampal.convert_pdb_to_ampal(test_path)
        self.assertEqual(len(isambard.ampal.interactions.find_salt_bridges(test_pdb)),29)

class MetPiInteractionFinderTestCase(unittest.TestCase):
    """Test Met-Pi interaction finding mechanism"""

    def test_3qy1(self):
        """Find number of Met-Pi interactions in 3qy1"""
        test_path = os.path.join(os.path.dirname(isambard.__file__), 'unit_tests','testing_files','3qy1.pdb')
        test_pdb = isambard.ampal.convert_pdb_to_ampal(test_path)
        polypeptide = test_pdb[0]
        self.assertEqual(len(isambard.ampal.interactions.find_Met_pi_interactions(polypeptide,acceptor_codes=['PHE','TYR','TRP'])),3)

class CationPiInteractionFinderTestCase(unittest.TestCase):
    """Test Cation-Pi finding mechanism"""

    def test_1ek9(self):
        test_path = os.path.join(os.path.dirname(isambard.__file__), 'unit_tests','testing_files','1ek9.pdb')
        test_pdb = isambard.ampal.convert_pdb_to_ampal(test_path)
        self.assertEqual(len(isambard.ampal.interactions.find_cation_pi_interactions(test_pdb)),4)

    def test_1gui(self):
        test_path = os.path.join(os.path.dirname(isambard.__file__),'unit_tests','testing_files','2ht0.pdb')
        test_pdb = isambard.ampal.convert_pdb_to_ampal(test_path)
        self.assertEqual(len(isambard.ampal.interactions.find_cation_pi_interactions(test_pdb)),1)

class PiPiInteractionFinderTestCase(unittest.TestCase):
    """Test PiPi interaction finding mechanism"""

    def test_3qy1(self):
        test_path=os.path.join(os.path.dirname(isambard.__file__), 'unit_tests','testing_files','3qy1.pdb')
        test_pdb = isambard.ampal.convert_pdb_to_ampal(test_path)
        face,edge = isambard.ampal.interactions.find_pi_pi_interactions(test_pdb)
        self.assertEqual(len(edge),2)

class C5HydrogenBondFinderTestCase(unittest.TestCase):

    def test_3qy1(self):
        test_path = os.path.join(os.path.dirname(isambard.__file__), 'unit_tests', 'testing_files', '3qy1.pdb')
        test_pdb = isambard.ampal.convert_pdb_to_ampal(test_path)
        c5hb = isambard.ampal.interactions.find_C5HydrogenBonds(test_pdb)
        self.assertEqual(len(c5hb),10)

    def test_1ek9(self):
        test_path = os.path.join(os.path.dirname(isambard.__file__), 'unit_tests', 'testing_files', '1ek9.pdb')
        test_pdb = isambard.ampal.convert_pdb_to_ampal(test_path)
        c5hb = isambard.ampal.interactions.find_C5HydrogenBonds(test_pdb)
        self.assertEqual(len(c5hb),44)
        
if __name__ == '__main__':
    unittest.main()
