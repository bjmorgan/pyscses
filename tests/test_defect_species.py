import unittest
from pyscses.defect_species import DefectSpecies

class TestDefectSpecies(unittest.TestCase):

    def test_init(self):
        defect_species = DefectSpecies(label='VO',
                                       valence=+2.0,
                                       mole_fraction=0.1,
                                       mobility=0.1,
                                       can_equilibrate=False)
        self.assertEqual(defect_species.label, 'VO')
        self.assertEqual(defect_species.valence, 2.0)
        self.assertEqual(defect_species.mole_fraction, 0.1)
        self.assertEqual(defect_species.mobility, 0.1)
        self.assertEqual(defect_species.can_equilibrate, False)

    def test_init_defaults(self):
        defect_species = DefectSpecies(label='VO',
                                       valence=+2.0,
                                       mole_fraction=0.1)
        self.assertEqual(defect_species.label, 'VO')
        self.assertEqual(defect_species.valence, 2.0)
        self.assertEqual(defect_species.mole_fraction, 0.1)
        self.assertEqual(defect_species.mobility, 0.0)
        self.assertEqual(defect_species.can_equilibrate, True)

    def test_init_raises_TypeError_if_label_is_incorrect_type(self):
        with self.assertRaises(TypeError):
            defect_species = DefectSpecies(label=3.0,
                                           valence=+2.0,
                                           mole_fraction=0.1)
    def test_init_raises_TypeError_if_valence_is_incorrect_type(self):
        with self.assertRaises(TypeError):
            defect_species = DefectSpecies(label='VO',
                                           valence='foo',
                                           mole_fraction=0.1)

    def test_init_raises_TypeError_if_mole_fraction_is_incorrect_type(self):
        with self.assertRaises(TypeError):
            defect_species = DefectSpecies(label='VO',
                                           valence=+2.0,
                                           mole_fraction='foo')


if __name__ == '__main__':
	unittest.main()
