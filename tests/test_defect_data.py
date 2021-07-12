import unittest
from pyscses.defect_data import DefectData

class TestDefectDataInit(unittest.TestCase):

    def test_DefectData_is_initialised(self):
        label = 'X'
        energy = -1.0
        defect_data = DefectData(label=label,
                                 energy=energy)
        self.assertEqual(defect_data.label, label)
        self.assertEqual(defect_data.energy, energy)

class TestDefectData(unittest.TestCase):

    def test_eq_returns_true(self):
        label = 'X'
        energy = -1.0
        defect_data_1 = DefectData(label=label,
                                   energy=energy)
        defect_data_2 = DefectData(label=label,
                                   energy=energy)
        self.assertTrue(defect_data_1 == defect_data_2)

    def test_eq_returns_false_on_unequal_labels(self):
        label1 = 'X'
        energy = -1.0
        label2 = 'Y'
        defect_data_1 = DefectData(label=label1,
                                   energy=energy)
        defect_data_2 = DefectData(label=label2,
                                   energy=energy)
        self.assertFalse(defect_data_1 == defect_data_2)

    def test_eq_returns_false_on_unequal_energies(self):
        label = 'X'
        energy1 = -1.0
        energy2 = +1.0
        defect_data_1 = DefectData(label=label,
                                   energy=energy1)
        defect_data_2 = DefectData(label=label,
                                   energy=energy2)
        self.assertFalse(defect_data_1 == defect_data_2)

    def test_eq_evalutates_as_False_for_invalid_comparison(self):
        label = 'X'
        energy = -1.0
        defect_data = DefectData(label=label,
                                 energy=energy)
        self.assertFalse(defect_data == 'foo')
