
import os
import unittest
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from scripts.soumission.submitting import Soumission
from home.scripts.utils.Constants import Constants


class TestSoumission(unittest.TestCase):
    def setUp(self):
        self.soumission = Soumission(Constants.BLAST_DB_PATH, 'XML')
        self.input_file_path = Constants.TEST_FILE_PATH

    def test_sequence_type(self):
        self.assertEqual(self.soumission.sequence_type('ACDEFGHIKLMNPQRSTVWY'), 'Protein')
        self.assertEqual(self.soumission.sequence_type('ACGT'), 'DNA')
        self.assertEqual(self.soumission.sequence_type('ACGU'), 'RNA')
        self.assertEqual(self.soumission.sequence_type('INV545ALID'), 'Invalid sequence format.')
        # you can add more testes here to test the function with different inputs

    def test_validate_sequences(self):
        try:
            self.soumission.validate_sequences(self.input_file_path)
        except ValueError:
            self.fail("validate_sequences raised ValueError unexpectedly!")
    def test_run_blast(self):
        sequence = SeqRecord(Seq("ACTG"), id="test")
        self.soumission.run_blast(sequence, 'test_output')
        self.assertTrue(os.path.exists('test_output.XML'))

    def test_parallel_blast(self):
        self.soumission.parallel_blast(self.input_file_path)

        self.assertTrue(os.path.exists('result_0.XML'))

if __name__ == '__main__':
    unittest.main()