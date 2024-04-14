import unittest
from blast.scripts.utils.Utils import Utils
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from io import StringIO

class TestUtils(unittest.TestCase):
    def setUp(self):
        self.utils = Utils()

    def test_is_valid_blast_program_with_valid_program(self):
        self.assertTrue(self.utils.is_valid_blast_program('blastp'))

    def test_is_valid_blast_program_with_invalid_program(self):
        self.assertFalse(self.utils.is_valid_blast_program('invalid'))

    def test_is_valid_blast_database_with_valid_database(self):
        self.assertTrue(self.utils.is_valid_blast_database('nr'))

    def test_is_valid_blast_database_with_invalid_database(self):
        self.assertFalse(self.utils.is_valid_blast_database('invalid'))

    def test_is_program_database_compatible_with_compatible_program_and_database(self):
        self.assertTrue(self.utils.is_program_database_compatible('blastp', 'nr'))

    def test_is_program_database_compatible_with_incompatible_program_and_database(self):
        self.assertFalse(self.utils.is_program_database_compatible('blastp', 'nt'))

    def test_get_databases_for_program_with_valid_program(self):
        self.assertEqual(self.utils.get_databases_for_program('blastp'), ['nr', 'refseq_protein', 'swissprot', 'pdb'])

    def test_get_databases_for_program_with_invalid_program(self):
        self.assertEqual(self.utils.get_databases_for_program('invalid'), [])

    def test_validate_and_parse_fasta_file_with_valid_fasta_file(self):
        fasta_file = StringIO(">seq1\nACGT")
        expected_result = [{"id": "seq1", "sequence": "ACGT", "type": "dna", "title": "seq1"}]
        self.assertEqual(self.utils.validate_and_parse_fasta_file(fasta_file), expected_result)

    def test_validate_and_parse_fasta_file_with_invalid_fasta_file(self):
        fasta_file = StringIO(">seq1\nXYZ")
        expected_result = [{"id": "seq1", "sequence": "XYZ", "type": "unknown", "title": "seq1"}]
        self.assertEqual(self.utils.validate_and_parse_fasta_file(fasta_file), expected_result)

    def test_validate_fasta_string_with_valid_fasta_string(self):
        fasta_string = ">seq1\nACGT"
        expected_result = [{"id": "seq1", "sequence": "ACGT", "type": "dna", "title": "seq1"}]
        self.assertEqual(self.utils.validate_fasta_string(fasta_string), expected_result)

    def test_validate_fasta_string_with_invalid_fasta_string(self):
        fasta_string = ">seq1\nXYZ"
        expected_result = [{"id": "seq1", "sequence": "XYZ", "type": "unknown", "title": "seq1"}]
        self.assertEqual(self.utils.validate_fasta_string(fasta_string), expected_result)

if __name__ == '__main__':
    unittest.main()