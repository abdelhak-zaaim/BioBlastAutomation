import ssl

from Bio.Blast import NCBIWWW

from blast.scripts.utils.BlastUtils import BlastUtils
from blast.scripts.utils.Constants import Constants


class Submission:

    def __init__(self, output_format, program):
        self.output_format = output_format
        self.program = program
        # check if program in programs
        if program in Constants.programs:
            self.program = program
        else:
            raise ValueError("Invalid program name , please use a valid program name.")

    def submit_blast_and_save(self, sequences, database):
        """
            this class is responsible for submitting the blast request to the NCBI server and saving the results to an xml file
            this fun will return the file name of the saved xml file
            """
        ssl._create_default_https_context = ssl._create_unverified_context

        result_handle = NCBIWWW.qblast(self.program, database, sequences[0].seq, format_type=self.output_format)

        blast_results = result_handle.read()

        return BlastUtils.save_blast_results_to_xml(blast_results)

    @staticmethod
    def get_databases_for_blast_program(program):
        nucleotide_databases = ["nt", "refseq_rna", "est", "gss", "sts", "pat", "dbsts", "htgs"]
        protein_databases = ["nr", "swissprot", "pdb", "refseq_protein"]
        translated_nucleotide_databases = ["nr", "swissprot", "pdb", "refseq_protein", "nt", "refseq_rna", "est", "gss",
                                           "sts", "pat", "dbsts", "htgs"]

        if program in ["blastn"]:
            return nucleotide_databases
        elif program in ["blastp"]:
            return protein_databases
        elif program in ["blastx", "tblastn", "tblastx"]:
            return translated_nucleotide_databases
        else:
            return []
