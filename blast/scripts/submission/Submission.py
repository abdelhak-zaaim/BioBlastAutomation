import datetime
import os
import ssl

from Bio import SeqIO
from Bio.Blast import NCBIWWW

from blast.database.DatabaseManager import DatabaseManager

from blast.scripts.utils import BlastUtils
from blast.scripts.utils.Constants import Constants
from pfe import settings


class Submission:

    def __init__(self, output_format, program):
        self.output_format = output_format
        self.program = program
        # check if program in programs
        if program in Constants.programs:
            self.program = program
        else:
            raise ValueError("Invalid program name , please use a valid program name.")

    def submit_blast_www(self, sequences, database):
        ssl._create_default_https_context = ssl._create_unverified_context

        result_handle = NCBIWWW.qblast(self.program, database, sequences[0].seq, format_type=self.output_format)

        blast_results = result_handle.read()
        BlastUtils.save_blast_results_to_xml(blast_results)
        return blast_results

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


# script to submit a sequence to the BLAST server for testing purposes
if __name__ == "__main__":
    soumission = Submission(output_format='XML', program="blastp")
    sequences = list(SeqIO.parse(os.path.join(settings.STATICFILES_DIRS[0], 'sequences.fasta'), "fasta"))
    for sequence in sequences:
        print("submitting sequence: ", sequence.id)
        result = soumission.submit_blast_www(sequence)

        # get current time to use it as a unique identifier for the output file
        time = datetime.datetime.now().strftime("%Y%m%d%H%M%S")

        BlastUtils.save_blast_results(result, f"result_{sequence.id}_{time}", soumission.output_format)
        print("done")
