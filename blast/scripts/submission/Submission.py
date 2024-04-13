import datetime
import os
import ssl

from Bio import SeqIO
from Bio.Blast import NCBIWWW

from blast.database.DatabaseManager import DatabaseManager
from blast.scripts.submission.utils.Utils import Utils
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

    def submit_blast_www(self, sequence):
        ssl._create_default_https_context = ssl._create_unverified_context

        if not Utils.is_program_compatible(self.program, sequence):
            raise ValueError("Incompatible program for the given sequence type.")

        sequence_type = Utils.sequence_type(sequence)
        # check if the sequence is valid
        if sequence_type == "invalid":
            raise ValueError("Invalid sequence format. Please use a valid sequence format.")

        # The database to search against
        database = "nt" if sequence_type == "DNA" or sequence_type == "RNA" else "swissprot"

        # Perform the BLAST search and specify the output format

        result_handle = NCBIWWW.qblast(self.program, database, sequence.seq, format_type=self.output_format)
        db_manager = DatabaseManager()
        db_manager.save_blast_results_to_db(result_handle, database)

        blast_results = result_handle.read()

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

    @staticmethod
    def get_database_description(database):
        descriptions = {
            "nr": "Non-redundant protein sequences",
            "nt": "Nucleotide collection (nt)",
            "swissprot": "Non-redundant sequences from Swiss-Prot",
            "pdb": "Protein Data Bank (3D)",
            "refseq_rna": "NCBI RefSeq RNA",
            "refseq_protein": "NCBI RefSeq protein",
            "est": "GenBank Expressed Sequence Tags",
            "gss": "Genomic Survey Sequences",
            "sts": "Sequence Tagged Sites",
            "pat": "Patents",
            "dbsts": "dbSTS",
            "htgs": "High Throughput Genomic Sequences"
        }

        return descriptions.get(database, "No description available")


# script to submit a sequence to the BLAST server for testing purposes
if __name__ == "__main__":
    soumission = Submission(output_format='XML', program="blastp")
    sequences = list(SeqIO.parse(os.path.join(settings.STATICFILES_DIRS[0], 'sequences.fasta'), "fasta"))
    for sequence in sequences:
        print("submitting sequence: ", sequence.id)
        result = soumission.submit_blast_www(sequence)

        # get current time to use it as a unique identifier for the output file
        time = datetime.datetime.now().strftime("%Y%m%d%H%M%S")

        Utils.save_blast_results(result, f"result_{sequence.id}_{time}", soumission.output_format)
        print("done")
