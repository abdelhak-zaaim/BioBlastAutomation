import ssl

from Bio.Blast import NCBIWWW

from blast.models.CustomBlastParameters import CustomBlastParameters
from blast.scripts.submission.blast.BlastSubmissionTask import BlastSubmissionTask



class Submission:

    def __init__(self, bast_parameters: CustomBlastParameters):

        self.bast_parameters = bast_parameters

    def submit_blast_and_save_async(self, sequences, database):
        # todo : this function need to be implemented and fixed to work with the new structure

        task = BlastSubmissionTask.submit_blast_and_save.delay(self.bast_parameters.program, database, sequences[0].seq,
                                                               self.bast_parameters.output_format)

        # Return the task ID so the status of the task can be checked later
        return task.id




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




