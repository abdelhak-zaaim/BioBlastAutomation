import datetime
import os
import ssl
from Bio import SeqIO
import subprocess
from Bio.Blast import NCBIWWW
from home.scripts.utils.Constants import Constants
from pfe import settings
from home.scripts.submission.utils.Utils import Utils


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

        # Read the results
        blast_results = result_handle.read()

        return blast_results

    def run_blast(self, sequence, output_file):
        # get file format number based on output_format
        if self.output_format == 'XML':
            outfmt = '5'
        elif self.output_format == 'JSON':
            outfmt = '15'
        elif self.output_format == 'CSV':
            outfmt = '10'

        else:
            raise ValueError("Invalid output format. Please use 'XML', 'JSON', or 'CSV'.")

        # Write the sequence to a temporary file
        with open("temp.fasta", "w") as temp_file:
            SeqIO.write(sequence, temp_file, "fasta")
        print("we strated the blast")
        command = [
            'blastn',
            '-db', self.db_path,
            '-query', "temp.fasta",
            '-out', output_file + '.' + self.output_format,
            '-outfmt', outfmt
        ]

        subprocess.run(command, check=True)


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
