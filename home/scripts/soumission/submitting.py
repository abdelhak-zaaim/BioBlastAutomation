from concurrent.futures import ProcessPoolExecutor
from Bio import SeqIO
import subprocess

from utils.constants import Constants


class Soumission:
    def __init__(self, db_path, output_format):
        self.db_path = db_path
        self.output_format = output_format

    def sequence_type(self, sequence):
        dna_bases = 'ACGT'
        rna_bases = 'ACGU'
        amino_acids = 'ACDEFGHIKLMNPQRSTVWY'

        if all(base in dna_bases for base in sequence):
            return "DNA"
        elif all(base in rna_bases for base in sequence):
            return "RNA"
        elif all(base in amino_acids for base in sequence):
            return "Protein"
        else:
            return "Invalid sequence format."

    def validate_sequences(self, query_file):
        # Read the sequences from the .fasta file
        sequences = list(SeqIO.parse(query_file, "fasta"))

        for sequence in sequences:
            sequence_type = self.sequence_type(sequence.seq)
            if sequence_type == "DNA" or sequence_type == "RNA":
                valid_bases = 'ACGTURYSWKMBDHVNZ'
            elif sequence_type == "Protein":
                valid_bases = 'ACDEFGHIKLMNPQRSTVWY'
            else:
                raise ValueError(f"Invalid sequence format for: {sequence.id}")

            for i, base in enumerate(sequence.seq):
                if base.upper() not in valid_bases:
                    print(f"Invalid residue '{base}' at position {i + 1} in sequence {sequence.id}")

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

    def parallel_blast(self, query_file):
        self.validate_sequences(query_file)
        # Read the sequences from the .fasta file
        sequences = list(SeqIO.parse(query_file, "fasta"))
        try:
            # Create a pool of worker processes
            with ProcessPoolExecutor() as executor:
                # Perform the BLAST search on all sequences in parallel

                executor.map(self.sequence_type, sequences)

                executor.map(self.run_blast, sequences, [f"result_{i}" for i in range(len(sequences))])
        finally:
            subprocess.run(['rm', 'temp.fasta'], check=True)



if __name__ == "__main__":
    blast_search = Soumission(Constants.BLAST_DB_PATH, 'XML')
    try:
        blast_search.validate_sequences(Constants.BLAST_QUERY_FILE)
        blast_search.parallel_blast(Constants.BLAST_QUERY_FILE)
    except ValueError as e:
        print(f"Error: {e}")
