import os

from Bio import SeqIO

from blast.scripts.utils.Constants import Constants
from pfe import settings


class Utils:
    @staticmethod
    def sequence_type(sequence):
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
            return "invalid"

    @staticmethod
    def validate_sequences(query_file):

        sequences = list(SeqIO.parse(query_file, "fasta"))

        for sequence in sequences:
            sequence_type = Utils.sequence_type(sequence.seq)
            if sequence_type == "DNA" or sequence_type == "RNA":
                valid_bases = 'ACGTURYSWKMBDHVNZ'
            elif sequence_type == "Protein":
                valid_bases = 'ACDEFGHIKLMNPQRSTVWY'
            else:
                raise ValueError(f"Invalid sequence format for: {sequence.id}")
                sequence_type = "Invalid sequence format."

            for i, base in enumerate(sequence.seq):
                if base.upper() not in valid_bases:
                    print(f"Invalid residue '{base}' at position {i + 1} in sequence {sequence.id}")

            return sequence_type

    @staticmethod
    def is_valid_output_format(output_format):
        return output_format in Constants.valid_output_formats

    @staticmethod
    def identify_query_type(query):
        nucleotide_bases = set('ATCG')
        amino_acids = set('ACDEFGHIKLMNPQRSTVWY')

        query_set = set(query.upper())

        if query_set.issubset(nucleotide_bases):
            return 'Nucleotide query'
        elif query_set.issubset(amino_acids):
            return 'Protein query'
        else:
            return 'Unknown query type'

    @staticmethod
    def is_program_compatible(program, sequence):
        sequence_type = Utils.sequence_type(sequence)

        if sequence_type == "DNA" and program not in ["blastn", "blastx", "tblastx"]:
            return False
        elif sequence_type == "RNA" and program not in ["blastn", "blastx", "tblastx"]:
            return False
        elif sequence_type == "Protein" and program not in ["blastp", "tblastn", "tblastx"]:
            return False
        else:
            return True

    @staticmethod
    def save_blast_results(blast_results, file_name, output_format):

        static_dir = settings.STATICFILES_DIRS[0]

        # Create the path to the output directory within the static directory
        output_dir = os.path.join(static_dir, 'blast_result')
        os.makedirs(output_dir, exist_ok=True)

        output_file_path = os.path.join(output_dir, f"{file_name}.{output_format}")
        with open(output_file_path, 'w') as output_file:
            if output_format == 'XML':
                output_file.write(blast_results)
            elif output_format == 'JSON':
                import json
                json.dump(blast_results, output_file)
            elif output_format == 'CSV':
                import csv
                writer = csv.writer(output_file)
                writer.writerows(blast_results)
            else:
                raise ValueError("Invalid output format. Please use 'XML', 'JSON', or 'CSV'.")
