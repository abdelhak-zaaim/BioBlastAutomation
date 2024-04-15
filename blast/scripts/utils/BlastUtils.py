import os
import random
import string
from io import StringIO

from Bio import SeqIO

from blast.scripts.utils.Constants import Constants
from pfe import settings


class BlastUtils:
    @staticmethod
    def is_valid_blast_program(blast_program):
        valid_blast_programs = ['blastp', 'blastn', 'blastx', 'tblastn', 'tblastx']
        return blast_program in valid_blast_programs

    @staticmethod
    def is_valid_blast_database(blast_database):
        valid_blast_databases = ['nr', 'nt', 'refseq_protein', 'refseq_rna', 'swissprot', 'pdb', 'env_nr', 'env_nt',
                                 'pat', 'refseq_genomic']
        return blast_database in valid_blast_databases

    @staticmethod
    def is_program_database_compatible(blast_program, blast_database):
        compatibility_map = {
            'blastp': ['nr', 'refseq_protein', 'swissprot', 'pdb'],
            'blastn': ['nt', 'refseq_rna', 'env_nt', 'refseq_genomic'],
            'blastx': ['nr', 'refseq_protein', 'swissprot', 'pdb'],
            'tblastn': ['nt', 'refseq_rna', 'env_nt', 'refseq_genomic'],
            'tblastx': ['nt', 'refseq_rna', 'env_nt', 'refseq_genomic']
        }
        return blast_database in compatibility_map.get(blast_program, [])

    @staticmethod
    def get_databases_for_program(blast_program):
        compatibility_map = {
            'blastp': ['nr', 'refseq_protein', 'swissprot', 'pdb'],
            'blastn': ['nt', 'refseq_rna', 'env_nt', 'refseq_genomic'],
            'blastx': ['nr', 'refseq_protein', 'swissprot', 'pdb'],
            'tblastn': ['nt', 'refseq_rna', 'env_nt', 'refseq_genomic'],
            'tblastx': ['nt', 'refseq_rna', 'env_nt', 'refseq_genomic']
        }
        return compatibility_map.get(blast_program, [])

    @staticmethod
    def validate_and_parse_fasta_file(file_fasta):
        sequence_list = []
        for record in SeqIO.parse(file_fasta, "fasta"):
            sequence = str(record.seq)
            title = record.description
            sequence_type = ""
            if set(sequence.upper()).issubset('ACGT'):
                sequence_type = "dna"
            elif set(sequence.upper()).issubset('ACGU'):
                sequence_type = "rna"
            elif set(sequence.upper()).issubset('ACDEFGHIKLMNPQRSTVWY'):
                sequence_type = "amino_acid"
            else:
                sequence_type = "unknown"
            sequence_list.append({"id": record.id, "sequence": sequence, "type": sequence_type, "title": title})
        return sequence_list

    @staticmethod
    def validate_fasta_string(fasta_string):
        sequence_list = []
        for record in SeqIO.parse(StringIO(fasta_string), "fasta"):
            sequence = str(record.seq)
            title = record.description
            sequence_type = ""
            if set(sequence.upper()).issubset('ACGT'):
                sequence_type = "dna"
            elif set(sequence.upper()).issubset('ACGU'):
                sequence_type = "rna"
            elif set(sequence.upper()).issubset('ACDEFGHIKLMNPQRSTVWY'):
                sequence_type = "amino_acid"
            else:
                sequence_type = "unknown"
            sequence_list.append({"id": record.id, "sequence": sequence, "type": sequence_type, "title": title})
        return sequence_list

    @staticmethod
    def is_sequence_type_compatible_with_program(sequence_type, blast_program):
        compatibility_map = {
            'blastp': ['amino_acid'],
            'blastn': ['dna'],
            'blastx': ['dna'],
            'tblastn': ['amino_acid'],
            'tblastx': ['dna']
        }
        return sequence_type in compatibility_map.get(blast_program, [])

    @staticmethod
    def are_sequences_valid(sequences):
        for sequence in sequences:
            if sequence['type'] == 'unknown':
                return False
        return True

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

    @staticmethod
    def get_sequences_from_fast_string(text_fasta):
        sequence_list = []

        for record in SeqIO.parse(StringIO(text_fasta), "fasta"):
            sequence_list.append(record)
        return sequence_list

    @staticmethod
    def get_sequence_type(sequence):
        if set(sequence.upper()).issubset('ACGT'):
            return "dna"
        elif set(sequence.upper()).issubset('ACGU'):
            return "rna"
        elif set(sequence.upper()).issubset('ACDEFGHIKLMNPQRSTVWY'):
            return "amino_acid"
        else:
            return "unknown"

    @staticmethod
    def save_blast_results_to_xml(blast_results):
        # Generate a unique file name using a random string of length 14
        random_string = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(14))

        file_name = f"BLAST_RESULT_{random_string}.xml"
        path = os.path.join(settings.STATIC_BLAST_RESULTS, file_name)
        with open(path, 'w') as file:
            file.write(blast_results)

        return file_name

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
            sequence_type = BlastUtils.sequence_type(sequence.seq)
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
        sequence_type = BlastUtils.sequence_type(sequence)

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

    @staticmethod
    def get_output_xml_file_path(file_name):
        output_dir = settings.STATIC_BLAST_RESULTS

        return os.path.join(output_dir, file_name)
