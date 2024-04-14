from io import StringIO

from Bio import SeqIO


class Utils:
    @staticmethod
    def is_valid_blast_program(blast_program):
        valid_blast_programs = ['blastp', 'blastn', 'blastx', 'tblastn', 'tblastx']
        return blast_program in valid_blast_programs

    @staticmethod
    def is_valid_blast_database(blast_database):
        valid_blast_databases = ['nr', 'nt', 'refseq_protein', 'refseq_rna', 'swissprot', 'pdb', 'env_nr', 'env_nt', 'pat', 'refseq_genomic']
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
