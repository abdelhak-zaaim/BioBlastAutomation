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
        sequence_map = {}
        for record in SeqIO.parse(file_fasta, "fasta"):
            sequence = str(record.seq)
            if set(sequence.upper()).issubset('ACGT'):
                sequence_map[record.id] = {"sequence": sequence, "type": "dna"}
            elif set(sequence.upper()).issubset('ACGU'):
                sequence_map[record.id] = {"sequence": sequence, "type": "rna"}
            elif set(sequence.upper()).issubset('ACDEFGHIKLMNPQRSTVWY'):
                sequence_map[record.id] = {"sequence": sequence, "type": "amino_acid"}
            else:
                sequence_map[record.id] = {"sequence": sequence, "type": "unknown"}
        return sequence_map
