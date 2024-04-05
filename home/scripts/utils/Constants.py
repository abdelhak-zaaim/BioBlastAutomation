class Constants:
    # *-------------------------------------------- todo: change the paths of the databases to the correct ones ***
    #  you can download it from https://ftp.ncbi.nlm.nih.gov/blast/db/ do not forget to extract the files in case you
    #  are using blast + you can use the command `makeblastdb -in your_file.fasta -dbtype nucl -parse_seqids -out
    #  your_file` *--------------------------------------------
    #11111111

    DATABASES = [{"name": "nt", "path": "/Users/macbook/Desktop/SMI/S6/PFE/blastdatabase/nt.001"},
                 {"name": "nr", "path": "/Users/macbook/Desktop/SMI/S6/PFE/blastdatabase/nr.08"}]

    # specify which database you are using (nt, nr, etc.)
    # NOTE : nt is recommended for nucleotide sequences and nr --
    # for protein sequences, however you can use any database you want

    BLAST_QUERY_FILE = "../../../sequences.fasta"
    # --------------------------------------------
    # todo: change the path to the correct one ***
    # you can download an example file from https://zaaim.me/src/sequences.fasta
    # --------------------------------------------
    TEST_FILE_PATH = "/Users/macbook/PycharmProjects/pfe-blast/sequences.fasta"

    # valid output formats for the BLAST search results (XML, JSON, CSV)
    valid_output_formats = ["XML", "JSON", "CSV"]

    # valid BLAST programs to use (blastn, blastp, blastx, tblastn, tblastx)
    programs = ["blastn", "blastp", "blastx", "tblastn", "tblastx"]


    OUTPUT_FILE_PATH = "/statics/blast_results"
