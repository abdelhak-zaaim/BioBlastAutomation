import ssl

from Bio.Blast import NCBIWWW

from blast.scripts.utils.BlastUtils import BlastUtils


class BlastSearchExtension:
    def __init__(self, custom_blast_parameters):
        self.custom_blast_parameters = custom_blast_parameters

    def submit_blast_search(self):
        print("the submission of the blast request is started ...  sssss")
        print(self.custom_blast_parameters.to_dict())
        ssl._create_default_https_context = ssl._create_unverified_context


        try:
            fasta_str = ""
            for record in self.custom_blast_parameters.sequences:
                fasta_str += f">{record.id} {record.description}\n{record.seq}\n"

            result_handle = NCBIWWW.qblast(self.custom_blast_parameters.program,
                                           self.custom_blast_parameters.database,
                                           fasta_str,
                                           format_type=self.custom_blast_parameters.output_format)
        except Exception as e:
            print(f"Error in NCBIWWW.qblast: {e}")
            return None



        print("the blast results are ready ...")

        blast_results = result_handle.read()

        return blast_results