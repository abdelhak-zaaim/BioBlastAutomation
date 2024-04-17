import ssl
import subprocess

from blast.models.CustomBlastParameters import CustomBlastParameters
from blast.scripts.utils.BlastPlusUtils import BlastPlusUtils


class BlastPlusSearchExtension:
    def __init__(self, custom_blast_parameters: CustomBlastParameters):
        self.custom_blast_parameters = custom_blast_parameters

    def submit_blast_search(self):


        ssl._create_default_https_context = ssl._create_unverified_context

        try:
            fasta_str = ""
            for record in self.custom_blast_parameters.sequences:
                fasta_str += f">{record.id} {record.description}\n{record.seq}\n"

            # use a blast plus version of the search
            result = subprocess.run(
                [self.custom_blast_parameters.program, '-query', '-', '-db', BlastPlusUtils.get_database_path_by_name(self.custom_blast_parameters.database), '-outfmt', '5'],
                input=fasta_str, text=True, stdout=subprocess.PIPE)
            blast_results = result.stdout.decode('utf-8')

        except Exception as e:
            print(f"Error in NCBIWWW.qblast: {e}")
            return None

        print("the blast results are ready ...")

        return blast_results
