from celery import Celery
from celery import group
from Bio.Blast import NCBIWWW
from blast.scripts.utils.BlastUtils import BlastUtils

# Initialize Celery
app = Celery('submit_blast_search', broker='submit_blast_search://abdelhak@localhost//')


class BlastSubmissionTask:
    """
        this class is responsible for submitting the blast request to the NCBI
        by separating the request into multiple requests and submitting them in parallel
        also this class is responsible for saving the results of the blast request to an xml file
        we are using cerely to submit the requestes and track the progress of the requests
    """
    def __init__(self, output_format, program, database, sequences):
        self.output_format = output_format
        self.program = program
        self.database = database

    @app.task
    def submit_unique_sequence_and_save(self, sequence):
        """
            this function is responsible for submitting the blast request to the NCBI server and saving the results to an xml file
        :param sequence:
        :return: result of the blast request as a file name of the saved xml file
        """
        print("submitting blast request to NCBI server is started ...")

        sequence = str(sequence)

        result_handle = NCBIWWW.qblast(self.program, self.database, sequence.seq, format_type=self.output_format)
        print("blast request is submitted successfully ...")
        blast_results = result_handle.read()

        return BlastUtils.save_blast_results_to_xml(blast_results)

    @app.task
    def submit_multiple_blast_and_save(self, sequences):
        """
            this class is responsible for submitting the blast request to the NCBI server and saving the results to an xml file
            this fun will return the file name of the saved xml file
        :param sequences: is a list of sequences
        :param database: is a string representing the database name
        :return: an xml file name of the saved xml file
        """
        print("submitting blast request to NCBI server is started ...")

        fasta_str = ""
        for record in sequences:
            fasta_str += f">{record.id} {record.description}\n{record.seq}\n"

        result_handle = NCBIWWW.qblast(self.program, self.database, fasta_str, format_type=self.output_format)
        print("blast requests are submitted successfully ...")
        blast_results = result_handle.read()

        return BlastUtils.save_blast_results_to_xml(blast_results)

    @app.task
    def submit_multiple_blast_parallel_and_save(self, sequences, database):
        """
            this class is responsible for submitting the blast request to the NCBI server and saving the results to an xml file
            this fun will return a list of file names of the saved xml files
        :param sequences:
        :param database:
        :return: a list of file names of the saved xml files
        """
        print("submitting multiple blast requests to NCBI server in parallel ...")

        # Create a group of tasks
        tasks = group(self.submit_unique_sequence_and_save.s(sequence, self.database) for sequence in sequences)

        result = tasks.apply_async()
        print("blast requests are submitted successfully ...")

        results = result.get()

        # Save all results to XML and collect file names
        file_names = []
        for i, blast_results in enumerate(results):
            file_name = BlastUtils.save_blast_results_to_xml(blast_results)
            file_names.append(file_name)

        return file_names
