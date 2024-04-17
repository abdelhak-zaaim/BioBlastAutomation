import queue
import threading
import uuid

from blast.models.CustomBlastParameters import CustomBlastParameters
from blast.scripts.submission.blast.BlastSearchExtension import BlastSearchExtension
from blast.scripts.utils.BlastUtils import BlastUtils
from blast.scripts.utils.Constants import Constants


# Initialize Celery


class BlastSubmissionTask:
    """
        this class is responsible for submitting the blast request to the NCBI
        by separating the request into multiple requests and submitting them in parallel
        also this class is responsible for saving the results of the blast request to an xml file
        we are using cerely to submit the requestes and track the progress of the requests
    """

    def __init__(self, custom_blast_parameters):
        self.custom_blast_parameters = custom_blast_parameters
        self.task_id = str(uuid.uuid4())

    def submit(self):

        """
        this function should supmit a blast request to the NCBI server and return the task id
        its using the BlastSubmissionTask class to submit the request
        """
        return self.submit_unique()
        # todo : this code
        if custom_blast_parameters.is_async and len(custom_blast_parameters.sequences) > 1:
            return BlastSubmissionTask.run_submit_async(custom_blast_parameters)
        else:
            return BlastSubmissionTask.run_submit_unique(custom_blast_parameters)

    def submit_unique(self):

        """
            this class is responsible for submitting the blast request to the NCBI server and saving the results to an xml file
            this fun will return the file name of the saved xml file

        """

        print("submitting blast request to NCBI server is started ...")

        blast_results = BlastSearchExtension(self.custom_blast_parameters).submit_blast_search()
        print("the blast results are ready ...")

        file_name = BlastUtils.write_blast_results_to_file(blast_results)


        return file_name

    def submit_async(self):
        """
            this class is responsible for submitting the blast request to the NCBI server and saving the results to an xml file
            this fun will return a list of task IDs
        :param sequences:
        :param database:
        :return: a list of task IDs
        """

        print("submitting multiple blast requests to NCBI server in parallel ...")

        # todo : this is a dummy implementation, we need to implement the real implementation
        return "task_ids"

    def run_submit_unique(self):
        try:
            result_queue = queue.Queue()
            task = threading.Thread(target=self.submit_unique, args=(result_queue,))
            task.start()
            return result_queue.get()
        except Exception as e:
            print("error in submitting the task ... :" + str(e))
            print(e)

        @staticmethod
        def is_task_alive(task_id, tasks):
            task = tasks.get(task_id)
            if task is None:
                return False
            return task.is_alive()

    @staticmethod
    def run_submit_async(custom_blast_parameters: CustomBlastParameters):

        task = BlastSubmissionTask.submit_async(custom_blast_parameters)

        return task


# a main function to test the code
if __name__ == '__main__':
    # test the code
    custom_blast_parameters = CustomBlastParameters(program="blastn", database="nt", sequences=[],
                                                    output_format="XML", is_async=True)
    BlastSubmissionTask.submit(custom_blast_parameters)
