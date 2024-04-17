import json
import os
import ssl

from django.http import HttpResponse
from django.shortcuts import render
from blast.scripts.utils.SecurityUtils import SecurityUtils
from blast.scripts.submission.Submission import Submission
from blast.scripts.submission.blast.BlastSubmissionTask import BlastSubmissionTask
from blast.scripts.utils.BlastUtils import BlastUtils
from blast.scripts.visualisation.AlignmentViewer import AlignmentViewer
from blast.scripts.visualisation.SequenceDataVisualize import SequenceDataVisualize
from pfe import settings

ssl._create_default_https_context = ssl._create_unverified_context


def blast(request):
    return SequenceDataVisualize.visualise_from_xml_file(request,
                                                         os.path.join(settings.STATICFILES_DIRS[0], 'test3.xml'))


def documentation(request):
    return render(request, 'documentation/index.html')


def alignement_viewer(request):
    data = json.loads(request.body)

    alignment_viewer = AlignmentViewer(data.get('query_seq', ""), data.get('midline_seq', ""),
                                       data.get('subject_seq', ""))
    resources = alignment_viewer.view_alignments()

    return HttpResponse(json.dumps(resources))


def submit_sequence(request):
    return render(request, 'submiting/index.html')


def submit_sequence_query(request):

    try:
        custom_blast_parameters = request.custom_blast_parameters
        task_id = BlastSubmissionTask(custom_blast_parameters).submit()
        print("the task id is : ",task_id)
        return HttpResponse(json.dumps({'req_id': str(task_id), 'status': 'submitted'}))
    except Exception as e:
        return HttpResponse(json.dumps({'error': str(e)}), status=400)

def result_viewer(request):
    # aaccess to an parameter from the request get
    req_id= request.GET.get('req_id')
    print("the requested id is : ",req_id)
    # check req_id is not null
    if req_id is not None:

        # check if the file is exist
        if BlastUtils.check_blast_results_file(req_id):
            # check req_id if contain any path traversaal attack

            if SecurityUtils.check_path_traversal(req_id):
                return HttpResponse("the parameter is not valid")
            # if the file exist return the result
            return SequenceDataVisualize.visualise_from_xml_file(request,BlastUtils.get_output_xml_file_path(req_id))
        else:
            # if the file is not exist return error
            return HttpResponse("The request result is not found.")

    else:
        return HttpResponse("please include a req_id parameter.")



