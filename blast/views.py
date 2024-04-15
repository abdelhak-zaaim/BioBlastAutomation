import json
import os
import ssl

from django.http import HttpResponse
from django.shortcuts import render

from blast.scripts.submission.Submission import Submission
from blast.scripts.utils import BlastUtils
from blast.scripts.visualisation.AlignmentViewer import AlignmentViewer
from blast.scripts.visualisation.SequenceDataVisualize import SequenceData
from pfe import settings

ssl._create_default_https_context = ssl._create_unverified_context


def blast(request):
    return SequenceData.visualise_from_xml_file(request, os.path.join(settings.STATICFILES_DIRS[0], 'test3.xml'))


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
    print(request.data.get('blast_program'))
    print(request.data.get('sequences')[0].seq)
    print(request.data.get('blast_database'))

    submision = Submission(program=request.data.get('blast_program'), output_format="XML")

    xml_file_name = submision.submit_blast_and_save(request.data.get('sequences'), request.data.get('blast_database'))

    return SequenceData.visualise_from_xml_file(request, BlastUtils.get_output_xml_file_path(xml_file_name))
