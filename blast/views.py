import json

from django.shortcuts import render
from django.http import HttpResponse

from blast.scripts.visualisation.AlignmentViewer import AlignmentViewer
from blast.scripts.visualisation.SequenceData import SequenceData
# import visualisation file

from Bio.Blast import NCBIWWW
import ssl

ssl._create_default_https_context = ssl._create_unverified_context



def blast(request):
    return SequenceData.home(request)


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
    # return request.data as a json
    return HttpResponse(json.dumps(request.data))