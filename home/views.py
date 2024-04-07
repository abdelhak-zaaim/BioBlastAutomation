import json

from django.shortcuts import render
from django.http import HttpResponse

from home.scripts.visualisation.AlignmentViewer import AlignmentViewer
# import visualisation file
from home.scripts.visualisation.Visualisation import Visualisation

from Bio.Blast import NCBIWWW
import ssl

ssl._create_default_https_context = ssl._create_unverified_context



def blast(request):
    return Visualisation.home(request)


def documentation(request):
    return render(request, 'documentation/index.html')


def alignement_viewer(request):
    data = json.loads(request.body)

    alignment_viewer = AlignmentViewer(data.get('query_seq', ""), data.get('midline_seq', ""),
                                       data.get('subject_seq', ""))
    resources = alignment_viewer.view_alignments()

    return HttpResponse(json.dumps(resources))
