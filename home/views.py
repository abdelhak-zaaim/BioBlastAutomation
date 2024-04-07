import json

from django.shortcuts import render
from django.http import HttpResponse

from home.scripts.visualisation.AlignmentViewer import AlignmentViewer
# import visualisation file
from home.scripts.visualisation.Visualisation import Visualisation

from Bio.Blast import NCBIWWW
import ssl

ssl._create_default_https_context = ssl._create_unverified_context


def blast_search(sequence, blast_program='blastn', database='nr', hitlist_size=10):
    try:
        result_handle = NCBIWWW.qblast(program=blast_program,
                                       database=database,
                                       format_type='HTML',
                                       sequence=sequence,
                                       hitlist_size=hitlist_size)
        return result_handle
    except Exception as e:
        print(f"An error occurred during BLAST search: {e}")
        return None


# Create your views here.
def home(request):
    sequence_to_search = "ACTGATCGATCGATCGATCGATCGATCGATCG"
    blast_result_handle = blast_search(sequence_to_search)
    if blast_result_handle:
        blast_result = blast_result_handle.read()
        return HttpResponse(blast_result)
    else:
        return HttpResponse("An error occurred during BLAST search.")


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
