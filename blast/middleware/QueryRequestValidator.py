from django.http import JsonResponse

from blast.scripts.utils.Utils import Utils


class QueryRequestValidator:
    def __init__(self, get_response):
        self.get_response = get_response

    def __call__(self, request):
        if request.path == '/submit_sequence_query' and request.method == 'POST':
            sequence_string = request.POST.get('sequence_string')
            from_value = request.POST.get('from')
            to_value = request.POST.get('to')
            file_fasta = request.FILES.get('file_fasta')
            jtitle = request.POST.get('jtitle')
            blast_program = request.POST.get('blast_program')
            blast_database = request.POST.get('blast_database')

            if not (sequence_string or file_fasta):
                return JsonResponse({'error': 'Invalid request'}, status=400)

            sequences = []
            if file_fasta:
                try:
                    file_fasta = file_fasta.read().decode('utf-8')
                    sequences = Utils.validate_fasta_file(file_fasta)
                except Exception as e:
                    return JsonResponse({'error': str(e)}, status=400)
            else:
                try:
                    sequences = Utils.validate_fasta_string(sequence_string)
                except Exception as e:
                    return JsonResponse({'error': str(e)}, status=400)

            request.data = {
                'sequences': sequences,
                'from': from_value,
                'to': to_value,
                'jtitle': jtitle,
                'blast_program': blast_program,
                'blast_database': blast_database
            }

            if not from_value or not to_value or not blast_program or not blast_database:
                return JsonResponse({'error': 'Invalid request'}, status=400)

            if not Utils.is_valid_blast_program(blast_program):
                return JsonResponse({'error': 'Invalid blast program'}, status=400)

        response = self.get_response(request)
        return response
