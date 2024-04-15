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


            # verify query_from and query_to are integers
            if from_value:
                try:
                    from_value = int(from_value)
                except ValueError:
                    return JsonResponse({'error': 'Invalid from value'}, status=400)

            if to_value:
                try:
                    to_value = int(to_value)
                except ValueError:
                    return JsonResponse({'error': 'Invalid to value'}, status=400)

            # verify from_value is less than to_value and both are positive
            if from_value and to_value:
                if from_value < 0 or to_value < 0:
                    return JsonResponse({'error': 'from and to values must be positive'}, status=400)
                if from_value > to_value:
                    return JsonResponse({'error': 'from value must be less than to value'}, status=400)


            if (from_value and not to_value) or (to_value and not from_value):
                return JsonResponse({'error': 'from and to values must be both exist or both not exist'}, status=400)

            if not (sequence_string or file_fasta):
                return JsonResponse({'error': 'Invalid request'}, status=400)

            if not blast_program or not blast_database:
                return JsonResponse({'error': 'Invalid request'}, status=400)

            if not Utils.is_valid_blast_program(blast_program):
                return JsonResponse({'error': 'Invalid blast program'}, status=400)

            if not Utils.is_valid_blast_database(blast_database):
                return JsonResponse({'error': 'Invalid blast database'}, status=400)

            if not Utils.is_program_database_compatible(blast_program, blast_database):
                return JsonResponse({'error': 'Incompatible blast program and database'}, status=400)

            sequences = []
            if file_fasta:
                try:
                    fasta_string = file_fasta.read().decode('utf-8')
                    sequences_file = Utils.get_sequences_from_fast_string(fasta_string)
                    sequences.extend(sequences_file)
                except Exception as e:
                    return JsonResponse({'error': 'Error reading file: ' + str(e)}, status=400)
            if sequence_string:
                try:
                    sequences_string = Utils.get_sequences_from_fast_string(sequence_string)
                    sequences.extend(sequences_string)
                except Exception as e:
                    return JsonResponse({'error': 'Error validating FASTA string: ' + str(e)}, status=400)

            if not sequences:
                return JsonResponse({'error': 'No valid sequences found'}, status=400)


            # Check if the sequences are compatible with the blast program
            for sequence in sequences:
                if not Utils.is_sequence_type_compatible_with_program(Utils.get_sequence_type(sequence), blast_program):
                    return JsonResponse(
                        {'error': 'Incompatible sequence type with blast program , sequence: ' + sequence.id},
                        status=400)

            request.data = {
                'sequences': sequences,
                'from': from_value,
                'to': to_value,
                'jtitle': jtitle,
                'blast_program': blast_program,
                'blast_database': blast_database
            }

        response = self.get_response(request)
        return response
