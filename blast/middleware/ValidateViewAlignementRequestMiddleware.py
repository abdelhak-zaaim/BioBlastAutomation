import json

from django.http import JsonResponse


class ValidateViewAlignementRequestMiddleware:
    def __init__(self, get_response):
        self.get_response = get_response

    def __call__(self, request):

        # this middleware is used to validate the request of the alignment viewer
        # the request must be a post request and must contain the query_seq, midline_seq and subject_seq
        # the sequences must be strings
        # the sequences must be in uppercase

        if request.path == '/alignement_viewer':
            if request.method == 'POST':
                data = json.loads(request.body)

                # check the parameters of the request
                if 'query_seq' not in data:
                    return JsonResponse({"error": "The 'query_seq' parameter is required."}, status=400)
                if 'midline_seq' not in data:
                    return JsonResponse({"error": "The 'midline_seq' parameter is required."}, status=400)
                if 'subject_seq' not in data:
                    return JsonResponse({"error": "The 'subject_seq' parameter is required."}, status=400)
                # check the type of the parameters
                if not isinstance(data['query_seq'], str):
                    return JsonResponse({"error": "The 'query_seq' parameter must be a string."}, status=400)
                if not isinstance(data['midline_seq'], str):
                    return JsonResponse({"error": "The 'midline_seq' parameter must be a string."}, status=400)
                if not isinstance(data['subject_seq'], str):
                    return JsonResponse({"error": "The 'subject_seq' parameter must be a string."}, status=400)

            else:
                return JsonResponse({"error": "Invalid request method. Only POST requests are allowed."}, status=405)
            # in case of success, we will convert the sequences to uppercase and add them to the request object
            data['query_seq'] = data['query_seq'].upper()
            data['midline_seq'] = data['midline_seq'].upper()
            data['subject_seq'] = data['subject_seq'].upper()

            request.data = data
        response = self.get_response(request)

        return response
