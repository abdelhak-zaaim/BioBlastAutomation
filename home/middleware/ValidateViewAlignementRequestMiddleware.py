import json
from django.http import JsonResponse

class ValidateViewAlignementRequestMiddleware:
    def __init__(self, get_response):
        self.get_response = get_response

    def __call__(self, request):

        if request.path == '/alignement_viewer':
            if request.method == 'POST':
                data = json.loads(request.body)

                # check the parameters
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

        response = self.get_response(request)

        return response