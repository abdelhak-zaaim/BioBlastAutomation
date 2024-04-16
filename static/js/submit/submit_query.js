$(document).ready(function () {
    $('form').on('submit', function (event) {


        $('#loadingModal').modal({
            backdrop: 'static',
            keyboard: false
        });

        event.preventDefault();  // Prevent the form from being submitted normally

        var formData = new FormData(this);

        $.ajax({
            url: '/submit_sequence_query',
            type: 'POST',
            data: formData,
            success: function (response) {
                $('#loadingModal').modal('hide');
                // get req_id from response and redirect to result_viewer page the response is a jon format
                var data = JSON.parse(response);
                window.location.href = "/result_viewer?req_id=" + data.req_id;
            },
            error: function (jqXHR, textStatus, errorThrown) {
                // Handle any errors
                $('#loadingModal').modal('hide');
                showAlert({
                    type: "error",
                    // the message should be from the response from the server
                    message: jqXHR.responseJSON.error
                })
            },
            cache: false,
            contentType: false,
            processData: false
        });
    });
});