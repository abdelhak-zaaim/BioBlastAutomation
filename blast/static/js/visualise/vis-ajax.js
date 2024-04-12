    $('#sequenceTable tbody').on('click', 'tr', function () {
        $('#loadingModal').modal({
            backdrop: 'static',
            keyboard: false
        });

        var subject_title = $(this).attr('data-subject-title');
        let csrftoken = document.querySelector('[name=csrfmiddlewaretoken]').value;

        let url = './alignement_viewer';
        let data = {
            query_seq: $(this).attr('data-query-seq'),
            midline_seq: $(this).attr('data-midline-seq'),
            subject_seq: $(this).attr('data-subject-seq'),
        };
        fetch(url, {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json',
                'X-CSRFToken': csrftoken,
            },
            body: JSON.stringify(data),

        })
            .then(response => response.json()
            )
            .then(data => {

                $('#loadingModal').modal('hide');
                $('#query_seq').html(data.query_seq);
                $('#midline_seq').html(data.midline_seq);
                $('#subject-title').text(subject_title);
                $('#subject_seq').html(data.subject_seq);
                $('#myModal').modal('show');
            })
            .catch((error) => {
                $('#loadingModal').modal('hide');

            });


    });