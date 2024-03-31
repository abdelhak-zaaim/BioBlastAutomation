$(document).ready(async function () {

    $('#sequenceTable').DataTable();


  //  var table = $('#sequenceTable').DataTable({
  //  "paging": false,
  //      "columns": [
  //          {
  //              "className": 'details-control',
  //              "orderable": false,
  //              "data": null,
  //              "defaultContent": ''
  //          },
  //          // Other columns...
  //      ],
  //
  //      "order": [[1, 'asc']]
  //  });

    // Add event listener for opening and closing details
    $('#sequenceTable tbody').on('click', 'td.details-control', function () {
        var tr = $(this).closest('tr');
        var row = table.row(tr);

        if (row.child.isShown()) {
            // This row is already open - close it
            row.child.hide();
            tr.removeClass('shown');
        } else {
            // Open this row
            row.child(format(row.data())).show();
            tr.addClass('shown');
        }
    });
});
$(document).ready(function () {
    $('body').tooltip({
        selector: '[data-toggle="tooltip"]'
    });
});

function format(data) {
    // `data` is the original data object for the row
    return '<table cellpadding="5" cellspacing="0" border="0" style="padding-left:50px;">' +
        '<tr>' +
        '<td>A bit Score:</td>' +
        '<td>' + data[3] + '</td>' +
        '</tr>' +
        '<tr>' +
        '<td>E-value:</td>' +
        '<td>' + data[4] + '</td>' +
        '</tr>' +
        '<tr>' +
        '<td>Identity:</td>' +
        '<td>' + data[5] + '</td>' +
        '</tr>' +
        // Add more rows as needed...
        '</table>';
}

