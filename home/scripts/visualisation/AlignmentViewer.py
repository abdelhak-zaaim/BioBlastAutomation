class AlignmentViewer:

    def _colour_midline(self, midline):
        return ''.join(['<span class="' + (
            'ml-diff' if chr == ' ' else 'ml-similar' if chr == '+' else 'ml-match') + '">&nbsp;</span>' for chr in
                        midline])

    def _color_seq(self, seq, seq_type):
        prefix = 'na-' if seq_type == 'nucleic_acid' else 'aa-'

        def _color(letter, position):
            html = '<span'
            if letter != '-':
                html += ' data-idx="' + str(position) + '"'
                html += ' class="' + prefix + letter.lower() + '"'
            else:
                html += ' class="gap"'
            html += '>' + letter + '</span>'
            return html

        letters = list(seq)
        coloured = [_color(letter, idx + 1) for idx, letter in enumerate(letters)]
        return ''.join(coloured)

    def view_alignments(self, query_seq, midline_seq, query_seq_type, subject_seq, query_def, query_id,
                        subject_seq_type, subject_def, subject_id):
        alignments = []

        alignment = {'query_seq': '  Query: ' + self._color_seq(query_seq, query_seq_type),
                     'midline_seq': '         ' + self._colour_midline(midline_seq),
                     'subject_seq': 'Subject: ' + self._color_seq(subject_seq, subject_seq_type)}

        alignments.append(alignment)

        return alignments


def main():
    # Create an instance of the AlignmentViewer class
    viewer = AlignmentViewer()

    # Define the hsps, query_seq_type, query_def, query_id, subject_seq_type, subject_def, and subject_id

    query_seq = 'ACTG',
    midline_seq = '||||',
    subject_seq = 'ACTG'

    query_seq_type = 'nucleic_acid'
    query_def = 'Query Definition'
    query_id = 'Query ID'

    subject_seq_type = 'nucleic_acid'
    subject_def = 'Subject Definition'
    subject_id = 'Subject ID'

    # Call the view_alignments method
    alignments = viewer.view_alignments(query_seq, midline_seq, subject_seq, query_seq_type, query_def, query_id,
                                        subject_seq_type, subject_def, subject_id)

    # Print the alignments
    for alignment in alignments:
        print(alignment)


if __name__ == "__main__":
    main()
