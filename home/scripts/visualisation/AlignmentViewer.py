class AlignmentViewer:
    def __init__(self, query_seq, midline_seq, subject_seq):
        self.query_seq = query_seq
        self.midline_seq = midline_seq
        self.subject_seq = subject_seq

    def colour_midline(self, midline):
        return ''.join(['<span class="' + (
            'ml-diff' if chr == ' ' else 'ml-similar' if chr == '+' else 'ml-match') + '">&nbsp;</span>' for chr in
                        midline])

    def is_nucleic_acid(self, seq):
        nucleic_acids = ['A', 'T', 'C', 'G', 'U']

        return all(char in nucleic_acids for char in seq)

    def color_seq(self, seq, seq_type):
        prefix = 'na-' if seq_type == 'nucleic_acid' else 'aa-'

        def color(letter, position):
            html = '<span'
            if letter != '-':
                html += ' data-idx="' + str(position) + '"'
                html += ' class="' + prefix + letter.lower() + '"'
            else:
                html += ' class="gap"'
            html += '>' + letter + '</span>'
            return html

        letters = list(seq)
        coloured = [color(letter, idx + 1) for idx, letter in enumerate(letters)]
        return ''.join(coloured)

    def view_alignments(self,  query_def, query_id,
                        subject_def, subject_id):

        subject_seq_type = 'nucleic_acid' if self.is_nucleic_acid(self.subject_seq) else 'amino_acid'
        query_seq_type = 'nucleic_acid' if self.is_nucleic_acid(self.query_seq) else 'amino_acid'

        alignment = {'query_seq': '  Query: ' + self.color_seq(self.query_seq, query_seq_type),
                     'midline_seq': '         ' + self.colour_midline(self.midline_seq),
                     'subject_seq': 'Subject: ' + self.color_seq(self.subject_seq, subject_seq_type)}

        return alignment


def main():
    viewer = AlignmentViewer('ACTG', '||||', 'ACTG')



    query_def = 'Query Definition'
    query_id = 'Query ID'

    subject_def = 'Subject Definition'
    subject_id = 'Subject ID'

    # Call the view_alignments method
    alignment = viewer.view_alignments(  query_def, query_id, subject_def,
                                        subject_id)



    print(alignment)


if __name__ == "__main__":
    main()
