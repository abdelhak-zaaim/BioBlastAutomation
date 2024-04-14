# this class is used to view the alignment of the sequences , it is used to color the sequences and the midline
# the midline is the line that shows the similarity between the two sequences
# the color of the midline is different from the color of the sequences
# the color of the sequences is different from the color of the midline


class AlignmentViewer:
    def __init__(self, query_seq, midline_seq, subject_seq):
        self.query_seq = query_seq
        self.midline_seq = midline_seq
        self.subject_seq = subject_seq

    # this function is used to color the midline
    def colour_midline(self, midline):
        return ''.join(['<span class="' + (
            'ml-diff' if chr == ' ' else 'ml-similar' if chr == '+' else 'ml-match') + '">&nbsp;</span>' for chr in
                        midline])

    # this function is used to check if the sequence is a nucleic acid or an amino acid
    def is_nucleic_acid(self, seq):
        nucleic_acids = ['A', 'T', 'C', 'G', 'U']

        return all(char in nucleic_acids for char in seq)

    # this function is used to color the sequences
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

    def view_alignments(self):
        # check if  the nucleic acid or amina acide
        subject_seq_type = 'nucleic_acid' if self.is_nucleic_acid(self.subject_seq) else 'amino_acid'
        # same thing for the query sequence
        query_seq_type = 'nucleic_acid' if self.is_nucleic_acid(self.query_seq) else 'amino_acid'
        # generate the response
        alignment = {'query_seq': '  Query: ' + self.color_seq(self.query_seq, query_seq_type),
                     'midline_seq': '         ' + self.colour_midline(self.midline_seq),
                     'subject_seq': 'Subject: ' + self.color_seq(self.subject_seq, subject_seq_type)}

        return alignment


# this is the main function that is used to test the class
def main():
    viewer = AlignmentViewer('ACTG', '||||', 'ACTG')

    alignment = viewer.view_alignments()

    print(alignment)


if __name__ == "__main__":
    main()
