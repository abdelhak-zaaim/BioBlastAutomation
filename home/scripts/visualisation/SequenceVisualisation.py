from Bio import pairwise2


class SequenceVisualisation:
    @staticmethod
    def get_color(score):
        if score < 40:
            return 'black'
        elif 40 <= score < 50:
            return 'blue'
        elif 50 <= score < 80:
            return 'green'
        elif 80 <= score < 200:
            return 'orange'
        else:
            return 'red'


    @staticmethod
    def perform_global_alignment(seq1, seq2):
        # Perform the global alignment
        alignments = pairwise2.align.globalxx(seq1, seq2)

        # If there are no alignments, return an empty string
        if not alignments:
            return ""

        # Get the first alignment
        first_alignment = alignments[0]

        # Format the first alignment
        formatted_alignment = SequenceVisualisation.format_custom_alignment(*first_alignment)

        # Return the formatted alignment
        return formatted_alignment

    @staticmethod
    def format_custom_alignment(align1, align2, score, begin, end):
        formatted_alignment = ""

        chunks1 = [align1[i:i + 120] for i in range(0, len(align1), 120)]
        chunks2 = [align2[i:i + 120] for i in range(0, len(align2), 120)]

        for i in range(len(chunks1)):
            start = i * 120 + 1
            end = start + len(chunks1[i]) - 1

            formatted_alignment += f"Query  {start:4}   {chunks1[i]}  {end:4}\n"
            formatted_alignment += " " * 14 + chunks1[i].replace('-', ' ') + "\n"
            formatted_alignment += f"Sbjct  {start:4}   {chunks2[i]}  {end:4}\n\n"

        # Return the formatted alignment
        return formatted_alignment

