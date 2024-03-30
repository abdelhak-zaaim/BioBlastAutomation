from Bio import pairwise2
from Bio.Align import substitution_matrices as sm
import matplotlib.pyplot as plt

def compare_sequences(seq1, seq2):
    # Get the scoring matrix
    matrix = sm.load("BLOSUM62")

    # Perform global alignment
    alignments = pairwise2.align.globalds(seq1, seq2, matrix, -10, -0.5)

    # Get the best alignment (first one in the list)
    best_alignment = alignments[0]

    # Return the alignment score
    aligned_seq1, aligned_seq2, score, begin, end = best_alignment
    return score

# Test the function with two sequences
if __name__ == "__main__":
    sequences = [("ACTTT", "ACTTT"), ("ACGTT", "ACTTT"), ("ACTTT", "ACGTT")]
    scores = [compare_sequences(seq1, seq2) for seq1, seq2 in sequences]

    # Create a bar chart of the scores
    plt.bar(range(len(scores)), scores)
    plt.xlabel('Sequence pair')
    plt.ylabel('Alignment score')
    plt.show()