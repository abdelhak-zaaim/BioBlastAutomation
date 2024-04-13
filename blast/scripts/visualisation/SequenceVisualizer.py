import base64
import io
import ssl
from collections import Counter

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from Bio import motifs
from Bio.Seq import Seq

ssl._create_default_https_context = ssl._create_unverified_context


class SequenceVisualizer:

    @staticmethod
    def create_sequence_logo(self, sequences, output_file):
        # Create a list of Seq objects
        instances = [Seq(seq) for seq in sequences]

        # Create a motif from these sequences
        m = motifs.create(instances)

        # Create sequence logo
        m.weblogo(output_file)

    @staticmethod
    def plot_kmer_frequency(sequence, k):
        # Generate all k-mers
        kmers = [sequence[i:i + k] for i in range(len(sequence) - k + 1)]

        # Count the frequency of each k-mer
        kmer_counts = Counter(kmers)

        # Create a bar plot of k-mer frequencies
        plt.bar(kmer_counts.keys(), kmer_counts.values())
        return plt

    @staticmethod
    def dotplot(seq1, seq2):
        # Create a 2D array filled with zeros
        matrix = np.zeros((len(seq1), len(seq2)))

        # Fill in 1s where the sequences match
        for i in range(len(seq1)):
            for j in range(len(seq2)):
                if seq1[i] == seq2[j]:
                    matrix[i][j] = 1
        plt.figure(dpi=300, figsize=(10, 6))
        # Plot the matrix
        plt.imshow(matrix)
        plt.xticks(np.arange(len(list(seq2))), list(seq2))
        plt.yticks(np.arange(len(list(seq1))), list(seq1))
        return plt

    @staticmethod
    def heatmap(seq1, seq2):
        # Create a 2D array filled with zeros
        matrix = np.zeros((len(seq1), len(seq2)))

        # Fill in the scores where the sequences match
        for i in range(len(seq1)):
            for j in range(len(seq2)):
                if seq1[i] == seq2[j]:
                    matrix[i][j] = 1  # or any score

        # Create a heatmap
        sns.heatmap(matrix, annot=True, fmt=".2f")
        return plt

    @staticmethod
    def convert_plot_to_base64(plt):
        # Save the plot to a BytesIO object
        buf = io.BytesIO()
        plt.savefig(buf, format='png')
        plt.close()

        # Get the value of the BytesIO object
        buf.seek(0)
        image = buf.read()

        # Encode the bytes to base64 and decode to string
        base64_string = base64.b64encode(image).decode()

        return base64_string



