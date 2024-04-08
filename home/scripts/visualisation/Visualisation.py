import os

import plotly.graph_objects as go
from Bio import pairwise2
from django.conf import settings
from django.shortcuts import render

from home.models import Sequence


class Visualisation:

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

    def home(self):
        sequences = Sequence.from_XML(os.path.join(settings.STATICFILES_DIRS[0], 'test3.xml'))

        matches = [sequence.score for sequence in sequences]
        colors = [Visualisation.get_color(match) for match in matches]
        sequence_info = [f"Name: {sequence.description}<br>Other info: {sequence.access}" for sequence in sequences]

        hits_per_sequence = {sequence.description: sequences.count(sequence) for sequence in sequences}
        categories = [sequence.access for sequence in sequences]
        hits_per_category = {category: categories.count(category) for category in categories}

        fig_hits_per_category = go.Figure(
            data=[go.Pie(values=list(hits_per_category.values()),
                         meta=[str(key) for key in hits_per_category.keys()],
                         hoverinfo='none',
                         hovertemplate='Séquence: %{meta}<br>Étiquette: %{label}<br> Pourcentage: %{percent}',
                         textfont_size=1)],
            layout=go.Layout(
                title_text='Distribution des coups directs parmi différentes catégories',
                autosize=True,
                paper_bgcolor='rgba(0,0,0,0)',  # Set the paper (entire figure) background to transparent
            ),
        )

        fig = go.Figure(
            data=[go.Bar(y=matches, text=matches, textposition='auto', orientation='v', hovertext=sequence_info, )],
            layout=go.Layout(
                title_text='Résultats de recherche BLAST',
                xaxis_title='Séquence',
                yaxis_title='Nombre de correspondances',
                yaxis_categoryorder='sum descending',
                autosize=True,
                paper_bgcolor='rgba(0,0,0,0)',  # Set the paper (entire figure) background to transparent
            ),
            frames=[
                go.Frame(
                    data=[go.Bar(y=matches, orientation='v', hovertext=sequence_info, marker_color=colors, )], )]
        )

        sequence_info = []
        for sequence in sequences:
            info = {
                "Name": sequence.description,
                "Other_info": sequence.access,
                "Score": sequence.score,
                "Bit_Score": sequence.bit_score,
                "E_value": sequence.e_value,
                "Identity": sequence.identity,
                "Length": sequence.length,
                "Per": sequence.per,
                "Query_Sequence": sequence.query_sequence,
                "Hit_Sequence": sequence.hit_sequence,
                "Midline": sequence.midline
            }
            sequence_info.append(info)

        fig_html = fig.to_html(full_html=False,
                               config={'displayModeBar': False, 'scrollZoom': False, 'displaylogo': False})
        fig_hits_per_category_html = fig_hits_per_category.to_html(full_html=False,
                                                                   config={'displayModeBar': False,
                                                                           'scrollZoom': False,
                                                                           'displaylogo': False})

        return render(self, "visualise/index.html", {
            'fig_html': fig_html,
            'fig_hits_per_sequence_html': fig_hits_per_category_html,
            'sequence_info': sequence_info
        })

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
        formatted_alignment = Visualisation.format_custom_alignment(*first_alignment)

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
