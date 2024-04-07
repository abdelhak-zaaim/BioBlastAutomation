import os
import xml.etree.ElementTree as ET

import plotly.graph_objects as go
from Bio import pairwise2
from django.conf import settings
from django.shortcuts import render


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
        tree = ET.parse(os.path.join(settings.STATICFILES_DIRS[0], 'test3.xml'))
        root = tree.getroot()

        matches = [int(hsp.find('Hsp_score').text) for hsp in root.findall('.//Hsp')]

        colors = [Visualisation.get_color(match) for match in matches]

        sequence_info = [f"Name: {hit.find('Hit_def').text}<br>Other info: {hit.find('Hit_accession').text}" for hit
                         in
                         root.findall('.//Hit')]

        sequences = [hit.find('Hit_def').text for hit in root.findall('.//Hit')]
        matches = [int(hsp.find('Hsp_score').text) for hsp in root.findall('.//Hsp')]

        hits_per_sequence = {sequence: sequences.count(sequence) for sequence in sequences}
        categories = [hit.find('Hit_accession').text for hit in root.findall('.//Hit')]

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
            data=[go.Bar(y=matches, text=matches, textposition='auto', orientation='v', hovertext=sequence_info,  )],
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
        # Extract additional information about the sequences
        sequence_info = []
        subject = root.find('.//BlastOutput_query-def').text

        for hit in root.findall('.//Hit'):
            per = int(hit.find('.//Hsp_identity').text) / int(hit.find('.//Hsp_align-len').text) * 100

            per = round(per, 2)
            info = {
                "Name": hit.find('Hit_def').text,
                "Other_info": hit.find('Hit_accession').text,
                "Score": hit.find('.//Hsp_score').text,
                "Bit_Score": hit.find('.//Hsp_bit-score').text,
                "E_value": hit.find('.//Hsp_evalue').text,
                "Identity": hit.find('.//Hsp_identity').text,
                "Gaps": hit.find('.//Hsp_gaps').text,
                "Length": hit.find('.//Hsp_align-len').text,
                "Query_Sequence": hit.find('.//Hsp_qseq').text,
                "Hit_Sequence": hit.find('.//Hsp_hseq').text,
                "Midline": hit.find('.//Hsp_midline').text,
                "Num": hit.find('.//Hsp_num').text,
                "Per": per,

            }

            sequence_info.append(info)

        fig_html = fig.to_html(full_html=False,
                               config={'displayModeBar': False, 'scrollZoom': False, 'displaylogo': False})
        fig_hits_per_category_html = fig_hits_per_category.to_html(full_html=False,
                                                                   config={'displayModeBar': False,
                                                                           'scrollZoom': False,
                                                                           'displaylogo': False})

        # Render the figure in a template
        return render(self, "visualise/index.html", {
            'fig_html': fig_html,
            'fig_hits_per_sequence_html': fig_hits_per_category_html,
            'sequence_info': sequence_info, 'subject': subject
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
