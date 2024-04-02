import os
from datetime import datetime
from urllib import request

from django.http import HttpResponse
from django.shortcuts import render
from flask import Flask, render_template_string, jsonify, render_template, current_app
import plotly.graph_objects as go
import xml.etree.ElementTree as ET
from django.conf import settings

from home.scripts.export_data import Export


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

        # Extract the data you're interested in
        sequences = [hit.find('Hit_def').text for hit in root.findall('.//Hit')]
        matches = [int(hsp.find('Hsp_score').text) for hsp in root.findall('.//Hsp')]

        # colors = ['red' if match >= high_similarity_threshold else 'blue' for match in matches]
        colors = [Visualisation.get_color(match) for match in matches]

        # Extract additional information about the sequences
        sequence_info = [f"Name: {hit.find('Hit_def').text}<br>Other info: {hit.find('Hit_accession').text}" for hit
                         in
                         root.findall('.//Hit')]

        sequences = [hit.find('Hit_def').text for hit in root.findall('.//Hit')]
        matches = [int(hsp.find('Hsp_score').text) for hsp in root.findall('.//Hsp')]

        # Extract the number of hits per sequence
        hits_per_sequence = {sequence: sequences.count(sequence) for sequence in sequences}
        categories = [hit.find('Hit_accession').text for hit in root.findall('.//Hit')]

        hits_per_category = {category: categories.count(category) for category in categories}

        fig_hits_per_category = go.Figure(
            data=[go.Pie(values=list(hits_per_category.values()),
                         meta=[str(key) for key in hits_per_category.keys()],
                         hoverinfo='none',

                         hovertemplate='Séquence: %{meta}<br>Étiquette: %{label}<br>Pourcentage: %{percent}',
                         textfont_size=1)],
            layout=go.Layout(
                title_text='Distribution des coups directs parmi différentes catégories',
                autosize=True,

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
            # per should be rounded to 2 decimal places
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
           # sequence_info.append(subject)
        # Convert the figure to HTML and remove the surrounding <html> and <body> tags
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
            'sequence_info': sequence_info , 'subject': subject
        })
