import os

import plotly.graph_objects as go
from django.conf import settings
from django.shortcuts import render

from blast.models.Query import Query
from blast.models.Sequence import Sequence
from blast.scripts.visualisation.SequenceVisualisation import SequenceVisualisation


class SequenceData:

    def home(self):
        sequences = Sequence.from_XML_File(os.path.join(settings.STATICFILES_DIRS[0], 'test3.xml'))

        matches = [sequence.score for sequence in sequences]
        colors = [SequenceVisualisation.get_color(match) for match in matches]
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
            sequence_info.append(sequence.get_sequence_info())

        fig_html = fig.to_html(full_html=False,
                               config={'displayModeBar': False, 'scrollZoom': False, 'displaylogo': False})
        fig_hits_per_category_html = fig_hits_per_category.to_html(full_html=False,
                                                                   config={'displayModeBar': False,
                                                                           'scrollZoom': False,
                                                                           'displaylogo': False})

        return render(self, "visualise/index.html", {
            'fig_html': fig_html,
            'fig_hits_per_sequence_html': fig_hits_per_category_html,
            'sequence_info': sequence_info, 'subject': "subject", "query_info": Query.from_XML_File(
                os.path.join(settings.STATICFILES_DIRS[0], 'test3.xml')).get_query_info()
        })
