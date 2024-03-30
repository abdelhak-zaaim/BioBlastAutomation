from Bio import AlignIO
from flask import Flask, render_template_string
import plotly.graph_objects as go

app = Flask(__name__)

@app.route('/')
def home():
    # Parse the alignment file
    alignment = AlignIO.read('../test.xml', 'fasta')

    # Extract the sequences and their lengths
    sequences = [record.id for record in alignment]
    lengths = [len(record.seq) for record in alignment]

    # Create a bar chart
    fig = go.Figure(data=go.Bar(x=sequences, y=lengths))

    # Convert the figure to HTML and remove the surrounding <html> and <body> tags
    fig_html = fig.to_html(full_html=False)

    # Render the figure in a template
    return render_template_string("""
        <html>
        <head>
            <title>Sequence Alignment</title>
        </head>
        <body>
            {{ fig_html | safe }}
        </body>
        </html>
    """, fig_html=fig_html)

if __name__ == '__main__':
    app.run(debug=True)