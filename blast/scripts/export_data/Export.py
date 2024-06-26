import csv
import xml.etree.ElementTree as ET
from xml.etree import ElementTree

from xhtml2pdf import pisa


class BlastResultExporter:
    def __init__(self, xml_file, html_template, output_pdf):
        self.xml_file = xml_file
        self.html_template = html_template
        self.output_pdf = output_pdf

    def blast_xml_to_csv(self, xml_file, csv_file):
        tree = ET.parse(xml_file)
        root = tree.getroot()

        data = []
        for iteration in root.find('BlastOutput_iterations'):
            for hit in iteration.find('Iteration_hits'):
                hit_data = {
                    'id': hit.find('Hit_id').text,
                    'def': hit.find('Hit_def').text,
                    'accession': hit.find('Hit_accession').text,
                    # Add more fields as needed
                }
                data.append(hit_data)

        # Write data to CSV
        with open(csv_file, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=data[0].keys())
            writer.writeheader()
            writer.writerows(data)

    def parse_xml(self):
        tree = ElementTree.parse(self.xml_file)
        root = tree.getroot()
        # Parse the XML data as needed
        # This will depend on the structure of your XML file
        data = {}
        return data

    def render_html(self, data):
        with open(self.html_template, 'r') as file:
            html = file.read().format(**data)
        return html

    def export_to_pdf(html):
        with open("tessttt.pdf", 'wb') as output_file:
            pisa_status = pisa.CreatePDF(html, dest=output_file)

        return not pisa_status.err

    def export(self):
        data = self.parse_xml()
        html = self.render_html(data)
        success = self.export_to_pdf(html)
        return success
