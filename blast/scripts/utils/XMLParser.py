import xml.etree.ElementTree as ET

from blast.models.Query import Query
from blast.models.Sequence import Sequence


class XMLParser:
    @staticmethod
    def parse_xml_file(xml_file_path):
        tree = ET.parse(xml_file_path)
        root = tree.getroot()
        return root

    @staticmethod
    def parse_xml_string(xml_string):
        root = ET.fromstring(xml_string)
        return root

    @staticmethod
    def create_query_from_xml(root, query_sequence):
        return Query(
            query_id=root.find('.//BlastOutput_query-ID').text,
            query_def=root.find('.//BlastOutput_query-def').text,
            query_len=int(root.find('.//BlastOutput_query-len').text),
            program=root.find('.//BlastOutput_program').text,
            version=root.find('.//BlastOutput_version').text,
            db=root.find('.//BlastOutput_db').text,
            query_sequence=query_sequence
        )

    @staticmethod
    def create_sequences_from_xml(root):
        sequences = []
        for hit in root.findall('.//Hit'):
            per = int(hit.find('.//Hsp_identity').text) / int(hit.find('.//Hsp_align-len').text) * 100
            per = round(per, 2)
            sequence = Sequence(
                description=hit.find('Hit_def').text,
                access=hit.find('Hit_accession').text,
                score=int(hit.find('.//Hsp_score').text),
                bit_score=float(hit.find('.//Hsp_bit-score').text),
                e_value="{:0.2f}".format(float(hit.find('.//Hsp_evalue').text)),  # round to 3 decimal places
                identity=int(hit.find('.//Hsp_identity').text),
                length=int(hit.find('.//Hsp_align-len').text),
                per=per,
                query_sequence=hit.find('.//Hsp_qseq').text,  # New field
                hit_sequence=hit.find('.//Hsp_hseq').text,  # New field
                midline=hit.find('.//Hsp_midline').text  # New field
            )
            sequences.append(sequence)
        return sequences

