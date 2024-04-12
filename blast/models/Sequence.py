import xml.etree.ElementTree as ET

from django.db import models


class Sequence(models.Model):
    description = models.CharField(max_length=2000)
    access = models.CharField(max_length=2000)
    score = models.IntegerField()
    bit_score = models.FloatField()
    e_value = models.FloatField()
    identity = models.IntegerField()
    length = models.IntegerField()
    per = models.FloatField()
    query_sequence = models.TextField()
    hit_sequence = models.TextField()
    midline = models.TextField()
    query_id = models.CharField(max_length=2000)

    def __str__(self):
        return self.description

    @classmethod
    def from_XML_File(cls, xml_file_path):
        tree = ET.parse(xml_file_path)
        root = tree.getroot()
        query_id = root.find('.//BlastOutput_query-ID').text,
        sequences = []

        for hit in root.findall('.//Hit'):
            per = int(hit.find('.//Hsp_identity').text) / int(hit.find('.//Hsp_align-len').text) * 100
            per = round(per, 2)
            sequence = cls(
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
                midline=hit.find('.//Hsp_midline').text,
                query_id=query_id

            )
            sequences.append(sequence)
        # get query_id

        return sequences

    @classmethod
    def from_XML_string(cls, xml_string):
        root = ET.fromstring(xml_string)
        sequences = []
        query_id = root.find('.//BlastOutput_query-ID').text,
        for hit in root.findall('.//Hit'):
            per = int(hit.find('.//Hsp_identity').text) / int(hit.find('.//Hsp_align-len').text) * 100
            per = round(per, 2)
            sequence = cls(
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
                midline=hit.find('.//Hsp_midline').text,
                query_id=query_id
            )
            sequences.append(sequence)

        return sequences
