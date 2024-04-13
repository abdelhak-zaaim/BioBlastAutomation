import xml.etree.ElementTree as ET

from django.db import models


class Sequence(models.Model):
    id = models.AutoField(primary_key=True)
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
    scientific_name = models.CharField(max_length=2000)

    def __str__(self):
        return self.description

    @classmethod
    def from_XML_File(cls, xml_file_path):
        tree = ET.parse(xml_file_path)
        root = tree.getroot()
        return cls._create_sequences_from_xml(root)

    @classmethod
    def from_XML_string(cls, xml_string):
        root = ET.fromstring(xml_string)
        return cls._create_sequences_from_xml(root)

    @classmethod
    def _create_sequences_from_xml(cls, root):
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
                query_id=query_id,
                scientific_name=cls.get_scientific_name_from_sequence_def(hit.find('Hit_def'))
            )
            sequences.append(sequence)

        return sequences

    def get_sequence_info(self):
        return {
            "Name": self.description,
            "Other_info": self.access,
            "Score": self.score,
            "Bit_Score": self.bit_score,
            "E_value": self.e_value,
            "Identity": self.identity,
            "Length": self.length,
            "Per": self.per,
            "Query_Sequence": self.query_sequence,
            "Hit_Sequence": self.hit_sequence,
            "Midline": self.midline,
            "Scientific_Name": self.scientific_name
        }

    @classmethod
    def get_scientific_name_from_sequence_def(cls, hit_def):

        if hit_def is not None:
            scientific_name = hit_def.text.split('[')[1].split(']')[0]
            return scientific_name
        else:
            return ""
