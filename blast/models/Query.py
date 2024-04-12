import xml.etree.ElementTree as ET

from blast import models
from blast.scripts.utils.XMLParser import XMLParser


class Query(models.Model):
    query_id = models.CharField(max_length=200)
    query_def = models.CharField(max_length=2000)
    query_len = models.IntegerField()
    program = models.CharField(max_length=200)
    version = models.CharField(max_length=200)
    db = models.CharField(max_length=200)
    query_sequence = models.CharField(max_length=5000)

    def __str__(self):
        return self.query_id

    @classmethod
    def from_XML_File(cls, xml_file_path, query_sequence):
        root = XMLParser.parse_xml_file(xml_file_path)
        return XMLParser.create_query_from_xml(root, query_sequence)

    @classmethod
    def from_XML_string(cls, xml_string, query_sequence):
        root = XMLParser.parse_xml_string(xml_string)
        return XMLParser.create_query_from_xml(root, query_sequence)

    @classmethod
    def to_string(cls, query):
        root = ET.Element("BlastOutput")
        query_id = ET.SubElement(root, "BlastOutput_query-ID")
        query_id.text = query.query_id
        query_def = ET.SubElement(root, "BlastOutput_query-def")
        query_def.text = query.query_def
        query_len = ET.SubElement(root, "BlastOutput_query-len")
        query_len.text = str(query.query_len)
        program = ET.SubElement(root, "BlastOutput_program")
        program.text = query.program
        version = ET.SubElement(root, "BlastOutput_version")
        version.text = query.version
        db = ET.SubElement(root, "BlastOutput_db")
        db.text = query.db
        query_sequence = ET.SubElement(root, "BlastOutput_query-seq")
        query_sequence.text = query.query_sequence
        return ET.tostring(root)
