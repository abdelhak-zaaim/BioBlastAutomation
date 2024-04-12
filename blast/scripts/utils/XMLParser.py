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



