import xml.etree.ElementTree as ET


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

    # this function used to convert float to scientific notation
    @staticmethod
    def float_to_evalue(num):
        return "{:.2e}".format(num)

