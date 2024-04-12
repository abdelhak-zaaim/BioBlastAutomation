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
        return cls._create_query_from_xml(root, query_sequence)

    @classmethod
    def from_XML_string(cls, xml_string, query_sequence):
        root = XMLParser.parse_xml_string(xml_string)
        return cls._create_query_from_xml(root, query_sequence)

    @classmethod
    def _create_query_from_xml(cls, root, query_sequence):
        return Query(
            query_id=root.find('.//BlastOutput_query-ID').text,
            query_def=root.find('.//BlastOutput_query-def').text,
            query_len=int(root.find('.//BlastOutput_query-len').text),
            program=root.find('.//BlastOutput_program').text,
            version=root.find('.//BlastOutput_version').text,
            db=root.find('.//BlastOutput_db').text,
            query_sequence=query_sequence
        )

    def get_query_info(self):
        return {
            "query-id": self.query_id,
            "query_definition": self.query_def,
            "query_length": self.query_len,
            "program": self.program,
            "version": self.version,
            "database": self.db,
            "query_sequence": self.query_sequence
        }
