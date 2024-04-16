from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


class CustomBlastParameters:
    """
    Class to store the parameters for a BLAST search.
    for now we are using this subpraameters:
    blast supports a lot of parameters that can be used to customize the search
    but for now we are using the default values for the parameters
    we can add more parameters to this class in the future
    its ok for now :)
    """

    def __init__(self, program: str, database: str, sequences: list[SeqRecord], output_format="XML", query_from=None, query_to=None,
                 expect_low=None, expect_high=None, is_async=False):
        self.program = program
        self.database = database
        self.sequences = sequences
        self.output_format = output_format
        self.query_from = query_from
        self.query_to = query_to
        self.expect_low = expect_low
        self.expect_high = expect_high
        self.is_async = is_async

    def to_dict(self):
        return {
            "program": self.program,
            "database": self.database,
            "sequences": [{"id": seq.id, "description": seq.description, "sequence": str(seq.seq)} if isinstance(seq, SeqRecord) else seq for seq in self.sequences],

            "output_format": self.output_format,
            "query_from": self.query_from,
            "query_to": self.query_to,
            "expect_low": self.expect_low,
            "expect_high": self.expect_high,
            "is_async": self.is_async
        }

    @staticmethod
    def from_dict(data):
        sequences = [SeqRecord(Seq(seq["sequence"]), id=seq["id"], description=seq["description"]) for seq in
                     data["sequences"]]
        return CustomBlastParameters(
            program=data["program"],
            database=data["database"],
            sequences=sequences,
            output_format=data["output_format"],
            query_from=data["query_from"],
            query_to=data["query_to"],
            expect_low=data["expect_low"],
            expect_high=data["expect_high"],
            is_async=data["is_async"]
        )