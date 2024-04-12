from Bio.Blast import NCBIXML
from BioSQL import BioSeqDatabase

from blast.scripts.utils.Constants import Constants


class DatabaseManager:
    def __init__(self):
        self.server = BioSeqDatabase.open_database(driver=Constants.DATABASE_DRIVER, user=Constants.DATABASE_USER, passwd=Constants.DATABASE_PASSWD, host=Constants.DATABASE_HOST, db=Constants.DATABASE_NAME)

    def save_blast_results_to_db(self, result_handle):
        # Parse the BLAST results into a SeqRecord object
        blast_record = NCBIXML.read(result_handle)

        # Load the SeqRecord object into the BioSQL database
        db = self.server[Constants.DATABASE_NAME]
        db.load(blast_record)