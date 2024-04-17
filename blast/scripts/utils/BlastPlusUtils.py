import os

from pfe import settings


class BlastPlusUtils:
    @staticmethod
    def get_database_path_by_name(database_name):
        """
        Get the path of the database by its name
        :param database_name: the name of the database
        :return: the path of the database
        """
        # Get the path of the database by its name
        return os.path.join(settings.BLAST_DATABASES_PATH, database_name)