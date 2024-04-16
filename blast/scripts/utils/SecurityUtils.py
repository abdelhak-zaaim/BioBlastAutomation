import os


class SecurityUtils:
    @staticmethod
    def check_path_traversal(filename):
        '''
        Check if the filename contains '/', '\', or '..'
        :param filename:
        :return: true if the filename contains '/', '\', or '..'
        '''
        # Check if the filename contains '/', '\', or '..'
        return not any(x in filename for x in ['/', '\\', '..'])