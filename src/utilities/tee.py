"""
File Name:          tee.py
Project:            bioseq-learning

File Description:

"""
import os
import sys


class Tee(object):
    """Tee class for storing terminal output to files.
    This class implements a tee class that flush std terminal output to a
    file for logging purpose.
    """

    def __init__(self, log_name, mode='a'):

        self.__stdout = sys.stdout

        self.__log_name = log_name
        self.__mode = mode

        try:
            os.makedirs(os.path.dirname(log_name))
        except FileExistsError:
            pass

    def __del__(self):
        sys.stdout = self.__stdout

    def write(self, data):

        with open(self.__log_name, self.__mode) as file:
            file.write(data)

        self.__stdout.write(data)

    def flush(self):
        self.__stdout.flush()

    def default_stdout(self):
        return self.__stdout