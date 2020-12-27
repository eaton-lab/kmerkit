#!/usr/bin/env python

"""
General utilities
"""

import sys
from loguru import logger



class KmpyError(Exception):
    """
    Exception handler that does clean exit for CLI, but also prints
    the traceback and cleaner message for API.
    """
    def __init__(self, *args, **kwargs):
        # raise the exception with this string message and a traceback
        Exception.__init__(self, *args, **kwargs)




class Logger():
    """
    Start a loguru logger to stdout and logfile
    """

    def __init__(self, logfile=None):
        self.logfile = (logfile if logfile else "/tmp/kmpy-log.txt")


    def start(self):
        """
        Config and start the logger
        """
        logger.configure(
            **{
            "handlers": [
                {    
                    "sink": sys.stdout, 
                    "format": (
                        "{time:YYYY-MM-DD-hh:mm} | "
                        "<cyan>{function}</cyan> | "
                        "<level>{message}</level>"
                    ),
                    "level": "DEBUG",
                },

                {
                    "sink": self.logfile,
                    "format": "{time:YYYY-MM-DD} | {function} | {message}",
                    "level": "INFO",
                }
            ],
            "extra": {"user": "deren"},
            }
        )
        logger.enable("Kmpy")
