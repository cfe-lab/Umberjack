"""
Stolen from http://streamhacker.com/2010/04/08/python-logging-filters/
Since logging.fileConfig is unable to handle filters, we need to override the logging handlers to do this.
"""
import logging
import logging.handlers
import filters

class HostnameStreamHandler(logging.StreamHandler):
    """
    Overrides StreamHandler and tacks on a HostName filter so that we log the hostname in the log messages
    """
    def __init__(self, *args, **kwargs):
        logging.StreamHandler.__init__(self, *args, **kwargs)
        self.addFilter(filters.HostnameFilter())


class HostnameRotatingFileHandler(logging.handlers.RotatingFileHandler):
    """
    Overrides RotatingFileHandler and tacks on a HostName filter so that we log the hostname in the log messages
    """
    def __init__(self, *args, **kwargs):
        logging.handlers.RotatingFileHandler.__init__(self, *args, **kwargs)
        self.addFilter(filters.HostnameFilter())
