import logging
import socket

class HostnameFilter(logging.Filter):
    def filter(self, record):
        record.hostname = socket.gethostname()
        return True  # a return of non-zero indicates the record should be logged
