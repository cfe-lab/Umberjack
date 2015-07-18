"""
Copied code from http://stackoverflow.com/questions/6728236/exception-thrown-in-multiprocessing-pool-not-detected
This entire script is a workaround for python multiprocessing module inability to bubble exception stack trace
    to the parent process from a child process.
Explicitly sets multiprocessing log level to WARN.
"""
import traceback
from multiprocessing.pool import Pool
import logging

LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.WARN)


# Shortcut to multiprocessing's logger
def error(msg, *args):
    """
    use multiprocessing's logger to log an ERROR
    """
    return LOGGER.error(msg, *args)


class LogExceptions(object):
    """
    Bubble exception stack trace from child process to parent
    """
    def __init__(self, callable):
        self.__callable = callable
        return

    def __call__(self, *args, **kwargs):
        try:
            result = self.__callable(*args, **kwargs)

        except Exception as e:
            # Here we add some debugging help. If multiprocessing's
            # debugging is on, it will arrange to log the traceback
            error(traceback.format_exc())
            # Re-raise the original exception so the Pool worker can
            # clean up
            raise Exception

        # It was fine, give a normal answer
        return result
    pass


class LoggingPool(Pool):
    """
    Wrapper class for multiprocessing.Pool that logs that full exception stack trace from a child process.
    """

    def apply_async(self, func, args=(), kwds={}, callback=None):
        """
        Override multiprocessing.Pool.apply_async() method such that it logs full exception stack trace from child process.
        """
        return Pool.apply_async(self, LogExceptions(func), args, kwds, callback)

    def imap_unordered(self, func, iterable, chunksize=1):
        """
        Override multiprocessing.Pool.imap_unordered() method such that it logs full exception stack trace from child process.
        """
        return Pool.imap_unordered(self, LogExceptions(func), iterable, chunksize)

    def imap(self, func, iterable, chunksize=1):
        """
        Override multiprocessing.Pool.imap() method such that it logs full exception stack trace from child process.
        """
        return Pool.imap(self, LogExceptions(func), iterable, chunksize)