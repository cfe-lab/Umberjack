"""
Base Pool Class that abstracts the implementations of process pools to child classes.
Master-replica relationship where master decides what tasks should be given to which replicas.
Replicas that finish first are expected to report back to the master for more work.
"""

class BasePool(object):
    def __init__(self):
        pass

    def is_parent(self):
        pass

    def spread_work(self, func, work_args_iter, callback=None):
        """
        Callback is only executed by the master process
        :param func:
        :param work_args_iter:
        :param callback:
        :return:
        """
        pass

    def terminate(self):
        """
        Shutdown master and child processes
        :return:
        """
        pass