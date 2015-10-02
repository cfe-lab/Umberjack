"""
Launches a pool of processes on a single node.
Useful when there is no MPI installed.
"""



import  logging
import pool.pool_traceback
import config.settings
from pool.basepool import  BasePool

config.settings.setup_logging()

LOGGER = logging.getLogger(__name__)


class OneNodePool(BasePool):
    """
    Uses python multiprocessing built in process pool.
    """

    def __init__(self, concurrent_windows, **kwargs):
        super(BasePool, self).__init__(**kwargs)
        self.concurrent_windows = concurrent_windows


    def is_parent(self):
        return True

    def spread_work(self, func, work_args_iter, callback=None):
        nodepool = pool.pool_traceback.LoggingPool(processes=self.concurrent_windows)
        process_results = []
        for work_args in work_args_iter:
            process_result = nodepool.apply_async(func=func, args=(), kwds=work_args, callback=callback)
            process_results.append(process_result)

        nodepool.close()  # no more work in the queue allowed
        LOGGER.debug("Done launching pool queue of " + str(len(process_results)) + " processes.  Wait for them to finish.")
        nodepool.join()

        done_results = []
        while len(done_results) != len(process_results):
            for process_result in process_results:
                try:
                    if process_result.ready():
                        done_results.append(process_result)
                        result = process_result.get()
                except Exception, e:
                    LOGGER.error("Error in one of replica processes:\n" + e.message)

        LOGGER.debug("Done pool queue of " + str(len(process_results)) + " processes")