"""
Launches a pool of processes on a single node.
Useful when there is no MPI installed.
"""



import  logging
import pool.pool_traceback
import config.settings
from UmberjackPool import UmberjackPool

config.settings.setup_logging()

LOGGER = logging.getLogger(__name__)


class OneNodeUmberjackPool(UmberjackPool):

    def __init__(self, **kwargs):
        super(OneNodeUmberjackPool, self).__init__(**kwargs)


    def is_parent(self):
        return True

    def spread_work(self, func, work_args_iter):
        nodepool = pool.pool_traceback.LoggingPool(processes=self.concurrent_windows)
        process_results = []
        for work_args in work_args_iter:
            process_result = nodepool.apply_async(func=func, args=(), kwds=work_args)
            process_results.append(process_result)

        nodepool.close()  # no more work in the queue allowed
        LOGGER.debug("Done launching pool queue of " + str(len(process_results)) + " processes.  Wait for them to finish.")


        done_results = []
        while len(done_results) != len(process_results):
            for process_result in process_results:
                try:
                    if process_result.ready():
                        done_results.append(process_result)
                        process_result.get()
                except Exception, e:
                    LOGGER.error("Error in one of replica processes:\n" + e.message)

        LOGGER.debug("Done pool queue of " + str(len(process_results)) + " processes")