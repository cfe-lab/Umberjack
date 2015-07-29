"""
MPI Pool implementation.
"""
from mpi4py import MPI
from UmberjackPool import UmberjackPool
import  logging
import traceback
from collections import deque
import config.settings

config.settings.setup_logging()

LOGGER = logging.getLogger(__name__)

# For MPI


class ReplicaInfo:
    """
    Keeps track of the replica information
    """
    def __init__(self, replica_rank, work_arguments, mpi_send_request, mpi_rcv_request):
        """
        :param replica_rank: integer replica rank (starts from 1)
        :param dict work_arguments:  dict of arguments sent to the replica to do work
        :param mpi_send_request: mpi request sent to replica
        :param mpi_rcv_request: mpi request received from replica
        """
        self.replica_rank = replica_rank
        self.work_arguments = work_arguments
        self.mpi_send_request = mpi_send_request
        self.mpi_rcv_request = mpi_rcv_request


class MPIUmberjackPool(UmberjackPool):
    """
    MPI implementation of work pool
    """
    PRIMARY_RANK = 0
    TAG_WORK = 1
    TAG_TERMINATE = 2

    def __init__(self, **kwargs):
        super(MPIUmberjackPool, self).__init__(**kwargs)
        self.rank = None
        try:
            self.rank = MPI.COMM_WORLD.Get_rank()
            LOGGER.debug("I am rank=" + str(self.rank))
            LOGGER.debug("I am on machine=" + str(MPI.Get_processor_name()))
        except Exception:
            LOGGER.exception("Uncaught Exception.  Aborting")
            MPI.COMM_WORLD.Abort()



    def is_parent(self):
        if self.rank == MPIUmberjackPool.PRIMARY_RANK:
            return True


    def launch_wait_replicas(self, work_args_iter):
        """
        Launch replica processes
        :return:
        """

        try:

            if self.rank != MPIUmberjackPool.PRIMARY_RANK:
                raise ValueError("Only the Parent can launch child replicas")

            pool_size = MPI.COMM_WORLD.Get_size()
            LOGGER.debug("Pool size = " + str(pool_size))

            available_replicas = deque(range(1, pool_size), maxlen=pool_size)
            busy_replica_2_request = {}

            for work_args in work_args_iter:
                # Wait until a replica is available
                while not available_replicas:
                    requests = [window_replica_info.mpi_rcv_request for window_replica_info in busy_replica_2_request.values()]
                    mpi_status = MPI.Status()
                    idx, err_msg = MPI.Request.waitany(requests, mpi_status)
                    done_replica_rank = mpi_status.Get_source()
                    available_replicas.append(done_replica_rank)
                    del busy_replica_2_request[done_replica_rank]
                    if err_msg:
                        LOGGER.error("Received error from replica=" + str(done_replica_rank) + " err_msg=" + str(err_msg))
                    else:
                        LOGGER.debug("Received success from replica=" + str(done_replica_rank))

                replica_rank = available_replicas.popleft()  # first in first out queue

                str_work_args = ', '.join('{}:{}'.format(key, val) for key, val in work_args.items())
                LOGGER.debug("Sending work_args to replica=" + str(replica_rank) + " args=" + str_work_args)

                send_request = MPI.COMM_WORLD.isend(obj=work_args, dest=replica_rank, tag=MPIUmberjackPool.TAG_WORK)  # non-blocking
                rcv_request = MPI.COMM_WORLD.irecv(dest=replica_rank, tag=MPI.ANY_TAG)  # non-blocking

                # MPI.MPI.COMM_WORLD.isend() stores a copy of the pickled (i.e. serialized) windows_args dict
                # in the buffer of a new MPI.Request object  (send_request)
                # The replica accesses the buffer in the send_request MPI.Request object when it calls MPI.MPI.COMM_WORLD.recv()
                # to retrieve work.
                # If the Primary does not keep send_request in scope, the send_request buffer can be overwitten by a random process by the time
                # the replica gets to it.
                # When the primary does a non-blocking call to retrieve an response form the replica in MPI.MPI.COMM_WORLD.irecv,
                # it gets another MPI.Request object (rcv_request) whose buffer will contain the response from the replica when it's finished.
                # The buffers within send_request and rcv_request at separate mem addresses.
                busy_replica_2_request[replica_rank] = ReplicaInfo(replica_rank=replica_rank,
                                                                         work_arguments=work_args,
                                                                         mpi_send_request=send_request,
                                                                         mpi_rcv_request=rcv_request)

            # Sent out all work requests.  Wait until all work is done
            while busy_replica_2_request:
                requests = [window_replica_info.mpi_rcv_request for window_replica_info in busy_replica_2_request.values()]
                mpi_status = MPI.Status()
                idx, err_msg = MPI.Request.waitany(requests, mpi_status)
                done_replica_rank = mpi_status.Get_source()
                available_replicas.append(done_replica_rank)
                del busy_replica_2_request[done_replica_rank]
                if err_msg:
                    LOGGER.error("Received error from replica=" + str(done_replica_rank) + " err_msg=" + str(err_msg))
                else:
                    LOGGER.debug("Received success from replica=" + str(done_replica_rank))
                

            # All the work is done.
            LOGGER.debug("Terminating replicas...")
            for replica_rank in range(1, pool_size):
                MPI.COMM_WORLD.isend(obj=None, dest=replica_rank, tag=MPIUmberjackPool.TAG_TERMINATE)
            LOGGER.debug("Done terminating replicas.")

        except Exception:
            LOGGER.exception("Uncaught Exception.  Aborting")
            MPI.COMM_WORLD.Abort()


    def replica_work(self, func):
        """
        This is a replica.  Waits for work from parent and does work.
        Tells parent when done.
        :return:
        """
        

        try:
            LOGGER.debug("Waiting for work")
            is_terminated = False
            while not is_terminated:
                try:
                    mpi_status = MPI.Status()
                    # block till the primary tells me to do something
                    work_args = MPI.COMM_WORLD.recv(source=MPIUmberjackPool.PRIMARY_RANK, tag=MPI.ANY_TAG, status=mpi_status)

                    if mpi_status.Get_tag() == MPIUmberjackPool.TAG_TERMINATE:
                        is_terminated = True
                        LOGGER.debug("Replica of rank %d directed to terminate by primary" % self.rank)
                    else:
                        str_work_args = ', '.join('{}:{}'.format(key, val) for key, val in work_args.items())
                        LOGGER.debug("Received work_args=" + str_work_args)
                        func(**work_args)
                        MPI.COMM_WORLD.send(obj=None, dest=MPIUmberjackPool.PRIMARY_RANK, tag=MPIUmberjackPool.TAG_WORK)
                except Exception:
                    LOGGER.exception("Failure in replica=" + str(self.rank))
                    err_msg = traceback.format_exc()
                    MPI.COMM_WORLD.send(obj=err_msg, dest=MPIUmberjackPool.PRIMARY_RANK, tag=MPIUmberjackPool.TAG_WORK)
        except Exception:
            LOGGER.exception("Uncaught Exception.  Aborting")
            MPI.COMM_WORLD.Abort()


    def spread_work(self, func, work_args_iter):
        """
        If parent process, then launches child replicas.
        If child replicas, then does work.
        :return:
        """
        try:
            if self.rank == MPIUmberjackPool.PRIMARY_RANK:  # parent
                self.launch_wait_replicas(work_args_iter)
            else:  # replica
                self.replica_work(func)
        except Exception:
            LOGGER.exception("Uncaught Exception.  Aborting")
            MPI.COMM_WORLD.Abort()

