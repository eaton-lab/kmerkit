#!/usr/bin/env python

"""
Starts an ipcluster (IPython cluster) and writes the pid/ipcluster.pid 
into a specified directory for killing cluster later.
"""

import os
import sys
import time
import platform
import tempfile
import shutil
import signal
import subprocess
from threading import Thread
from contextlib import AbstractContextManager
import ipyparallel as ipp
from loguru import logger


class Cluster(AbstractContextManager):
    """
    Stores cluster information.
    """
    def __init__(self, cores=4, tmpdir=None, ipyclient=None):

        # get random name for this cluster to avoid conflicts
        self.name = str(os.getpid())
        self.tmpdir = (tmpdir if tmpdir is not None else tempfile.gettempdir())

        # cluster config for starting an ipcluster
        self.cores = cores
        self.ipyclient = ipyclient

        # to be filled: client, and pids for shutting down
        self.auto_started = False
        self.engine_pids = {}

        # cluster config files directory path
        self.clusterdir = os.path.join(
            self.tmpdir,
            self.name + "_tmpfiles",
            ".kmerkit-cluster",
        )

        # paths to controller files created in a hidden dir in workdir
        self.hub_config_file = os.path.join(
            self.clusterdir, "security", 
            f"ipcontroller-{self.name}-client.json")

        # futures should be stored in this dict
        self.rasyncs = {}


    def __enter__(self):
        """
        Connect to a IPython cluster that either already exists or 
        is automatically started. If already exists then user should
        enter an ipyclient object. If auto-starting user only provides
        ncores argument (0 means all available cores).
        """
        self.auto_connect()
        self.store_pids()
        return self


    def store_pids(self):
        """
        Store process ids of ipcluster and engine processes
        """
        # store pids for interrupting engine jobs
        self.engine_pids = self.ipyclient[:].apply(os.getpid).get_dict()


    def auto_connect(self):
        """
        Start an ipcluster instance, or connect to a running instance.
        """
        # the user provided an already active ipclient through the API
        if hasattr(self.ipyclient, "ids"):
            return

        # start an ipcluster with a local PID stored in workdir/pids
        self.auto_started = True
        self.start_ipcluster()

        # It can sometimes take a long time for the ipcluster to startup
        # engines on some machines, especially HPC, so we check the 
        # number of engines it finds until all requested engines started.
        # Loop until all engines are found
        while 1:

            # if ipcluster file is not made yet this will raise an
            # OSError, in which case just want to wait a tiny bit.
            try: 
                # connect to client
                self.ipyclient = ipp.Client(self.hub_config_file)

                # if all engines are found then break
                if len(self.ipyclient) == self.cores:
                    logger.info(f"cluster established: {len(self.ipyclient)} engines")
                    break

                # close it again (prevents too many open files)
                self.ipyclient.close()

            except OSError:
                pass


    def start_ipcluster(self):
        """
        Start ipcluster by calling it from subprocess. We give it the
        argument profile_dir so that it will write the pid of the job
        to a local file so we can easily shut it down later. We also
        """
        # avoid starting a null cluster
        assert self.cores > 0, (
            'start ipcluster using the .run function, or set .cores > 0')

        # existing .clustdir will raise an exception
        if os.path.isdir(self.clusterdir):
            shutil.rmtree(self.clusterdir)

        # path to the ipcluster binary
        ipcluster_bin = os.path.join(sys.prefix, "bin", "ipcluster")

        # build the command
        cmd = [
            ipcluster_bin,
            "start",
            "--cluster-id={}".format(self.name),
            "--profile-dir={}".format(self.clusterdir),
            "--n={}".format(self.cores),
            "--quiet",
            "--daemonize",            
        ]

        # start binary running on fork
        subprocess.run(cmd, check=True)
        logger.debug(" ".join(cmd))
        time.sleep(1)


    def shutdown(self):
        """
        Calls stop from the ipcluster binary 
        This is really really really the only reliable way we could 
        find to get everything to stop without zombies.
        """
        # path to the ipcluster binary
        ipcluster_bin = os.path.join(sys.prefix, "bin", "ipcluster")

        # build the command
        cmd = [
            ipcluster_bin,
            "stop",
            "--cluster-id={}".format(self.name),
            "--profile-dir={}".format(self.clusterdir),
            "--quiet",
        ]

        # run command and block until finished
        logger.debug(" ".join(cmd))
        subprocess.run(cmd, check=True)
        time.sleep(1)


    def cleanup(self, rasyncs=None):
        """
        If ipcluster was auto-started then it is shutdown, otherwise
        we simply cancell all jobs.
        """
        # protect from keyboard interrup while cleaning up
        try:
            # auto-started: stop ipcluster
            if self.auto_started:
                # stop future jobs
                self.ipyclient.abort()
                # stop current jobs
                if rasyncs:
                    for job in rasyncs:
                        if rasyncs[job].stdout:
                            pid = int(rasyncs[job].stdout.strip())
                            os.kill(pid, signal.SIGINT)
                # give it a second
                time.sleep(1)
                # send shutdown to cluster and controller
                self.shutdown()
                logger.info('ipcluster stopped')

            # user-entered: leave ipcluster alive, just stop jobs
            else:
                self.ipyclient.abort()
                for pid in self.engine_pids.values():
                    os.kill(pid, signal.SIGINT)
                time.sleep(1)
                if not self.ipyclient.outstanding:
                    self.ipyclient.purge_everything()

        except KeyboardInterrupt:
            self.ipyclient.close()
            logger.warning("cleaning up...")


    def __exit__(self, exc_type, exc, exc_tb):
        """
        Start cleanup on a separate thread that cannot be interrupted
        by a keyboard interrupt.
        """
        # log any errors
        if exc is not None:
            logger.exception("An exception occurred at the following:")
            # logger.error(f"{str(exc_type)}: {exc}")

        # send shutdown job
        waitjob = Thread(target=self.cleanup, args=(self.rasyncs,))
        waitjob.start()
        waitjob.join()



def get_num_cpus():
    """
    Return the effective number of CPUs in the system as an integer for
    either Unix or MacOSX (Darwin). This code is mostly copied from a
    similar implementation in IPython.
    If it can't find a sensible answer, it returns 1.
    """
    if platform.system() == "Linux":
        ncpus = os.sysconf("SC_NPROCESSORS_ONLN")
    else:
        proc = subprocess.run(
            ['sysctl', '-n', 'hw.ncpu'], check=True, capture_output=True)
        ncpus = proc.stdout.decode().strip()
    try:
        ncpus = max(1, int(ncpus))
    except:
        ncpus = 1
    return ncpus


if __name__ == "__main__":

    import kmerkit
    kmerkit.set_loglevel("DEBUG")

    with Cluster(cores=2) as c:
        # raise kmerkit.utils.KmerkitError("aaah")
        print(c.ipyclient.ids, c.rasyncs)

    # import superbpp
    # superbpp.set_loglevel("DEBUG")
    # PROJECT = superbpp.load_project("/tmp/test1.json")

    # CLUSTER = Cluster(PROJECT)
    # CLUSTER.start(cores=2)
    # time.sleep(1)
    # CLUSTER.cleanup()



# def run(self, cores=4, ipyclient=None):
#     """
#     Distribute jobs on a local or distributed cluster. Almost all 
#     users will prefer the auto-start option (ipyclient=None), but 
#     in some cases such as large HPC clusters the user could start
#     an ipcluster with MPI over many machines and distribute to 
#     those cores here.
#     """
#     # wrap tracker for cluster cleanup -------------------------
#     try:
#         # start or connect to a running ipcluster instance. 
#         self.cluster.start(cores, ipyclient)

#         # distribute initial set of jobs, this will be called again
#         # inside of the tracker to start new nodes as others finish
#         sub = NodeSubmitter(self.project, self.cluster.ipyclient, self.rasyncs)
#         self.rasyncs = sub.rasyncs

#         # start tracker which will re-call NodeWriter and NodeSubmitter
#         # until all jobs are finished, or interrupted.
#         self.iptracker()

#     except KeyboardInterrupt:
#         logger.warning("keyboard interrupt by user, cleaning up.")

#     except SuperbppError:
#         logger.error("Error occurred, cleaning up.")
#         raise

#     except Exception:
#         logger.error("An unexpected error occurred, cleaning up.")
#         raise

#     finally:
#         self.cluster.cleanup_safely(self.rasyncs)
