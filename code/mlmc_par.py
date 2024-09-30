# MLMC_PORESCALE
# python3 software to generate random packings, CFD and statistical analysis
# -------------------------------------------------------------------------------------------
# AUTHOR: Matteo Icardi
# November 2014
#--------------------------------------------------------------------------------------------
# module for defining parallel processes

import multiprocessing
import logging
#import collections
multiprocessing.log_to_stderr(logging.WARNING)
from general import memory

taskQueue = multiprocessing.JoinableQueue
resultQueue = multiprocessing.Queue

# --- Parallel workers
class mlmc_worker(multiprocessing.Process):
    def __init__(self, task_queue, result_queue,ramlimit):
        multiprocessing.Process.__init__(self)
        self.task_queue   = task_queue
        self.result_queue = result_queue
        self.ramlimit     = ramlimit
        #result_queue.cancel_join_thread()

    def run(self):
        import random
        import time
        random.seed(hash(time.time())+hash(self.name))
        while True:
            if memory()>self.ramlimit:
                print("RAM limit: Pausing worker ", self.name)
                time.sleep(60)
                continue
            next_task = self.task_queue.get()
            if next_task is None:
                # Poison pill means shutdown
                print("Shutting down worker ",self.name)
                self.task_queue.task_done()
                break
            print(self.name,": job ",next_task," started")
            answer = next_task()
            print(self.name,": job ",next_task," done\n result = ", answer)
            self.task_queue.task_done()
            self.result_queue.put(answer)
        self.result_queue.close()
        #self.terminate()
        return
  
