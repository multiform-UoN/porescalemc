# MLMC_PORESCALE
# python3 software to generate random packings, CFD and statistical analysis
# -------------------------------------------------------------------------------------------
# AUTHOR: Matteo Icardi
# November 2014
#--------------------------------------------------------------------------------------------
# module for defining serial process

import queue

taskQueue = queue.Queue
resultQueue = queue.Queue


# --- Serial (single) worker
class mlmc_worker(object):
    def __init__(self, task_queue, result_queue, ramlimit):
        self.task_queue   = task_queue
        self.result_queue = result_queue
        self.ramlimit     = ramlimit

    def start(self):
        while True:
            next_task = self.task_queue.get()
            if next_task is None:
                # Poison pill means shutdown
                print("End of serial jobs")
                self.task_queue.task_done()
                break
            print("Serial job started ",next_task)
            answer = next_task()
            print("Serial job done ", next_task, "\n result = ", answer)
            self.task_queue.task_done()
            self.result_queue.put(answer)
        return
    
    def join(self):
        return
