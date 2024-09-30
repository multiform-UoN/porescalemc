# MLMC_PORESCALE
# python3 software to generate random packings, CFD and statistical analysis
# -------------------------------------------------------------------------------------------
# AUTHOR: Matteo Icardi
# November 2014
#--------------------------------------------------------------------------------------------
# module for defining parallel distributed processes

from mlmc_par import *
mlmc_worker_base=mlmc_worker

import os
import sys
import rpyc
rpyc.core.protocol.DEFAULT_CONFIG['allow_pickle'] = True

class MLMCService(rpyc.Service):
    def on_connect(self):
        # code that runs when a connection is created
        # (to init the serivce, if needed)
        pass

    def on_disconnect(self):
        # code that runs when the connection has already closed
        # (to finalize the service, if needed)
        pass

    def exposed_runsample(self,mlmcsample): # this is an exposed method
        # it runs the MLMC sample (object defined in mlmc.py)
        return mlmcsample()
    
    def exposed_stress_test(self):
        for i in range(10000):
            for j in range(10000):
                a=i**1.123+j**3.323
    
    def exposed_kill(self): # kill server
        os._exit(0)
        

# --- Parallel workers
class mlmc_worker(mlmc_worker_base):
    def __init__(self, task_queue, result_queue):
        mlmc_worker_base.__init__(self,task_queue,result_queue)
        # find available host
        for i in range(len(hostnames)):
            if hostavail[i]<i[1]:
                self.port=18861+hostavail[i]
                self.c=rpyc.connect(hostnames[i],self.port)
                self.hostname=hostnames[i]
                hostavail[i]=False
                break

    def run(self):
        proc_name = self.name
        while True:
            next_task = self.task_queue.get()
            if next_task is None:
                # Poison pill means shutdown
                print("Shutting down worker ",proc_name)
                self.task_queue.task_done()
                try:
                    self.c.root.exposed_kill() # kill the server
                except:
                    self.c.close()
                break
            print(proc_name,": job ",next_task," started in host ",self.hostname)
            answer = self.c.root.runsample(next_task)
            print(proc_name,": job ",next_task," done in host ",self.hostname, answer)
            self.task_queue.task_done()
            self.result_queue.put(answer)
        return

def run_rpyc_server(hosts):
    from paramiko import SSHClient, SSHConfig
    
    # ssh config file
    config = SSHConfig()
    home = os.getenv('HOME')
    mlmcdir = os.getenv('MLMC_PORESCALE')
    config.parse(open(home+'/.ssh/config'))
    
    ssh=[]
    port=18861
    for i in hosts:
        o = config.lookup(i)
        # ssh client
        ssh_client = SSHClient()
        ssh_client.load_system_host_keys()
        ssh_client.connect(o[i[0]])#, username=o['user'], key_filename=o['identityfile'])
        ncores = 1  # one single rpyc server and different processes or threads (to do)
        # run a command
        print("SSH connected to ",i[0])
        ## ------ every line is a host with a specified number of cores. port is sequentially incremented in each host
        #for j in range(i[1]):
            #cmd = mlmcdir+"/script/run_rpyc_server.sh "+str(j)+" "+str(ncores)
            #cmd_bg="nohup "+cmd+" 0<&- &> run_rpyc_server.log &"
##           cmd = "screen -t "+i+" -d -m -L "+mlmcdir+"/script/run_rpyc_server.sh"  # "screen -X quit" to kill
        ## ------ every line is a core in a host. port is sequentially incremented among all cores
        cmd = mlmcdir+"/script/run_rpyc_server.sh "+str(port)+" "+str(ncores)
        cmd_bg="nohup "+cmd+" 0<&- &> run_rpyc_server.log &"
#        cmd = "screen -t "+i+" -d -m -L "+mlmcdir+"/script/run_rpyc_server.sh"  # "screen -X quit" to kill
        stdin, stdout, stderr = ssh_client.exec_command(cmd_bg)
        port+=1
        ssh.append([ssh_client, stdin, stdout, stderr])

    return ssh


if __name__ == "__main__": # file loaded from run_rpyc_server.sh script
    from rpyc.utils.server import ThreadedServer
    t = ThreadedServer(MLMCService, port = int(sys.argv[1]))
    try:
        t.start()
    except:
        print("rpyc already running or other error preventinv rpyc to run")
else:   # file imported as a module in mlmc.py
    taskQueue = multiprocessing.Queue
    resultQueue = multiprocessing.Queue

    # Open hostname file list and save it
    import os
    hostfile = os.getenv('HOSTFILE')
    f = open(hostfile, 'r')
    # Loop over lines and extract all hosts
    hostnames=[]
    for line in f:
        line = line.strip()
        columns = line.split()
        hostnames.append([columns[0], int(columns[1])])
    f.close()
    hostavail=[True for i in range(len(hostnames))]

