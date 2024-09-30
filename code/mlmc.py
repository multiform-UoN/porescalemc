# MLMC_PORESCALE
# python3 software to generate random packings, CFD and statistical analysis
# -------------------------------------------------------------------------------------------
# AUTHOR: Matteo Icardi
# November 2013
#--------------------------------------------------------------------------------------------
# module for MLMC estimators

try:
    random
except:
    from general import *

from mlmcDefault import *

try:
    exec(open("mlmcDict.py").read(),globals())  # moved in the run.py script
except IOError:
    print("Warning! file mlmcDict.py not found or not readable. Loaded default values")

try:
    exec(open("runDict.py").read(),globals())  # moved in the run.py script
except IOError:
    print("Warning! file runDict.py not found or not readable. Using module specific options")

# this works only if mlmc.py is imported as a module
def reread_input(globals_=globals()):
    exec(open("mlmcDict.py").read(),globals_)
    
# compute confidence constant Calpha
from scipy.stats import norm as norm_scipy
calpha= norm_scipy.ppf(confidence)

# definitions for serial/parallel/distributed computing
if (distributed): # or for distributed
    from mlmc_hybrid import *
    print("launching all python rpyc servers")
    sshclients=run_rpyc_server()
elif (nthreads=="adaptive" or nthreads>1):
    from mlmc_par import *
else:
    #from mlmc_ser import *
    from mlmc_par import * # use parallel library also for 1 process


# --- class definining a sample
class mlmc_sample(object):
    def __init__(self, l, name, targetf):
        self.l = l
        self.name = name
        self.targetf = targetf
        
    def __call__(self):
        t1=time.time()
        res=self.targetf(self.l,self.name)
        return self.l,res,time.time()-t1
    
    def __str__(self):
        return "sample "+self.name
# ---

# --- class defining a MLMC or MC estimator
class estimator:
    def __init__(self, name , sol, geo, tol=final_tol, split=error_split):
        self.name=name
        self.geo=geo
        self.sol=sol
        # setup the first guess
        self.tot_level=tot_level
        self.m0=int(m0)
        self.ml=[int(self.m0)]
        self.final_tol=tol
        self.error_split=split
        for i in range(tot_level-1):
            self.ml.append(int(self.ml[i]/mlmc_mratio))
        self.q=None
        self.w=None
        self.cprod=array([[1.] for i in range(tot_level)])
        self.csum=array([[0.] for i in range(tot_level)])

    # ---------------------------------
    def solve(self):
        if (mlmc_algorithm=="giles"):
            self.mlmc_adaptive()
        elif (mlmc_algorithm=="none"):
            self.mlmc_run()
        elif (mlmc_algorithm=="mesh_convergence"):
            self.mesh_convergence()
            return self.q,self.w
        self.save_qw()
        self.compute_errors()
        if (scale_qoi in ("prod","sum")):
            print("Rescaled results:")
            self.scale_mlmc(scale_qoi)
        return final_postprocess(self.QQ, self.stat_err, self.bias)

    # generate all cases without running
    def setup(self):
        #TODO
        return 1

    def mlmc_adaptive(self):
        self.converged=False
        self.ml_relax=update_ml_relaxation
        i=0
        # adaptive loop
        print("--------- MLMC SIMULATION STARTED")
        while (not self.converged):
            i+=1
            print("loop with ",self.tot_level, " levels", self.ml)
            self.mlmc_run()
            self.ml_relax=1-(1-update_ml_relaxation)**i
            self.update_ml()
            if (intermediate_output):
                self.save_qw()

            self.mlmc_run()

            self.check_convergence()
            if (intermediate_output):
                self.save_qw()

        print("--------- MLMC SIMULATION ENDED")

    def scale_mlmc(self,scale="sumprod"):
        if (scale in ("sumprod","prodsum")):
            self.scale_mlmc("default")
            self.compute_errors()
            self.scale_mlmc("sum")
            return
        q=self.q
        nvar=int(len(self.QQ)/2)
        r=[]
        r.append(zeros(nvar))
        rsum=[]
        rsum.append(zeros(nvar))
        for i in range(1,len(q)):
            dql=self.dq[i][:nvar]
            ql=self.dq[i][nvar:]
            if (scale=="prod"):
                r.append(ql/(ql-dql)) # only dependent samples
#                               r.append(ql/self.dq[i-1][nvar:]) # only independent samples
#                               r.append(ql/(self.ml[i-1]*self.dq[i-1][nvar:]-self.ml[i]*(dql-ql))*(self.ml[i]+self.ml[i-1]))  # both
                for j in range(len(q[i])):
                    q[i][j][:nvar]=(1-r[i])*q[i][j][nvar:]+r[i]*q[i][j][:nvar] #prod
            elif (scale=="sum"):
                r.append(dql)     # only dependent samples
#                               for j in range(len(q[i])):
#                                       q[i][j][:nvar]-=r[i]  #sum1
            else:
                vdq=array(self.dqvar[i][:nvar])
                vl=array(self.qvar[i])
                vlm1=array(self.qm1var[i])
                cov=array(self.covar[i])
                r.append(-cov/(2*vlm1))
                for j in range(len(q[i])):
                    q[i][j][:nvar]=(1-r[i])*q[i][j][nvar:]+r[i]*q[i][j][:nvar] #prod
        if (scale=="sum"):
            c=zeros_like(self.dq) # sum
            prodr=r[-1]*self.cprod[-1][-2:]
        else:
            c=ones_like(self.dq) # prod
            prodr=r[-1]
        for i in range(len(q)-2,-1,-1):
            c[i]=concatenate([prodr,prodr])
            if (scale=="sum"):
                for j in range(len(q[i+1])):
                    q[i+1][j][:nvar]=q[i+1][j][:nvar]-prodr  #sum2
                prodr=prodr+r[i]*self.cprod[i][-2:] # sum2
#                               prodr=r[i] # sum1
            else:
                prodr=prodr*r[i] # prod

        print("Scaling variables: ",r)
        print("Correctors: ",scale,c)

        self.q=q
        if (scale=="sum"):
            self.csum=array(self.csum)+array(c)  # sum
        else:
            self.cprod=array(self.cprod)*array(c) # prod
        self.r=r


    def update_ml(self):
        print("------- update number of samples")

        self.compute_errors()
        w=array(self.w)/array([len(x) for x in self.q])
        print('mean work ',w)
        y=min_samples*ones_like(self.ml)

        # select variable to compute the number of samples
        if (active_var<0 or active_var>len(self.QQ)):
            vrange=range(int(len(self.QQ)/2))
        else:
            vrange=[active_var]

        tol=array(final_tol)*ones(int(len(self.QQ)/2))

        print("update samples relaxation ",self.ml_relax)
        # compute the number of samples
        for i in vrange:
            y=maximum(y,(calpha/((1-self.error_split)*tol[i]*self.QQ[i]))**2*sqrt(transpose(self.dqvar)/array(w))[i,:]*sum(sqrt(array(w)*array(self.dqvar)[:,i])))
        print("------- proposed number of samples", y)
        y=maximum(array(self.ml),self.ml_relax*y+(1.-self.ml_relax)*array(self.ml))
        #y[1:]=minimum(y[1:],maximum(min_samples,y[:-1]/mlmc_mratio))
        y=maximum(y,min_samples)
        y=ceil(y).astype(int)
        print("------- new number of samples", y)
        print("------- additional  samples", array(y)-array(self.ml))
        self.ml=y
        return None

    def compute_errors(self):
        q=self.q
        nvar=int(len(q[0][0])/2)
        est=sum([mean(x,0) for x in q],0)
        modq=qoi_twopass(est[:nvar])
        q=[[concatenate([modq(x[nvar:])-modq(x[nvar:]-x[:nvar]),modq(x[nvar:])]) for x in y] for y in q] # works only for unweighted pair #TODO extend for weighted sums with cprod and csum
        nvar=int(len(q[0][0])/2)
        q[0]=[concatenate([x[nvar:],x[nvar:]]) for x in q[0]] # fix first level where second part can be screwed up by the nonlinear transformation
        w=array(self.w)
        dq=[]
        dq_abs=[]
        dqvar=[]
        dqvar_err=[]
        qvar=[]
        covar=[]
        qm1var=[]
        Ml_eff=[]
        cumsum=[]
        for i in range(len(q)):
            ci=array(self.cprod[i])*ones_like(q[0][0])
            ki=array(self.csum[i])*ones_like(q[0][0])
            dq.append(ci*mean(q[i],0)+ki)                                                 # E(Ql-Ql-1)
            dq_abs.append(mean(abs(ci*array(q[i])+ki),0))
            dqvar.append(samplevar(ci*q[i]+ki,0))                                         # Var(Ql -Ql-1)
            qvar.append(dqvar[-1][nvar:])                                                 # Var(Ql)
            qm1=array([[z[x+nvar]-z[x] for x in range(nvar)] for z in q[i] ])             # Var(Ql-1)
            qm1var.append(samplevar(ci[-nvar:]*qm1+ki[-nvar:],0)) # Var(Ql-1)
            covar.append(dqvar[-1][:nvar]-dqvar[-1][nvar:]-samplevar(ci[-nvar:]*qm1+ki[-nvar:],0))        # -2Cov(Ql,Ql-1)
            dqvar_err.append(calpha*samplevar(ci*q[i]+ki,0)*sqrt(2.0/(len(q[i]))))        # error using CLT
            Ml_eff.append(len(q[i]))
            cumsum.append(sum(dq,0))
        varnorm_eff=transpose(transpose(dqvar)/array(Ml_eff))
        stat_err=calpha*sqrt(sum((varnorm_eff),0))
        QQ=sum(dq,0)

        self.dq=array(dq)
        self.dq_err=array(calpha*sqrt(varnorm_eff))
        self.dqvar=array(dqvar)
        self.dqvar_err=array(dqvar_err)
        self.QQ=array(QQ)
        self.covar=covar
        self.qvar=qvar
        self.qm1var=qm1var
        self.stat_err=array(stat_err)
        self.bias=array(dq[-1])
        print("Mean: ", self.dq)
        print("Variances: ", self.dqvar)
        print("Correctors: ", self.csum,self.cprod)
        print("Estimator: ", self.QQ[:nvar])
        print("Statistical error (relative): ", self.stat_err[:nvar]/self.QQ[:nvar])
        print("Bias error (relative): ", self.bias[:nvar]/self.QQ[:nvar])



    def check_convergence(self):
        print("--------- check convergence")
        self.compute_errors()
        self.converged=True

        # select variable to compute the bias
        if (active_var<0 or active_var>len(self.QQ)):
            vrange=range(int(len(self.QQ)/2))
        else:
            vrange=[active_var]

        tol=array(final_tol)*ones(int(len(self.QQ)/2))

        for i in vrange:
            if ((abs(self.bias[i]/self.QQ[i]))>self.error_split*tol[i] and self.tot_level<max_level):
                print("BIAS ABOVE THE TOLERANCE, INCREASE LEVELS")
                if (self.tot_level>=max_level):
                    print("REACHED MAX_LEVEL")
                else:
                    self.converged=False
                    self.tot_level+=1
                    self.ml=append(self.ml,max(5,int(min(self.ml[-1]**2/self.ml[-2],self.ml[-1]/mlmc_mratio))))
                    self.q.append([])
                    self.cprod=append(self.cprod,ones_like(self.cprod[0]))
                    self.csum=append(self.csum,zeros_like(self.cprod[0]))
                break
            if ((abs(self.stat_err[i]/self.QQ[i]))>(1-self.error_split)*tol[i]):
                print("STATISTICAL ERROR ABOVE THE TOLERANCE")
                self.ml=self.ml+2
                self.converged=False
                break

    def makeplot(self,fileformat):
        makeplot(self.q,self.w,self.name,fileformat)

    def save_qw(self):
        save_qw(self.q,self.w,self.name)

    ## -------- OVERLOADED OPERATORS
    #def __call__(self, level): # (re)set level and compute QoI
        #self.geom.set_input(level)
        #self.solver.set_input(level
        #[res,log]=solver(self,level)
        #return res

    #def __call__(self):       # compute QoI
        #[res,log]=solver(self)
        #return res

    #def __add__(self, other):
        #return self()+other()

    #def __sub__(self, other):
        #return self()-other()


    #----------------- RUN A HIERARCHY CONVERGENCE STUDY ON DIFFERENT LEVELS FOR A FIXED REALIZATION
    def mesh_convergence(self):
        name=self.name
        l=0
        work=[]
        q=[]
        t=time.time()
        g=geo_interface(l,name,hierarchy)
        s=sol_interface(g,hierarchy)
        log=s.setup(g)
        [res,log1]=s.solve()
        res=array(flatten(final_qoi(res,g,s)))
        log+=log1
        q.append(res)
        levtime=(time.time()-t)
        print(res,levtime,s.name)
        work.append(levtime) # real CPU time cost
        for l in range(1,tot_level):
            t=time.time()
            g.set_input(l)
            log+=s.setup(g)
            [res,log1]=s.solve()
            res=array(flatten(final_qoi(res,g,s)))
            log+=log1
            q.append(res)
            levtime=(time.time()-t)
            print(res,levtime,s.name)
            work.append(levtime) # real CPU time cost
        s.close(cleardat)
        self.q=q
        self.w=work
        #return q,work
    #--------------------------

    #----------------- RUN A SINGLE SAMPLE
    def single(self,l,g,s):
        t=time.time()
        name=self.name
        g=geo_interface(l,name,hierarchy)
        s=sol_interface(g,hierarchy)
        log=s.setup(g)
        [res,log1]=s.solve()
        res=array(flatten(final_qoi(res,g,s)))
        log+=log1
        s.close(cleardat)
        self.q=res
        self.w=time.time()-t
        self.g=g
        self.s=s
        #return res
    #--------------------------

    #--------------------- RUN THE MLMC CODE
    # ------Unsuccesful samples are simply dropped
    def mlmc_run(self):

        # retrieve available results from the queue
        def try_getting_results():
            #ires=0
            #itime=0
            #out=False
            while True:
                #itime+=100
                #if (intermediate_output and itime%intermediate_output_t==0):
                    #save_qw(qoi,[sum(i) for i in t],self.name)
                try:
                    val=results.get(False)
                except queue.Empty:
                    break
                #val=results.get()
                #results.task_done()
                if (isnan(val[1]).any()):
                    continue
                l=val[0]
                qoi[l].append(val[1])
                t[l].append(val[2])
                #out=True
                #if timeout==True:
                #    if tasks.empty():
                #        break
                #ires+=1
                #if (intermediate_output and ires%intermediate_output_n==0):
                    #save_qw(qoi,[sum(i) for i in t],self.name)
            return #out

        basename=self.name

        #  tasks and results lists/queues
        tasks = taskQueue()
        results = resultQueue()

        # --- choice of the difference operator
        if (mlmc_estimator=='triplet'):
            targetf=triplet
        elif (mlmc_estimator=='test'):
            targetf=test_function
        else:
            targetf=pair

        qoi=[]
        work=[]
        t=[]

        #TODO check what to do for nprocs=max or nthreads=adaptive
        #for l in range(len(self.ml)):
            ## ----- find how many cores we can use per each sample
            #if (nprocs=="max"):
                #np=int(get_cpu()/pow(refratio,len(self.ml)-level-1))
            #else:
                #np=int(nprocs/pow(refratio,len(self.ml)-level-1))
            #np=max(1,np)

            ## -----  create dummy sample to check how many cores and threads we can use
            #runcase=geo_interface(level,"test_sample_l"+str(level),hierarchy)
            #s=sol_interface(runcase,hierarchy)
            #s.nprocs=min(np,s.max_procs())
            #if (nthreads=="adaptive"):
                #nt = max(1,int(get_cpu()/s.nprocs))
            #else:
                #nt = max(1,int(nthreads/nthreads_scale**level))
        nt=nthreads  # now we simply assume that user chose well
        np=ones_like(self.ml)*1 # and we assume same number of cores (1)for everything
        #nt=np/sum(np)     # nthreads per level

        # --- list of working units (defined in mlmc_par.py or mlmc_ser.py)
        workers=[]
        #workers=[mlmc_worker(tasks,results,ramlimit)]
        ## serial job cannot be run async
        #if (nthreads>1):
            #workers[0].start()
        
        for i in range(len(workers),nt):
            worker=mlmc_worker(tasks,results,ramlimit)
            worker.start()
            workers.append(worker)


        # -----  loop on levels
        ntot=1
        submitted=zeros_like(self.ml)

        # ---- check how many samples we need to compute and how to split them
        tosubmit=[]
        for level in range(len(self.ml)):
            qoi.append([])
            t.append([])
            if (self.q is None) or (not reuse_sample) :
                tosubmit.append(max(0,self.ml[level]))
            else:
                tosubmit.append(max(0,self.ml[level] - len(self.q[level])))
        tosubmit=array(tosubmit)
        if self.w is not None:
            cost=array(self.w)/array([len(x) for x in self.q[0:len(self.w)]])
            for missing_level in range(len(self.ml)-len(cost)):
                cost=concatenate([cost,[cost[-1]**2/cost[-2]]])
        else:
            cost=array([mlmc_mratio**i for i in range(len(self.ml))])
        nthreads_per_level=(1./cost)/sum(1./cost)*nt
        #nthreads_per_level=(tosubmit/sum(tosubmit)*nt).astype(int)
        #nthreads_per_level=maximum(1,nthreads_per_level)
        #nthreads_per_level[0]=nt-sum(nthreads_per_level[1:])
        nthreads_comulative=zeros_like(nthreads_per_level)

        if (intermediate_output):
            t0=time.time()

        while (submitted<tosubmit).any():
            nthreads_comulative+=nthreads_per_level
            # --- loop on levels
            for level in range(len(tosubmit)):

            ## -----  create additional workers for distributed computing
            #for i in range(len(workers),nt):
                #worker=mlmc_worker(tasks,results,ramlimit)
                #worker.start()
                #workers.append(worker)
            ## -----  send poison pill to extra workers
            #for i in range(nt,len(workers)):
                #tasks.put(None)


                # ----- loop on number of samples
                n=0
                while submitted[level]<min(nthreads_comulative[level],tosubmit[level]):
                    name=basename+"_"+"l%01d"%level+"_"+"%06d"%submitted[level]
                    ntot+=1
                    tasks.put(mlmc_sample(level, name, targetf))
                    submitted[level]+=1
                    #print("job ",submitted[level],"of",tosubmit[level],"level",level,"added to queue")
            #if try_getting_results(cost[0]/2):
                #if (intermediate_output):
                    #t1=time.time()
                    #if t1>t0+intermediate_output_t or ntot%intermediate_output_n==0:
                        #save_qw(qoi,[sum(i) for i in t],self.name)
                        #t0=time.time()
            # --- end loop on levels

        ## serial job cannot be run async
        #if nthreads==1:
            #workers[0].start()

        
        print("--- All simulations queued...waiting")
        # --- send poison pill to terminate all the processes
        for i in workers:
            tasks.put(None)
        #tasks.close()
        #tasks.join_thread()
        #for i in workers:
        #    i.join()
            #if i.is_alive():
                #i.terminate()
        tasks.join()

        try_getting_results()

        #results.close()
        #results.join_thread()

        #for i in workers:
            #i.join()
            ##i.terminate()
##        tasks.join()

        # compute total work per level
        work = [sum(i) for i in t]

        # merge samples and store results
        if (reuse_sample):
            self.q,self.w=merge_data(qoi,work,self.q,self.w)
        else:
            self.q=qoi
            self.w=work
        return
    #--------------------- END MLMC_RUN
    
# --- end estimator class definitions



# --- OTHER FUNCTIONS

#----------------- RUN A MONTECARLO PAIR/TRIPLET
def pair(l,name):
    #s.nprocs=min(np,s.max_procs())  #TODO override solver nprocs based on s.max_procs()
    res0=0
    log0=''.encode('ascii')
    if (l>0):
        g=geo_interface(l-1,name,hierarchy,tot_level)
        s=sol_interface(g,hierarchy)
        log0+=s.setup(g)
        [res0,log1]=s.solve()
        log0+=log1
        if accept(res0,g,s): # MCMC-type rejection #TODO
            res0=array(flatten(final_qoi(res0,g,s)))
        else:
            print("sample rejected")
            return nan
        if isnan(res0).any():
            print("-----> error in sample "+s.name,res0)
            return nan
        g.set_input(l)
    else:
        g=geo_interface(l,name,hierarchy,tot_level)
        s=sol_interface(g,hierarchy)
    log0=log0+s.setup(g)
    [res,log1]=s.solve()
    res=array(flatten(final_qoi(res,g,s)))
    log=log0+log1
    if isnan(res).any():
        print("-----> error in sample "+s.name,res)
        return nan
    res = concatenate([res-res0,res])
    s.close(cleardat)
    return res

def triplet(l,name):
    resm1=0
    logm1=''.encode('ascii')
    resm2=0
    logm2=''.encode('ascii')
    if (l>1):
        g=geo_interface(l-2,name,hierarchy,tot_level)
        s=sol_interface(g,hierarchy)
        logm2+=s.setup(g)
        [resm2,log1]=s.solve()
        resm2=array(flatten(final_qoi(resm2,g,s)))
        logm2+=log1
        if isnan(resm2).any():
            print("-----> error in sample "+s.name,resm2)
            return nan
        g.set_input(l-1)
    elif l==1:
        g=geo_interface(l-1,name,hierarchy,tot_level)
        s=sol_interface(g,hierarchy)
    if (level>0):
        logm1+=s.setup(g)
        [resm1,log1]=s.solve()
        resm1=array(flatten(final_qoi(resm1,g,s)))
        logm1+=log1
        if isnan(resm1).any():
            print("-----> error in sample "+s.name,resm1)
            return nan
        g.set_input(l)
    else:
        g=geo_interface(l,name,hierarchy,tot_level)
        s=sol_interface(g,hierarchy)
    log=s.setup(g)
    [res,log1]=s.solve()
    res=array(flatten(final_qoi(res,g,s)))
    log+=log1
    log=logm2+logm1+log
    if isnan(res).any():
        print("-----> error in sample "+s.name,res)
        return nan
    if (level>1):
        kappa= 2.0 ** -beta  # given by weak convergence rate
#               kappa= 1.0          # linear Taylor expansion
#               kappa= minimum(1.0,abs((res-resm1)/(resm1-resm2)))  # tuned on the single realization
#               print('kappa',kappa)
        res=res-resm1-kappa*(resm1-resm2)
    else:
        res-=resm1
    s.close(cleardat)
    return res
#--------------------------


#-------------------- DUMMY FUNC FOR TESTING PARALLELIZATION
def test_function(l,name):
    print("test func "+name)
    r=[random.random()-.5 for i in range(3)]
    if l==0:
        x= array([y/(l+1)+l/(l+1) for y in r])
    else:
        x= array([y/(l+1)+l/(l+1) for y in r])
        x-= array([y/(l)+(l-1)/(l) for y in r])
    return x
#--------------------

# -------------- FUNCTIONS TO LOAD AND SAVE OUTPUT DATA q and w
# save q and w after a simulation
def save_qw(qq=None,ww=None,basename=""):
    qfile=basename+"_Q"
    wfile=basename+"_W"
    if (qq is None):
        pickle.dump(q,open(qfile,"wb")) # save the quantities of interest for each level and sample
    else:
        pickle.dump(qq,open(qfile,"wb")) # save the quantities of interest for each level and sample
    if (ww is None):
        pickle.dump(w,open(wfile,"wb")) # save the work loads for each level
    else:
        pickle.dump(ww,open(wfile,"wb")) # save the work loads for each level


# load q and w from file
def load_qw(basename=""):
    try:
        q=pickle.load( open( basename+"_Q", "rb" ) )
    except IOError:
        print("file Q not found")
#               name=input("load unfinished results? enter the name or empty to exit\n")
        q=load_unfinished(basename)
    try:
        w=pickle.load( open( basename+"_W", "rb" ) )
    except IOError:
        print("file W not found")
        w=[]
#       q=[[y for y in x if y is not None] for x in q] #remove None values
#       q=[[y for y in x if isnan(y).any()==False] for x in q]  #remove NaN values
    return q,w

# remove the outliers from 2d array q
def remove_outliers(q,nsigma=4):
    q2=[]
    qrem=[]
    for l in range(len(q)):
        m=mean(q[l],0)
        sd=samplestd(q[l],0)
        q2.append([x for x in q[l] if ((abs(x-m)<nsigma*sd).all() and (abs(x)<1.e5).all())])
        qrem.append([x for x in q[l] if not ((abs(x-m)<nsigma*sd).all() and (abs(x)<1.e5).all())])
        print("at level", l, "removed", len(qrem[l]), "samples")
    return q2

# remove NaN from 2d array q
def remove_nan(q):
    q=[[y for y in x if y is not nan] for x in q]
    return q

# merge data from 2 simulations
# number of levels of q1 should not be smaller than q2
def merge_data(q1,w1,q2,w2):
    if (q2 is None or w2 is None):
        return q1,w1
    for l in range(len(q2)):
        if len(q2[l])==0:
            continue
        for n in range(len(q2[l])):
            q1[l].append(q2[l][n])
        w1[l]+=w2[l]
    return q1,w1


# load unfinished simulations using log files (only vector q)
def load_unfinished(name,loadpair=True,all_iterations=False):
    import re
    import fnmatch
    q1=[]
    q2=[]
    q=[]
#       q3=[]
#       q4=[]
    covq=[]
    maxl=0
    maxn=[]
    for file in os.listdir('.'):
        if fnmatch.fnmatch(file, name+'*.qoi'):
            prop=re.findall('\d+',file[len(name):])
            l=int(prop[1])
            n=int(prop[0])
            pair=(l==int(prop[2][0]))
            while (l>maxl-1):
                q1.append([])
                q2.append([])
                maxn.append(0)
                maxl+=1
            while (n>maxn[l]-1):
                q1[l].append(nan)
                q2[l].append(nan)
                maxn[l]+=1
            res=array(flatten(final_qoi(qoi(l,file[:-8],all_iterations)) )) # TO DO change to manage vector valued qoi
            if (not isnan(res).any()):
                if (pair):
                    q1[l][n]=res
                else:
                    q2[l][n]=res
    for l in range(maxl):
        q.append([])
#               q3.append([])
#               q4.append([])
        for n in range(len(q1[l])):
            if (q1[l][n] is not nan):
                if (l>0):
                    if (q2[l-1][n] is not nan):
                        if (loadpair):
                            q[l].append(q1[l][n]-q2[l-1][n])
                        else:
                            q[l].append(q1[l][n])
                            q[l-1].append(q2[l-1][n])
#                                                       q3[l].append(q1[l][n])
#                                                       q4[l-1].append(q2[l-1][n])
                else:
                    q[l].append(q1[l][n])
#               if (not loadpair and l>0):
#                       covq.append(cov(array([array(q4[l-1]),array(q3[l])])))
#                       cova=sqrt((covq[l-1]))
#                       print('mean=',mean(q[l]),' var=',samplevar(q[l]),' det=',linalg.det(cova),' diff_std=',(cova[0][0]-cova[1][1])**2,' corr=',cova[1][0]*cova[0][1]/cova[0][0]/cova[1][1])
    return q #, covq
# -----------------------


# -------------- FUNCTIONS TO COMPUTE AND PLOT MLMC STATISTICS
def makeplot(q1=None,w1=None,basename="",fileformat='png',totwork=True):
    # totwork is for retro-compatibility only, when work is loaded as average value per realization
    plt.close('all')
    if (q1 is None and 'q' not in globals()):
        (q1,w1)=load_qw(basename)
    q=q1
    nvar=int(len(q[0][0])/2)
    est=sum([mean(x,0) for x in q],0)
    modq=qoi_twopass(est[:nvar])
    q=[[concatenate([modq(x[nvar:])-modq(x[nvar:]-x[:nvar]),modq(x[nvar:])]) for x in y] for y in q] # works only for unweighted pair #TODO extend for weighted sums with cprod and csum
    nvar=int(len(q[0][0])/2)
    q[0]=[concatenate([x[nvar:],x[nvar:]]) for x in q[0]] # fix first level where second part can be screwed up by the nonlinear transformation
    w=array(w1)
    if totwork:
        w=w/array([len(x) for x in q])
    x=range(len(q))
    dq=[]
    dq2=[]
    dq_abs=[]
    dqvar=[]
    dqvar_err=[]
    Ml=[]
    Ml_eff=[]
    cumsum=[]
    totwork=0
    for i in range(len(q)):
        dq.append(mean(q[i],0))
        dq2.append(mean(array(q[i])**2,0))
        dq_abs.append(mean(abs(array(q[i])),0))
        dqvar.append(samplevar(q[i],0))
        dqvar_err.append(calpha*dqvar[-1]*sqrt(2.0/(len(q[i])-1)))
        if ('w' in vars()):
            totwork+=len(q[i])*w[i]
        Ml.append(m0/(mlmc_mratio**i))
        Ml_eff.append(len(q[i]))
        cumsum.append(sum(dq,0))
    varnorm=(transpose(dqvar)/array(Ml))
    varnorm_eff=(transpose(dqvar)/array(Ml_eff))
    stat_err=calpha*sqrt(sum(transpose(varnorm_eff),0))
    if ('w' is not None):
        stat_errSL=calpha*sqrt(dqvar[0]*totwork/w[-1])
    else:
        stat_errSL=stat_err
    QQ=sum(dq,0)

    print("Estimator: "+ str(QQ))
    print("Statistical error: ", str(array(stat_err)), "Percentage:", str(array(stat_err)/QQ))
    print("Bias: ", str(dq[-1]), "Percentage:", str(dq[-1]/QQ))
    print("Error single level: ", str(stat_errSL), "Improvement:", str(stat_errSL/stat_err))
    filename=basename+'_stats.dat'
    f=open(filename,'w')
    f.write("\nEstimator "+ str(QQ))
    f.write("\nStatistical error "+ str(array(stat_err)))
    f.write("\nStatistical error percentage "+ str(array(stat_err)/QQ))
    f.write("\nBias error"+ str(dq[-1]))
    f.write("\nBias error percentage"+ str(dq[-1]/QQ))
    f.close()

    for nv in range(len(q[0][0])):
        # ---- PLOT SAMPLES OF THE PAIRS PER EACH LEVEL
        fig=plt.figure()
        plt.hold(True)
    #       plt.title('Individual samples')
        plt.xlabel(r'level $\ell$')
        plt.ylabel(r'$Q_\ell-Q_{\ell-1}$')
        plt.grid(True)
        #plt.xscale('log')
        plt.yscale('log')
        plt.xticks(x)
        for i in range(len(q)):
            y=abs(array(q[i]))[:,nv]
            plt.plot(ones_like(y)*i,y,'x',color='k')
            plt.text(i, max(y)*1.2, str(len(y)))
        plt.savefig(basename+'_'+str(nv)+'_samples.'+fileformat,transparent=True, bbox_inches='tight', pad_inches=0.2)
        plt.close(fig)
        # ---- PLOT HISTOGRAMS OF THE PAIRS PER EACH LEVEL
        for i in range(len(q)):
            if (len(q[i])<20):
                continue
            fig=plt.figure()
            plt.hold(True)
    #               plt.title('Histogram of pair realizations in level '+str(i))
            plt.xlabel(r'$Q_\ell-Q_{\ell-1}$')
            plt.ylabel('PDF with '+str(len(q[i]))+' samples')
            plt.grid(True)
            plt.gca().xaxis.set_major_formatter(ffloat)
            plt.hist(array(q[i])[:,nv], min(50,int(len(q[i])/5)), normed=1, facecolor='g', alpha=0.75, histtype='bar')#,log=True)
            plt.savefig(basename+'_'+str(nv)+'_samples'+str(i)+'.'+fileformat,transparent=True, bbox_inches='tight', pad_inches=0.2)
            plt.close(fig)


        # ---- BEGIN MULTIPLE FIGURE
        fig=plt.figure()
        plt.hold(True)
        plt.xticks(x)
        plt.xlabel(r'level $\ell$')
        plt.grid(True)
        #plt.xscale('log')
        plt.yscale('log')
        #plt.title('Sample mean')
        plt.hold(True)
        plt.gca().xaxis.set_major_locator(plt.MaxNLocator(integer=True))
        # ---- PLOT EXPECTED VALUE OF THE PAIRS PER EACH LEVEL
        y=abs(array(dq))[:,nv]
        #y=y/y[0]
        plt.axis([min(x)-.1,max(x)+.1,min(y)*0.8,max(y)*1.2])
        y_error=calpha*sqrt(varnorm_eff[nv,:])
        y_error0=minimum(y_error,y-1e-13)
        special.errorfill(x, y, yerr=array([y_error0,y_error]),alpha_fill=0.1, ls='dashed', lw=3,color='k')
        plt.plot(x, y, ls='None', marker='o', markersize=10, label=r'$E[Q_{\ell}-Q_{\ell-1}]$',color='k')
        rate=-log(y[:-1]/y[1:])/log(2.0)
        for i,v in enumerate(rate):
            if i==0:
                continue
            plt.annotate("%.1f"%v,((x[i]+x[i+1])/2,sqrt(y[i]*y[i+1])),color='k')
        # ---- PLOT SECOND MOMENT OF THE PAIRS PER EACH LEVEL
        y=abs(array(dq2))[:,nv]
        #y=y/y[0]
        ax=plt.axis()
        plt.axis([min(x)-.1,max(x)+.1,min(min(y)*0.8,ax[2]),max(max(y)*1.2,ax[3])])
        p2, = plt.plot(x, y, clip_on=False, ls='dashed', lw=3, marker='o', markersize=10, label=r'$E[(Q_{\ell}-Q_{\ell-1})^2]$',color='b')
        #rate=-log(y[:-1]/y[1:])/log(2.0)
        #for i,v in enumerate(rate):
            #if i==0:
                #continue
            #plt.annotate("%.1f"%v,((x[i]+x[i+1])/2,sqrt(y[i]*y[i+1])))
        # ---- PLOT MEAN SQUARE OF THE PAIRS PER EACH LEVEL
        y=abs(array(dq)**2)[:,nv]
        #y=y/y[0]
        ax=plt.axis()
        plt.axis([min(x)-.1,max(x)+.1,min(min(y)*0.8,ax[2]),max(max(y)*1.2,ax[3])])
        p2, = plt.plot(x, y, clip_on=False, ls='dashed', lw=3, marker='o', markersize=10, label=r'$E[Q_{\ell}-Q_{\ell-1}]^2$',color='c')
        #rate=-log(y[:-1]/y[1:])/log(2.0)
        #for i,v in enumerate(rate):
            #if i==0:
                #continue
            #plt.annotate("%.1f"%v,((x[i]+x[i+1])/2,sqrt(y[i]*y[i+1])))
        # ----- PLOT VARIANCES PER EACH LEVEL
        y=array(dqvar)[:,nv]
        #y=y/y[0]
        ax=plt.axis()
        plt.axis([min(x)-.1,max(x)+.1,min(min(y)*0.8,ax[2]),max(max(y)*1.2,ax[3])])
        y_error=array(dqvar_err)[:,nv]
        y_error0=minimum(y_error,y-1e-13)
        special.errorfill(x, y, yerr=array([y_error0,y_error]),alpha_fill=0.1, ls='dashed', lw=3,color='r')
        plt.plot(x, y, ls='None', marker='o', markersize=10, label=r'$V[Q_\ell-Q_{\ell-1}]$',color='r')
        rate=-log(y[:-1]/y[1:])/log(2.0)
        for i,v in enumerate(rate):
            if i==0:
                continue
            plt.annotate("%.1f"%v,((x[i]+x[i+1])/2,sqrt(y[i]*y[i+1])),color='r')
        # ---- PLOT NORMALIZED VARIANCES PER EACH LEVEL
        #y=sqrt(array(varnorm[nv,:]))
        #y=y/y[0]
        #plt.plot(x,y,clip_on=False,label="proposed")
        #rate=-log(y[:-1]/y[1:])/log(2.0)
        #for i,v in enumerate(rate):
            #if i==0:
                #continue
        #       plt.annotate("%.1f"%v,((x[i]+x[i+1])/2,sqrt(y[i]*y[i+1])))
        y=sqrt(array(varnorm_eff[nv,:]))
        #y=y/y[0]
        ax=plt.axis()
        plt.axis([min(x)-.1,max(x)+.1,min(min(y)*0.8,ax[2]),max(max(y)*1.2,ax[3])])
        #plt.plot(x,y,clip_on=False,label=r'$\sqrt{V[Q_\ell-Q_{\ell-1}]}/\sqrt{M_\ell}$')
        plt.plot(x,y,clip_on=False,label=r'$Stat. err.}$',color='m')
        #rate=-log(y[:-1]/y[1:])/log(2.0)
        #for i,v in enumerate(rate):
            #if i==0:
                #continue
        #       plt.annotate("%.1f"%v,((x[i]+x[i+1])/2,sqrt(y[i]*y[i+1])))

        if ('w' in vars() and len(w)==len(q)):
            # ---- PLOT AVERAGE WORK LOAD PER LEVEL PER SINGLE RUN
            y=array(w)
            y=y/y[0]
            ax=plt.axis()
            plt.axis([min(x)-.1,max(x)+.1,min(min(y)*0.8,ax[2]),max(max(y)*1.2,ax[3])])
            plt.plot(x,y,clip_on=False,label=r'$W_{\ell}$',color='g')
            rate=-log(y[:-1]/y[1:])/log(2.0)
            for i,v in enumerate(rate):
                plt.annotate("%.1f"%v,((x[i]+x[i+1])/2,sqrt(y[i]*y[i+1])),color='g')
            # ---- PLOT AVERAGE WORK TIMES VARIANCE
            #y=sqrt(array(w)*array(dqvar)[:,nv])
            ##y=y/y[0]
            #plt.plot(x,y,clip_on=False,label=r'$\sqrt{W_{\ell}V_{\ell}}$')
            #rate=-log(y[:-1]/y[1:])/log(2.0)
            #for i,v in enumerate(rate):
                #plt.annotate("%.1f"%v,((x[i]+x[i+1])/2,sqrt(y[i]*y[i+1])))
            ## ---- PLOT OPTIMAL Ml
            #y=(calpha/(error_split*final_tol*QQ[nv]))**2*sqrt(transpose(dqvar)/array(w))[nv,:]*sum(sqrt(array(w)*array(dqvar)[:,nv]))
            #y=y/y[0]
            #plt.plot(x,y,clip_on=False,label=r'optimal $M_\ell$')
            #rate=-log(y[:-1]/y[1:])/log(2.0)
            #for i,v in enumerate(rate):
                #if i==0:
                    #continue
                #plt.annotate("%.1f"%v,((x[i]+x[i+1])/2,sqrt(y[i]*y[i+1])))
            ##y=array(Ml_eff)
            ##plt.plot(x,y,clip_on=False,label="effective")
            ##rate=-log(y[:-1]/y[1:])/log(2.0)
            #for i,v in enumerate(rate):
                ##plt.annotate("%.1f"%v,(sqrt(x[i]+x[i+1])/2,sqrt(y[i]*y[i+1])))

        #plt.legend([p1,p2,p3,p4],[r'$E[Q_{\ell}-Q_{\ell-1}]$',r'$E[(Q_{\ell}-Q_{\ell-1})^2]$',r'$E[Q_\ell]$',r'$V[Q_\ell-Q_{\ell-1}]$'])
        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        plt.savefig(basename+'_'+str(nv)+'.'+fileformat,transparent=True, bbox_inches='tight', pad_inches=0.2)
        plt.close(fig)

        # ---- SINGLE FIGURES
        # ---- PLOT EXPECTED VALUE OF THE PAIRS PER EACH LEVEL
        fig=plt.figure()
        plt.gca().xaxis.set_major_locator(plt.MaxNLocator(integer=True))
        plt.hold(True)
        y=abs(array(dq))[:,nv]
        plt.axis([min(x)-.1,max(x)+.1,min(y)*0.8,max(y)*1.2])
        y_error=calpha*sqrt(varnorm_eff[nv,:])
        y_error0=minimum(y_error,y-1e-13)
        special.errorfill(x, y, yerr=array([y_error0,y_error]),alpha_fill=0.1, ls='dashed', lw=3,color='k')
        plt.plot(x, y, ls='None', marker='o', markersize=10, label=r'$E[Q_{\ell}-Q_{\ell-1}]$',color='k')
        plt.xticks(x)
        rate=-log(y[:-1]/y[1:])/log(2.0)
        for i,v in enumerate(rate):
            if i==0:
                continue
            plt.annotate("%.1f"%v,((x[i]+x[i+1])/2,sqrt(y[i]*y[i+1])),color='k')
    #       plt.title('Sample mean')
        plt.xlabel(r'level $\ell$')
        plt.ylabel(r'$E[Q_{\ell}-Q_{\ell-1}]$')
        plt.grid(True)
        plt.yscale('log')
        plt.savefig(basename+'_'+str(nv)+'_mean.'+fileformat,transparent=True, bbox_inches='tight', pad_inches=0.2)
        plt.close(fig)

        # ---- PLOT SECOND MOMENT OF THE PAIRS PER EACH LEVEL
        fig=plt.figure()
        plt.gca().xaxis.set_major_locator(plt.MaxNLocator(integer=True))
        plt.hold(True)
        y=abs(array(dq2))[:,nv]
        p2, = plt.plot(x, y, clip_on=False, ls='dashed', lw=3, marker='o', markersize=10, label=r'$E[(Q_{\ell}-Q_{\ell-1})^2]$',color='b')
        plt.xticks(x)
        rate=-log(y[:-1]/y[1:])/log(2.0)
        for i,v in enumerate(rate):
            if i==0:
                continue
            plt.annotate("%.1f"%v,((x[i]+x[i+1])/2,sqrt(y[i]*y[i+1])),color='b')
    #       plt.title('Sample mean')
        plt.xlabel(r'level $\ell$')
        plt.ylabel(r'$E[(Q_{\ell}-Q_{\ell-1})^2]$')
        plt.grid(True)
        plt.yscale('log')
        plt.savefig(basename+'_'+str(nv)+'_mean2.'+fileformat,transparent=True, bbox_inches='tight', pad_inches=0.2)
        plt.close(fig)

        # ---- PLOT MEAN SQUARE OF THE PAIRS PER EACH LEVEL
        fig=plt.figure()
        plt.gca().xaxis.set_major_locator(plt.MaxNLocator(integer=True))
        plt.hold(True)
        y=abs(array(dq)**2)[:,nv]
        p2, = plt.plot(x, y, clip_on=False, ls='dashed', lw=3, marker='o', markersize=10, label=r'$E[Q_{\ell}-Q_{\ell-1}]^2$',color='c')
        plt.xticks(x)
        rate=-log(y[:-1]/y[1:])/log(2.0)
        for i,v in enumerate(rate):
            if i==0:
                continue
            plt.annotate("%.1f"%v,((x[i]+x[i+1])/2,sqrt(y[i]*y[i+1])),color='c')
    #       plt.title('Sample mean')
        plt.xlabel(r'level $\ell$')
        plt.ylabel(r'$E[(Q_{\ell}-Q_{\ell-1})^2]$')
        plt.grid(True)
        plt.yscale('log')
        plt.savefig(basename+'_'+str(nv)+'_mean2e.'+fileformat,transparent=True, bbox_inches='tight', pad_inches=0.2)
        plt.close(fig)

        # ---- PLOT EXPECTED VALUE OF THE  ABS OF THE PAIRS PER EACH LEVEL
        fig=plt.figure()
        plt.gca().xaxis.set_major_locator(plt.MaxNLocator(integer=True))
        plt.hold(True)
        y=abs(array(dq_abs))[:,nv]
        plt.plot(x,y,clip_on=False,color='k')
        plt.xticks(x)
        rate=-log(y[:-1]/y[1:])/log(2.0)
        for i,v in enumerate(rate):
            if i==0:
                continue
            plt.annotate("%.1f"%v,((x[i]+x[i+1])/2,sqrt(y[i]*y[i+1])),color='k')
    ##      plt.title('Sample mean ABS')
        plt.xlabel(r'$\ell$')
        plt.ylabel(r'$E[Q_\ell-Q_{\ell-1}]$')
        plt.grid(True)
        plt.yscale('log')
        plt.savefig(basename+'_'+str(nv)+'_mean_abs.'+fileformat,transparent=True, bbox_inches='tight', pad_inches=0.2)
        plt.close(fig)

        # ---- PLOT EXPECTED VALUE PER EACH LEVEL
        fig=plt.figure()
        plt.gca().xaxis.set_major_locator(plt.MaxNLocator(integer=True))
        plt.hold(True)
    #       plt.title('Cumulative sample mean')
        plt.xlabel(r'level $\ell$')
        plt.ylabel(r'$E[Q_\ell]$')
        plt.grid(True)
        plt.xticks(x)
        y=array(cumsum)[:,nv]
        plt.yticks(y)
        #plt.yscale('log')
        p3, = plt.plot(x,y,clip_on=False, label=r'$E[Q_\ell]$',color='y')
        plt.savefig(basename+'_'+str(nv)+'_cumsum.'+fileformat,transparent=True, bbox_inches='tight', pad_inches=0.2)
        plt.close(fig)

        # ----- PLOT VARIANCES PER EACH LEVEL
        fig=plt.figure()
        plt.gca().xaxis.set_major_locator(plt.MaxNLocator(integer=True))
        plt.hold(True)
    #       plt.title('Sample variance')
        y=array(dqvar)[:,nv]
        plt.axis([min(x)-.1,max(x)+.1,min(y)*0.8,max(y)*1.2])
        y_error=array(dqvar_err)[:,nv]
        y_error0=minimum(y_error,y-1e-13)
        special.errorfill(x, y, yerr=array([y_error0,y_error]),alpha_fill=0.1, ls='dashed', lw=3,color='r')
        plt.plot(x, y, ls='None', marker='o', markersize=10, label=r'$V[Q_\ell-Q_{\ell-1}]$',color='r')
        rate=-log(y[:-1]/y[1:])/log(2.0)
        for i,v in enumerate(rate):
            if i==0:
                continue
            plt.annotate("%.1f"%v,((x[i]+x[i+1])/2,sqrt(y[i]*y[i+1])),color='r')
        plt.xlabel(r'level $\ell$')
        plt.ylabel(r'$V[Q_\ell-Q_{\ell-1}]$')
        plt.grid(True)
        plt.yscale('log')
        plt.savefig(basename+'_'+str(nv)+'_var.'+fileformat,transparent=True, bbox_inches='tight', pad_inches=0.2)
        plt.close(fig)

        ## ---- PLOT NORMALIZED VARIANCES PER EACH LEVEL
        fig=plt.figure()
        plt.gca().xaxis.set_major_locator(plt.MaxNLocator(integer=True))
        plt.hold(True)
        y=sqrt(array(varnorm[nv,:]))
        #plt.title('Normalized sample variance')
        plt.xticks(x)
        #plt.plot(x,y,clip_on=False,label="proposed")
        #rate=-log(y[:-1]/y[1:])/log(2.0)
        #for i,v in enumerate(rate):
            #if i==0:
                #continue
        #       plt.annotate("%.1f"%v,((x[i]+x[i+1])/2,sqrt(y[i]*y[i+1])))
        y=sqrt(array(varnorm_eff[nv,:]))
        plt.plot(x,y,clip_on=False,color='m')
        plt.legend(loc=3,fontsize=10)
        rate=-log(y[:-1]/y[1:])/log(2.0)
        for i,v in enumerate(rate):
            if i==0:
                continue
            plt.annotate("%.1f"%v,((x[i]+x[i+1])/2,sqrt(y[i]*y[i+1])),color='m')
        plt.xlabel(r'level $\ell$')
        plt.ylabel(r'$\sqrt{V[Q_\ell-Q_{\ell-1}]}/\sqrt{M_\ell}$')
        plt.grid(True)
        plt.yscale('log')
        plt.savefig(basename+'_'+str(nv)+'_varnorm.'+fileformat,transparent=True, bbox_inches='tight', pad_inches=0.2)
        plt.close(fig)

        if ('w' in vars() and len(w)==len(q)):
            # ---- PLOT AVERAGE WORK LOAD PER LEVEL PER SINGLE RUN
            fig=plt.figure()
            plt.hold(True)
            y=array(w)
            plt.plot(x,y,clip_on=False,color='g')
            plt.xticks(x)
            plt.gca().xaxis.set_major_locator(plt.MaxNLocator(integer=True))
    #               plt.title('Mean work per level')
            rate=-log(y[:-1]/y[1:])/log(2.0)
            for i,v in enumerate(rate):
                plt.annotate("%.1f"%v,((x[i]+x[i+1])/2,sqrt(y[i]*y[i+1])),color='g')
            plt.xlabel(r'level $\ell$')
            plt.ylabel(r'$Work$')
            plt.grid(True)
            plt.yscale('log')
            plt.savefig(basename+'_work.'+fileformat,transparent=True, bbox_inches='tight', pad_inches=0.2)
            plt.close(fig)

            # ---- PLOT AVERAGE WORK TIMES VARIANCE
            fig=plt.figure()
            plt.gca().xaxis.set_major_locator(plt.MaxNLocator(integer=True))
            plt.hold(True)
            y=sqrt(array(w)*array(dqvar)[:,nv])
            plt.plot(x,y,clip_on=False,color='g')
            rate=-log(y[:-1]/y[1:])/log(2.0)
            for i,v in enumerate(rate):
                plt.annotate("%.1f"%v,((x[i]+x[i+1])/2,sqrt(y[i]*y[i+1])),color='g')
            plt.xticks(x)
##                      plt.title('Mean work times variance per level')
            plt.xlabel(r'level $\ell$')
            plt.ylabel(r'$\sqrt{W_{\ell}V_{\ell}}$')
            plt.grid(True)
            plt.yscale('log')
            plt.savefig(basename+'_'+str(nv)+'_workvariance.'+fileformat,transparent=True, bbox_inches='tight', pad_inches=0.2)
            plt.close(fig)

            ## ---- PLOT OPTIMAL Ml
            #fig=plt.figure()
            #plt.gca().xaxis.set_major_locator(plt.MaxNLocator(integer=True))
            #plt.hold(True)
            #y=(calpha/(error_split*final_tol*QQ[nv]))**2*sqrt(transpose(dqvar)/array(w))[nv,:]*sum(sqrt(array(w)*array(dqvar)[:,nv]))
            #plt.plot(x,y,clip_on=False,color='k')
            #rate=-log(y[:-1]/y[1:])/log(2.0)
            #for i,v in enumerate(rate):
                #plt.annotate("%.1f"%v,((x[i]+x[i+1])/2,sqrt(y[i]*y[i+1])),color='k')
            ##y=array(Ml_eff)
            ##plt.plot(x,y,clip_on=False,label="effective")
            ##rate=-log(y[:-1]/y[1:])/log(2.0)
            ##for i,v in enumerate(rate):
                ##plt.annotate("%.1f"%v,(sqrt(x[i]+x[i+1])/2,sqrt(y[i]*y[i+1])))
            #plt.xticks(x)
            #plt.xlabel(r'level $\ell$')
            #plt.ylabel(r'$M_\ell$')
            #plt.grid(True)
            ##plt.legend(loc=3,fontsize=10)
            #plt.yscale('log')
            #plt.savefig(basename+'_'+str(nv)+'_optimalMl.'+fileformat,transparent=True, bbox_inches='tight', pad_inches=0.2)
            #plt.close(fig)

    #plt.show()
    plt.close('all')
# -----------------------

## ------- PLOT VARIANCES FOR FIXED LEVEL AND DIFFERENT ITERATIONS
#def makeplot_iterations(q,fileformat='png'):
    #plt.close('all')
    #niter=[max([len(x) for x in y]) for y in q]
    #fig=plt.figure()
    #plt.grid(True)
    #plt.hold(True)
    #for l in range(len(q)):
        #mm=array([mean([x[int(y*1.2)]-x[int(y)] for x in q[l] if y<int(len(x)/1.2)]) for y in range(int(niter[l]/1.2)-1)])
        #mm=abs(mm)
        #x=array(range(niter[l]))[-len(mm):]
        #rate=-log(mm[:-5]/mm[5:])/log(x[:-5]/x[5:])
        #for i,v in enumerate(rate[::5]):
            #plt.annotate("%.1f"%v,(x[::5][i],mm[::5][i]),fontsize=10)
        #plt.plot(x,mm,'--', label="Mean, l="+str(l))
        #vv=array([samplevar([x[int(y*1.2)]-x[int(y)] for x in q[l] if y<int(len(x)/1.2)]) for y in range(int(niter[l]/1.2)-1)])
        #x=array(range(niter[l]))[-len(vv):]
        #rate=-log(vv[:-5]/vv[5:])/log(x[:-5]/x[5:])
        #for i,v in enumerate(rate[::3]):
            #plt.annotate("%.1f"%v,(x[::3][i],vv[::3][i]),fontsize=10)
        #plt.plot(x,vv,'-o', label="Var, l="+str(l))
    #plt.yscale('log')
    #plt.xscale('log')
    #plt.xlabel('iteration i')
    #plt.ylabel(r'$V[Q_{2i}-Q_{i}],E[Q_{2i}-Q_{i}]$')
    #plt.xlim(left=20,right=120)
    #plt.ylim(bottom=1e-2,top=1e4)
    #plt.legend(loc=3,fontsize=10)
    #plt.savefig(basename+'_var_iterations_x2.'+fileformat,transparent=True, pad_inches=0.4)
    #plt.close(fig)
    #fig=plt.figure()
    #plt.grid(True)
    #plt.hold(True)
    #for l in range(len(q)):
        #mm=array([mean([x[y]-x[y-1] for x in q[l] if y<len(x)-2]) for y in range(niter[l]-2)])
        #mm=abs(mm)
        #x=array(range(niter[l]))[-len(mm):]
        #rate=-log(mm[:-3]/mm[3:])/log(x[:-3]/x[3:])
        #for i,v in enumerate(rate[::3]):
            #plt.annotate("%.1f"%v,(x[::3][i],mm[::3][i]),fontsize=10)
        #plt.plot(x,mm,'--', label="Mean, l="+str(l))
        #vv=array([samplevar([x[y]-x[y-1] for x in q[l] if y<len(x)-2]) for y in range(niter[l]-2)])
        #x=array(range(niter[l]))[-len(vv):]
        #rate=-log(vv[:-5]/vv[5:])/log(x[:-5]/x[5:])
        #for i,v in enumerate(rate[::5]):
            #plt.annotate("%.1f"%v,(x[::5][i],vv[::5][i]),fontsize=10)
        #plt.plot(x,vv,'-o', label="Var, l="+str(l))
    #plt.yscale('log')
    #plt.xscale('log')
    #plt.xlabel('iteration i')
    #plt.ylabel(r'$V[Q_{i}-Q_{i-1}],E[Q_{i}-Q_{i-1}]$')
    #plt.xlim(left=20,right=120)
    #plt.ylim(bottom=1e-4,top=1e2)
    #plt.legend(loc=3,fontsize=10)
    #plt.savefig(basename+'_var_iterations_1.'+fileformat,transparent=True, pad_inches=0.4)
    #plt.close(fig)
    ##plt.show()
    #plt.close('all')
## -----------------------

class postprocess(estimator):
    def __init__(self, name):
        self.name=name
        # setup the first guess
        self.final_tol=final_tol
        self.error_split=error_split
        (q,w)=load_qw(name)
        self.q=q
        self.w=w
        self.tot_level=len(w)
        self.ml=[len(y) for y in q]
        self.cprod=array([[1.] for i in range(tot_level)])
        self.csum=array([[0.] for i in range(tot_level)])
        self.compute_errors()






# ------------------------------------------------- END FILE
#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.
