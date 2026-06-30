import multiprocessing as mp
from bsand import *
import time

def fun(id):
	t0=time.time()
	g=geom(name="bsand_only"+str(id))
	g.sample()
	g.write()
	print(g.porosity_out)
	print(time.time()-t0,len(g.g))

pp=[]
for i in range(5):
	p=mp.Process(target=fun,args=(i,))
	p.start()
	pp.append(p)

x=[p.join() for p in pp]
