from bsand import *
g=geom(name="cyl")
#g.sample()
from openfoam import *
s=solver(g)
s.setup(g)
s.solve()
