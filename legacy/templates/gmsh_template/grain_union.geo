Include 'mesh_size.geo';

pores = { };

// Instead of using included files, we now use a user-defined function
// in order to carve some holes in the cube:
Function CheeseHole 

  // In the following commands we use the reserved variable name
  // `newp', which automatically selects a new point number. This
  // number is chosen as the highest current point number, plus
  // one. (Note that, analogously to `newp', the variables `newc',
  // `news', `newv' and `newreg' select the highest number amongst
  // currently defined curves, surfaces, volumes and `any entities
  // other than points', respectively.)

  p1 = newp; Point(p1) = {x,  y,  z,  ref3} ;
  p2 = newp; Point(p2) = {x+r1,y,  z,  ref3} ;
  p3 = newp; Point(p3) = {x,  y+r2,z,  ref3} ;
  p4 = newp; Point(p4) = {x,  y,  z+r3,ref3} ;
  p5 = newp; Point(p5) = {x-r1,y,  z,  ref3} ;
  p6 = newp; Point(p6) = {x,  y-r2,z,  ref3} ;
  p7 = newp; Point(p7) = {x,  y,  z-r3,ref3} ;
  p8 = newp; Point(p8) = {x+r1/2, y, z, ref3} ;
  p9 = newp; Point(p9) = {x, y+r2/2, z, ref3} ;
  p10 = newp; Point(p10) = {x, y, z+r3/2, ref3} ;

  c1 = newreg; Ellipse(c1) = {p2,p1,p8,p7};
  c2 = newreg; Ellipse(c2) = {p7,p1,p8,p5};
  c3 = newreg; Ellipse(c3) = {p5,p1,p8,p4};
  c4 = newreg; Ellipse(c4) = {p4,p1,p8,p2};
  c5 = newreg; Ellipse(c5) = {p2,p1,p8,p3};
  c6 = newreg; Ellipse(c6) = {p3,p1,p8,p5};
  c7 = newreg; Ellipse(c7) = {p5,p1,p8,p6};
  c8 = newreg; Ellipse(c8) = {p6,p1,p8,p2};
  c9 = newreg; Ellipse(c9) = {p7,p1,p9,p3};
  c10 = newreg; Ellipse(c10) = {p3,p1,p9,p4};
  c11 = newreg; Ellipse(c11) = {p4,p1,p9,p6};
  c12 = newreg; Ellipse(c12) = {p6,p1,p9,p7};

  // We need non-plane surfaces to define the spherical holes. Here we
  // use ruled surfaces, which can have 3 or 4 sides:

  l1 = newreg; Line Loop(l1) = {c5,c10,c4};   Ruled Surface(newreg) = {l1};
  l2 = newreg; Line Loop(l2) = {c9,-c5,c1};   Ruled Surface(newreg) = {l2};
  l3 = newreg; Line Loop(l3) = {c12,-c8,-c1}; Ruled Surface(newreg) = {l3};
  l4 = newreg; Line Loop(l4) = {c8,-c4,c11};  Ruled Surface(newreg) = {l4};
  l5 = newreg; Line Loop(l5) = {-c10,c6,c3};  Ruled Surface(newreg) = {l5};
  l6 = newreg; Line Loop(l6) = {-c11,-c3,c7}; Ruled Surface(newreg) = {l6};
  l7 = newreg; Line Loop(l7) = {-c2,-c7,-c12};Ruled Surface(newreg) = {l7};
  l8 = newreg; Line Loop(l8) = {-c6,-c9,c2};  Ruled Surface(newreg) = {l8};

  // We then store the surface loops identification numbers in a list
  // for later reference (we will need these to define the final
  // volume):

  theloops[t] = newreg ; 

  Rotate {  {a1, a2, a3}, {x, y, z}, th}  {Surface{l8+1,l5+1,l1+1,l2+1,l3+1,l7+1,l6+1,l4+1};}

  Surface Loop(theloops[t]) = {l8+1,l5+1,l1+1,l2+1,l3+1,l7+1,l6+1,l4+1};

  Physical Surface(t) = {l8+1,l5+1,l1+1,l2+1,l3+1,l7+1,l6+1,l4+1};
  thehole = newreg ; 
  Volume(thehole) = theloops[t];
  Physical Volume(t) = thehole;

Return

// We can use a `For' loop to generate five holes in the cube:

Include 'grains.geo';

For t In {1:ngrains}

  x = xx[t-1] ;
  y = yy[t-1] ;
  z = zz[t-1] ;
  r1 = rr1[t-1] ;
  r2 = rr2[t-1] ;
  r3 = rr3[t-1] ;
  a1 = ax1[t-1] ;
  a2 = ax2[t-1] ;
  a3 = ax3[t-1] ;
  th = tx[t-1] ;
  Call CheeseHole ;

  //Physical Surface(t) = theloops[t];
  //Physical Volume(t) = thehole ;

EndFor

// We finally define a physical volumes and surfaces

//Physical Surface(t+1) = {pores[]};
