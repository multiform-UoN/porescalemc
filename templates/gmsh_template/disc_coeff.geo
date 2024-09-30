Include 'mesh_size.geo';
Include 'box.geo';

Point(1) = {x2, y2, z2, ref1};
Point(2) = {x1, y1, z1, ref1};
Point(3) = {x2, y1, z1, ref1};
Point(4) = {x2, y2, z1, ref1};
Point(5) = {x1, y2, z1, ref1};
Point(6) = {x2, y1, z2, ref1};
Point(7) = {x1, y2, z2, ref1};
Point(8) = {x1, y1, z2, ref1};

Line(4) = {3, 2};
Line(5) = {2, 8};
Line(6) = {8, 6};
Line(7) = {6, 3};

Line(8) = {1, 4};
Line(9) = {4, 5};
Line(10) = {5, 7};
Line(11) = {7, 1};

Line Loop(17) = {4,5,6,7};
Line Loop(19) = {8,9,10,11};
Plane Surface(18) = {17};
Plane Surface(21) = {19};


Line(22) = {3, 4};
Line(23) = {6, 1};
Line(24) = {8, 7};
Line(25) = {2, 5};
Line Loop(26) = {7, 22, -8, -23};
Plane Surface(27) = {26};
Line Loop(28) = {23, -11, -24, 6};
Plane Surface(29) = {28};
Line Loop(30) = {24, -10, -25, 5};
Plane Surface(31) = {30};
Line Loop(32) = {22, 9, -25, -4};
Plane Surface(33) = {32};

//Surface Loop(34) = {33, 27, 18, 31, 29, 21};

//Periodic Surface 27 {26} = 31 {30};
theloops[0] = newreg ;

//Volume(35) = {34, 2};
//Physical Volume(43) = {35};

Surface Loop(theloops[0]) = {33, 27, 18, 31, 29, 21} ;

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

  Rotate {  {Rand(1), Rand(1), Rand(1)}, {x, y, z}, Rand(1)}  {Surface{l8+1,l5+1,l1+1,l2+1,l3+1,l7+1,l6+1,l4+1};}

  Surface Loop(theloops[t]) = {l8+1,l5+1,l1+1,l2+1,l3+1,l7+1,l6+1,l4+1};

  pores[] += {l8+1,l5+1,l1+1,l2+1,l3+1,l7+1,l6+1,l4+1};
  thehole[t-1] = newreg ; 
  Volume(thehole[t-1]) = theloops[t] ;

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
  Call CheeseHole ;

  // We define a physical volume for each hole:

  //Physical Volume(t) = thehole ;

  Printf("Hole %g (center = {%g,%g,%g,%g,%g}, radius = %g) has number %g!",
	 t, x, y, z, r1, r2, r3, thehole) ;

EndFor

// The volume of the cube, without the 5 holes, is now defined by 6
// surface loops (the exterior surface and the five interior loops).
// To reference an array of variables, its identifier is followed by
// '[]':

thevol = newreg;
Volume(thevol) = {theloops[]} ;

// We finally define a physical volumes and surfaces

Physical Surface(1) = 31;
Physical Surface(2) = 27;
Physical Surface(3) = 29;
Physical Surface(4) = 18;
Physical Surface(5) = 21;
Physical Surface(6) = 33;
Physical Surface(7) = {pores[]};
Physical Volume (8) = thevol ;
Physical Volume (9) = {thehole[]} ;
