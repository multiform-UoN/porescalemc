
P2 = 0; // set to 1 for second order elements
V2 = 1; // set to 1 for second order elements

Group {
  Omega = Region[ 8 ];
  Inlet = Region[ 1 ];
  Outlet = Region[ 2 ];
  InletOutlet = Region[{1 2}];
  Lateral = Region[ {3 4 5 6} ];
  Grains = Region[ 7 ];
  DomainV = Region[ { Omega, Inlet, Outlet } ];
  DomainP = Region[ { Omega, Lateral, Grains, Inlet, Outlet } ];
  AllDomain = Region[ {1 2 3 4 5 6 7 8} ];
  AllBoundaries = Region[{1 2 3 4 5 6 7}];
}

Function {
  g[] = 0;
  rho[] = 1;
  visc[] = 1;
  // f[] = X[] > 0 ? Sin[X[]] : 0;
}

Include "DirBCs"

Constraint {
  { Name pres; Type Assign;
    Case {
      { Region Inlet; Value inletp; } // PRESSURE INFLOW
      { Region Outlet; Value outletp; } // PRESSURE OUTFLOW
    }
  }
  { Name velx; Type Assign;
    Case {
//      { Region AllBoundaries; Value 0; }
      { Region Lateral; Value 0; }
      { Region Grains; Value 0; }
//      { Region Inlet; Value 0; }  // VELOCITY INFLOW
//      { Region Outlet; Value 0; }  // VELOCITY INFLOW
    }
  }
  { Name vely; Type Assign;
    Case {
//      { Region AllBoundaries; Value 0; }
      { Region Lateral; Value 0; }
      { Region Grains; Value 0; }
      { Region Inlet; Value 0; }  // VELOCITY INFLOW
      { Region Outlet; Value 0; }  // VELOCITY INFLOW
    }
  }
  { Name velz; Type Assign;
    Case {
//      { Region AllBoundaries; Value 0; }
      { Region Lateral; Value 0; }
      { Region Grains; Value 0; }
      { Region Inlet; Value 0; }  // VELOCITY INFLOW
      { Region Outlet; Value 0; }  // VELOCITY INFLOW
    }
  }
}

Jacobian {
  { Name JVol;
    Case { { Region AllDomain; Jacobian Vol; } }
  }
  { Name JSur;
    Case { { Region AllBoundaries; Jacobian Sur; } }
  }
}

Integration {
  { Name I1;
    Case { { Type Gauss; Case {
          { GeoElement Point ; NumberOfPoints  1 ; }
          { GeoElement Line ; NumberOfPoints  4 ; }
          { GeoElement Triangle ; NumberOfPoints  6 ; }
          { GeoElement Quadrangle ; NumberOfPoints 7 ; }
          { GeoElement Tetrahedron ; NumberOfPoints 15 ; }
          { GeoElement Hexahedron ; NumberOfPoints 34 ; }
          { GeoElement Prism       ; NumberOfPoints 9 ; }
  } } }
  }
}

FunctionSpace {
  { Name H_u ; Type Vector ;
    BasisFunction {
      { Name sxn ; NameOfCoef uxn ; Function BF_NodeX ;
        Support DomainV ; Entity NodesOf[ All ] ; }
      { Name syn ; NameOfCoef uyn ; Function BF_NodeY ;
        Support DomainV ; Entity NodesOf[ All ] ; }
      { Name szn ; NameOfCoef uzn ; Function BF_NodeZ ;
        Support DomainV ; Entity NodesOf[ All ] ; }
      If(V2)
      { Name sx2e ; NameOfCoef ux2e ; Function BF_NodeX_2E ;
        Support DomainV ; Entity NodesOf[ All ] ; }
      { Name sy2e ; NameOfCoef uy2e ; Function BF_NodeY_2E ;
        Support DomainV ; Entity NodesOf[ All ] ; }
      { Name sz2e ; NameOfCoef uz2e ; Function BF_NodeZ_2E ;
        Support DomainV ; Entity NodesOf[ All ] ; }
      EndIf
    }
    Constraint {
      { NameOfCoef uxn ;
        EntityType NodesOf ; NameOfConstraint velx ; }
      { NameOfCoef uyn ;
        EntityType NodesOf ; NameOfConstraint vely ; }
      { NameOfCoef uzn ;
        EntityType NodesOf ; NameOfConstraint velz ; }
      If(V2)
      { NameOfCoef ux2e ;
        EntityType NodesOf ; NameOfConstraint velx ; }
      { NameOfCoef uy2e ;
        EntityType NodesOf ; NameOfConstraint vely ; }
      { NameOfCoef uz2e ;
        EntityType NodesOf ; NameOfConstraint velz ; }
      EndIf
    }
  }
  { Name H_p; Type Form0; 
    BasisFunction {
      { Name sn1n; NameOfCoef wn1n; Function BF_Node; Support DomainP; Entity NodesOf[All]; }
      If(P2)
      { Name sn2e; NameOfCoef wn2e; Function BF_Node_2E; Support DomainP; Entity EdgesOf[All]; }
      EndIf
    }
    Constraint {
      { NameOfCoef wn1n; EntityType NodesOf; NameOfConstraint pres; }
      If(P2)
      { NameOfCoef wn2e; EntityType EdgesOf; NameOfConstraint pres; }
       EndIf
    }
  }
}

Formulation {
  { Name Pde; Type FemEquation; 
    Quantity { 
      { Name p; Type Local; NameOfSpace H_p; }
      { Name u; Type Local; NameOfSpace H_u; }
    }
    Equation {
       Galerkin { [ visc[] * Dof{D1 u} , {D1 u} ] ;
                  In Omega ; Jacobian JVol ; Integration I1 ; } 
       Galerkin { [ 0.5*visc[] * Dof{D2 u} , {D2 u} ] ;
                  In Omega ; Jacobian JVol ; Integration I1 ; } 
       Galerkin { [ 2*visc[] * CompX[Dof{D1 u}] , CompX[{u}] ] ;
                  In Inlet ; Jacobian JSur ; Integration I1 ; } 
//       Galerkin { [ visc[] * CompX[Dof{D2 u}] , CompY[{u}] ] ;
//                  In Inlet ; Jacobian JSur ; Integration I1 ; } 
//       Galerkin { [ visc[] * CompZ[Dof{D2 u}] , CompZ[{u}] ] ;
//                  In Inlet ; Jacobian JSur ; Integration I1 ; } 
       Galerkin { [ -2*visc[] * Dof{D1 u} , {u} ] ;
                  In Outlet ; Jacobian JSur ; Integration I1 ; } 
//       Galerkin { [ -visc[] * CompX[Dof{D2 u}] , CompY[{u}] ] ;
//                  In Outlet ; Jacobian JSur ; Integration I1 ; } 
//       Galerkin { [ -visc[] * CompZ[Dof{D2 u}] , CompZ[{u}] ] ;
//                  In Outlet ; Jacobian JSur ; Integration I1 ; } 
// PRESSURE TERM
      Galerkin { [ -Dof{p} , {Div u} ]; 
                 In Omega; Integration I1; Jacobian JVol;  }
      Galerkin { [ -Dof{p} , CompX[{u}] ]; 
                 In Inlet; Integration I1; Jacobian JSur;  }
      Galerkin { [ Dof{p} , CompX[{u}] ]; 
                 In Outlet; Integration I1; Jacobian JSur;  }
// STABILIZATION
      Galerkin { [ 0.00001*Dof{ Grad p} , {Grad p} ]; 
                 In Omega; Integration I1; Jacobian JVol;  }
// DIVERGENCE FREE CONDITION
      Galerkin { [ {p} , Dof{Div u} ]; 
                 In Omega; Integration I1; Jacobian JVol;  }
//      Galerkin { [ {Grad p} , Dof{u} ]; 
//                 In Omega; Integration I1; Jacobian JVol;  }
//      Galerkin { [ {p} ,  CompX[Dof{u}] ]; 
//                 In Inlet; Integration I1; Jacobian JSur;  }
//      Galerkin { [ -{p} ,  CompX[Dof{u}] ]; 
//                 In Outlet; Integration I1; Jacobian JSur;  }
//      Galerkin { [ Vector[1.0,0.0,0.0] , {u} ];     // FIXED PRESSURE GRADIENT
//                 In Omega; Integration I1; Jacobian JVol;  }
// SOURCE TERM (GRAVITY)
      Galerkin { [ -Vector[1.0,0.0,0.0]*rho[]*g[] , {u} ]; 
                 In Omega; Integration I1; Jacobian JVol;  }
    }
  }
}

Resolution {
  { Name Pde;
    System {
      { Name A; NameOfFormulation Pde; }
    }
    Operation { Generate[A]; Solve[A]; SaveSolution[A]; }
  }
}

PostProcessing {
  { Name Pde; NameOfFormulation Pde;
    Quantity {
      { Name p; Value{ Local{ [ {p} ]; In AllDomain; Jacobian JVol; } } }
      { Name U; Value{ Local{ [ CompX[{ u }]  ]; Integration I1; In AllDomain; Jacobian JVol;} } }
      { Name V; Value{ Local{ [ CompY[{ u }]  ]; Integration I1; In AllDomain; Jacobian JVol;} } }
      { Name W; Value{ Local{ [ CompZ[{ u }]  ]; Integration I1; In AllDomain; Jacobian JVol;} } }
      { Name velmag; Value{ Local{ [ Norm[{ u }]  ]; Integration I1; In AllDomain; Jacobian JVol;} } }
      { Name divu; Value{ Local{ [ { Div u }  ]; Integration I1; In AllDomain; Jacobian JVol;} } }
      { Name IntFlux; Value{ Integral{ [ CompX[{ u }]  ]; Integration I1; In AllDomain; Jacobian JVol;} } }
    }
  }
}

PostOperation {
  { Name Pde; NameOfPostProcessing Pde;
    Operation {
      Print[ p , OnElementsOf Omega , File "p.pos", Depth 1];
      Print[ U , OnElementsOf Omega , File "U.pos", Depth 1];
      Print[ V , OnElementsOf Omega , File "V.pos", Depth 1];
      Print[ W , OnElementsOf Omega , File "W.pos", Depth 1];
      Print[ velmag , OnElementsOf Omega , File "velmag.pos", Depth 1];
      Print[ divu , OnElementsOf Omega , File "divu.pos", Depth 1];
      Print[ IntFlux[Omega] , OnGlobal , File "IntFlux.dat", Depth 0, Format Table];
    }
  }
}
