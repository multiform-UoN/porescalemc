
P2 = 0; // set to 1 for second order elements

Group {
  Omega = Region[ {8 9} ];
  Gamma0 = Region[ {1 2 3 4 5 6 }];
  Spheres = Region[ 7 ];
  AllDomain = Region[ {7 8 9} ];
  Void     = Region[ 8 ];
  Balls     = Region[ 9 ];
}

Include "DirBCs"

Function {
f[AllDomain] = force;//*($X*$Y*$Z-1); //Norm[$1];
k[Void] = permeability;
k[Balls] = permeability2;
g[AllDomain] = flux;
}

Constraint {
  { Name u0;
    Case {
      { Region Gamma0; Value outletp; }
      { Region Spheres; Value inletp; }
    }
  }
}

Jacobian {
  { Name JVol;
    Case { { Region All; Jacobian Vol; } }
  }
  { Name JSur;
    Case { { Region All; Jacobian Sur; } }
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
  { Name Hgrad_u; Type Form0; 
    BasisFunction {
      { Name sn1n; NameOfCoef wn1n; Function BF_Node; Support Omega; Entity NodesOf[All]; }
      If(P2)
      { Name sn2e; NameOfCoef wn2e; Function BF_Node_2E; Support Omega; Entity EdgesOf[All]; }
      EndIf
    }
    Constraint {
      { NameOfCoef wn1n; EntityType NodesOf; NameOfConstraint u0; }
      If(P2)
      { NameOfCoef wn2e; EntityType EdgesOf; NameOfConstraint u0; }
      EndIf
    }
  }
}

Formulation {
  { Name Pde; Type FemEquation; 
    Quantity { 
      { Name u; Type Local; NameOfSpace Hgrad_u; }
    }
    Equation {
      Galerkin { [ k[]*Dof{Grad u} , {Grad u} ]; 
                 In Omega; Integration I1; Jacobian JVol;  }
      Galerkin { [ - f[] , {u} ]; 
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
      { Name u; Value{ Local{ [ {u} ]; In Omega; Jacobian JVol; } } }
      //{ Name f; Value{ Local{ [ f[] ] ; Integration I1; In Omega; Jacobian JVol;} } }
      //{ Name vel; Value{ Local{ [ Norm[{ Grad u }]  ]; Integration I1; In Omega; Jacobian JVol;} } }
      { Name IntFlux; Value{ Integral{ [ CompX[{ Grad u }]  ]; Integration I1; In Omega; Jacobian JVol;} } }
      { Name Int; Value{ Integral{ [ { u }  ]; Integration I1; In Omega; Jacobian JVol;} } }
      { Name IntFlux2; Value{ Integral{ [ CompX[{ Grad u }]^2  ]; Integration I1; In Omega; Jacobian JVol;} } }
      { Name Int2; Value{ Integral{ [ { u }^2  ]; Integration I1; In Omega; Jacobian JVol;} } }
      { Name Void; Value{ Integral{ [ 1  ]; Integration I1; In Omega; Jacobian JVol;} } }
    }
  }
}

PostOperation {
  { Name Pde; NameOfPostProcessing Pde;
    Operation {
      Print[ u , OnElementsOf Omega , File "u.pos", Depth P2*2+1];
      //Print[ f , OnElementsOf Omega , File "f.pos", Depth P2*2+1];
      //Print[ vel , OnElementsOf Omega , File "vel.pos", Depth P2*2+1];
      Print[ IntFlux2[Omega] , OnGlobal , File "IntFlux2.dat", Depth 0, Format Table];
      Print[ Int2[Omega] , OnGlobal , File "Int2.dat", Depth 0, Format Table];
      Print[ IntFlux[Omega] , OnGlobal , File "IntFlux.dat", Depth 0, Format Table];
      Print[ Int[Omega] , OnGlobal , File "Int.dat", Depth 0, Format Table];
      Print[ Void[Omega] , OnGlobal , File "Por.dat", Depth 0, Format Table];
    }
  }
}
