(* ::Package:: *)

BeginPackage["Emphasis`SpecialTensors`"];
(*Introduced in version 1.0*) 
Needs["Emphasis`IndexGymnastics`"];

Kummer::usage = "calculate the Kummer tensor from a (4,0)-tensor whose first two and second two indices are antisymmetric. 
The Kummer tensor is given by equations (45)-(46) in P. Baekler et al., Annals of Physics (NY), 349 (2014) 297\[Dash]324.
Example: 
(*Tensor based on Hodge star of Lorentzian metric*)
h=DiagonalMatrix[{-1,1,1,1}];
w=TensorProduct[h,h];
x=Transpose[w,{1,3,2,4}];
y=2*Normal[Symmetrize[x,Antisymmetric[{1,2}]]];
r=f*y/Sqrt[-Det[h]];
(*Ensuing Kummer tensor*)
t=Kummer[r]; 
(*Check*) 
u=-(f^3)*(7*w-Transpose[w,{1,4,3,2}])/Sqrt[-Det[h]]; 
Union[Flatten[FullSimplify[t-u]]]"

TammRubilar::usage = "calculate the Tamm-Rubilar tensor from a (4,0)-tensor whose first two and second two indices are antisymmetric. 
The Tamm-Rubilar tensor is given by equations (50)-(51) in P. Baekler et al., Annals of Physics (NY), 349 (2014) 297\[Dash]324.
Example:
(*Tensor based on Hodge star of Lorentzian metric*)
h=DiagonalMatrix[{-1,1,1,1}];
w=TensorProduct[h,h];
x=Transpose[w,{1,3,2,4}];
y=2*Normal[Symmetrize[x,Antisymmetric[{1,2}]]];
r=f*y/Sqrt[-Det[h]];
(*Ensuing Tamm-Rubilar tensor*)
t=TammRubilar[r]; 
(*Check*)
u=-(f^3)*Normal[Symmetrize[w]]/Sqrt[-Det[h]];
Union[Flatten[FullSimplify[t-u]]]"

Bel::usage ="Bel[r,g,h] computes the Bel tensor from the Riemann curvature r of a manifold with n dimensions, using the Lorentzian metric g and its inverse h. 
In particular, the tensors r and g are assumed to have type (0,4) and (0,2) respectively. 
Moreover, the Bel tensor is defined as in equation (30) of J.M.M. Senovilla, Class. Quantum Grav., 17 (2000) 2799.
Example:
(*Build Riemann tensor of general relativity*)
w=Table[Subscript[ww,i,j,k,l],{i,4},{j,4},{k,4},{l,4}];
(*Demand first two and second two indices are antisymmetric*)
x=Symmetrize[w,{Cycles[{{1,2}}],-1}]; 
y=Symmetrize[x,{Cycles[{{3,4}}],-1}];  
(*Demand symmetry under exchange of antisymmetric index pairs*)
z=Symmetrize[y,{Cycles[{{1,3},{2,4}}],1}]; 
(*Demand totally antisymmetric part is zero*)
r=z-Symmetrize[z,Antisymmetric[{1,2,3,4}]];
(*Metric and its inverse*)
g=DiagonalMatrix[{-1,1,1,1}];
h=Inverse[g];
(*Ensuing Bel tensor*)
t=Bel[r,g,h];
(*Check symmetries*)
TensorSymmetry[t]"

Begin["Private`"]	
	(*Kummer tensor*)
	Kummer[r_]:=Module[{p,pr,prp,q,t},
	p=LeviCivitaTensor[4];
	pr=TwoContract[p,r,{{3,1},{4,2}}]/2;
	prp=TwoContract[pr,p,{{3,1},{4,2}}]/2;
	q=TwoContract[r,prp,{{1,1},{3,3}}];
	t=TwoContract[q,r,{{3,1},{4,3}}];
	t]

	(*Tamm-Rubilar tensor*)
	TammRubilar[r_]:=Symmetrize[Kummer[r]]/6 	
	
	(*Bel tensor*)
	Bel[r_,g_,h_]:=Module[{r24,r124,r234,r1234,w,wa,wb,x,y,z,t},
	r24=RaiseLower[r,{2,4},h];
	r124=RaiseLower[r24,{1},h];
	r234=RaiseLower[r24,{3},h];
	r1234=RaiseLower[r234,{1},h];
	w=TwoContract[r,r24,{{2,2},{4,4}}]; 
	wa=Transpose[w,{1,3,2,4}]; 
	wb=Transpose[w,{1,4,2,3}];
	x=-TensorProduct[g,TwoContract[r,r124,{{1,1},{2,2},{4,4}}]]/2;
	y=-TensorProduct[TwoContract[r,r234,{{2,2},{3,3},{4,4}}],g]/2;
	z=TwoContract[r,r1234,{{1,1},{2,2},{3,3},{4,4}}]*TensorProduct[g,g]/8; 
	t=wa+wb+x+y+z;
	t]
End[]
EndPackage[]
