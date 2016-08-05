(* ::Package:: *)

BeginPackage["Emphasis`TensorDecompositions`"];
(*Introduced in version 1.1*)

ShapeOfTableau::usage="ShapeOfTableau[tab] returns the shape of Young tableau tab as a list.
Example:
(*Step 1: Create tableau with shape {2,1}*)
tab = {{1,3},{2}};
(*Step 2: Recover shape of tableau*)
shp = ShapeOfTableau[tab];
(*Step 3: Check result*)
shp == {2,1}"

IrreduciblePart::usage="IrreduciblePart[ten,tab,opt] returns the irreducible part of the tensor ten as specified by the Young tableau tab.
-If the argument opt is omitted, or if opt is equal to 0, the irreducible part is calculated by 1) antisymmetrising the indices of ten according to the columns of tab, and 2) symmetrizing the indices of ten according to the rows of of tab. 
-If the argument opt is equal to 1, the irreducible part is calculated by applying 1) and 2) in reverse order.
Example: 
(*Find a specific irreducible part of the general (3,0)-tensor*)
tab ={{1,2},{3}};
gnr30 = Table[gnr[i,j,k],{i,1,5},{j,1,5},{k,1,5}];
(*****Convention: First antisymmetrize, then symmetrize*****)
(*Step 1: Make use of IrreduciblePart*)
ipt30a =IrreduciblePart[gnr30,tab,0];
(*Step 2: Do the calculation by hand*)
ref30a = ConstantArray[0,{5,5,5}];
Do[ref30a[[i,j,k]]=(gnr30[[i,j,k]]+gnr30[[j,i,k]]
-gnr30[[k,j,i]]-gnr30[[k,i,j]])/3;
,{i,1,5},{j,1,5},{k,1,5}];
(*Step 3: Compare results*)
Union[Flatten[FullSimplify[ipt30a-ref30a]]]
(*****Convention: First symmetrize, then antisymmetrize*****)
(*Step 1: Make use of IrreduciblePart*)
ipt30b =IrreduciblePart[gnr30,tab,1];
(*Step 2: Do the calculation by hand*)
ref30b = ConstantArray[0,{5,5,5}];
Do[ref30b[[i,j,k]]=(gnr30[[i,j,k]]+gnr30[[j,i,k]]
-gnr30[[k,j,i]]-gnr30[[j,k,i]])/3;
,{i,1,5},{j,1,5},{k,1,5}];
(*Step 3: Compare results*)
Union[Flatten[FullSimplify[ipt30b-ref30b]]]"

CanonicalPart::usage="CanonicalPart[ten,pat] returns the canonical part of the tensor ten as defined by the integer partition (Young diagram shape) pat.
Example: 
(*Decompose an electromagnetic medium tensor canonically*)
(*See F.W. Hehl and Y.N. Obukhov, Foundations of Classical Electrodynamics, Birkh\"auser, Boston, 2003*)
chimat = Table[chi[i,j],{i,1,6},{j,1,6}]; 
chi40 = MatrixToTensorA[chimat];
(*Step 1: Decompose using CanonicalPart*)
ord = ArrayDepth[chi40];
par = IntegerPartitions[ord];
chi40p1=FullSimplify[CanonicalPart[chi40,par[[1]]]]; 
chi40p2=FullSimplify[CanonicalPart[chi40,par[[2]]]]; 
chi40p3=FullSimplify[CanonicalPart[chi40,par[[3]]]]; 
chi40p4=FullSimplify[CanonicalPart[chi40,par[[4]]]]; 
chi40p5=FullSimplify[CanonicalPart[chi40,par[[5]]]]; 
(*Step 2: Decompose by hand*)
epsmat = 
{
{0,0,0,1,0,0},
{0,0,0,0,1,0},
{0,0,0,0,0,1},
{1,0,0,0,0,0},
{0,1,0,0,0,0},
{0,0,1,0,0,0}
};
chimatax = Tr[epsmat.chimat]*epsmat/6;
chimatsk= Symmetrize[chimat, Antisymmetric[{1,2}]];
chimatpr = Symmetrize[chimat, Symmetric[{1,2}]]-chimatax;
chi40ax = MatrixToTensorA[chimatax];
chi40sk = MatrixToTensorA[chimatsk];
chi40pr = MatrixToTensorA[chimatpr];
(*Step 3: Compare results*)
Union[Flatten[FullSimplify[chi40p3-chi40pr]]]
Union[Flatten[FullSimplify[chi40p4-chi40sk]]]
Union[Flatten[FullSimplify[chi40p5-chi40ax]]]
(*Step 4: Vanishing canonical parts*)
Union[Flatten[FullSimplify[chi40p1]]]
Union[Flatten[FullSimplify[chi40p2]]]"

Begin["Private`"]
Quiet@Block[{$ContextPath},Needs["Combinatorica`"]]

(*Symmetrise without normalisation*)
SymmetrizeNoNormalize[ten_,sym_]:=Module[{nrm,out},
nrm =Factorial[Length[sym[[1]]]];
out = nrm*Symmetrize[ten,sym];
out]

(*Shape of Young tableau*)
ShapeOfTableau[tab_]:=Map[Length,tab]

(*Irreducible part of tensor corresponding to given tableau*)
IrreduciblePart[ten_,tab_,opt_:0]:=Module[{SymFirst,SymSecond,tab1st,tab2nd,ten1st,ten2nd,ntx,nom},
If[opt==0,
(*Default convention*)
(*1st: Antisymmetrize over columns of tableau*)
SymFirst = Antisymmetric;
tab1st = Combinatorica`TransposeTableau[tab];
(*2nd: Symmetrize over rows of tableau*)
SymSecond = Symmetric;
tab2nd = tab;
,
(*Alternative convention*)
(*1st: Symmetrize over rows of tableau*)
SymFirst = Symmetric;
tab1st = tab;
(*2nd: Antisymmetrize over columns of tableau*)
SymSecond = Antisymmetric;
tab2nd = Combinatorica`TransposeTableau[tab];
];
(*Apply 1st*)
ten1st = ten;
Do[
ten1st  = SymmetrizeNoNormalize[ten1st,SymFirst[tab1st[[i]]]];
,{i,1,Length[tab1st]}];
(*Apply 2nd*)
ten2nd = ten1st;
Do[
ten2nd = SymmetrizeNoNormalize[ten2nd,SymSecond[tab2nd[[i]]]];
,{i,1,Length[tab2nd]}];
(*Normalize action of tableau*)
ntx = Combinatorica`NumberOfTableaux[ShapeOfTableau[tab]];
nom = ntx/Factorial[ArrayDepth[ten]];
nom*ten2nd]

(*Canonical part of tensor corresponding to given integer partition*)
CanonicalPart[ten_,pat_]:=Module[{tax,dim,ctn},
tax = Combinatorica`Tableaux[pat];
dim = Dimensions[ten];
ctn = ConstantArray[0,dim];
(*Sum over irreducible parts, that is, over Young Tableaux*)
Do[ctn = ctn + IrreduciblePart[ten,tax[[cc]]],{cc,1,Length[tax]}];
ctn]

End[]
EndPackage[]
