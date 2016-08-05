(* ::Package:: *)

BeginPackage["Emphasis`IndexGymnastics`"]; 
(*Introduced in version 1.0*)

TwoContract::usage="TwoContract[x,y,{{a,i},{b,j},...}] contracts the indices {a,b,...} of a tensor x with the indices {i,j,...} of a tensor y.
Example:
(*Define tensors x and y*)
x=Table[Subscript[xx,i,j,k],{i,4},{j,4},{k,4}] ; 
y=Table[Subscript[yy,i,j],{i,4},{j,4}] ;
(*Contract 1st and 3rd indices of x with indices of y*)
t=TwoContract[x,y,{{1,1},{3,2}}];
(*Check*)
u=ConstantArray[0,4];
Do[u[[i]]=
    Sum[x[[m,i,n]]*y[[m,n]]
     ,{m,4},{n,4}];
  ,{i,4}];
Union[Flatten[FullSimplify[t-u]]]"

RaiseLower::usage="RaiseLower[x,{i,j,...},h] raises the indices {i,j,...} of the tensor x using the inverse metric h. 
RaiseLower[x,{i,j,...},g] lowers the {i,j,...} indices of the tensor x using the metric g.
Example:
(*Define (3,0)-tensor x*)
x=Table[Subscript[xx,i,j,k],{i,4},{j,4},{k,4}] ;
(*Define metric g*)
g=DiagonalMatrix[{a,b,1,1}];
(*Lower 1st and 3rd indices of x*)
t=RaiseLower[x,{1,3},g];
(*Check*)
u=ConstantArray[0,{4,4,4}];
Do[
  u[[i,j,k]]=
    Sum[x[[p,j,q]]*g[[p,i]]*g[[q,k]]
     ,{p,4},{q,4}];
  ,{i,4},{j,4},{k,4}];
Union[Flatten[FullSimplify[t-u]]]
(*Raise 1st and 3rd indices of t*)
h=Inverse[g];
w=RaiseLower[t,{1,3},h];
(*Check*)
Union[Flatten[FullSimplify[w-x]]]"

MoveUp::usage="MoveUp[x,{i,j,...}] raises the indices {i,j,...} of the tensor x, if the metric is Lorentian and the basis is orthonormal.
MoveUp[x,{i,j,...},s] adopts the convention 'mostly plus' (s=1) or 'mostly minus' (s=-1) for the metric signature. The default option is 'mostly plus'.
Example:
(*Define (0,3)-tensor x*)
x=Table[Subscript[xx,i,j,k],{i,4},{j,4},{k,4}] ;
(*Raise 1st and 3rd indices of x using the inverse metric diag(-1,1,1,1)*)
t=MoveUp[x,{1,3}] ;
(*Check*)
h=DiagonalMatrix[{-1,1,1,1}];
u=RaiseLower[x,{1,3},h];
Union[Flatten[FullSimplify[t-u]]]
(*Raise 1st and 3rd indices of x using the inverse metric diag(1,-1,-1,-1)*)
t=MoveUp[x,{1,3},-1] ;
(*Check*)
h=DiagonalMatrix[{1,-1,-1,-1}];
u=RaiseLower[x,{1,3},h];
Union[Flatten[FullSimplify[t-u]]]"

MoveDown::usage="MoveDown[x,{i,j,...}] lowers the indices {i,j,...} of the tensor x, if the metric is Lorentian and the basis is orthonormal. 
MoveDown is operationally the same as MoveUp. Having two functions makes it easier to keep track of covariant and contravariant indices.  
Example:
(*Define (3,0)-tensor x*)
x=Table[Subscript[xx,i,j,k],{i,4},{j,4},{k,4}] ;
(*Lower 1st and 3rd indices of x using the metric diag(-1,1,1,1)*)
t=MoveDown[x,{1,3}] ;
(*Check*)
u=MoveUp[t,{1,3}];
Union[Flatten[FullSimplify[x-u]]]
(*Lower 1st and 3rd indices of x using the metric diag(1,-1,-1,-1)*)
t=MoveDown[x,{1,3},-1] ;
(*Check*)
u=MoveUp[t,{1,3},-1];
Union[Flatten[FullSimplify[x-u]]]"

Begin["Private`"]
	(*Contract two tensors*)
	TwoContract[x_,y_,i_]:=Module[{xi,yi,xd,yd,xicp,yicp,xnst,ynst,xflt,yflt,t},
	(*Indices of x and y to be contracted*)
	xi=i[[All,1]]; 
	yi=i[[All,2]]; 
	(*Remaining indices of x and y*)
	xd=ArrayDepth[x]; 
	yd=ArrayDepth[y];
	xicp=Complement[Table[a,{a,xd}],xi];
	yicp=Complement[Table[a,{a,yd}],yi];
	(*Contract x and y using Flatten and Dot*)
	xnst=Append[Transpose[{xicp}],xi];
	ynst=Prepend[Transpose[{yicp}],yi];
	t=Flatten[x,xnst].Flatten[y,ynst];
	t]

	(*Raising or lowering of indices*)
	RaiseLower[x_,i_,g_]:=Module[{t,n,p},
	t=x;
	n=ArrayDepth[x];
	(*Raise or lower index j*)
	Do[
	(*Swap j and n*)
	p=Table[k,{k,n}];
	p[[{j,n}]]=p[[{n,j}]];
	t=Transpose[t,p];
	(*Contract with metric or its inverse*)
	t=t.g;
	(*Undo swap of j and n*)
	t=Transpose[t,p];
	,{j,i}];
	t]

	(*Raise indices with Lorentzian metric in orthonormal basis*)
	MoveUp[x_,i_,s_:1]:=Module[{t,n,m,a,u,p,b,v,q},
	t=x;
	n=ArrayDepth[x];
	m=Length[i];
	(*Loop with counter j over the odd numbers between 1 and m*)
	a=Table[c,{c,1,m,2}];
	u=Table[2;;All,{c,m}];
	Do[
	(*Build vector p with p[[1;;j]] and p[[j;;m-j]] equal to 1 and 2;;All*)
	p=u;
	p[[1;;j]]=1;
	(*Loop over permutations of the entries in p*)
	b=Permutations[p];
	v=Table[All,{c,n}];
	Do[
	(*Build vector q that specifies the relevant entries of t*)
	q=v;
	q[[i]]=b[[k]];
	(*Change the sign of the entries in t specified by q*)
	t[[Sequence @@ q]]=-t[[Sequence @@ q]];
    ,{k,Length[b]}];
    ,{j,a}];
	(*Allow for 'mostly minus' signature convention*)
     (s^m)*t]

	(*Lower indices with Lorentzian metric in orthonormal basis*)
	MoveDown[x_,i_,s_:1]:=MoveUp[x,i,s]
End[]
EndPackage[]



