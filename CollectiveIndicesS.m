(* ::Package:: *)

BeginPackage["Emphasis`CollectiveIndicesS`"];
(*Introduced in version 1.0*)

SymmetricToCollective::usage = "the function {1,1}->1, {1,2}->2, {1,3}->3, {1,4}->4, {2,2}->5, {3,3}->6, {4,4}->7, {3,4}->8, {4,2}->9, {2,3}->10. 
It is used to map pairs of symmetric indices into collective indices. 
Example:
a={4,2};
SymmetricToCollective[a]"

CollectiveToSymmetric::usage = "the function 1->{1,1}, 2->{1,2}, 3->{1,3}, 4->{1,4}, 5->{2,2}, 6->{3,3}, 7->{4,4}, 8->{3,4}, 9->{4,2}, 10->{2,3}.
It is used to map collective indices into pairs of symmetric indices.
Example:
c=9;
CollectiveToSymmetric[c]"

TensorToMatrixS::usage="map 4th-order tensors in 4 dimensions whose 1st and 2nd pairs of indices are symmetric to 10x10 matrices. 
Example:
(*Define 4th-order tensor with symmetries*)
gs=DiagonalMatrix[{-1,1,1,1}];
t=TensorProduct[gs,gs];

(*Map to 10x10 matrix*)
m=TensorToMatrixS[t];
MatrixForm[m]

(*Check result*)
gc={-1,0,0,0,1,1,1,0,0,0};
Union[Flatten[FullSimplify[m-TensorProduct[gc,gc]]]]"

MatrixToTensorS::usage="map 10x10 matrices to 4th-order tensors in 4 dimensions whose 1st and 2nd pairs of indices are symmetric.
Example: 
(*Define 10x10 matrix*)
gc={-1,0,0,0,1,1,1,0,0,0};
m=TensorProduct[gc,gc];
 
(*Map to 4th-order tensor with symmetries*)
t=MatrixToTensorS[m];
MatrixForm[t]

(*Check result*)
gs=DiagonalMatrix[{-1,1,1,1}];
Union[Flatten[FullSimplify[t-TensorProduct[gs,gs]]]]"

Begin["Private`"]
	(*From symmetric indices to collective indices*)
	SymmetricToCollective[x_]:=Module[{a,y},
	a={{1,1},{1,2},{1,3},{1,4},{2,2},{3,3},{4,4},{3,4},{4,2},{2,3}};
	y=Position[a,x][[1,1]];
	y]

	(*From collective indices to symmetric indices*)
	CollectiveToSymmetric[x_] := Module[{a,y},
	a={{1,1},{1,2},{1,3},{1,4},{2,2},{3,3},{4,4},{3,4},{4,2},{2,3}};
	y=a[[x,All]];
	y]

	(*From tensors with two pairs of symmetric indices to 10x10 matrices*)
	TensorToMatrixS[x_]:=Module[{y,m,n,p,q,r,s},
		y=ConstantArray[0,{10,10}];
		Do[
			m=CollectiveToSymmetric[i];
			p=m[[1]]; 
			q=m[[2]];	
			n=CollectiveToSymmetric[j];
			r=n[[1]]; 
			s=n[[2]];
			y[[i,j]]=x[[p,q,r,s]]
		,{i,1,10},{j,1,10}];
		y]

	(*From 10x10 matrices to tensors with two pairs of symmetric indices*)
	MatrixToTensorS[x_]:=Module[{y,m,n,p,q,r,s},
	y=ConstantArray[0,{4,4,4,4}];
	Do[
	m=CollectiveToSymmetric[i];
	p=m[[1]]; 
	q=m[[2]];	
	n=CollectiveToSymmetric[j];
	r=n[[1]]; 
	s=n[[2]];
	y[[p,q,r,s]]=x[[i,j]];
	y[[q,p,r,s]]=x[[i,j]];
	y[[p,q,s,r]]=x[[i,j]];
	y[[q,p,s,r]]=x[[i,j]];
	,{i,1,10},{j,1,10}];
	y]
End[]
EndPackage[]



