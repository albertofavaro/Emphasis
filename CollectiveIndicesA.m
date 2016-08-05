(* ::Package:: *)

BeginPackage["Emphasis`CollectiveIndicesA`"];
(*Introduced in version 1.0*)

AntisymmetricToCollective::usage = "the function {1,2}->1, {1,3}->2, {1,4}->3, {3,4}->4, {4,2}->5, {2,3}->6. 
It is used to map pairs of antisymmetric indices into collective indices. 
Example:
a={4,2};
AntisymmetricToCollective[a]"

CollectiveToAntisymmetric::usage = "the function 1->{1,2}, 2->{1,3}, 3->{1,4}, 4->{3,4}, 5->{4,2}, 6->{2,3}. 
It is used to map collective indices into pairs of antisymmetric indices.
Example:
c=5;
CollectiveToAntisymmetric[c]"

TensorToMatrixA::usage="map 4th-order tensors in 4 dimensions whose 1st and 2nd pairs of indices are antisymmetric to 6x6 matrices. 
Example: 
t=LeviCivitaTensor[4];
m=TensorToMatrixA[t];
MatrixForm[m]"

MatrixToTensorA::usage="map 6x6 matrices to 4th-order tensors in 4 dimensions whose 1st and 2nd pairs of indices are antisymmetric.
Example: 
m={{0,0,0,1,0,0},{0,0,0,0,1,0},{0,0,0,0,0,1},{1,0,0,0,0,0},{0,1,0,0,0,0},{0,0,1,0,0,0}};
t=MatrixToTensorA[m]; 
MatrixForm[t]
Union[Flatten[FullSimplify[t-LeviCivitaTensor[4]]]]"

Begin["Private`"]

	(*From antisymmetric indices to collective indices*)
	AntisymmetricToCollective[x_]:=Module[{a,y},
	a={{1,2},{1,3},{1,4},{3,4},{4,2},{2,3}};
	y=Position[a,x][[1,1]];
	y]

	(*From collective indices to antisymmetric indices*)
	CollectiveToAntisymmetric[x_] := Module[{a,y},
	a={{1,2},{1,3},{1,4},{3,4},{4,2},{2,3}};
	y=a[[x,All]];
	y]

	(*From tensors with two pairs of antisymmetric indices to 6x6 matrices*)
	TensorToMatrixA[x_]:=Module[{y,m,n,p,q,r,s},
		y=ConstantArray[0,{6,6}];
		Do[
			m=CollectiveToAntisymmetric[i];
			p=m[[1]]; 
			q=m[[2]];	
			n=CollectiveToAntisymmetric[j];
			r=n[[1]]; 
			s=n[[2]];
			y[[i,j]]=x[[p,q,r,s]]
		,{i,1,6},{j,1,6}];
		y]

	(*From 6x6 matrices to tensors with two pairs of antisymmetric indices*)
	MatrixToTensorA[x_]:=Module[{y,m,n,p,q,r,s},
	y=ConstantArray[0,{4,4,4,4}];
	Do[
	m=CollectiveToAntisymmetric[i];
	p=m[[1]]; 
	q=m[[2]];	
	n=CollectiveToAntisymmetric[j];
	r=n[[1]]; 
	s=n[[2]];
	y[[p,q,r,s]]=x[[i,j]];
	y[[q,p,r,s]]=-x[[i,j]];
	y[[p,q,s,r]]=-x[[i,j]];
	y[[q,p,s,r]]=x[[i,j]];
	,{i,1,6},{j,1,6}];
	y]
End[]

EndPackage[]



