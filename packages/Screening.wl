(* ::Package:: *)

BeginPackage["Screening`"];


(* ::Text:: *)
(*Functions used for computing Screening (this is a modularized version of Ina's screening code)*)


(* ::Subsubsection:: *)
(*Public Declarations*)


(* ::Text:: *)
(*Source Term Functions for use in GetScreened\[Phi]*)


derivF::usage = "";
complexF::usage = "";
rootF::usage = "";


(* ::Text:: *)
(*Get Screening*)


GetScreened\[Phi]::usage = "Compute the screened potential (assuming neutral plasma) given a value of the potential at Earth's boundary";


(* ::Text:: *)
(*Get captured number density*)


betaSolver::usage = "";
Get\[Beta]rfunc::usage = "";


(* ::Section:: *)
(*Private*)


Begin["Private`"];


(* ::Subsubsection:: *)
(*Source Term (F) Functions*)


(* ::Input::Initialization:: *)
derivF[i_,RMax_,r_,phi_?NumericQ,n_,vInfList_,RMaxList_,MBfactor_]:=
(*Used to find a good endpoint for solving the F[phi]=0 equation (there can be a superfluous root at higher phi)*)Re[N[(-3E^(-3phi)/2-3E^(3phi/2)/2)/n+((MBfactor/(2*\[Sqrt](vInfList[[i]]^2+phi-RMax^2/r^2 (vInfList[[i]]^2+RMaxList[[i,2]]))vInfList[[i]]))),20]];

complexF[i_,RMax_,r_,vInfList_,RMaxList_]:=
(*find where the function F[phi] (where we are solving F[phi]=0 to find phi) goes complex, which will cause the root finding protocol to not work well*)
N[RMax^2/r^2 (vInfList[[i]]^2+RMaxList[[i,2]])-vInfList[[i]]^2,50];

rootF[i_,RMax_,r_,phi_?NumericQ,n_,vInfList_,RMaxList_,MBfactor_]:=
Re[N[(E^(-3phi/2)-E^(3phi/2))/n+MBfactor/vInfList[[i]] Sqrt[vInfList[[i]]^2+phi-RMax^2/r^2 (vInfList[[i]]^2+RMaxList[[i,2]])],50]];


(* ::Subsubsection:: *)
(*Get Screened Potential Outside the Earth*)


(* ::Input::Initialization:: *)
GetScreened\[Phi][phiEarth_:0.05,vInfList_:PowerRange[1/20,4,1+2/10](*,pointNumber_:100*)]:=Module[{MBNorm,MBfactor,pointNumber=100,phiList={{1,0}},RMaxList={{1,0}},minPoint =SetPrecision[ 0.5,20],scanPoints,earthFlag,zoomFlag,phiListWorking,delNum,zoomCount,phiMin,phiMax,res,successplot,phiNew,phiFunc,ckplot,insertPoints,f,MBfactorList,\[Phi]Dict,\[Beta]sol,\[Beta]max},

(*pointNumber=100;*)(*keep pointnumber low for fast run times, ~100 works well*)
(*
vInfList - [Subscript[v, \[Infinity]]/Subscript[v, 0]] velocities at infinity being sampled
phiEarth - [(3 Subscript[v, esc]^2)/(2 Subscript[v, 0]^2)] Electrostatic potential
scanpoints - []
RMaxList - []
*)

MBNorm=SetPrecision[Sum[vInfList[[i]]^2 Exp[-(3/2) vInfList[[i]]^2],{i,1,Length[vInfList]}],20];
MBfactor[i_]:=SetPrecision[ vInfList[[i]]^2 Exp[-(3/2) vInfList[[i]]^2]/MBNorm,20];

(* start of recursion to scan over different r's to find escape potential *)
(*phiList={{1,0}};*)
(*RMaxList={{1,0}};*)
(*minPoint =SetPrecision[ 0.5,20];*)(*CHANGED*)
Clear[scanPoints];
scanPoints=SetPrecision[Table[1-i (1-minPoint)/pointNumber,{i,pointNumber}],20];
earthFlag=False (*triggers True when escape potential matches Earth potential*);
zoomFlag=False (*triggers True when the recursion zooms in to finer value of the radius*);
(*recursion over different values of v_inf*)
Do[
phiListWorking = {};
delNum=0;
zoomCount=0;
If[earthFlag,Break[]];
If[zoomFlag,Break[]];
Print[j];
(*recursion over different values of r*);
Do[
If[zoomCount>12,zoomFlag=True;Print["You've zoomed in more than 12 times, problem?"];
Break[]];

phiMin=Max[Table[complexF[l,RMaxList[[l,1]],scanPoints[[k]],vInfList,RMaxList],{l,j}]]+10^(-15);
Do[If[Sqrt[vInfList[[l]]^2+phiMin-RMaxList[[l,1]]^2/scanPoints[[k]]^2(vInfList[[l]]^2+RMaxList[[l,2]])]<=10^-8,
phiMin+=10^-12],
{l,j}];

phiMax=phi/.FindRoot[Sum[derivF[i,RMaxList[[i,1]],scanPoints[[k]],phi,j,vInfList,RMaxList,MBfactor[i]],{i,1,j}],{phi,phiMin,0.1},Method->"Brent"];

res=Check[FindRoot[Sum[rootF[i,RMaxList[[i,1]],scanPoints[[k]],phi,j,vInfList,RMaxList,MBfactor[i]],{i,1,j}],
{phi,phiMin-10^(-12),phiMax},
Method->"Brent",WorkingPrecision->20],$Failed[Evaluate@$MessageList],
FindRoot::bbrac];


If[
FreeQ[res,$Failed](*check whether FindRoot was succesful*),
successplot=Plot[Sum[rootF[i,RMaxList[[i,1]],scanPoints[[k]],phi,j,vInfList,RMaxList,MBfactor[i]],{i,1,j}],{phi,phiMin-10^(-12),phiMax}];(*ADDED*)
phiNew=phi/.res;
If[phiNew>phiEarth,(*CHANGED*)
earthFlag=True;phiList=Join[phiList,phiListWorking];Print["Reached Earth phi"];
Break[]];

(*make updated potential function so far*);
AppendTo[phiListWorking,{scanPoints[[k]],phiNew}];
phiFunc=Interpolation[Join[phiList,phiListWorking]];
delNum+=1;

(*check if we have hit RMax(v_inf)*);
If[phiFunc'[scanPoints[[k]]]*scanPoints[[k]]+2phiFunc[scanPoints[[k]]]+2(vInfList[[j+1]])^2<=0,AppendTo[RMaxList,{scanPoints[[k]],phiNew}];
phiList=Join[phiList,phiListWorking];
scanPoints=Drop[scanPoints,delNum];
pointNumber-=k;
Return[phiList[[-1]]];];,

(*if FindRoot was not succesful, go back a step in r and insert extra points to "fine-grain" your scan over r*)
Print["Zooming In"];
ckplot=Plot[Sum[rootF[i,RMaxList[[i,1]],scanPoints[[k]],phi,j,vInfList,RMaxList,MBfactor[i]],{i,1,j}],{phi,phiMin-10^(-12),phiMax}];(*ADDED*)
zoomCount+=1;
(*pick extra points to insert*);
If[k==1,insertPoints=Table[RMaxList[[j,1]]-(RMaxList[[j,1]]-scanPoints[[k]])/10 p,{p,0,9}],insertPoints=Table[scanPoints[[k-1]]+(scanPoints[[k]]-scanPoints[[k-1]])/10 p,{p,0,9}]];

scanPoints=Join[scanPoints[[1;;k-1]],insertPoints,scanPoints[[k;;]]];
pointNumber+=10;
delNum+=1
],

{k,10000}],
{j,1,Length[vInfList]-1}];

MBfactorList = Table[MBfactor[i],{i,Length[vInfList]}];

f=Interpolation[phiList];

\[Phi]Dict=<|"\[Phi](r)"->f,"phiList"->phiList,"RMaxList"->RMaxList,"\[Phi]E"->phiEarth,"v\[Infinity]List"->vInfList,"MBfactorList"->MBfactorList,"MBNorm"->MBNorm,"rE"->phiList[[-1,1]]|>;

\[Beta]sol = Get\[Beta]rfunc[\[Phi]Dict];

\[Beta]max = Max[\[Beta]sol["ValuesOnGrid"]];

Append[\[Phi]Dict,{"\[Beta](r)"->\[Beta]sol,"\[Beta]max"->\[Beta]max}]
]


(* ::Subsubsection:: *)
(*Get \[Beta](r) from \[Phi] solution*)


betaSolver[\[Phi]Dict_,r_]:=Module[{j,RMaxList,phiE,vInfList,MBfactorList},
(*Use to solve for overall charge density in the Earth as a function of radius. Subscript[n, H] - Subscript[n, e] = Subscript[\[Alpha], H ]* \[Beta](r) * (Subscript[n, H]^F + Subscript[n, e]^F)*)

{RMaxList,phiE,vInfList,MBfactorList}= \[Phi]Dict[#]&/@{"RMaxList","\[Phi]E","v\[Infinity]List","MBfactorList"};

j = Length[RMaxList];

-Sum[rootF[i,RMaxList[[i,1]],r,phiE,j,vInfList,RMaxList,MBfactorList[[i]]],{i,1,j}]
]


Get\[Beta]rfunc[\[Phi]Dict_,N_:100]:=Module[{rE,interTable},
rE = \[Phi]Dict["rE"];
interTable=Table[{r,betaSolver[\[Phi]Dict,r]},{r,rE/#,rE(1-1/#),rE/#}]&[N];
Interpolation[interTable]
]


(* ::Section::Closed:: *)
(*End*)


End[];


EndPackage[];