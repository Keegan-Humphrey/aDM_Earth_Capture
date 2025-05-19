(* ::Package:: *)

Needs["Constants`",NotebookDirectory[]<>"/Constants.wl"]
Needs["Capture`",NotebookDirectory[]<>"/Capture.wl"]


BeginPackage["Screening`"];


(* ::Text:: *)
(*Functions used for computing Screening (this is a modularized version of Ina's screening code)*)


(* ::Subsection:: *)
(*Public Declarations*)


(* ::Subsubsection:: *)
(*Ina's code*)


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


(* ::Text:: *)
(*Get total captured number*)


GetTotalCapturedN::usage = "";


(* ::Subsubsection:: *)
(*flat potential code*)


getncap::usage = "";


getrmax::usage = "";


Ncapof\[Delta]::usage = "";


\[Delta]ofNcap::usage = "";


get\[Delta]fromPDict::usage = "";


Compute\[Delta]onScan::usage = "";


Get\[Delta]Truth::usage = "";


(* ::Subsubsection:: *)
(*Potential outside the flat region*)


get\[Phi]Dnoncap::usage = "";


get\[Phi]DAnalyticEstimates::usage = "";


(* ::Section:: *)
(*Private*)


Begin["Private`"];


(* ::Subsection:: *)
(*Ina's code and developments therein*)


(* ::Subsubsection::Closed:: *)
(*Source Term (F) Functions*)


(* ::Input::Initialization:: *)
derivF[i_,RMax_,r_,phi_?NumericQ,n_,vInfList_,RMaxList_,MBfactor_]:=
(*Used to find a good endpoint for solving the F[phi]=0 equation (there can be a superfluous root at higher phi)*)Re[N[(-3E^(-3phi)/2-3E^(3phi/2)/2)/n+((MBfactor/(2*\[Sqrt](vInfList[[i]]^2+phi-RMax^2/r^2 (vInfList[[i]]^2+RMaxList[[i,2]]))vInfList[[i]]))),20]];

complexF[i_,RMax_,r_,vInfList_,RMaxList_]:=
(*find where the function F[phi] (where we are solving F[phi]=0 to find phi) goes complex, which will cause the root finding protocol to not work well*)
N[RMax^2/r^2 (vInfList[[i]]^2+RMaxList[[i,2]])-vInfList[[i]]^2,50];

rootF[i_,RMax_,r_,phi_?NumericQ,n_,vInfList_,RMaxList_,MBfactor_]:=
Re[N[(E^(-3phi/2)-E^(3phi/2))/n+MBfactor/vInfList[[i]] Sqrt[vInfList[[i]]^2+phi-RMax^2/r^2 (vInfList[[i]]^2+RMaxList[[i,2]])],50]];


(* ::Subsubsection::Closed:: *)
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


(* ::Subsubsection::Closed:: *)
(*Get \[Beta](r) from \[Phi] solution*)


betaSolver[\[Phi]Dict_,r_]:=Module[{j,RMaxList,phiE,vInfList,MBfactorList},
(*Use to solve for overall charge density in the Earth as a function of radius. Subscript[n, H] - Subscript[n, e] = \[Beta](r) * Subscript[n, H]^F *)

{RMaxList,phiE,vInfList,MBfactorList}= \[Phi]Dict[#]&/@{"RMaxList","\[Phi]E","v\[Infinity]List","MBfactorList"};

j = Length[RMaxList];

-Sum[rootF[i,RMaxList[[i,1]],r,phiE,j,vInfList,RMaxList,MBfactorList[[i]]],{i,1,j}]
]


Get\[Beta]rfunc[\[Phi]Dict_,N_:100]:=Module[{rE,interTable},
rE = \[Phi]Dict["rE"];
interTable=Table[{r,betaSolver[\[Phi]Dict,r]},{r,rE/#,rE(1-1/#),rE/#}]&[N];
Interpolation[interTable]
]


(* ::Subsubsection::Closed:: *)
(*Get Total captured number*)


GetTotalCapturedN[\[Phi]Dict_,mDMtot_:"m"2 \!\(\*SuperscriptBox[\(10\), \("\<logm\>" - 6\)]\)/.Constants`SIConstRepl,fD_:0.05]:=Module[{n\[Infinity]aDM,\[Beta]ofr,\[Beta]dom},

(*logm - [eV] Subscript[Log, 10]Subscript[m, Subscript[p, D]]*)

n\[Infinity]aDM=(Capture`naDM[fD,mDMtot]m^-3);

\[Beta]ofr = \[Phi]Dict["\[Beta](r)"];
\[Beta]dom = \[Beta]ofr["Domain"];
n\[Infinity]aDM(4 \[Pi] (("rE")/\[Beta]dom[[1,2]])^3 m^3/.Constants`EarthRepl)NIntegrate[r^2 \[Beta]ofr[r],{r,\[Beta]dom[[1,1]],\[Beta]dom[[1,2]]}](*Subscript[n, \[Infinity]]/V \[Integral] \[DifferentialD]^3r \[Beta](r)*)
]


(* ::Subsection:: *)
(*Consistent Flat potential *)


(* ::Subsubsection:: *)
(*get Subscript[n, C](r)*)


getncap[mpD_,meD_,\[Beta]D_,nF_,\[Alpha]Dby\[Alpha]_]:=Module[{niofr,npofr,neofr,nc,\[Phi]g,\[Phi]gE},
\[Phi]g = Function[{r},Capture`Vgrav[r]];
\[Phi]gE = \[Phi]g["rE"/.Constants`EarthRepl];
niofr=1/(mi^2 Sqrt[2 \[Pi]] \[Beta]D^2) nF (mi \[Beta]D)^(3/2) (2 mi vesc \[Beta]D+E^(1/2 mi vesc^2 \[Beta]D) Sqrt[2 \[Pi]] Sqrt[mi \[Beta]D]-E^(1/2 mi vesc^2 \[Beta]D) Sqrt[mi] Sqrt[2 \[Pi]] Sqrt[\[Beta]D] Erf[(Sqrt[mi] vesc Sqrt[\[Beta]D])/Sqrt[2]]);(*via: niofr=nF (mi \[Beta]D)^(3/2)Sqrt[2/\[Pi]]Integrate[v^2E^(- \[Beta]D mi/2(v^2-vesc^2)),{v,vesc,\[Infinity]}]//Normal*)
npofr = niofr/.mi->mpD/.vesc->Sqrt[(-2 ("\[Delta]" mpD \[Phi]gE))/mpD];
neofr = niofr/.mi->meD/.vesc->Sqrt[(-2((mpD+meD)\[Phi]g[r]-"\[Delta]" mpD \[Phi]gE))/meD];
nc=(mpD "G" "ME" )/(("e")^2/("\[Epsilon]0") \[Alpha]Dby\[Alpha]((4 \[Pi])/3 ("rE")^3)) Boole["rE" -r>0] - (npofr-neofr)/.Constants`SIConstRepl/.Constants`EarthRepl;
nc
]


(* ::Subsubsection:: *)
(*get Subscript[r, max](\[Delta])*)


getrmax[ncof\[Delta]andr_]:=rmax/.FindRoot[Re@(ncof\[Delta]andr/."\[Delta]"->#/.r->rmax),{rmax,Evaluate[("rE")/(2 #)/.Constants`EarthRepl]}]&


(* ::Subsubsection:: *)
(*get Subscript[N, C](\[Delta])*)


Ncapof\[Delta][ncof\[Delta]andr_,\[Delta]t_?NumberQ]:= 4 \[Pi] NIntegrate[r^2 ncof\[Delta]andr/."\[Delta]"->\[Delta]t,{r,0,getrmax[ncof\[Delta]andr][\[Delta]t]}]


(* ::Subsubsection:: *)
(*get \[Delta](Subscript[N, C])*)


\[Delta]ofNcap[ncof\[Delta]andr_,NC_]:= \[Delta]t/.FindRoot[Evaluate@(Ncapof\[Delta][ncof\[Delta]andr,\[Delta]t]==NC),{\[Delta]t,(*10^-2*)10^-7,10^-10,1}]


(* ::Subsubsection:: *)
(*get \[Delta] from a PDict*)


get\[Delta]fromPDict[pdict_,\[Alpha]Dby\[Alpha]_:1,fD_:0.05]:=Module[{nD,nCap,NC,Ncapof\[Delta]table,\[Delta]t,rmax},
(*The case of \[Delta]>1 still needs to be sorted out*)

nD =Capture`naDM[fD,"meD" + "mpD"]/.pdict;
nCap=getncap["mpD" ,"meD" ,"\[Beta]D",nD,\[Alpha]Dby\[Alpha]]/.pdict;

NC = "Nc"/.pdict/."nD"->nD;

\[Delta]t=\[Delta]ofNcap[nCap,NC];
rmax=getrmax[nCap][\[Delta]t];

<|"\[Delta]"->\[Delta]t,"rmax"->rmax,"NC"->NC|>
]


(* ::Subsubsection:: *)
(*get \[Delta] from a scan *)


Clear[Compute\[Delta]onScan]
Compute\[Delta]onScan[\[Delta]List_,v0ind_:1]:= Module[{return\[Delta]List,\[Delta]Listcurrent,PDict,v0,mratio,\[Delta]table},
(*
\[Delta]List - of the form: {<|"\[Delta]"->\[Delta],"file"->"\\path\\to\\file"|>}
*)
return\[Delta]List ={};

(*Print[\[Delta]List];
Print[Length[\[Delta]List]];*)
(*Print[Utilities`ReadIt[\[Delta]List[[1]]["file"]]];*)
(*Print[Dimensions@Utilities`ReadIt[\[Delta]List[[1]]["file"]]];
Print[\[Delta]List];
Print[Length[\[Delta]List]];*)

Monitor[
Do[
(*Table[*)
(*Print["test 1"];
Print[\[Delta]];*)
\[Delta]Listcurrent = \[Delta]List[[\[Delta]]];
(*Print[\[Delta]Listcurrent];*)
(*Print[Utilities`ReadIt[\[Delta]Listcurrent["file"]][[1]]];

Print[\[Delta]Listcurrent["\[Delta]"]];

Print[\[Delta]Listcurrent];
Print[\[Delta]Listcurrent["file"]];*)

PDict= Utilities`ReadIt[\[Delta]Listcurrent["file"]][[;;9,v0ind,;;]];
(*PDict = PDict[[;;9,v0ind,;;]];*)

(*Print[PDict];*)
(*Print[{"v0","meD"/"mpD"}/.PDict];*)


(*v0, mratio = First@First[{"v0","meD"/"mpD"}/.PDict];*)
{v0, mratio} = {"v0","meD"/"mpD"}/.First@First[PDict];

(*Print[v0,"\n",mratio];*)

\[Delta]table=Flatten[Table[Table[{Log10[("mpD" ("c")^2)/("JpereV")]/.PDict[[i,j]]/.Constants`SIConstRepl,Log10["\[Kappa]"]/.PDict[[i,j]],Log10@"\[Delta]"/.get\[Delta]fromPDict[PDict[[i,j]]]},{i,Dimensions[PDict][[1]]}],{j,Dimensions[PDict][[2]]}],1];

(*PDict["\[Delta]table"]=\[Delta]table;*)
(*AppendTo[\[Delta]Listcurrent,<|"\[Delta]table"->\[Delta]table|>];*)
AppendTo[\[Delta]Listcurrent,<|"\[Delta]table"->\[Delta]table,"v0"->v0,"meDbympD"->mratio|>];
AppendTo[return\[Delta]List,\[Delta]Listcurrent];

,{\[Delta],Length[\[Delta]List]}];
,{\[Delta],j,i}];

Print["Finished running the PDicts for: \n{v0,meD/mpD}=",First@First[{"v0","meD"/"mpD"}/.PDict]];

return\[Delta]List
]


(* ::Subsubsection:: *)
(*Get \[Delta] Truth*)


Get\[Delta]Truth[\[Delta]Scan_]:=Module[{\[Delta]truthtable,\[Delta]outof\[Delta]intable,\[Delta]outof\[Delta]in,\[Delta]truth},

\[Delta]truthtable ={};

Do[
(*Print@Union[\[Delta]Scan[[;;]]["\[Delta]table"][[i,;;2]]//N//Round];*)
\[Delta]outof\[Delta]intable =Table[{Log10@"\[Delta]"/.\[Delta]Scan[[j]],\[Delta]Scan[[j]]["\[Delta]table"][[i,3]]},{j,Length[\[Delta]Scan]}];
\[Delta]outof\[Delta]in = Interpolation[\[Delta]outof\[Delta]intable,InterpolationOrder->1];
\[Delta]truth = \[Delta]/.FindRoot[10^\[Delta]outof\[Delta]in[Log10@\[Delta]]-\[Delta],{\[Delta],0.01}];
AppendTo[\[Delta]truthtable,{\[Delta]Scan[[1]]["\[Delta]table"][[i,1]],\[Delta]Scan[[1]]["\[Delta]table"][[i,2]],Log10@\[Delta]truth}];
,{i,Length[\[Delta]Scan[[1]]["\[Delta]table"]]}];(*assuming all \[Delta]tables have the same ordering size*)

\[Delta]truthtable
]


(* ::Subsection:: *)
(*Potential outside the flat region*)


(* ::Subsubsection:: *)
(*Get EOM*)


Clear[get\[Phi]Dnoncap]
get\[Phi]Dnoncap[mpD_,meD_,\[Beta]D_,nF_,\[Alpha]Dby\[Alpha]_]:=Module[{niofr,npofr,neofr,nc,\[Psi]g,EOM},
\[Psi]g = Function[{r},- mpD \[Beta]D Capture`Vgrav[r]];
niofr= nF (-E^(((mi \[Psi]g[r])/mpD)) (-1+Erf[(Sqrt[mi] Sqrt[\[Psi]g[r]])/Sqrt[mpD]])+E^((mi \[Psi]g[r])/mpD) ("e")Sqrt[\[Alpha]Dby\[Alpha]] qD \[Beta]D (-1+Erf[(Sqrt[mi] Sqrt[\[Psi]g[r]])/Sqrt[mpD]]) \[Phi]D[r]+(2 Sqrt[mi] Sqrt[\[Psi]g[r]])/(Sqrt[mpD] Sqrt[\[Pi]]));(*linearized assuming eD \[Phi]D \[Beta]D << 1 (an ansatz that must be checked self consistently)*)
(*niofr=nF-("e")Sqrt[\[Alpha]Dby\[Alpha]]nF qD \[Beta]D \[Phi]D[r]+(mi nF \[Psi]g[r])/mpD;*)
npofr = niofr/.mi->mpD/.qD->1;
neofr = niofr/.mi->meD/.qD->-1;
EOM = -1/r^2 D[r^2 \[Phi]D'[r],r]==("e")/("\[Epsilon]0") Sqrt[\[Alpha]Dby\[Alpha]](npofr-neofr)/.Constants`SIConstRepl/.Constants`EarthRepl
]


(* ::Subsubsection:: *)
(*get analytic estimate for \[Phi]D in the small charge regime*)


get\[Phi]DAnalyticEstimates[mpD_,meD_,\[Beta]D_,nF_,\[Alpha]Dby\[Alpha]_,\[Phi]g_:Capture`Vgrav] := Module[{niofr,npD,neD,\[Rho]S,\[Lambda]D},
niofr = 1/(mi^2 Sqrt[2 \[Pi]] \[Beta]D^2) nF (mi \[Beta]D)^(3/2) (2 mi vesc \[Beta]D+E^(1/2 mi vesc^2 \[Beta]D) Sqrt[2 \[Pi]] Sqrt[mi \[Beta]D]-E^(1/2 mi vesc^2 \[Beta]D) Sqrt[mi] Sqrt[2 \[Pi]] Sqrt[\[Beta]D] Erf[(Sqrt[mi] vesc Sqrt[\[Beta]D])/Sqrt[2]]);
npD=niofr/.vesc->Sqrt[-((2 Vgp)/mi)]/.mi->mpD;
neD=niofr/.vesc->Sqrt[-((2 Vge)/mi)]/.mi->meD;
\[Rho]S = ("e")/("\[Epsilon]0") Sqrt[\[Alpha]Dby\[Alpha]](npD-neD)/.{Vgp->mpD \[Phi]g[r],Vge->meD \[Phi]g[r]};
(*\[Lambda]D = (("e")^2/("\[Epsilon]0")\[Alpha]Dby\[Alpha](D[npD,Vgp]+D[neD,Vge]))^(-1/2)/.{Vgp->mpD \[Phi]g[r],Vge->mpD \[Phi]g[r]};*)
\[Lambda]D =(-(("e")^2/("\[Epsilon]0"))\[Alpha]Dby\[Alpha](D[npD,Vgp]+D[neD,Vge])/.{Vgp->mpD \[Phi]g[r],Vge->meD \[Phi]g[r]})^(-1/2);
<|"\[Phi]D"->\[Rho]S \[Lambda]D^2,"\[Rho]S"->\[Rho]S,"\[Lambda]D"->\[Lambda]D|>
]


(* ::Section::Closed:: *)
(*End*)


End[];


EndPackage[];
