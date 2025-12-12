(* ::Package:: *)

Needs["Constants`",NotebookDirectory[]<>"/Constants.wl"]
Needs["Capture`",NotebookDirectory[]<>"/Capture.wl"]
(*Needs["ThermTimescales`",NotebookDirectory[]<>"/ThermTimescales.wl"]*)


Needs["CollisionRate`",NotebookDirectory[]<>"../packages/CollisionRate.wl"]
Needs["Utilities`",NotebookDirectory[]<>"../packages/Utilities.wl"]


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


Compute\[Delta]DictonScan::usage="";


Get\[Delta]DictTruth::usage="";


ScanOverUnfixedParameters::usage="";


checkforVeffMaxima::usage="";


(* ::Subsubsection:: *)
(*Potential outside the flat region*)


get\[Phi]Dnoncap::usage = "";


get\[Phi]DAnalyticEstimates::usage = "";


get\[Phi]D::usage = "";


binsearch::usage = "";


(* ::Subsubsection:: *)
(*Plasma Interactions*)


(* ::Text:: *)
(*Modules for Capture*)


Get\[Tau]A::usage = "";
Get\[Tau]C::usage = "";


(* ::Text:: *)
(*Timescale Modules*)


Get\[Tau]AA::usage = "";
Get\[Tau]CC::usage = "";
Get\[Tau]CE::usage = "";
Get\[Tau]AE::usage = "";
Get\[Tau]CA::usage = "";


(* ::Text:: *)
(*Get timescale table*)


Get\[Tau]sFromPDict::usage = "";


(* ::Text:: *)
(*Get rate table*)


GetPlasmaRates::usage = "";


(* ::Chapter:: *)
(*Private*)


Begin["Private`"];


(* ::Subsection::Closed:: *)
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


(* ::Section:: *)
(*Screening*)


(* ::Subsection:: *)
(*Consistent Flat potential *)


(* ::Subsubsection::Closed:: *)
(*get Subscript[n, C](r)*)


(*getncap[mpD_,meD_,\[Beta]D_,nF_,\[Alpha]Dby\[Alpha]_]:=Module[{niofr,npofr,neofr,nc,\[Phi]g,\[Phi]gE},
\[Phi]g = Function[{r},Capture`Vgrav[r]];
\[Phi]gE = \[Phi]g["rE"/.Constants`EarthRepl];
niofr=1/(mi^2 Sqrt[2 \[Pi]] \[Beta]D^2) nF (mi \[Beta]D)^(3/2) (2 mi vesc \[Beta]D+E^(1/2 mi vesc^2 \[Beta]D) Sqrt[2 \[Pi]] Sqrt[mi \[Beta]D]-E^(1/2 mi vesc^2 \[Beta]D) Sqrt[mi] Sqrt[2 \[Pi]] Sqrt[\[Beta]D] Erf[(Sqrt[mi] vesc Sqrt[\[Beta]D])/Sqrt[2]]);(*via: niofr=nF (mi \[Beta]D)^(3/2)Sqrt[2/\[Pi]]Integrate[v^2E^(- \[Beta]D mi/2(v^2-vesc^2)),{v,vesc,\[Infinity]}]//Normal*)
npofr = niofr/.mi->mpD/.vesc->Sqrt[(-2 ("\[Delta]" mpD \[Phi]gE))/mpD];
neofr = niofr/.mi->meD/.vesc->Sqrt[(-2((mpD+meD)\[Phi]g[r]-"\[Delta]" mpD \[Phi]gE))/meD];
nc=(mpD "G" "ME" )/(("e")^2/("\[Epsilon]0") \[Alpha]Dby\[Alpha]((4 \[Pi])/3 ("rE")^3)) Boole["rE" -r>0] - (npofr-neofr)/.Constants`SIConstRepl/.Constants`EarthRepl;
nc
]*)


(* ::Text:: *)
(*The reason this is screwing up is that the vescs are wrong for r>rmax, it assumes rmax < 0*)


Clear[getncap]
getncap[mpD_,meD_,\[Beta]D_,nF_,\[Alpha]Dby\[Alpha]_]:=Module[{niofr,npofr,neofr,nG,nc,\[Phi]g,\[Phi]gE},
\[Phi]g = Function[{r},Capture`Vgrav[r]];
\[Phi]gE = \[Phi]g["rE"/.Constants`EarthRepl];
niofr=1/(mi^2 Sqrt[2 \[Pi]] \[Beta]D^2) nF (mi \[Beta]D)^(3/2) (2 mi vesc \[Beta]D+E^(1/2 mi vesc^2 \[Beta]D) Sqrt[2 \[Pi]] Sqrt[mi \[Beta]D]-E^(1/2 mi vesc^2 \[Beta]D) Sqrt[mi] Sqrt[2 \[Pi]] Sqrt[\[Beta]D] Erf[(Sqrt[mi] vesc Sqrt[\[Beta]D])/Sqrt[2]]);(*via: niofr=nF (mi \[Beta]D)^(3/2)Sqrt[2/\[Pi]]Integrate[v^2E^(- \[Beta]D mi/2(v^2-vesc^2)),{v,vesc,\[Infinity]}]//Normal*)
npofr = niofr/.mi->mpD/.vesc->Sqrt[(-2 ("\[Delta]" mpD \[Phi]gE))/mpD];
(*neofr = niofr/.mi->meD/.vesc->(*Re@*)Sqrt[(-2((mpD+meD)\[Phi]g[r]-"\[Delta]" mpD \[Phi]gE))/meD];*)
neofr = niofr/.mi->meD/.vesc->Sqrt[If[#>0,#,0]&@(-2((mpD+meD)\[Phi]g[r]-"\[Delta]" mpD \[Phi]gE)/meD)];
nG = (mpD "G" "ME" )/(("e")^2/("\[Epsilon]0") \[Alpha]Dby\[Alpha]((4 \[Pi])/3 ("rE")^3)) Boole["rE" -r>0]/.Constants`SIConstRepl/.Constants`EarthRepl;
nc = nG - (npofr - neofr)/.Constants`SIConstRepl/.Constants`EarthRepl;
(*Re@*)nc (*nc can develop an imaginary part from vesc used in neofr if r>rmax.  This restrict nc to be real (as needed by the root finding functions), in reality it is 0 beyond there but this will give it a non-zero (but negative) value. So nc is only valid for r<=rmax. The restriction must be placed explicitly while evaluating the function. *)
(*<|"npofr"->npofr,"neofr"->neofr,"nG"->nG,"nc"->nc|>*)
]


(* ::Subsubsection::Closed:: *)
(*get Subscript[r, max](\[Delta])*)


(*getrmax[ncof\[Delta]andr_]:=rmax/.FindRoot[Re@(ncof\[Delta]andr/."\[Delta]"->#/.r->rmax),{rmax,Evaluate[("rE")/(2 #)/.Constants`EarthRepl]},Method->"Secant"]&*)


 Clear[binsearch]
binsearch[f_,var_,domain_:{1,11},n_:20,DEBUG_:False(*,ofmmin_:1,ofmmax_:13*)(*,LogQ_:True*)]:=Module[{midlist={},interval,fposQ,fmingrfmaxQ,mid,root(*,DEBUG=False*)},
		interval=domain;
		
		 (*Monitor[*)
		 Do[mid = Mean[interval];
		 
		 fposQ = (f/.var->10^mid)>0; (*if positive, decrease, if negative increase*)
		 fmingrfmaxQ = (f/.var->10^interval[[1]])>(f/.var->10^interval[[2]]);(*if need to decrease, take smaller, if need to increase take larger*)
			interval =If[fposQ,
			(*need to decrease*)
			If[fmingrfmaxQ,
				(*f(xmin) is larger, so take interval with xmax*)
				{mid,interval[[2]]},
				(*f(xmin) is smaller, so take interval with xmin*)
				{interval[[1]],mid}
				],
			(*need to increase*)
			If[fmingrfmaxQ,
				(*f(xmin) is larger, so take interval with xmin*)
				{interval[[1]],mid},
				(*f(xmin) is smaller, so take interval with xmax*)
				{mid,interval[[2]]}
				]
			];
			
			If[DEBUG,AppendTo[midlist,{mid,(f/.var->10^mid)}]];
			
			,{k,n}];
			
			(*,k];*)
			
		If[DEBUG,Print[ListPlot[midlist,PlotRange->All]]];
		
		root = 10^mid;
		
		(*If[#<0,Print[#]]&@(f/.var->root);*)
		
		(*Print@(f/.var->root);
		Print[root];*)
		
		N[{root,f/.var->root}]
	]


(*getrmax[ncof\[Delta]andr_]:=binsearch[Abs@(*Re@*)(ncof\[Delta]andr/."\[Delta]"->#),r,{Log10@Evaluate[("rE")/(2 #)/.Constants`EarthRepl]-1,Log10@Evaluate[("rE")/(2 #)/.Constants`EarthRepl]+1},30][[1]]&*)


(*getrmax[ncof\[Delta]andr_]:=binsearch[Abs@(*Re@*)(ncof\[Delta]andr/."\[Delta]"->#),r,{Log10@Evaluate[("rE")/(2 #)/.Constants`EarthRepl]-2,Log10@Evaluate[("rE")/(2 #)/.Constants`EarthRepl]+1},30][[1]]&*)


getrmax[ncof\[Delta]andr_]:=binsearch[(*Abs@*)(*Re@*)(ncof\[Delta]andr/."\[Delta]"->#),r,{Log10@Evaluate[("rE")/(2 #)/.Constants`EarthRepl]-2,Log10@Evaluate[("rE")/(2 #)/.Constants`EarthRepl]+1},30][[1]]&


(*getrmax[ncof\[Delta]andr_]:=10^Log10rmax/.FindRoot[Abs@(*Re@*)(ncof\[Delta]andr/."\[Delta]"->#/.r->10^Log10rmax),{Log10rmax,Log10@Evaluate[("rE")/(10 #)/.Constants`EarthRepl]},Method->"Secant"(*,MaxIterations->50,AccuracyGoal->10,PrecisionGoal->10*)]&*)


(*getrmax[ncof\[Delta]andr_]:=binsearch[Abs@(*Re@*)(ncof\[Delta]andr/."\[Delta]"->#),r,{1,1+Log10@Evaluate[("rE")/(10 #)/.Constants`EarthRepl]},70][[1]]&*)


(*getrmax[ncof\[Delta]andr_]:=rmax/.FindRoot[Re@(ncof\[Delta]andr/."\[Delta]"->#/.r->rmax),{rmax,Evaluate[("rE")/(100 #)/.Constants`EarthRepl]}]&*)


(*getrmax[ncof\[Delta]andr_]:=Module[{rmaxfunc,x0fators={2,100}},
(*rmaxfunc=rmax/.FindRoot[Re@(ncof\[Delta]andr/."\[Delta]"->#/.r->rmax),{rmax,Evaluate[("rE")/(x0fators[[1]] #)/.Constants`EarthRepl]}]&;*)
If[Evaluate[ncof\[Delta]andr/."\[Delta]"->#/.r->rmax/.FindRoot[Re@(ncof\[Delta]andr/."\[Delta]"->#/.r->rmax),{rmax,Evaluate[("rE")/(x0fators[[1]] #)/.Constants`EarthRepl]}]>0],
Evaluate[rmax/.FindRoot[Re@(ncof\[Delta]andr/."\[Delta]"->#/.r->rmax),{rmax,Evaluate[("rE")/(x0fators[[1]] #)/.Constants`EarthRepl]}]],
Evaluate[rmax/.FindRoot[Re@(ncof\[Delta]andr/."\[Delta]"->#/.r->rmax),{rmax,Evaluate[("rE")/(x0fators[[2]] #)/.Constants`EarthRepl]}]]
]&
]
*)


(*getrmax[ncof\[Delta]andr_]:=With[{result=Quiet@Check[FindRoot[Re@(ncof\[Delta]andr/."\[Delta]"->#/.r->10^Log10rmax),{Log10rmax,Log10@Evaluate[("rE")/(2 #)/.Constants`EarthRepl]}],"Failed"](*x0fators={2,100}*)},
(*rmaxfunc=rmax/.FindRoot[Re@(ncof\[Delta]andr/."\[Delta]"->#/.r->rmax),{rmax,Evaluate[("rE")/(x0fators[[1]] #)/.Constants`EarthRepl]}]&;*)


Print[result];

If[!StringQ[result],rmax/.result,
Evaluate[10^Log10rmax/.FindRoot[Re@(ncof\[Delta]andr/."\[Delta]"->#/.r->10^Log10rmax),{Log10rmax,Log10@Evaluate[("rE")/(100 #)/.Constants`EarthRepl]}]]
]
]&

*)


(* ::Subsubsection:: *)
(*get Subscript[N, C](\[Delta])*)


Ncapof\[Delta][ncof\[Delta]andr_,\[Delta]t_?NumberQ]:= 4 \[Pi] NIntegrate[r^2 ncof\[Delta]andr/."\[Delta]"->\[Delta]t,{r,0,getrmax[ncof\[Delta]andr][\[Delta]t]}]


(* ::Subsubsection::Closed:: *)
(*get \[Delta](Subscript[N, C])*)


(*\[Delta]ofNcap[ncof\[Delta]andr_,NC_]:= \[Delta]t/.FindRoot[Evaluate@(Ncapof\[Delta][ncof\[Delta]andr,\[Delta]t]==NC),{\[Delta]t,(*10^-2*)10^-7,10^-10,1}]*)


(*\[Delta]ofNcap[ncof\[Delta]andr_,NC_]:= binsearch[(*Evaluate@*)(Ncapof\[Delta][ncof\[Delta]andr,\[Delta]t]-NC),\[Delta]t,{-7,Log10[0.75]},30][[1]]*)


(*\[Delta]ofNcap[ncof\[Delta]andr_,NC_]:= Module[{\[Delta]dom={-7,Log10[0.74]},N\[Delta]s=10,rmaxintrp,\[Delta]min,NCof\[Delta],NCof\[Delta]intrp},
(*speed up \[Delta] root finding byinterpolating over rmax first*)

Print["test1"];
rmaxintrp = Interpolation[Table[{\[Delta]t,Log10@getrmax[ncof\[Delta]andr][10^\[Delta]t]},{\[Delta]t,\[Delta]dom[[1]],\[Delta]dom[[2]],(\[Delta]dom[[2]]-\[Delta]dom[[1]])/N\[Delta]s}],InterpolationOrder->1];

Print["test2"];
(*Print[Plot[rmaxintrp[\[Delta]t],{\[Delta]t,\[Delta]dom[[1]],\[Delta]dom[[2]]}]];
*)
(*\[Delta]min = binsearch[10^rmaxintrp[\[Delta]t],\[Delta]t,{-1,\[Delta]dom[[2]]},30];

Print[\[Delta]min];*)

NCof\[Delta] = 4 \[Pi] NIntegrate[r^2 ncof\[Delta]andr/."\[Delta]"->#,{r,0,Evaluate[10^rmaxintrp[Log10[#]]]}]&;

Print["test3"];

(*Print@Table[Log10@NCof\[Delta][10^\[Delta]t],{\[Delta]t,\[Delta]dom[[1]],\[Delta]dom[[2]]}];

Print[((NCof\[Delta][#]-NC)&[\[Delta]t]/.\[Delta]t->10^(\[Delta]dom[[2]]))];

(*Print[ListLinePlot[Table[Log10@NCof\[Delta][\[Delta]t],{\[Delta]t,\[Delta]dom[[1]],\[Delta]dom[[2]]}]]];*)

Print[NC];*)

NCof\[Delta]intrp = Interpolation[Table[{\[Delta]t,Log10@NCof\[Delta][10^\[Delta]t]},{\[Delta]t,\[Delta]dom[[1]],\[Delta]dom[[2]],Evaluate[(\[Delta]dom[[2]]-\[Delta]dom[[1]])/N\[Delta]s]}],InterpolationOrder->1];
Print["test4"];
(*Print[LogPlot[10^NCof\[Delta]intrp[\[Delta]t]-NC,{\[Delta]t,\[Delta]dom[[1]],\[Delta]dom[[2]]}(*,PlotRange->{0,12}*),PlotRange->All]];
*)
(*Print@Table[{\[Delta]t,10^NCof\[Delta]intrp[10^\[Delta]t]},{\[Delta]t,\[Delta]dom[[1]],\[Delta]dom[[2]],(\[Delta]dom[[2]]-\[Delta]dom[[1]])/N\[Delta]s}];*)

(*binsearch[(NCof\[Delta][\[Delta]t]-NC),\[Delta]t,\[Delta]dom,30][[1]]*)
binsearch[10^NCof\[Delta]intrp[Log10@\[Delta]t]-NC,\[Delta]t,\[Delta]dom,30][[1]]
]*)


\[Delta]ofNcap[ncof\[Delta]andr_,NC_]:= Module[{\[Delta]dom={-7,(*Log10[1.5]*)(*Log10[0.75]*)Log10[0.7494]},N\[Delta]s=10,rmaxintrp,\[Delta]min,NCof\[Delta],\[Delta]minNC,NCof\[Delta]intrp},
(*speed up \[Delta] root finding byinterpolating over rmax first*)

(*Print["test1"];*)
rmaxintrp = Interpolation[Table[{\[Delta]t,Log10@getrmax[ncof\[Delta]andr][10^\[Delta]t]},{\[Delta]t,\[Delta]dom[[1]],\[Delta]dom[[2]],(\[Delta]dom[[2]]-\[Delta]dom[[1]])/N\[Delta]s}],InterpolationOrder->1];

(*Print["test2"];*)
(*Print[Plot[rmaxintrp[\[Delta]t],{\[Delta]t,\[Delta]dom[[1]],\[Delta]dom[[2]]}]];
*)
(*\[Delta]min = binsearch[10^rmaxintrp[\[Delta]t],\[Delta]t,{-1,\[Delta]dom[[2]]},30];

Print[\[Delta]min];*)

NCof\[Delta] = 4 \[Pi] NIntegrate[r^2 ncof\[Delta]andr/."\[Delta]"->#,{r,0,Evaluate[10^rmaxintrp[Log10[#]]]}]&;

(*Print["test3"];*)

(*Print@Table[Log10@NCof\[Delta][10^\[Delta]t],{\[Delta]t,\[Delta]dom[[1]],\[Delta]dom[[2]]}];

Print[((NCof\[Delta][#]-NC)&[\[Delta]t]/.\[Delta]t->10^(\[Delta]dom[[2]]))];

(*Print[ListLinePlot[Table[Log10@NCof\[Delta][\[Delta]t],{\[Delta]t,\[Delta]dom[[1]],\[Delta]dom[[2]]}]]];*)

Print[NC];*)

NCof\[Delta]intrp = Interpolation[Table[{\[Delta]t,NCof\[Delta][10^\[Delta]t]},{\[Delta]t,\[Delta]dom[[1]],\[Delta]dom[[2]],Evaluate[(\[Delta]dom[[2]]-\[Delta]dom[[1]])/N\[Delta]s]}],InterpolationOrder->1];
(*Print["test4"];*)
(*
Print[LogPlot[NCof\[Delta]intrp[\[Delta]t] - NC,{\[Delta]t,\[Delta]dom[[1]],\[Delta]dom[[2]]}(*,PlotRange->{0,12}*),PlotRange->All]];
Print[LogPlot[NCof\[Delta]intrp[\[Delta]t],{\[Delta]t,\[Delta]dom[[1]],\[Delta]dom[[2]]}(*,PlotRange->{0,12}*),PlotRange->All]];

\[Delta]minNC = binsearch[NCof\[Delta]intrp[Log10@\[Delta]t],\[Delta]t,{-2,\[Delta]dom[[2]]},100,True][[1]];
Print[\[Delta]minNC];
Print[NCof\[Delta]intrp[\[Delta]dom[[2]]]];

Print[{rmaxintrp[Log10@\[Delta]minNC],NCof\[Delta]intrp[Log10@\[Delta]minNC]}];*)

(*Print@Table[{\[Delta]t,10^NCof\[Delta]intrp[10^\[Delta]t]},{\[Delta]t,\[Delta]dom[[1]],\[Delta]dom[[2]],(\[Delta]dom[[2]]-\[Delta]dom[[1]])/N\[Delta]s}];*)

(*binsearch[(NCof\[Delta][\[Delta]t]-NC),\[Delta]t,\[Delta]dom,30][[1]]*)
binsearch[NCof\[Delta]intrp[Log10@\[Delta]t]-NC,\[Delta]t,\[Delta]dom,30][[1]]
]


(* ::Subsubsection:: *)
(*get \[Delta] from a PDict*)


(*get\[Delta]fromPDict[pdict_,\[Alpha]Dby\[Alpha]_:1,fD_:0.05]:=Module[{meD,mpD,\[Beta]D,v0,\[Kappa],nD,nCap,NC,Ncapof\[Delta]table,\[Delta]t,rmax},
(*The case of \[Delta]>1 still needs to be sorted out*)

(*meD = "meD"/.pdict;
mpD = "mpD"/.pdict;
\[Beta]D = "\[Beta]D"/.pdict;
v0 = "v0"/.pdict;
\[Kappa] = "\[Kappa]"/.pdict;*)
{meD, mpD, \[Beta]D, v0, \[Kappa]}={"meD","mpD","\[Beta]D","v0","\[Kappa]"}/.pdict;

(*nD =Capture`naDM[fD,"meD" + "mpD"]/.pdict;
nCap=getncap["mpD" ,"meD" ,"\[Beta]D",nD,\[Alpha]Dby\[Alpha]]/.pdict;
*)
nD =Capture`naDM[fD,meD+mpD];
nCap=getncap[mpD, meD, \[Beta]D, nD, \[Alpha]Dby\[Alpha]];

NC = "Nc"/.pdict/."nD"->nD;

\[Delta]t=\[Delta]ofNcap[nCap,NC];
rmax=getrmax[nCap][\[Delta]t];

<|"\[Delta]"->\[Delta]t,"rmax"->rmax,"NC"->NC,"meD"->meD,"mpD"->mpD,"\[Beta]D"->\[Beta]D,"nD"->nD,"v0"->v0,"\[Kappa]"->\[Kappa],"fD"->fD,"\[Alpha]Dby\[Alpha]"->\[Alpha]Dby\[Alpha]|>
]*)


InterpolationFunctionQ[expr_]:=Head[expr]===InterpolatingFunction


get\[Delta]fromPDict[pdict_,\[Alpha]Dby\[Alpha]_:1,fD_:0.05,\[Delta]ofNcapintrp_:Null,\[Delta]dom_:{-5,Log10[1.5]}]:=Module[{meD,mpD,\[Beta]D,v0,\[Kappa],nD,nCap,NCE,NC,\[CapitalGamma]cap,\[CapitalGamma]evap,\[CapitalGamma]capE,\[CapitalGamma]evapE,log\[Tau]s,Ncapof\[Delta]table,\[Delta]t,rmax},

(*Print["in get \[Delta]"];*)
{meD, mpD, \[Beta]D, v0, \[Kappa]}={"meD","mpD","\[Beta]D","v0","\[Kappa]"}/.pdict;

nD =Capture`naDM[fD,meD+mpD];
nCap=getncap[mpD, meD, \[Beta]D, nD, \[Alpha]Dby\[Alpha]];

NCE = "Nc"/.pdict/."nD"->nD;
NC = "Nctot"/.pdict/."nD"->nD;

(*Print[NC];
Print["LeafCount[pdictwtherm] = ", LeafCount[pdict]];
Print["ByteCount[pdictwtherm]  = ", ByteCount[pdict]];
Print[pdict];*)

(*\[CapitalGamma]cap = "dNcMaxdt"/.pdict/."nD"->nD; (*capture rate averaged over Earth*)
\[CapitalGamma]evap = "\[CapitalGamma]total"/.pdict;(*per particle evaporation rate averaged over Earth*)*)

\[CapitalGamma]cap = "\[CapitalGamma]captot"/.pdict/."nD"->nD; (*total capture rate - including plasma*)
\[CapitalGamma]evap = "\[CapitalGamma]evaptot"/.pdict;(*per particle evaporation rate - including plasma*)
\[CapitalGamma]capE = "\[CapitalGamma]capE"/.pdict/."nD"->nD; (*total capture rate - Earth only*)
\[CapitalGamma]evapE = "\[CapitalGamma]evapE"/.pdict;(*per particle evaporation rate - Earth only*)

(*log\[Tau]s = "log\[Tau]s"/.pdict;*)

(*Print["before \[Delta]t line"];*)

(*\[Delta]t=\[Delta]ofNcap[nCap,NC];*)
\[Delta]t=If[!InterpolationFunctionQ[\[Delta]ofNcapintrp],\[Delta]ofNcap[nCap,NC],binsearch[\[Delta]ofNcapintrp[Log10@\[Delta]t]-NC,\[Delta]t,\[Delta]dom][[1]]];
rmax=getrmax[nCap][\[Delta]t];
(*\[Delta]t={Print@#[[1]],#[[2]]}[[2]]&@AbsoluteTiming[\[Delta]ofNcap[nCap,NC]];
rmax={Print@#[[1]],#[[2]]}[[2]]&@AbsoluteTiming[getrmax[nCap][\[Delta]t]];*)

(*Print["at the end"];*)

<|"\[Delta]"->\[Delta]t,"rmax"->rmax,"NCE"->NCE,"NC"->NC,"\[CapitalGamma]cap"->\[CapitalGamma]cap,"\[CapitalGamma]evap"->\[CapitalGamma]evap,"\[CapitalGamma]capE"->\[CapitalGamma]capE,"\[CapitalGamma]evapE"->\[CapitalGamma]evapE,"log\[Tau]s"->log\[Tau]s,"meD"->meD,"mpD"->mpD,"\[Beta]D"->\[Beta]D,"nD"->nD,"v0"->v0,"\[Kappa]"->\[Kappa],"fD"->fD,"\[Alpha]Dby\[Alpha]"->\[Alpha]Dby\[Alpha]|>
]


(* ::Subsubsection::Closed:: *)
(*get \[Delta] from a scan  - deprecated*)


Clear[Compute\[Delta]onScan]
Compute\[Delta]onScan[\[Delta]List_,v0ind_:1]:= Module[{return\[Delta]List,\[Delta]Listcurrent,PDict,v0,mratio,\[Delta]table},
(*
\[Delta]List - of the form: {<|"\[Delta]"->\[Delta],"file"->"\\path\\to\\file"|>}
*)
return\[Delta]List ={};

Monitor[
Do[

\[Delta]Listcurrent = \[Delta]List[[\[Delta]]];

PDict= Utilities`ReadIt[\[Delta]Listcurrent["file"]][[;;9,v0ind,;;]];

{v0, mratio} = {"v0","meD"/"mpD"}/.First@First[PDict];

(*\[Delta]table=Flatten[Table[Table[{Log10[("mpD" ("c")^2)/("JpereV")]/.PDict[[i,j]]/.Constants`SIConstRepl,Log10["\[Kappa]"]/.PDict[[i,j]],Log10@"\[Delta]"/.get\[Delta]fromPDict[PDict[[i,j]]]},{i,Dimensions[PDict][[1]]}],{j,Dimensions[PDict][[2]]}],1];*)
\[Delta]table=Flatten[Table[Table[{Log10[("mpD" ("c")^2)/("JpereV")]/.PDict[[i,j]]/.Constants`SIConstRepl,Log10["\[Kappa]"]/.PDict[[i,j]],Log10@"\[Delta]"/.get\[Delta]fromPDict[PDict[[i,j]]]},{i,Dimensions[PDict][[1]]}],{j,Dimensions[PDict][[2]]}],1];

AppendTo[\[Delta]Listcurrent,<|"\[Delta]table"->\[Delta]table,"v0"->v0,"meDbympD"->mratio|>];
AppendTo[return\[Delta]List,\[Delta]Listcurrent];

,{\[Delta],Length[\[Delta]List]}];
,{\[Delta],j,i}];

Print["Finished running the PDicts for: \n{v0,meD/mpD}=",First@First[{"v0","meD"/"mpD"}/.PDict]];

return\[Delta]List
]


(* ::Subsubsection::Closed:: *)
(*Get \[Delta] Truth  - deprecated*)


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


(* ::Subsubsection:: *)
(*get \[Delta] Dict from a scan *)


Clear[GetNcapof\[Delta]Interpolation]
GetNcapof\[Delta]Interpolation[mpD_, meD_, \[Beta]D_,fD_,\[Alpha]Dby\[Alpha]_,\[Delta]dom_:{-5,(*Log10[0.7494]*)Log10[1.5]}]:=Module[{nD,ncof\[Delta]andr,N\[Delta]=(*10*)5,\[Delta]points},

nD =Capture`naDM[fD,meD+mpD];
ncof\[Delta]andr=getncap[mpD, meD, \[Beta]D, nD, \[Alpha]Dby\[Alpha]];

(*\[Delta]points = Subdivide[\[Delta]dom[[1]],\[Delta]dom[[2]],N\[Delta]];*)
\[Delta]points = Table[Sign[i]i^2,{i,Subdivide[Sign[\[Delta]dom[[1]]]Sqrt[Abs@\[Delta]dom[[1]]],Sign[\[Delta]dom[[2]]]Sqrt[Abs@\[Delta]dom[[2]]],N\[Delta]]}]; (*weight the larger values of \[Delta] more heavily, since the \[Delta] is linear in Subscript[N, C] for small \[Delta]*)

(*Interpolation[Table[{log\[Delta],Ncapof\[Delta][ncof\[Delta]andr,10^log\[Delta]]},{log\[Delta],Subdivide[\[Delta]dom[[1]],\[Delta]dom[[2]],N\[Delta]]}],InterpolationOrder->1]*)
Interpolation[Table[{log\[Delta],Ncapof\[Delta][ncof\[Delta]andr,10^log\[Delta]]},{log\[Delta],\[Delta]points}],InterpolationOrder->1]
]


Clear[Compute\[Delta]DictonScan]
Compute\[Delta]DictonScan[\[Delta]List_,(*\[Sigma]dicte_,\[Sigma]dictNuc_,*)v0ind_:1,\[Alpha]Dby\[Alpha]_:1,fD_:0.05]:= Module[{return\[Delta]List,\[Delta]dom,\[Delta]Listcurrent,PDict,pdictwtherm,v0,mratio,\[Delta]table,\[Delta]ofNcapintrp},
(*
\[Delta]List - of the form: {<|"\[Delta]"->\[Delta],"file"->"\\path\\to\\file"|>,...}
*)
return\[Delta]List ={};

(*need to get \[Delta]dom and pass it*)

\[Delta]dom = {Log10@Min[#],Log10@Max[#]}&@("\[Delta]"/.\[Delta]List);

(*Print["test"];*)

Monitor[
Do[

\[Delta]Listcurrent = \[Delta]List[[\[Delta]]];

(*PDict= Utilities`ReadIt[\[Delta]Listcurrent["file"]][[;;9,v0ind,;;]];*)
PDict= Utilities`ReadIt[\[Delta]Listcurrent["file"]][[;;,v0ind,;;]];

(*Print["test test"];

Print[Head[PDict]];*)
(*Print[Utilities`ReadIt];*)

(*pdictwtherm = ThermTimescales`Get\[Tau]sFromPDict[PDict,"\[Delta]"/.\[Delta]Listcurrent,\[Sigma]dicte,\[Sigma]dictNuc,\[Alpha]Dby\[Alpha],fD]; *)
(*Print["bouta run: ",PDict];*)

pdictwtherm = Get\[Tau]sFromPDict[PDict,"\[Delta]"/.\[Delta]Listcurrent,(*\[Sigma]dicte,\[Sigma]dictNuc,*)\[Alpha]Dby\[Alpha],fD];

(*Print["test test"];

(*pdictwtherm = PDict;*)
(*Return[Print["Test over"],Module];*)
Print["LeafCount[PDict] = ", LeafCount[PDict]];
Print["ByteCount[PDict]  = ", ByteCount[PDict]];*)

{v0, mratio} = {"v0","meD"/"mpD"}/.First@First[PDict];

(*Print[pdictwtherm];*)
(*Print["LeafCount[pdictwtherm] = ", LeafCount[pdictwtherm]];
Print["ByteCount[pdictwtherm]  = ", ByteCount[pdictwtherm]];
Print[pdictwtherm];*)
(*\[Delta]table=Flatten[Table[Table[{Log10[("mpD" ("c")^2)/("JpereV")]/.PDict[[i,j]]/.Constants`SIConstRepl,Log10["\[Kappa]"]/.PDict[[i,j]],Log10@"\[Delta]"/.get\[Delta]fromPDict[PDict[[i,j]]]},{i,Dimensions[PDict][[1]]}],{j,Dimensions[PDict][[2]]}],1];*)
(*\[Delta]table=Flatten[Table[Table[get\[Delta]fromPDict[PDict[[i,j]]],{i,Dimensions[PDict][[1]]}],{j,Dimensions[PDict][[2]]}],\[Alpha]Dby\[Alpha],fD];*)

(*\[Delta]table=Flatten[Table[Table[get\[Delta]fromPDict[PDict[[i,j]],\[Alpha]Dby\[Alpha],fD],{i,Dimensions[PDict][[1]]}],{j,Dimensions[PDict][[2]]}],\[Alpha]Dby\[Alpha],fD];*)
(*\[Delta]table=Flatten[Table[Table[get\[Delta]fromPDict[PDict[[i,j]],\[Alpha]Dby\[Alpha],fD],{i,Dimensions[PDict][[1]]}],{j,Dimensions[PDict][[2]]}]];*)
(*Print["entering \[Delta]table"];*)
\[Delta]table=Monitor[Flatten[Table[	
			\[Delta]ofNcapintrp = GetNcapof\[Delta]Interpolation["mpD"/.#,"meD"/.#,"\[Beta]D"/.#,fD,\[Alpha]Dby\[Alpha],\[Delta]dom]&[pdictwtherm[[i,1]]];
			(*\[Delta]ofNcapintrp = {Print[#[[1]]],#[[2]]}[[2]]&@AbsoluteTiming@GetNcapof\[Delta]Interpolation["mpD"/.#,"meD"/.#,"\[Beta]D"/.#,fD,\[Alpha]Dby\[Alpha],\[Delta]dom]&[PDict[[i,1]]];*)
			(*Print[\[Delta]ofNcapintrp];*)
			Table[(*If[Log10@("\[Kappa]"/.pdictwtherm[[i,j]])>-10,Print["stuff"];Continue[]];*)
				get\[Delta]fromPDict[pdictwtherm[[i,j]],\[Alpha]Dby\[Alpha],fD,\[Delta]ofNcapintrp,\[Delta]dom]
				,{j,(*1*)Dimensions[pdictwtherm][[2]](*\[Kappa]*)}],
			{i,(*1*)Dimensions[pdictwtherm][[1]](*mpD*)}]]
			(*Can interpolate of getNcap in this \[Delta]table loop. We just need to put it in the mass loop (and make sure it is the outer one so we can use it for all \[Kappa])*)
			,{i,j}];

(*Print[\[Delta]table];*)

(*AppendTo[\[Delta]Listcurrent,<|"\[Delta]table"->\[Delta]table,"v0"->v0,"meDbympD"->mratio|>];
AppendTo[return\[Delta]List,\[Delta]Listcurrent];*)

\[Delta]Listcurrent=Append[\[Delta]Listcurrent,<|"\[Delta]table"->\[Delta]table,"v0"->v0,"meDbympD"->mratio|>];
return\[Delta]List=Append[return\[Delta]List,\[Delta]Listcurrent];

,{\[Delta],Length[\[Delta]List]}];

(*Return[Print["Test over"],Module];*)

,{\[Delta]}];

(*Return[Print["Test over"],Module];*)

(*Print["Finished running the PDicts for: \n{v0,meD/mpD}=",{"v0","meD"/"mpD"}/.First@First[pdictwtherm]];*)
Print["Finished running the PDicts for: \n{v0,meD/mpD}=",{"v0","meD"/"mpD"}/.return\[Delta]List[[1]]["\[Delta]table"]];

(*Print[return\[Delta]List];*)

return\[Delta]List
]


(* ::Subsubsection::Closed:: *)
(*Compute\[Delta]dictonscan - old*)


(*Clear[Compute\[Delta]DictonScan]
Compute\[Delta]DictonScan[\[Delta]List_,\[Sigma]dicte_,\[Sigma]dictNuc_,v0ind_:1,\[Alpha]Dby\[Alpha]_:1,fD_:0.05]:= Module[{return\[Delta]List,\[Delta]dom,\[Delta]Listcurrent,PDict,v0,mratio,\[Delta]table,\[Delta]ofNcapintrp},
(*
\[Delta]List - of the form: {<|"\[Delta]"->\[Delta],"file"->"\\path\\to\\file"|>,...}
*)
return\[Delta]List ={};

(*need to get \[Delta]dom and pass it*)

\[Delta]dom = {Log10@Min[#],Log10@Max[#]}&@("\[Delta]"/.\[Delta]List);

Monitor[
Do[

\[Delta]Listcurrent = \[Delta]List[[\[Delta]]];

(*PDict= Utilities`ReadIt[\[Delta]Listcurrent["file"]][[;;9,v0ind,;;]];*)
PDict= Utilities`ReadIt[\[Delta]Listcurrent["file"]][[;;,v0ind,;;]];

Print["test test"];

pdictwtherm = ThermTimescales`Get\[Tau]sFromPDict[PDict,"\[Delta]"/.\[Delta]Listcurrent,\[Sigma]dicte,\[Sigma]dictNuc,\[Alpha]Dby\[Alpha],fD];

Return[Print["Test over"],Module];

{v0, mratio} = {"v0","meD"/"mpD"}/.First@First[PDict];

(*\[Delta]table=Flatten[Table[Table[{Log10[("mpD" ("c")^2)/("JpereV")]/.PDict[[i,j]]/.Constants`SIConstRepl,Log10["\[Kappa]"]/.PDict[[i,j]],Log10@"\[Delta]"/.get\[Delta]fromPDict[PDict[[i,j]]]},{i,Dimensions[PDict][[1]]}],{j,Dimensions[PDict][[2]]}],1];*)
(*\[Delta]table=Flatten[Table[Table[get\[Delta]fromPDict[PDict[[i,j]]],{i,Dimensions[PDict][[1]]}],{j,Dimensions[PDict][[2]]}],\[Alpha]Dby\[Alpha],fD];*)

(*\[Delta]table=Flatten[Table[Table[get\[Delta]fromPDict[PDict[[i,j]],\[Alpha]Dby\[Alpha],fD],{i,Dimensions[PDict][[1]]}],{j,Dimensions[PDict][[2]]}],\[Alpha]Dby\[Alpha],fD];*)
(*\[Delta]table=Flatten[Table[Table[get\[Delta]fromPDict[PDict[[i,j]],\[Alpha]Dby\[Alpha],fD],{i,Dimensions[PDict][[1]]}],{j,Dimensions[PDict][[2]]}]];*)
\[Delta]table=Flatten[Table[
			\[Delta]ofNcapintrp = GetNcapof\[Delta]Interpolation["mpD"/.#,"meD"/.#,"\[Beta]D"/.#,fD,\[Alpha]Dby\[Alpha],\[Delta]dom]&[PDict[[i,1]]];
			(*\[Delta]ofNcapintrp = {Print[#[[1]]],#[[2]]}[[2]]&@AbsoluteTiming@GetNcapof\[Delta]Interpolation["mpD"/.#,"meD"/.#,"\[Beta]D"/.#,fD,\[Alpha]Dby\[Alpha],\[Delta]dom]&[PDict[[i,1]]];*)
			Table[get\[Delta]fromPDict[PDict[[i,j]],\[Alpha]Dby\[Alpha],fD,\[Delta]ofNcapintrp,\[Delta]dom]
				,{j,(*1*)Dimensions[PDict][[2]](*\[Kappa]*)}],
			{i,(*1*)Dimensions[PDict][[1]](*mpD*)}]];
			(*Can interpolate of getNcap in this \[Delta]table loop. We just need to put it in the mass loop (and make sure it is the outer one so we can use it for all \[Kappa])*)

AppendTo[\[Delta]Listcurrent,<|"\[Delta]table"->\[Delta]table,"v0"->v0,"meDbympD"->mratio|>];
AppendTo[return\[Delta]List,\[Delta]Listcurrent];

,{\[Delta],Length[\[Delta]List]}];

(*Return[Print["Test over"],Module];*)

,{\[Delta],i,j}];

(*Return[Print["Test over"],Module];*)

Print["Finished running the PDicts for: \n{v0,meD/mpD}=",First@First[{"v0","meD"/"mpD"}/.PDict]];

return\[Delta]List
]*)


(* ::Subsubsection:: *)
(*Scan over unfixed parameters*)


Clear[ScanOverUnfixedParameters]
ScanOverUnfixedParameters[\[Delta]ScanList_,\[Alpha]Dby\[Alpha]s_:{1,2},fDs_:{0.01,0.05}]:=Module[{\[Delta]Dict},
Do[
Print["Running: ",{v0ind,\[Alpha]Dby\[Alpha],fD}];
\[Delta]Dict=(*Quiet@*)Compute\[Delta]DictonScan[\[Delta]ScanList,v0ind,\[Alpha]Dby\[Alpha],fD];(*this was the culpable line (I think) it was expecting \[Sigma] dicts to be passed too.*)
(*Print[\[Delta]Dict];*)
(*SaveIt[NotebookDirectory[]<>StringForm["\[Delta]Dict_``_aDbya_``_fD_``",{v0ind,\[Alpha]Dby\[Alpha],fD}],\[Delta]Dict];*)
Capture`ExportDatFileToDir[\[Delta]Dict,ToString@StringForm["\[Delta]Dict_``_aDbya_``_fD_``",v0ind,\[Alpha]Dby\[Alpha],fD],"\[Delta]Dicts"];
,{v0ind,(*4*)1},{\[Alpha]Dby\[Alpha],\[Alpha]Dby\[Alpha]s},{fD,fDs}](*this system of having a hardcoded v0ind is stupid, this needs to be read in from the \[Delta]scan list some how or not at all. 
Ie. the accounting and printing can be done internally to Compute\[Delta]DictonScan*)
]


(* ::Subsubsection:: *)
(*Get \[Delta] Dict Truth*)


Clear[Get\[Delta]DictTruth]
Get\[Delta]DictTruth[\[Delta]Scan_]:=Module[{DEBUG=False,\[Delta]truthtable,\[Delta]outof\[Delta]intable,\[Delta]outof\[Delta]in,\[Delta]truth,\[Delta]Scantruth,keytable,keyintrp,keytruth},

(*take the \[Delta]Scans, and replace \[Delta] by the truth value*)
(**)

\[Delta]truthtable ={};

Monitor[Do[
(*Print[i];*)
(*Print["test 1"];*)
\[Delta]outof\[Delta]intable =Table[{Log10@"\[Delta]"/.\[Delta]Scan[[j]],Log10@"\[Delta]"/.Flatten[\[Delta]Scan[[j]]["\[Delta]table"]][[i]]},{j,Length[\[Delta]Scan]}];
(*Print[\[Delta]outof\[Delta]intable];*)
\[Delta]outof\[Delta]in = Interpolation[\[Delta]outof\[Delta]intable,InterpolationOrder->1];
(*\[Delta]truth = \[Delta]/.FindRoot[10^\[Delta]outof\[Delta]in[Log10@\[Delta]]-\[Delta],{\[Delta],0.01}]*)
\[Delta]truth = binsearch[10^\[Delta]outof\[Delta]in[Log10@\[Delta]]-\[Delta],\[Delta],{-5,(*Log10[1.5]*)(*Log10[0.75]*)Log10[0.7494]}][[1]];

(*Print["test 2"];*)

\[Delta]Scantruth = Flatten[\[Delta]Scan[[1]]["\[Delta]table"]][[i]];
\[Delta]Scantruth["\[Delta]"]=\[Delta]truth;

(*Print["test 3"];*)

(*Print[\[Delta]truth];
Print[Flatten[\[Delta]Scan[[1]]["\[Delta]table"]][[i]]];

Print["Before"];*)

(*Print[\[Delta]Scantruth];
Print[#&@(getncap["mpD","meD","\[Beta]D","nD","\[Alpha]Dby\[Alpha]"]/.\[Delta]Scantruth)];
(*Print[getrmax[getncap["mpD","meD","\[Beta]D","nD","\[Alpha]Dby\[Alpha]"]/.\[Delta]Scantruth][\[Delta]truth]]*)
Print@Ncapof\[Delta][#,\[Delta]truth]&@(getncap["mpD","meD","\[Beta]D","nD","\[Alpha]Dby\[Alpha]"]/.\[Delta]Scantruth(*//Quiet*));
Print@getrmax[#][\[Delta]truth]&@(getncap["mpD","meD","\[Beta]D","nD","\[Alpha]Dby\[Alpha]"]/.\[Delta]Scantruth(*//Quiet*));*)
(*{\[Delta]Scantruth["rmax"],\[Delta]Scantruth["NC"]}={getrmax[#]["\[Delta]"],Ncapof\[Delta][#,"\[Delta]"]}&@getncap["mpD","meD","\[Beta]D","nD","\[Alpha]Dby\[Alpha]"]/.\[Delta]Scantruth(*//Quiet*);*)
{\[Delta]Scantruth["rmax"],\[Delta]Scantruth["NC"]}={getrmax[#][\[Delta]truth],Ncapof\[Delta][#,\[Delta]truth]}&@(getncap["mpD","meD","\[Beta]D","nD","\[Alpha]Dby\[Alpha]"]/.\[Delta]Scantruth(*//Quiet*));



(*(\[Delta]Scantruth[#]&/@{"\[CapitalGamma]evap","\[CapitalGamma]cap"})*)\[Delta]Scantruth= ReplacePart[(*Flatten[\[Delta]Scan[[1]]["\[Delta]table"]][[i]]*)\[Delta]Scantruth,#->With[{key=#},
keytable =Table[{Log10@"\[Delta]"/.\[Delta]Scan[[j]],Log10@Max[(key/.Flatten[\[Delta]Scan[[j]]["\[Delta]table"]][[i]]),10^-100]},{j,Length[\[Delta]Scan]}];
keyintrp = Interpolation[keytable,InterpolationOrder->1];
keytruth = 10^keyintrp[Log10[\[Delta]truth]]]&/@{"\[CapitalGamma]evap","\[CapitalGamma]cap","\[CapitalGamma]evapE","\[CapitalGamma]capE","NCE"}];

(*If[i==1,Print@\[Delta]Scantruth];*)

(*Print["After"];*)

(*Print["test 4"];*)

AppendTo[\[Delta]Scantruth,get\[Phi]DAnalyticEstimates["mpD","meD","\[Beta]D","nD","\[Alpha]Dby\[Alpha]"]/.\[Delta]Scantruth];

If[DEBUG,If[(#[[1]]==13&&#[[2]]==-13)&@(Log10@{"mpD" ("c")^2/("JpereV"),"\[Kappa]"}/.\[Delta]Scan[[1]]["\[Delta]table"][[i]]/.Constants`SIConstRepl//N//Round),Print[i];Print[\[Delta]truth];Print[\[Delta]Scantruth]]];

AppendTo[\[Delta]truthtable,\[Delta]Scantruth];

,{i,(*1*)Length[Flatten[\[Delta]Scan[[1]]["\[Delta]table"]]]}],{i}];(*assuming all \[Delta]tables have the same ordering size*)

Gather[\[Delta]truthtable,("\[Kappa]"/.#1)==("\[Kappa]"/.#2)&]
]


(* ::Subsection:: *)
(*Potential outside the flat region*)


(* ::Subsubsection::Closed:: *)
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


(* ::Subsubsection::Closed:: *)
(*get analytic estimate for \[Phi]D *)


get\[Phi]DAnalyticEstimates[mpD_,meD_,\[Beta]D_,nF_,\[Alpha]Dby\[Alpha]_,\[Phi]g_:Capture`Vgrav] := Module[{niofr,npD,neD,\[Rho]S,\[Lambda]D,nC},
niofr = 1/(mi^2 Sqrt[2 \[Pi]] \[Beta]D^2) nF (mi \[Beta]D)^(3/2) (2 mi vesc \[Beta]D+E^(1/2 mi vesc^2 \[Beta]D) Sqrt[2 \[Pi]] Sqrt[mi \[Beta]D]-E^(1/2 mi vesc^2 \[Beta]D) Sqrt[mi] Sqrt[2 \[Pi]] Sqrt[\[Beta]D] Erf[(Sqrt[mi] vesc Sqrt[\[Beta]D])/Sqrt[2]]);
npD=niofr/.vesc->Sqrt[-((2 Vgp)/mi)]/.mi->mpD;
neD=niofr/.vesc->Sqrt[-((2 Vge)/mi)]/.mi->meD;
(*npD=niofr/.vesc->Sqrt[-((2 Vgp)/mi)]/.Vgp->mpD \[Phi]g[r]/.mi->mpD;*)
(*neD=niofr/.vesc->Sqrt[-((2 Vge)/mi)]/.Vge->meD \[Phi]g[r]/.mi->meD;*)
\[Rho]S = ("e")/("\[Epsilon]0") Sqrt[\[Alpha]Dby\[Alpha]](npD-neD)/.{Vgp->mpD \[Phi]g[r],Vge->meD \[Phi]g[r]};
(*\[Lambda]D = (("e")^2/("\[Epsilon]0")\[Alpha]Dby\[Alpha](D[npD,Vgp]+D[neD,Vge]))^(-1/2)/.{Vgp->mpD \[Phi]g[r],Vge->mpD \[Phi]g[r]};*)
\[Lambda]D =(-(("e")^2/("\[Epsilon]0"))\[Alpha]Dby\[Alpha](D[npD,Vgp]+D[neD,Vge])/.{Vgp->mpD \[Phi]g[r],Vge->meD \[Phi]g[r]})^(-1/2);
nC = (mpD "G" "ME" )/(("e")^2/("\[Epsilon]0") \[Alpha]Dby\[Alpha]((4 \[Pi])/3 ("rE")^3)) Boole["rE" - r > 0] - (npD-neD)/.Constants`SIConstRepl/.Constants`EarthRepl;

(*the potentials here need to include \[Delta] (for the densities) or you will get negative nC. It is fine for the \[Rho]S and \[Lambda]D*)

<|"\[Phi]D"->\[Rho]S \[Lambda]D^2,"\[Rho]S"->\[Rho]S,"\[Lambda]D"->\[Lambda]D,"npD"->npD,"neD"->neD,"nC"->nC|>
]


(* ::Subsubsection::Closed:: *)
(*get full \[Phi]D potential *)


Clear[get\[Phi]D]
get\[Phi]D[\[Phi]g_:Capture`Vgrav]:=Module[{\[Phi]S,\[Delta]\[Phi],\[Phi]Dh,\[Phi]D},
\[Phi]S = ("\[Lambda]D")^2 "\[Rho]S" Piecewise[{{0,r<"rmax"},{(1-E^(-(r/("\[Lambda]D"))) Cosh[("rmax")/("\[Lambda]D")]),r>="rmax">"rE"},{1+("\[Lambda]D" E^(-(("rE")/("\[Lambda]D"))-r/("\[Lambda]D")))/(2 Min["rE",r])-("\[Lambda]D" E^(-(Abs[-"rE"+r]/("\[Lambda]D"))))/(2 Min["rE",r])-("rmax" E^(-(r/("\[Lambda]D"))) Cosh[("rmax")/("\[Lambda]D")])/Min["rE",r]+("\[Lambda]D" E^(-(r/("\[Lambda]D"))) Sinh[("rmax")/("\[Lambda]D")])/Min["rE",r],"rE">="rmax"}}];
\[Delta]\[Phi] = ("mpD")/("e" Sqrt["\[Alpha]Dby\[Alpha]"]) (-"\[Delta]" ("vesc")^2/2-\[Phi]g[r]);
\[Phi]Dh = ("rmax")/r E^-((r-"rmax")/("\[Lambda]D")) (\[Delta]\[Phi]-\[Phi]S/.r->"rmax");

(*Print[(\[Delta]\[Phi]-\[Phi]S/.r->"rmax")];*)

\[Phi]D = Piecewise[{{\[Delta]\[Phi],r<"rmax"},{\[Phi]Dh+\[Phi]S,r>="rmax"}}];

<|"\[Phi]S"->\[Phi]S,"\[Delta]\[Phi]"->\[Delta]\[Phi],"\[Phi]Dh"->\[Phi]Dh,"\[Phi]D"->\[Phi]D|>
]


(* ::Subsubsection::Closed:: *)
(*Check for Maxima*)


(* ::Text:: *)
(*this doesn't work because of some scope issues with the variable r*)


Clear[checkforVeffMaxima]
checkforVeffMaxima[paramrepl_,\[Phi]g_:Capture`Vgrav,consts_:Join[Association[Constants`SIConstRepl],Constants`EarthRepl]]:= Module[{DEBUG=True,Pots,Lexsq,Lexsqder,rmax,\[Lambda]D,MaximaIntervalLength},

(*Pots = {"meD"\[Phi]g[r]-"e" "\[Phi]D","mpD"\[Phi]g[r]+"e" "\[Phi]D"}/.get\[Phi]D[]/.paramrepl/.Private`r->r/.consts;*)
Pots = {"meD"\[Phi]g[r]-"e" "\[Phi]D","mpD"\[Phi]g[r]+"e" "\[Phi]D"}/.get\[Phi]D[]/.paramrepl/.Global`r->r/.consts;
Lexsq = r^3 D[Pots,r]//Expand;
Lexsqder = D[Lexsq,r]//Expand;

rmax = "rmax"/.paramrepl/.consts;

(*Print["what's the symbol: ",Cases[Pots, s_Symbol :> Context[s] <> SymbolName[s], \[Infinity]]];*)
(*Print[Pots/.r->rmax];*)
(*Print[Pots/.Global`r->r/.r->rmax];*)

(*\[Lambda]D = Function[{r},("\[Lambda]D"/.paramrepl/.r->Private`r)](*/.Private`r->rmax/.r->rmax*);*)
(*\[Lambda]D = Simplify@Evaluate@PiecewiseExpand@("\[Lambda]D"/.paramrepl);
(*Print[Lexsq/.r->rmax+\[Lambda]D];
Print[Lexsqder/.r->rmax+\[Lambda]D];*)
Print[rmax];
Print[Simplify@Evaluate@(\[Lambda]D/.r->rmax)];
(*Print[\[Lambda]D/.Private`r->rmax];*)
Print[Cases[\[Lambda]D, r, \[Infinity]]];*)
\[Lambda]D = "\[Lambda]D"/.paramrepl/.Global`r->rmax/.consts;
(*Print["\[Lambda]D"/.paramrepl/.r->rmax];
Print["\[Lambda]D"/.paramrepl/.Global`r->rmax];*)
(*\[Lambda]D = Evaluate[PiecewiseExpand@("\[Lambda]D" /. paramrepl)];

Print["rmax = ", rmax];
Print["\[Lambda]D(rmax) = ", \[Lambda]D /. r -> rmax];  (* or \[Lambda]D[rmax] if it\[CloseCurlyQuote]s a Function *)
Print["r appears in \[Lambda]D: ", Cases[\[Lambda]D, r, \[Infinity]]];
Print["Private`r appears in \[Lambda]D: ", Cases[\[Lambda]D, Private`r, \[Infinity]]];
Print["what's the symbol: ",Cases["\[Lambda]D" /. paramrepl, s_Symbol :> Context[s] <> SymbolName[s], \[Infinity]]];
Print[Head["\[Lambda]D" /. paramrepl]];*)

(*MaximaIntervalLength = NIntegrate[Boole[Lexsq>0]Boole[Lexsqder<0],{r,rmax,rmax+10 \[Lambda]D}];*)
MaximaIntervalLength = NIntegrate[Boole[Lexsq[[#]]>0]Boole[Lexsqder[[#]]<0],{r,rmax,rmax+10 \[Lambda]D}]&/@{1,2};

(*Print["what's the symbol: ",Cases[Evaluate@Pots/.Global`r->r/.consts, s_Symbol :> Context[s] <> SymbolName[s], \[Infinity]]];*)

If[DEBUG,
Print@Plot[Evaluate@(Pots/.consts),{r,rmax,rmax+5 \[Lambda]D}];
Print@Plot[Evaluate@(r^-3 Lexsq/.consts),{r,rmax,rmax+5 \[Lambda]D}];
Print@Plot[Evaluate@(Lexsq/.consts),{r,rmax,rmax+5 \[Lambda]D}];
Print@Plot[Evaluate@(Lexsqder/.consts),{r,rmax,rmax+5 \[Lambda]D}];
Pots = {"meD"\[Phi]g[r]-"e" ("\[Phi]D"(*-"\[Rho]S" ("\[Lambda]D")^2*)),"mpD"\[Phi]g[r]+"e" ("\[Phi]D"(*-"\[Rho]S" ("\[Lambda]D")^2*))}/.get\[Phi]D[]/.paramrepl/.Global`r->r/.consts;
Lexsq = r^3 D[Pots,r]//Expand;
Lexsqder = D[Lexsq,r]//Expand;
Print[Pots/.r->rmax];
Print@Plot[Evaluate@((Sign[#]Abs@Log10@Abs@#)&@Expand@Simplify@(Pots/.consts)),{r,rmax,rmax+50 \[Lambda]D}];
Print@Plot[Evaluate@((Sign[#]Abs@Log10@Abs@#)&@Expand@Simplify@(r^-3 Lexsq/.consts)),{r,rmax,rmax+50 \[Lambda]D}];
Print@Plot[Evaluate@((Sign[#]Abs@Log10@Abs@#)&@Expand@Simplify@(Lexsq/.consts)),{r,rmax,rmax+50 \[Lambda]D}];
Print@Plot[Evaluate@((Sign[#]Abs@Log10@Abs@#)&@Expand@Simplify@(Lexsqder/.consts)),{r,rmax,rmax+50 \[Lambda]D}];
(*Pots = {(*"meD"\[Phi]g[r]*)-"e" ("\[Phi]D"(*-"\[Rho]S" ("\[Lambda]D")^2*)),(*"mpD"\[Phi]g[r]+*)"e" ("\[Phi]D"(*-"\[Rho]S" ("\[Lambda]D")^2*))}[[1]]/.get\[Phi]D[]/.paramrepl/.Global`r->r/.consts;
Lexsq = r^3 D[Pots,r]//Expand;
Lexsqder = D[Lexsq,r]//Expand;
Print[Pots/.r->rmax];
Print@Plot[Evaluate@(Pots/.consts),{r,rmax,rmax+5 \[Lambda]D}];
Print@Plot[Evaluate@(r^-3Lexsq/.consts),{r,rmax,rmax+5 \[Lambda]D}];
Print@Plot[Evaluate@(Lexsq/.consts),{r,rmax,rmax+5 \[Lambda]D}];
Print@Plot[Evaluate@(Lexsqder/.consts),{r,rmax,rmax+5 \[Lambda]D}];*)
(*Print@Plot[Evaluate@(Pots/.consts),{r,rmax-5\[Lambda]D,rmax+5 \[Lambda]D}];
Print@Plot[Evaluate@(Lexsq/.consts),{r,rmax-5\[Lambda]D,rmax+5 \[Lambda]D}];
Print@Plot[Evaluate@(Lexsqder/.consts),{r,rmax-5\[Lambda]D,rmax+5 \[Lambda]D}];
Pots = {(*"meD"\[Phi]g[r]*)-"e" ("\[Phi]D"(*-"\[Rho]S" ("\[Lambda]D")^2*)),(*"mpD"\[Phi]g[r]*)+"e" ("\[Phi]D"(*-"\[Rho]S" ("\[Lambda]D")^2*))}/.get\[Phi]D[]/.paramrepl/.Global`r->r/.consts;
Lexsq = r^3 D[Pots,r];
Lexsqder = D[Lexsq,r];
Print[Pots/.r->rmax];
Print@Plot[Evaluate@(Pots/.consts),{r,rmax-5\[Lambda]D,rmax+5 \[Lambda]D}];
Print@Plot[Evaluate@(Lexsq/.consts),{r,rmax-5\[Lambda]D,rmax+5 \[Lambda]D}];
Print@Plot[Evaluate@(Lexsqder/.consts),{r,rmax-5\[Lambda]D,rmax+5 \[Lambda]D}];*)
(*Print@LogLogPlot[Evaluate@Pots/.Global`r->r/.consts,{r,rmax-5\[Lambda]D,rmax+5 \[Lambda]D}];
Print@LogLogPlot[Evaluate@Lexsq/.Global`r->r/.consts,{r,rmax-5\[Lambda]D,rmax+5 \[Lambda]D}];
Print@LogLogPlot[Evaluate@Lexsqder/.Global`r->r/.consts,{r,rmax-5\[Lambda]D,rmax+5 \[Lambda]D}];
Pots = {"meD"\[Phi]g[r]-"e" ("\[Phi]D"-"\[Rho]S" ("\[Lambda]D")^2),"mpD"\[Phi]g[r]+"e" ("\[Phi]D"-"\[Rho]S" ("\[Lambda]D")^2)}/.get\[Phi]D[]/.paramrepl/.Private`r->r/.consts;
Lexsq = r^3 D[Pots,r];
Lexsqder = D[Lexsq,r];
Print@LogLogPlot[Evaluate@Pots/.Global`r->r/.consts,{r,rmax-5\[Lambda]D,rmax+5 \[Lambda]D}];
Print@LogLogPlot[Evaluate@Lexsq/.Global`r->r/.consts,{r,rmax-5\[Lambda]D,rmax+5 \[Lambda]D}];
Print@LogLogPlot[Evaluate@Lexsqder/.Global`r->r/.consts,{r,rmax-5\[Lambda]D,rmax+5 \[Lambda]D}];*)
];



(*Print[\[Lambda]D];
Print[rmax];
Print[Boole[r^3D[-1/r^3,r]>0]/.r->rmax+\[Lambda]D];
Print[NIntegrate[Boole[r^3D[-1/r^3,r]>0]Boole[D[r^3D[-1/r^3,r],r]<0],{r,rmax,rmax+10 \[Lambda]D}]&/@{1,2}];
*)
(*I NEED TO DEAL WITH SCOPE, UNTIL I DO, I WON'T TRUST THESE RESULTS*)

MaximaIntervalLength
]


(* ::Section:: *)
(*Plasma Interactions*)


(* ::Subsection::Closed:: *)
(*Define Base Functions (Simple relaxation times)*)


(* ::Text:: *)
(*Should change the inputs to relative and target speeds*)


(* ::Subsubsection::Closed:: *)
(*Module for Captured Relaxation Times*)


(* ::Input::Initialization:: *)
Clear[Get\[Tau]C]
Get\[Tau]C[Logm_,v0kms_?NumericQ,\[Delta]t_:10^-4]:=Module[{ncap,ncapofr,rmax,Ntot,(*mTs,vrels,*)rat,ravg,\[Tau]relCintegrand,\[Tau]relC,meDbympD=10^-4,fD=0.05,\[Alpha]Dby\[Alpha]=1,DEBUG=False},
(*ncap = Screening`getncap["mpD","meD","\[Beta]D","nD",\[Alpha]Dby\[Alpha]]/."nD"->Capture`naDM[fD,"mpD"+"meD"]/."\[Beta]D"->(*3/("mpD"("v0")^2)*)3 /("mpD" ("v0")^2)/.{Private`v->"v0"}/."v0"->v0kms 10^3/."meD"->"mpD"meDbympD/."mpD"->("JpereV")/("c")^2 10^Logm/.Constants`SIConstRepl/.Constants`EarthRepl;*)
ncap = getncap["mpD","meD","\[Beta]D","nD",\[Alpha]Dby\[Alpha]]/."nD"->Capture`naDM[fD,"mpD"+"meD"]/."\[Beta]D"->(*3/("mpD"("v0")^2)*)3 /("mpD" ("v0")^2)/.{Private`v->"v0"}/."v0"->v0kms 10^3/."meD"->"mpD"meDbympD/."mpD"->("JpereV")/("c")^2 10^Logm/.Constants`SIConstRepl/.Constants`EarthRepl;

(*Print[ncap];*)

(*rmax=Screening`getrmax[ncap][\[Delta]t];*)
rmax=getrmax[ncap][\[Delta]t];

ncapofr=Evaluate[ncap/.{"\[Delta]"->\[Delta]t,Private`r->#}]&;

(*mTs = {"mpD","meD"};
vrels = {v0kms 10^3,v0kms 10^3 Sqrt[("mpD")/("meD")]};*)

(*Ntot = Screening`Ncapof\[Delta][ncap,\[Delta]t];*)
Ntot = Ncapof\[Delta][ncap,\[Delta]t];

ravg =Re[ (4 \[Pi])/Ntot NIntegrate[r^3 ncap/."\[Delta]"->\[Delta]t/.Private`r->r,{r,0,rmax}]];

\[Tau]relCintegrand=(v /(("m\[Chi]")/2 v^2) CollisionRate`dEkdrcoul )^-1/."nF"->ncap/.{"m\[Chi]"->"mpD","mT"->"mpD"}/.{"\[Alpha]D"->"\[Alpha]","meD"->"mpD"meDbympD}/.{"nD"->Capture`naDM[fD,"mpD"]}/."\[Beta]D"->3/("mpD" ("v0")^2)/.Private`v->v/.v->"vesc" \[Delta]t /."v0"->v0kms 10^3/."mpD"->("JpereV")/("c")^2 10^Logm/.\[Delta]->\[Delta]t/.Constants`SIConstRepl/.Constants`EarthRepl/.{Private`r->#,"\[Delta]"->\[Delta]t}&;
(*Print@Plot[Log10@\[Tau]relCintegrand[10^logr],{logr,0,Log10@rmax},PlotRange->{0,2}];*)
\[Tau]relC  = 1/rmax NIntegrate[\[Tau]relCintegrand[r],{r,0,0.95 rmax}];(*averaged relaxation time through an orbit assuming right through Earth, so this is really an lower bound since some orbits will miss earth*)
(*Print[\[Tau]relC];*)


If[DEBUG,
Print[\[Delta]t];

(*Print[rmax];*)
Print["ncap ",ncap/.{"\[Delta]"->\[Delta]t,Private`r->0.95rmax}];

(*Print["ncap ",ncap/.{"\[Delta]"->\[Delta]t,Private`r->0.1rmax}];*)
Print[Log10@rmax];
(*ncapofr=ncap/.{"\[Delta]"->\[Delta]t,Private`r->#}&;*)
(*Print[Evaluate[Log10@rmax]];*)
(*Print[Evaluate[0.1rmax]];*)
(*Print@Plot[Capture`Vgrav[r],{r,0,Evaluate[10^2rmax]},PlotRange->All];*)
(*Print[ncapofr[4 10^7]];*)
Print@Plot[Log10@Abs@ncapofr[10^logr],{logr,0,Log10@Evaluate[10^2 rmax]},PlotRange->All,PlotLabel->"ncap"];
Print@Plot[ncapofr[10^logr],{logr,Log10@Evaluate[rmax]-1,Log10@Evaluate[rmax]+1},PlotRange->All,PlotLabel->"ncap"];
Print@Plot[ncapofr'[10^logr],{logr,Log10@Evaluate[rmax]-1,Log10@Evaluate[rmax]+1},PlotRange->All];
Print@Plot[Log10@Abs@\[Tau]relCintegrand[10^logr],{logr,0,Log10@Evaluate[10^2 rmax]},PlotRange->All,PlotLabel->"\[Tau]relc"];
(*Print@Print[Evaluate[D[ncapofr[r],r]]/.r->rmax];*)
];

<|"rmax"->rmax,"ravg"->ravg,"Ntot"->Ntot,"\[Tau]relC"->\[Tau]relC|>
]


(* ::Subsubsection::Closed:: *)
(*Module for Ambient Relaxation times*)


(* ::Input::Initialization:: *)
Clear[Get\[Tau]A]
Get\[Tau]A[Logm_,v0kms_?NumericQ,\[Delta]t_:10^-4(*,vrat_:10^-4*)]:=Module[{ncap,Ntot,rmax,rat,mTs,vrels,ravg,\[Tau]relCintegrand,\[Tau]relA,meDbympD=10^-4,fD=0.05,\[Alpha]Dby\[Alpha]=1},

(*ncap = Screening`getncap["mpD","meD","\[Beta]D","nD",\[Alpha]Dby\[Alpha]]/."nD"->Capture`naDM[fD,"mpD"]/."\[Beta]D"->3/("mpD" ("v0")^2)/.{Private`v->"v0",v->"v0"}/."v0"->v0kms 10^3/."meD"->"mpD"meDbympD/."mpD"->("JpereV")/("c")^2 10^Logm/.Constants`SIConstRepl/.Constants`EarthRepl;

rmax=Screening`getrmax[ncap][\[Delta]t];

Ntot = Screening`Ncapof\[Delta][ncap,\[Delta]t];*)
ncap = getncap["mpD","meD","\[Beta]D","nD",\[Alpha]Dby\[Alpha]]/."nD"->Capture`naDM[fD,"mpD"]/."\[Beta]D"->3/("mpD" ("v0")^2)/.{Private`v->"v0",v->"v0"}/."v0"->v0kms 10^3/."meD"->"mpD"meDbympD/."mpD"->("JpereV")/("c")^2 10^Logm/.Constants`SIConstRepl/.Constants`EarthRepl;

rmax=getrmax[ncap][\[Delta]t];

Ntot = Ncapof\[Delta][ncap,\[Delta]t];

ravg =Re[ (4 \[Pi])/Ntot NIntegrate[r^3 ncap/."\[Delta]"->\[Delta]t/.Private`r->r,{r,0,rmax}]];

mTs = {"mpD","meD"};
vrels = {v0kms 10^3,v0kms 10^3 Sqrt[("mpD")/("meD")]};

\[Tau]relA=Sum[(v /(("m\[Chi]")/2 v^2) CollisionRate`dEkdrcoul )^-1/.{"m\[Chi]"->"mpD","mT"->mTs[[Tind]]}/.Private`v->v/.v->vrels[[Tind]],{Tind,Length@mTs}]/.{"\[Alpha]D"->"\[Alpha]"\[Alpha]Dby\[Alpha],"meD"->"mpD"meDbympD}/."nF"->Capture`naDM[fD,"mpD"]/."mpD"->("JpereV")/("c")^2 10^Logm(*"vesc"vrat*)(* \[Delta]t *)/."v0"->v0kms 10^3/.Constants`SIConstRepl/.Constants`EarthRepl;


(*\[Tau]relC  = 1/rmax NIntegrate[\[Tau]relCintegrand[r],{r,0,0.99 rmax}];(*averaged relaxation time through an orbit*)*)

(*Print@Plot[Log10@\[Tau]relCintegrand[10^logr],{logr,0,Log10@rmax},PlotRange->{0,2}];*)

<|"rmax"->rmax,"ravg"->ravg,"Ntot"->Ntot,"\[Tau]relA"->\[Tau]relA|>
]


(* ::Subsection:: *)
(*Timescale Modules (Simple relaxation times)*)


(* ::Subsubsection::Closed:: *)
(*\[Tau]AA*)


Clear[Get\[Tau]AA]
Get\[Tau]AA[logm_,log\[Kappa]_,pdict_,\[Delta]in_]:=Module[{flatpdict,pdictind,\[Tau]relAdict,\[Tau]relAA,\[Tau]nearEarth},

flatpdict = Flatten[pdict,1];

pdictind=FirstPosition[{Round@Log10@("mpD"//N),Round@Log10@("\[Kappa]"//N)}/.flatpdict,{Round@Log10@(10^logm ("JpereV")/("c")^2/.Constants`SIConstRepl),log\[Kappa]}][[1]];
(*\[Delta]temp="\[Delta]"/.pdict[[\[Delta]dictind]];*)

\[Tau]relAdict=Get\[Tau]A[logm,10^-3 "v0"/.flatpdict[[pdictind]],\[Delta]in(*,("v0")/("vesc")/.flatpdict[[pdictind]]/.Constants`EarthRepl*)];
\[Tau]relAA = "\[Tau]relA"/.\[Tau]relAdict;

\[Tau]nearEarth = 10 ("rE")/("v0")/.flatpdict[[pdictind]]/.Constants`EarthRepl;
<|"logm"->logm,"log\[Kappa]"->log\[Kappa],"logassrat"->Log10[\[Tau]relAA/\[Tau]nearEarth],"log\[Tau]rel"->Log10@\[Tau]relAA,"log\[Tau]comp"->Log10@\[Tau]nearEarth|>
]


(* ::Subsubsection::Closed:: *)
(*\[Tau]CC*)


Clear[Get\[Tau]CC]
Get\[Tau]CC[logm_,log\[Kappa]_,pdict_,\[Delta]in_]:=Module[{flatpdict,pdictind,\[Tau]relCdict,\[Tau]relCC,\[Tau]evap},

flatpdict = Flatten[pdict,1];

pdictind=FirstPosition[{Round@Log10@("mpD"//N),Round@Log10@("\[Kappa]"//N)}/.flatpdict,{Round@Log10@(10^logm ("JpereV")/("c")^2/.Constants`SIConstRepl),log\[Kappa]}][[1]];
(*\[Delta]temp="\[Delta]"/.pdict[[pdictind]];*)

\[Tau]relCdict=Get\[Tau]C[logm,10^-3 "v0"/.flatpdict[[pdictind]],\[Delta]in];
\[Tau]relCC = "\[Tau]relC"/.\[Tau]relCdict;

(*\[Tau]evap =("\[Kappa]")^-2 ("\[CapitalGamma]evap")^-1/.flatpdict[[pdictind]];*)
(*\[Tau]evap = (("\[Kappa]")^-2 /.#)If[AnyTrue[StringQ/@("\[CapitalGamma]evap"/.#),#&],"\[CapitalGamma]total"/.#,"\[CapitalGamma]evap"/.#]^-1&@flatpdict[[pdictind]];*)
\[Tau]evap = (("\[Kappa]")^-2 /.#)If[StringQ@("\[CapitalGamma]evap"/.#),"\[CapitalGamma]total"/.#,"\[CapitalGamma]evap"/.#]^-1&@flatpdict[[pdictind]];

<|"logm"->logm,"log\[Kappa]"->log\[Kappa],"logassrat"->Log10[\[Tau]relCC/\[Tau]evap],"log\[Tau]rel"->Log10@\[Tau]relCC,"log\[Tau]comp"->Log10@\[Tau]evap|>
]


(* ::Subsubsection::Closed:: *)
(*\[Tau]CE*)


Clear[Get\[Tau]CE]
Get\[Tau]CE[logm_,log\[Kappa]_,pdict_,\[Delta]in_,\[Sigma]dicte_,\[Sigma]dictNuc_]:=Module[{flatpdict,pdictind,eind,Nucind,\[Tau]relCdict,\[Tau]relCE,\[Tau]evap},

flatpdict = Flatten[pdict,1];

pdictind=FirstPosition[{Round@Log10@("mpD"//N),Round@Log10@("\[Kappa]"//N)}/.flatpdict,{Round@Log10@(10^logm ("JpereV")/("c")^2/.Constants`SIConstRepl),log\[Kappa]}][[1]];
(*\[Delta]temp="\[Delta]"/.pdict[[pdictind]];*)

eind = FirstPosition["m\[Chi]"/.\[Sigma]dicte,10^logm ("JpereV")/("c")^2/.Constants`SIConstRepl][[1]];
Nucind = FirstPosition["m\[Chi]"/.\[Sigma]dictNuc,10^logm ("JpereV")/("c")^2/.Constants`SIConstRepl][[1]];

\[Tau]relCdict=Get\[Tau]C[logm,10^-3 "v0"/.flatpdict[[pdictind]],\[Delta]in];

Print[\[Tau]relCdict];

(*\[Tau]relCE = "\[Tau]relC"/.\[Tau]relCdict;*)
\[Tau]relCE=Quiet[((("rE")/("ravg"))^-1 v\[Chi](("\[Kappa]")^2/.flatpdict[[pdictind]])Quiet@Max[10^-100,("dEdlofr"/.Capture`GetTotalELFunction[\[Sigma]dicte[[eind]],\[Sigma]dictNuc[[Nucind]]])[1 10^6,v\[Chi]]]/((10^logm ("JpereV")/("c")^2/.Constants`SIConstRepl) v\[Chi]^2/2))^-1/.v\[Chi]->"vesc"\[Delta]in/.flatpdict[[pdictind]]/.Constants`EarthRepl/.\[Tau]relCdict];
(*\[Tau]relCE=1;*)
(*Print[\[Tau]relCE];

Print[Max[10^-100,("dEdlofr"/.Quiet@Capture`GetTotalELFunction[\[Sigma]dicte[[eind]],\[Sigma]dictNuc[[Nucind]]])[1 10^6,"vesc"\[Delta]in/.flatpdict[[pdictind]]/.Constants`EarthRepl]]];
*)
(*\[Tau]evap =("\[Kappa]")^-2 ("\[CapitalGamma]evap")^-1/.flatpdict[[pdictind]];*)
(*\[Tau]evap =(("\[Kappa]")^-2 /.#)If[AnyTrue[StringQ/@("\[CapitalGamma]evap"/.#),#&],"\[CapitalGamma]total"/.#,"\[CapitalGamma]evap"/.#]^-1&@flatpdict[[pdictind]];*)
\[Tau]evap =(("\[Kappa]")^-2 /.#)If[StringQ@("\[CapitalGamma]evap"/.#),"\[CapitalGamma]total"/.#,"\[CapitalGamma]evap"/.#]^-1&@flatpdict[[pdictind]];

<|"logm"->logm,"log\[Kappa]"->log\[Kappa],"logassrat"->Log10[\[Tau]relCE/\[Tau]evap],"log\[Tau]rel"->Log10@\[Tau]relCE,"log\[Tau]comp"->Log10@\[Tau]evap|>
]


(* ::Subsubsection::Closed:: *)
(*\[Tau]AE*)


Clear[Get\[Tau]AE]
Get\[Tau]AE[logm_,log\[Kappa]_,pdict_,\[Delta]in_,\[Sigma]dicte_,\[Sigma]dictNuc_]:=Module[{flatpdict,pdictind,eind,Nucind,\[Tau]relAdict,\[Tau]relAE,\[Tau]EarthPass},

flatpdict = Flatten[pdict,1];

pdictind=FirstPosition[{Round@Log10@("mpD"//N),Round@Log10@("\[Kappa]"//N)}/.flatpdict,{Round@Log10@(10^logm ("JpereV")/("c")^2/.Constants`SIConstRepl),log\[Kappa]}][[1]];
(*\[Delta]temp="\[Delta]"/.pdict[[pdictind]];*)

eind = FirstPosition["m\[Chi]"/.\[Sigma]dicte,10^logm ("JpereV")/("c")^2/.Constants`SIConstRepl][[1]];
Nucind = FirstPosition["m\[Chi]"/.\[Sigma]dictNuc,10^logm ("JpereV")/("c")^2/.Constants`SIConstRepl][[1]];

(*\[Tau]relAdict=Quiet@Get\[Tau]C[logm,10^-3 "v0"/.flatpdict[[pdictind]],\[Delta]in];*)

\[Tau]relAE=Quiet[(v\[Chi](("\[Kappa]")^2/.flatpdict[[pdictind]])Quiet@Max[10^-100,("dEdlofr"/.Capture`GetTotalELFunction[\[Sigma]dicte[[eind]],\[Sigma]dictNuc[[Nucind]]])[1 10^6,v\[Chi]]]/((10^logm ("JpereV")/("c")^2/.Constants`SIConstRepl) v\[Chi]^2/2))^-1/.v\[Chi]->"v0"/.flatpdict[[pdictind]]/.Constants`EarthRepl(*/.\[Tau]relAdict*)];
(*\[Tau]relAE=1;*)
\[Tau]EarthPass = ("rE")/("v0")/.flatpdict[[pdictind]]/.Constants`EarthRepl;

<|"logm"->logm,"log\[Kappa]"->log\[Kappa],"logassrat"->Log10[\[Tau]relAE/\[Tau]EarthPass],"log\[Tau]rel"->Log10@\[Tau]relAE,"log\[Tau]comp"->Log10@\[Tau]EarthPass|>
]


(* ::Subsubsection::Closed:: *)
(*\[Tau]CA*)


Clear[Get\[Tau]CA]
Get\[Tau]CA[logm_,log\[Kappa]_,pdict_,\[Delta]in_]:=Module[{flatpdict,pdictind,\[Delta]temp,\[Tau]relCAdict,\[Tau]relCA,\[Tau]evap},

flatpdict = Flatten[pdict,1];

pdictind=FirstPosition[{Round@Log10@("mpD"//N),Round@Log10@("\[Kappa]"//N)}/.flatpdict,{Round@Log10@(10^logm ("JpereV")/("c")^2/.Constants`SIConstRepl),log\[Kappa]}][[1]];
(*\[Delta]temp="\[Delta]"/.pdict[[pdictind]];*)

(*\[Tau]relCdict=Get\[Tau]C[logm,10^-3"v0"/.pdict[[pdictind]],\[Delta]temp];
\[Tau]relCC = "\[Tau]relC"/.\[Tau]relCdict;*)

\[Tau]relCAdict=Get\[Tau]A[logm,10^-3 "v0"/.flatpdict[[pdictind]],\[Delta]in(*,("v0")/("vesc")/.flatpdict[[pdictind]]/.Constants`EarthRepl*)];
\[Tau]relCA = "\[Tau]relA"/.\[Tau]relCAdict;

(*\[Tau]evap =("\[Kappa]")^-2 ("\[CapitalGamma]evap")^-1/.flatpdict[[pdictind]];*)
\[Tau]evap = (("\[Kappa]")^-2 /.#)If[StringQ@("\[CapitalGamma]evap"/.#),"\[CapitalGamma]total"/.#,"\[CapitalGamma]evap"/.#]^-1&@flatpdict[[pdictind]];

<|"logm"->logm,"log\[Kappa]"->log\[Kappa],"logassrat"->Log10[\[Tau]relCA/\[Tau]evap],"log\[Tau]rel"->Log10@\[Tau]relCA,"log\[Tau]comp"->Log10@\[Tau]evap|>
]


(* ::Subsection:: *)
(*Rate Modules (Hardcore Capture and Evaporation)*)


(* ::Subsubsection::Closed:: *)
(*Get Probe Functions*)


Clear[GetProbeFunctions]
GetProbeFunctions[]:=Module[{Ficap,Fievap,FipList},

Ficap =ni Sqrt[(2 \[Beta]i mi)/\[Pi]] HeavisideTheta[1/4 (mj/mi z + Abs[u])^2-zesci^2](E^(-4 (\[Beta]i mi)/(\[Beta]j mj) zltcap^2)-E^(-4 (\[Beta]i mi)/(\[Beta]j mj) zgtcap^2))/.{zltcap->Sqrt[Max[zesci^2,mj/mi u z]],zgtcap ->Sqrt[Min[zesci^2 + mj/mi u z, 1/4 (mj/mi z + u)^2]]};
Fievap = ni Sqrt[(2 \[Beta]i mi)/\[Pi]] HeavisideTheta[1/4 (mj/mi z + Abs[u])^2-zesci^2](E^(-4 (\[Beta]i mi)/(\[Beta]j mj) zltcap^2)-E^(-4 (\[Beta]i mi)/(\[Beta]j mj) zgtcap^2))/.{zltcap->Sqrt[Max[0,zesci^2-mj/mi Abs[u] z]],zgtcap ->Sqrt[Min[zesci^2 , 1/4 (mj/mi z -Abs[u])^2]]};
FipList = <|"cap"->Ficap,"evap"->Fievap|>/.{\[Beta]i->"\[Beta]i",\[Beta]j->"\[Beta]j",ni->"ni",nj->"nj",zesci->"zesci",zescj->"zescj",zD->"zD"};

FipList
]


(* ::Subsubsection::Closed:: *)
(*Get Target Response Functions*)


Clear[GetTargetResponses]
GetTargetResponses[]:=Module[{V\[Chi]Therm,V\[Chi]Bnd,V\[Chi]Amb,V\[Chi]List,TargetjList},

V\[Chi]Therm = zD^2/z^2 1/(4 Sqrt[\[Pi]]z) ((#/.a->u+z)-(#/.a->u-z))&@ {2 Sqrt[\[Pi]]DawsonF[a],-\[Pi] E^-a^2}(*{Re,Im}*);
V\[Chi]Bnd =zD^2/z^2 1/(4 Sqrt[\[Pi]]z) ((#/.a->u+z)-(#/.a->u-z))&@ {2 zescj/a (1-E^-a^2)-(E^-zescj^2+E^-a^2)Log[Abs[(a-zescj)/(a+zescj)]],-\[Pi](E^-zescj^2+E^-a^2)HeavisideTheta[zescj - Abs[a]]};(*{Re,Im} - valid for zesc<<1*)
V\[Chi]Amb = V\[Chi]Therm-V\[Chi]Bnd//Simplify;(*{Re,Im} - valid for zesc<<1*)

V\[Chi]List = <|"therm"->V\[Chi]Therm,"bnd"->V\[Chi]Bnd,"amb"->V\[Chi]Amb|>;
TargetjList = Association@Table[key->("e")^2/("\[Epsilon]0") 2/(1-E^(-4 u z)) V\[Chi]List[key][[2]]/Total@(({1,0}+V\[Chi]List[key])^2),{key,{"therm","bnd","amb"}}]/.{\[Beta]i->"\[Beta]i",\[Beta]j->"\[Beta]j",ni->"ni",nj->"nj",zesci->"zesci",zescj->"zescj",zD->"zD"};
TargetjList
]


(* ::Subsubsection::Closed:: *)
(*Get Plasma Rates*)


Clear[GetPlasmaRates]
GetPlasmaRates[pdict_,\[Delta]_,rmax_,\[Alpha]Dby\[Alpha]_:1,fD_:0.01]:=Module[{FipList,TargetjList,j=0,ratetable,VpDT,VeDT,zD,zesci,zescj,\[Beta]i,\[Beta]j,npD,neD,nC,nD,ni,nj,plasmaparamlist,(*badreg,*)integrand,integrandnum,rate},

FipList = GetProbeFunctions[];
TargetjList = GetTargetResponses[];

(*It looks like npD, neD are not being loaded in.
There's also a Null["cap"] somewhere*)
(*what is the right prescription for npD and neD? Well, the consistent \[Delta] value it won't matter. We can use \[Delta] to get it*)

(*get\[Phi]DAnalyticEstimates[mpD,meD,\[Beta]D,nF,\[Alpha]Dby\[Alpha],\[Phi]g];*)
(*get\[Phi]DAnalyticEstimates["mpD","meD","\[Beta]D","nD","\[Alpha]Dby\[Alpha]"]/.pdict;*)

(*Print[Head["npD"/.pdict]];
Print["npD"/.pdict];*)

nD = If[Head["nD"/.pdict]==String,Capture`naDM[fD,"mpD"/.pdict],"nD"];

{npD,neD,nC} = If[Head["npD"/.pdict]==String||Head["neD"/.pdict]==String||Head["nC"/.pdict]==String,{"npD","neD","nC"}/.(get\[Phi]DAnalyticEstimates["mpD","meD","\[Beta]D","nD","\[Alpha]Dby\[Alpha]"]/.pdict),{"npD","neD","nC"}]/."nD"->nD/."\[Alpha]Dby\[Alpha]"->\[Alpha]Dby\[Alpha];
(*neD = If[Head["neD"/.pdict]==String,(get\[Phi]DAnalyticEstimates["mpD","meD","\[Beta]D","nD","\[Alpha]Dby\[Alpha]"]/.pdict)["npD"],"npD"];*)

(*Print[{npD,neD,nC}];*)

(*Print[{npD,neD}];*)
(*Print[nD];*)

(*Print["nD"/.pdict];*)
(*Print[Keys@pdict];*)

(*trackingind = 0;*)
ratetable = {};
Do[

If[jkey=="therm",Continue[]];
(*If[jkey=="amb",Continue[]];*)
(*If[jkey=="bnd",Continue[]];*)
(*If[{ikey,jkey,mT}=={"cap","amb","mpD"},Continue[]];*)

(*If[{ikey,jkey,mT}!={"evap","amb","meD"}&&{ikey,jkey,mT}!={"cap","amb","mpD"},Continue[]];*)
(*If[{ikey,jkey,mT}!={"evap","amb","mpD"}&&{ikey,jkey,mT}!={"cap","amb","meD"},Continue[]];*)

{VpDT,VeDT}= {-\[Delta] ("mpD")/2 "vesc","mpD"Capture`Vgrav[rmax/2]+\[Delta] ("mpD")/2 "vesc"};

If[ikey =="cap",
\[Beta]i = "\[Beta]D"/.pdict;
(*ni =  "npD"/.pdict/.Private`r->rmax/2/.{Private`Vgp->VpDT}/.pdict/.Constants`EarthRepl//Simplify;*)
ni =  npD/.pdict/.Private`r->rmax/2/.{Private`Vgp->VpDT}/.pdict/.Constants`EarthRepl//Simplify;
,(*evap*)
\[Beta]i = "\[Beta]core"/.Constants`EarthRepl;
(*ni= "nC"/.pdict/.Private`r->rmax/2/.{Private`Vgp->VpDT,Private`Vge->VeDT}/.pdict/.Constants`EarthRepl//Simplify;*)
ni= nC/.pdict/.Private`r->rmax/2/.{Private`Vgp->VpDT,Private`Vge->VeDT}/.pdict/.Constants`EarthRepl//Simplify;
];

zesci = Sqrt[\[Beta]i/4 \[Delta] ("mpD")/2 ("vesc")^2]/.Constants`EarthRepl/.pdict;(*assumes probe is Subscript[m, Subscript[p, D]]*)


If[jkey =="bnd",
\[Beta]j = "\[Beta]core"/.Constants`EarthRepl;
(*nj =  "nC"/.pdict/.Private`r->rmax/2/.{Private`Vgp->VpDT,Private`Vge->VeDT}/.pdict/.Constants`EarthRepl//Simplify;*)
nj =  nC/.pdict/.Private`r->rmax/2/.{Private`Vgp->VpDT,Private`Vge->VeDT}/.pdict/.Constants`EarthRepl//Simplify;
,(*amb*)
\[Beta]j = "\[Beta]D";
(*nj = {"npD","neD"}[[If[mT=="mpD",1,2]]]/.pdict/.Private`r->rmax/2/.{Private`Vgp->VpDT,Private`Vge->VeDT}/.pdict/.Constants`EarthRepl//Simplify;*)
nj = {npD,neD}[[If[mT=="mpD",1,2]]]/.pdict/.Private`r->rmax/2/.{Private`Vgp->VpDT,Private`Vge->VeDT}/.pdict/.Constants`EarthRepl//Simplify;
	];
zescj= Sqrt[-(\[Beta]j/4)If[mT=="mpD",VpDT,VeDT]]/.Constants`EarthRepl/.pdict;
zD = Sqrt[(\[Beta]j ("\[HBar]")^2)/(8 mT) (2 nj \[Beta]j ("e")^2/("\[Epsilon]0"))]/.Constants`SIConstRepl/.pdict;

(*Print[ni];
Print[nj];*)

plasmaparamlist = {"\[Beta]i"->\[Beta]i,"\[Beta]j"->\[Beta]j,"ni"->ni,"nj"->nj,"zesci"->zesci,"zescj"->zescj,"zD"->zD,"\[Delta]"->\[Delta],"rmax"->rmax}/.pdict;
(*Print[plasmaparamlist];*)


integrand=\[Alpha]Dby\[Alpha]/(\[Pi]^2 \[Beta]j ("\[HBar]")^2) FipList[ikey] TargetjList[jkey] (*HeavisideTheta[x-Abs[y]]/.{u->(x+y)/2,z->(x-y)/2}*)/.plasmaparamlist/.mj->mT/.mi->"mpD"/.pdict/.Constants`SIConstRepl;
integrandnum[ut_?NumericQ,zt_?NumericQ] := integrand/.{u->ut,z->zt};
(*integrandnum[ut_,zt_]/;VectorQ[{ut,zt},NumericQ]:= integrand/.{u->ut,z->zt};*)

(*Print[integrand];*)

(*Print[integrand/.{x->1,y->10^-3}];*)
(*If[(integrand/.{x->1,y->10^-3})<0,
Print[ikey];
Print[jkey];
Print[mT];
Print@DensityPlot[integrand,{x,0,1.1},{y,0,10^-2},PlotLegends->Automatic,PlotRange->All];
Print[FipList[ikey]/.plasmaparamlist/.mj->mT/.mi->"mpD"/.pdict/.Constants`SIConstRepl/.{u->(x+y)/2,z->(x-y)/2}/.{x->1,y->10^-3}];
Print[TargetjList[jkey]/.plasmaparamlist/.mj->mT/.mi->"mpD"/.pdict/.Constants`SIConstRepl/.{u->(x+y)/2,z->(x-y)/2}/.{x->1,y->10^-3}];
];*)

(*badreg = ImplicitRegion[u==z,{u,z}];*)

(*rate = 1/2NIntegrate[integrand(*[ikey,jkey]*),{x,0,\[Infinity]},{y,-#,#},Method->"LocalAdaptive"]&[If[jkey=="bnd",zescj,1]];*)
(*rate = NIntegrate[integrand(*[ikey,jkey]*),{u,z}\[Element] ImplicitRegion[Abs[u-z]<=#&&u>0&&z>0,{u,z}],(*Method->"LocalAdaptive",*)AccuracyGoal->3,PrecisionGoal->3]&[If[jkey=="bnd",zescj,2]];*)(*slow but most stable*)

(*rate = NIntegrate[integrand(*[ikey,jkey]*),{u,z}\[Element] ImplicitRegion[(*10^-1Min[zD,zescj,zesci]<=Abs[u-z]&&*)Abs[u-z]<=#&&u>0&&z>0,{u,z}],(*Method->"LocalAdaptive",*)AccuracyGoal->3,PrecisionGoal->3]&[If[jkey=="bnd",zescj,2]];
*)(*current*)
rate = NIntegrate[integrand(*[ikey,jkey]*),{z,0,100},{u,Max[z-#,0],z+#},(*Method->"LocalAdaptive",*)AccuracyGoal->3,PrecisionGoal->3]&[If[jkey=="bnd",zescj,2]];(*only one that seems to work (ie. is quick enough and accurate) for amb and bnd at small masses*)
(*z is IR dominated (ie. by z << 1) so 0 - 100 is fine. *)
(*rate = NIntegrate[integrand(*[ikey,jkey]*),{z,0,100},{u,Max[z-#,0],z+#},Method->"LocalAdaptive",AccuracyGoal->3,PrecisionGoal->3]&[If[jkey=="bnd",zescj,2]];*)
(*rate = NIntegrate[integrandnum[u,z](*[ikey,jkey]*),{z,0,100},{u,Max[z-#,0],z+#},Method->"LocalAdaptive",AccuracyGoal->3,PrecisionGoal->3]&[If[jkey=="bnd",zescj,2]];*)


(*rate = NIntegrate[integrandnum[u,z](*[ikey,jkey]*),{u,z}\[Element]RegionDifference[ ImplicitRegion[Abs[u-z]<=#&&u>0&&z>0,{u,z}],Thickening[badreg]](*,Method->"LocalAdaptive"*)(*,AccuracyGoal->3,PrecisionGoal->3*)(*Method->{(*"MonteCarlo","MaxPoints"->10^5*)"SymbolicProcessing"->0}*),WorkingPrecision->MachinePrecision,Exclusions->{u==z}]&[If[jkey=="bnd",zescj,2]];
*)

(*Print[zD/10];*)
(*rate = NIntegrate[integrandnum[u,z](*[ikey,jkey]*),{u,z}\[Element] ImplicitRegion[(*(zD/10)<*)Abs[u-z]<=#&&u>0&&z>0&&Min[zD,zescj,zesci]/10<Abs[u-z]&&(1-Min[zD,zescj,zesci]/10>u+z||u+z>1+Min[zD,zescj,zesci]/10),{u,z}](*,Method->"LocalAdaptive"*),AccuracyGoal->10,PrecisionGoal->10,WorkingPrecision->40(*Method->{(*"MonteCarlo","MaxPoints"->10^5*)"SymbolicProcessing"->0}*)(*,WorkingPrecision->MachinePrecision*)]&[If[jkey=="bnd",zescj,2]];*)(*current in progress, fast but unstabel*)
(*rate = NIntegrate[integrand(*[ikey,jkey]*),{u,z}\[Element] ImplicitRegion[Abs[u-z]<=#&&u>0&&z>0&&Min[zD,zescj,zesci]/10<Abs[u-z]&&(1-Min[zD,zescj,zesci]/10>u+z||u+z>1+Min[zD,zescj,zesci]/10),{u,z}],(*Method->"LocalAdaptive",*)AccuracyGoal->3,PrecisionGoal->3]&[If[jkey=="bnd",zescj,2]];*)

(*rate = NIntegrate[integrand(*[ikey,jkey]*),{u,0,\[Infinity]},{z,0,\[Infinity]},Method->"LocalAdaptive",AccuracyGoal->3,PrecisionGoal->3]&[If[jkey=="bnd",zescj,2]];*)

(*Print[{ikey,jkey,mT}];
Print@DensityPlot[Log10@integrand,{u,0,1.1},{z,0,1.2},PlotLegends->Automatic,PlotRange->All];
Print@DensityPlot[Log10@integrandnum,{u,0,1.1},{z,0,1.2},PlotLegends->Automatic,PlotRange->All];
Print@DensityPlot[Log10@FipList[ikey] /.plasmaparamlist/.mj->mT/.mi->"mpD"/.pdict/.Constants`SIConstRepl,{u,0,11},{z,0,12},PlotLegends->Automatic,PlotRange->All];
Print@DensityPlot[Log10@TargetjList[jkey]  /.plasmaparamlist/.mj->mT/.mi->"mpD"/.pdict/.Constants`SIConstRepl,{u,0,11},{z,0,12},PlotLegends->Automatic,PlotRange->All];
*)
(*Print[ nC/.pdict/.Private`r->r/.{Private`Vgp->VpDT,Private`Vge->VeDT}/.pdict/.Constants`EarthRepl/.r->(rmax/2)//Simplify];*)

If[False(*(rate)<0*),
Print[ikey];
Print[jkey];
Print[mT];
Print[{Log10@"\[Kappa]",Log10["mpD" ("c")^2/("JpereV")]}/.pdict/.Constants`SIConstRepl//N];
Print@DensityPlot[Log10@integrand,{u,0,100},{z,0,100},PlotLegends->Automatic,PlotRange->All];
Print@DensityPlot[Log10[-integrand],{u,0,100},{z,0,100},PlotLegends->Automatic,PlotRange->All];
Print[rmax/2/.pdict];
Print[LogLogPlot[nC/.pdict/.Private`r->r/.{Private`Vgp->VpDT,Private`Vge->VeDT}/.pdict/.Constants`EarthRepl//Simplify,{r,0,rmax}]];
Print[Plot[{-\[Delta] ("mpD")/2 "vesc","mpD"Capture`Vgrav[r]+\[Delta] ("mpD")/2 "vesc"}/.pdict/.Constants`EarthRepl//Simplify,{r,0,rmax}]];
(*Print[FipList[ikey]/.plasmaparamlist/.mj->mT/.mi->"mpD"/.pdict/.Constants`SIConstRepl/.{u->(x+y)/2,z->(x-y)/2}/.{x->1,y->10^-3}];
Print[TargetjList[jkey]/.plasmaparamlist/.mj->mT/.mi->"mpD"/.pdict/.Constants`SIConstRepl/.{u->(x+y)/2,z->(x-y)/2}/.{x->1,y->10^-3}];*)
];

AppendTo[ratetable,<|"rate"->Abs@rate,"keys"->{ikey,jkey,mT},"params"->plasmaparamlist(*,"sign"->Sign[rate]*)|>];
(*the Abs on rate is a band aid, I haven't been able to diagnose the occasional incorrect sign. Maybe related to low accuracy / precision goals?*)

j= j+1;

,{ikey,{"cap","evap"}},{jkey,{"therm","bnd","amb"}},{mT,{"mpD","meD"}}];

ratetable
]


(* ::Subsection:: *)
(*Pipeline Modules*)


(* ::Subsubsection::Closed:: *)
(*Thoughts and TODOs*)


(* ::Text:: *)
(*The place to include these timescales I think would be in Get\[Delta]dicttruth *)
(*	I think the right way to implement would be to create a module to recompute NC given the new therm timescales . For now leave them as is, we can leave them as is to get the pipeline sorted . *)
(*   		Later we' ll need to update the cross - section, also modify the rates to account for kinematics / distribution stuff etc . Maybe we can use the scaling estimate we were doing for eg the capture rate . *)
(*   	This will include the two new contributions to the evaporation rate from self scattering and ambient - captured scattering .*)
(*   	We should also as a side quest estimate the ambient self interaction rate and find a way to compare to the capture rate from the Earth . We would like to know where each dominates . *)


(* ::Text:: *)
(**)


(* ::Text:: *)
(*- Actually NO! we need do do this during the scanoverunfixedparameters step since that computes \[Delta]*)
(*	Ie. in get\[Delta]fromaPDict we compute \[Delta]out using NC *)


(* ::Subsubsection::Closed:: *)
(*Get therm timescale table - old*)


(* ::Text:: *)
(**** Decide where and how to include this in the pipeline then come back to this.*)


(* ::Text:: *)
(*I think I need to also pass \[Delta] since it hasn't been included in the pdict yet*)


(*Clear[Get\[Tau]sFromPDict]
Get\[Tau]sFromPDict[pdict_,\[Delta]_,\[Sigma]dicte_,\[Sigma]dictNuc_]:= Module[{logms,log\[Kappa]s,\[Tau]stable,\[Tau]AA,\[Tau]CC,\[Tau]CA,\[Tau]CE,\[Tau]AE,log\[Tau]sdict,\[Tau]evaptot},

{logms ,log\[Kappa]s}=Log10@{Union@Flatten@(("m\[Chi]" ("c")^2)/("JpereV")/.pdict/.Constants`SIConstRepl),Union@Flatten@("\[Kappa]"/.pdict)};

\[Tau]stable=Flatten[Table[

(*can probably do the index search globally and pass it as needed.*)
\[Tau]AA = Get\[Tau]AA[logm,log\[Kappa],pdict,\[Delta]];
\[Tau]CC = Get\[Tau]CC[logm,log\[Kappa],pdict,\[Delta]];
\[Tau]CA = Get\[Tau]CA[logm,log\[Kappa],pdict,\[Delta]];
\[Tau]CE = Get\[Tau]CE[logm,log\[Kappa],pdict,\[Delta],\[Sigma]dicte,\[Sigma]dictNuc];
\[Tau]AE = Get\[Tau]AE[logm,log\[Kappa],pdict,\[Delta],\[Sigma]dicte,\[Sigma]dictNuc];

(*add these to the pdict as an association*)

log\[Tau]sdict = <|"\[Tau]AA"->"log\[Tau]rel"/.\[Tau]AA,"\[Tau]CC"->"log\[Tau]rel"/.\[Tau]CC,"\[Tau]CA"->"log\[Tau]rel"/.\[Tau]CA,"\[Tau]CE"->"log\[Tau]rel"/.\[Tau]CE,"\[Tau]AE"->"log\[Tau]rel"/.\[Tau]AE,"\[Tau]Epass"->"log\[Tau]comp"/.\[Tau]AE,"\[Tau]evapE"->"log\[Tau]comp"/.\[Tau]CE|>;

Print[log\[Tau]sdict];

\[Tau]evaptot = Total[10^-#&/@({"\[Tau]CC","\[Tau]CA","\[Tau]evapE"}/.log\[Tau]sdict)];

Print[\[Tau]evaptot];

pdictind=FirstPosition[{Round@Log10@("mpD"//N),Round@Log10@("\[Kappa]"//N)}/.flatpdict,{Round@Log10@(10^logm ("JpereV")/("c")^2/.Constants`SIConstRepl),log\[Kappa]}][[1]];


(*recompute NC - maybe separate out the code segment we want in Capture and call it here.*)

,{logm,logms[[;;1]]},{log\[Kappa],log\[Kappa]s[[;;1]]}],1];

Print[\[Tau]stable];

Print[\[Tau]AA];
Print[\[Tau]CC];
Print[\[Tau]CA];

\[Tau]stable
]*)


(* ::Subsubsection:: *)
(*Get therm timescale table*)


Clear[Get\[Tau]sFromPDict]
Get\[Tau]sFromPDict[pdict_,\[Delta]_,(*\[Sigma]dicte_,\[Sigma]dictNuc_,*)\[Alpha]Dby\[Alpha]_:1,fD_:0.05]:= Module[{pdictnew,logm,log\[Kappa],ratetable,\[Tau]stable,\[Tau]AA,\[Tau]CC,\[Tau]CA,\[Tau]CE,\[Tau]AE,log\[Tau]sdict,rmax,\[CapitalGamma]pC,\[CapitalGamma]pA,NE,NA,NCfrom\[Delta],\[CapitalGamma]evapE,\[CapitalGamma]capE,\[CapitalGamma]evapP,\[CapitalGamma]capP,\[CapitalGamma]evaptot,\[CapitalGamma]captot,tE,timeenoughforequil,Nctot,newstuffassoc(*,timeenoughforequilE,NcE*)},

(*Print@Dimensions[pdict];*)

pdictnew = {};

(*Print@Dimensions[pdict];
Print["in get \[Tau]s"];*)

\[Tau]stable=Monitor[Flatten[Table[

AppendTo[pdictnew,{}];

Table[

(*Print["testy test 0"];

Print[Head[pdict]];
Print[Head[pdict[[mind,\[Kappa]ind]]]];*)

(*{logm ,log\[Kappa]}=Log10@{(("m\[Chi]" ("c")^2)/("JpereV")/.pdict[[mind,\[Kappa]ind]]/.Constants`SIConstRepl),("\[Kappa]"/.pdict[[mind,\[Kappa]ind]])};*)
{logm ,log\[Kappa]}=Log10@{(("mpD" ("c")^2)/("JpereV")/.pdict[[mind,\[Kappa]ind]]/.Constants`SIConstRepl),("\[Kappa]"/.pdict[[mind,\[Kappa]ind]])};

(*Print[{logm,log\[Kappa]}];*)
(*If[log\[Kappa]>-10, Continue[]];*)

(*Print["testy test 1"];

Print["Head[logm] = ", Head[logm]];
Print["NumericQ[logm] = ", NumericQ[logm]];
Print["LeafCount[logm] = ", LeafCount[logm]];
Print["ByteCount[logm]  = ", ByteCount[logm]];*)

(*
Print[Dimensions[logm]];
Print[logm[[1]]];
Print[logm[[2]]];
Print[Dimensions[log\[Kappa]]];
Print[log\[Kappa][[1]]];
Print[log\[Kappa][[2]]];*)

(*Print[pdict[[mind,\[Kappa]ind]]];*)

(*Print[{logm ,log\[Kappa]}];*)

(*Print["testy test"];*)

(*Print[Keys@pdict[[mind,\[Kappa]ind]]];*)

rmax = getrmax[getncap["mpD","meD","\[Beta]D",Capture`naDM[fD,"mpD"],\[Alpha]Dby\[Alpha]]/.pdict[[mind,\[Kappa]ind]]][\[Delta]];

(*Print["rmax is ", rmax];*)

ratetable = GetPlasmaRates[pdict[[mind,\[Kappa]ind]],\[Delta],rmax,\[Alpha]Dby\[Alpha],fD];

(*Print[ratetable];*)
(*Print@Gather[ratetable,#1["keys"][[1]]==#2["keys"][[1]]&];*)

(*Print[Dimensions@ratetable];*)

(*can probably do the index search globally and pass it as needed.*)
(*need to update these with the relative velocity rather than the velocity of the probe.*)
(*\[Tau]AA = Get\[Tau]AA[logm,log\[Kappa],pdict,\[Delta]];
\[Tau]CC = Get\[Tau]CC[logm,log\[Kappa],pdict,\[Delta]];
\[Tau]CA = Get\[Tau]CA[logm,log\[Kappa],pdict,\[Delta]];
\[Tau]CE = Get\[Tau]CE[logm,log\[Kappa],pdict,\[Delta],\[Sigma]dicte,\[Sigma]dictNuc];
\[Tau]AE = Get\[Tau]AE[logm,log\[Kappa],pdict,\[Delta],\[Sigma]dicte,\[Sigma]dictNuc];*)

(*add these to the pdict as an association*)

(*log\[Tau]sdict = <|"\[Tau]AA"->"log\[Tau]rel"/.\[Tau]AA,"\[Tau]CC"->"log\[Tau]rel"/.\[Tau]CC,"\[Tau]CA"->"log\[Tau]rel"/.\[Tau]CA,"\[Tau]CE"->"log\[Tau]rel"/.\[Tau]CE,"\[Tau]AE"->"log\[Tau]rel"/.\[Tau]AE,"\[Tau]Epass"->"log\[Tau]comp"/.\[Tau]AE,"\[Tau]evapE"->"log\[Tau]comp"/.\[Tau]CE|>;
*)
(*log\[Tau]sdict = <|"\[Tau]AA"->0.02,"\[Tau]CC"->0.04,"\[Tau]CA"->0.03,"\[Tau]CE"->0.01,"\[Tau]AE"->0.001,"\[Tau]Epass"->0.00001,"\[Tau]evapE"->0.004|>;

Print[log\[Tau]sdict];*)

(*rmax = Screening`getrmax[Screening`getncap["mpD","meD","\[Beta]D",Capture`naDM[fD,"mpD"],\[Alpha]Dby\[Alpha]]/.pdict[[mind,\[Kappa]ind]]][\[Delta]];*)

(*Print[Capture`naDM[fD,"mpD"]/.pdict[[mind,\[Kappa]ind]]];*)
(*
These should really be rates that are integrated over r 

WE NEED TO DO THIS BEFORE SCANNING OR WE WILL BE SIGNIFICANTLY OFF.
*)

(*the missing piece here is to find the total per particle evaporation rate. At the mo they are total rates. *)
(*to do that we just need to divide by the total expected number of captured charges*)
(*(again, point is that they will agree with the rate definition once a solution is found) - so this means just divide through by 4 \[Pi] \[Integral] \[DifferentialD]r r^2Subscript[n, C](r), ie. just get NC of \[Delta]*)
(*Ie. we are solving Subscript[\[CapitalGamma], evap]=Subscript[\[CapitalGamma], cap] <=> Subscript[N, C]= Subscript[\[CapitalGamma], cap]/(Subscript[\[CapitalGamma], evap]/Subscript[N, C])= Subscript[\[CapitalGamma], cap]/(\[CapitalGamma]evap,per-part)*)

(*\[CapitalGamma]pC = Total[10^-#&/@({"\[Tau]CC","\[Tau]CA"}/.log\[Tau]sdict)];
\[CapitalGamma]pA = Total[10^-#&/@({"\[Tau]AA","\[Tau]CA"}/.log\[Tau]sdict)];*)

(*NE = Capture`naDM[fD,"mpD"]((4\[Pi])/3 ("rE")^3)/.pdict[[mind,\[Kappa]ind]]/.Constants`EarthRepl;
NA = Capture`naDM[fD,"mpD"]((4\[Pi])/3 rmax^3)/.pdict[[mind,\[Kappa]ind]];*)

NCfrom\[Delta] = Ncapof\[Delta][getncap["mpD","meD","\[Beta]D",Capture`naDM[fD,"mpD"],\[Alpha]Dby\[Alpha]]/.pdict[[mind,\[Kappa]ind]],\[Delta]]; (*NC computed from the given value of \[Delta]*)


(*Capture rate for earth is total rate accounting for flux incident on Earth. Evap is per particle so needs to account for number
of captured charges in Earth's volume at a time (ie. a ratio of rmax volume to earth volume)*)
\[CapitalGamma]evapE = (("rE")/rmax)^3 "\[CapitalGamma]total"/.pdict[[mind,\[Kappa]ind]]/.Constants`EarthRepl;
\[CapitalGamma]capE=(*(("rE")/rmax)^3*)"dNcMaxdt"/.pdict[[mind,\[Kappa]ind]]/."nD"->Capture`naDM[fD,"mpD"]/.pdict[[mind,\[Kappa]ind]]/.Constants`EarthRepl; (*per particle rate*)

(*\[CapitalGamma]capP = ((4\[Pi])/3 rmax^3) Sum["rate"/.#[[i]],{i,Length[#]}]&@Gather[ratetable,#1["keys"][[1]]==#2["keys"][[1]]&][[1]];*)
\[CapitalGamma]capP = ((4\[Pi])/3 rmax^3) Sum[If[#["keys"][[1]]=="cap","rate"/.#[[i]],0],{i,Length[#]}]&@ratetable;

(*\[CapitalGamma]evapP = ((4\[Pi])/3 rmax^3)1/NCfrom\[Delta]Sum["rate"/.#[[i]],{i,Length[#]}]&@Gather[ratetable,#1["keys"][[1]]==#2["keys"][[1]]&][[2]]; (*per particle rate*)*)
\[CapitalGamma]evapP = ((4\[Pi])/3 rmax^3) 1/NCfrom\[Delta] Sum[If[#["keys"][[1]]=="evap","rate"/.#[[i]],0],{i,Length[#]}]&@ratetable; (*per particle rate*)

(*Nctot = (\[CapitalGamma]capE-\[CapitalGamma]evapE+ NA \[CapitalGamma]pA)/\[CapitalGamma]pC;*)

\[CapitalGamma]captot = \[CapitalGamma]capE+ \[CapitalGamma]capP;(*total cap rate*)
\[CapitalGamma]evaptot = \[CapitalGamma]evapE+ \[CapitalGamma]evapP;(*total per particle evap rate*)
(*ATTN: we should really do this more carefully - this is assuming a particular order for the gather which is precarious. *)
(*ATTN: double check that factors of \[Kappa] etc are accounted for in Earth rates*)
(*ah no there are volume factors for the Earth rates that need to be included.*)

(*Print[\[CapitalGamma]evapE];
Print[\[CapitalGamma]evapP];
Print[\[CapitalGamma]capE];
Print[\[CapitalGamma]capP];*)

(*am I missing volume factors for the plasma rates?*)

(*Return[Print["nah b"],Module];*)


tE = 4.543 10^9 (356  24 3600);(*[s] age of the Earth*)

(*these rates may need to be modified slightly to account for the varying rates with r.*)
(*\[CapitalGamma]evaptot = \[CapitalGamma]pC + \[CapitalGamma]evapE/Nctot;
\[CapitalGamma]captot = NA \[CapitalGamma]pA + \[CapitalGamma]capE;*)

timeenoughforequil = 1/tE<(\[CapitalGamma]evaptot);

Nctot = If[timeenoughforequil, \[CapitalGamma]captot/\[CapitalGamma]evaptot, \[CapitalGamma]captot tE];

(*Print[Nctot];*)

(*timeenoughforequilE = 1/tE<\[CapitalGamma]evapE/Nctot;

NcE = If[timeenoughforequilE,Nctot\[CapitalGamma]capE/\[CapitalGamma]evapE,\[CapitalGamma]capE tE];*)

(*Print[\[CapitalGamma]capE/\[CapitalGamma]evapE];*)


(*Print[timeenoughforequil];*)

newstuffassoc = <|"Nctot"->Nctot,(*"NcE"->NcE,*)"\[CapitalGamma]captot"->\[CapitalGamma]captot,"\[CapitalGamma]evaptot"->\[CapitalGamma]evaptot,"\[CapitalGamma]evapE"-> \[CapitalGamma]evapE,"\[CapitalGamma]capE"->\[CapitalGamma]capE,"\[CapitalGamma]evapP"-> \[CapitalGamma]evapP,"\[CapitalGamma]capP"->\[CapitalGamma]capP,"timeenoughforequil"->timeenoughforequil,"ratetable"->ratetable|>;

(*Print["made it herer"];*)

(*AppendTo[pdictnew[[-1]],pdict[[mind,\[Kappa]ind]]];*)
(*pdictnew = ReplacePart[pdictnew, -1 -> Append[pdictnew[[-1]], pdict[[mind,\[Kappa]ind]]]];*)

(*Print[Dimensions[pdictnew]];*)
(*AssociateTo[pdictnew[[-1,-1]],{"Nctot"->Nctot,(*"NcE"->NcE,*)"\[CapitalGamma]captot"->\[CapitalGamma]captot,"\[CapitalGamma]evaptot"->\[CapitalGamma]evaptot,"\[CapitalGamma]evapE"-> \[CapitalGamma]evapE,"\[CapitalGamma]capE"->\[CapitalGamma]capE,"\[CapitalGamma]evapP"-> \[CapitalGamma]evapP,"\[CapitalGamma]capP"->\[CapitalGamma]capP,"timeenoughforequil"->timeenoughforequil,"ratetable"->ratetable(*,"timeenoughforequilE"->timeenoughforequilE*)}];*)
(*AssociateTo[pdictnew[[-1,-1]],<|"Nctot"->Nctot,(*"NcE"->NcE,*)"\[CapitalGamma]captot"->\[CapitalGamma]captot,"\[CapitalGamma]evaptot"->\[CapitalGamma]evaptot,"\[CapitalGamma]evapE"-> \[CapitalGamma]evapE,"\[CapitalGamma]capE"->\[CapitalGamma]capE,"\[CapitalGamma]evapP"-> \[CapitalGamma]evapP,"\[CapitalGamma]capP"->\[CapitalGamma]capP,"timeenoughforequil"->timeenoughforequil,"ratetable"->ratetable(*,"timeenoughforequilE"->timeenoughforequilE*)|>];*)
(*pdictnew = ReplacePart[pdictnew, {-1, -1} -> Append[pdictnew[[-1, -1]], newstuffassoc]];*)
pdictnew = ReplacePart[pdictnew, {-1} -> Append[pdictnew[[-1]],Append[pdict[[mind,\[Kappa]ind]], newstuffassoc]]];

(*Print[Dimensions[pdictnew]];*)

(*Print["also made it here"];*)

,{\[Kappa]ind,(*7,7*)Dimensions[pdict][[2]]}]
,{mind,(*4,4*)Dimensions[pdict][[1]]}],1];
,{\[Kappa]ind,mind}];

pdictnew
]


(* ::Chapter::Closed:: *)
(*End*)


End[];


EndPackage[];
