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


Compute\[Delta]DictonScan::usage="";


Get\[Delta]DictTruth::usage="";


ScanOverUnfixedParameters::usage="";


checkforVeffMaxima::usage="";


(* ::Subsubsection:: *)
(*Potential outside the flat region*)


get\[Phi]Dnoncap::usage = "";


get\[Phi]DAnalyticEstimates::usage = "";


get\[Phi]D::usage = "";


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
]


(* ::Subsubsection:: *)
(*get Subscript[r, max](\[Delta])*)


(*getrmax[ncof\[Delta]andr_]:=rmax/.FindRoot[Re@(ncof\[Delta]andr/."\[Delta]"->#/.r->rmax),{rmax,Evaluate[("rE")/(2 #)/.Constants`EarthRepl]},Method->"Secant"]&*)


 Clear[binsearch]
binsearch[f_,var_,domain_:{1,11},n_:20,DEBUG_:False(*,ofmmin_:1,ofmmax_:13*)(*,LogQ_:True*)]:=Module[{midlist={},interval,fposQ,fmingrfmaxQ,mid,root(*,DEBUG=True*)},
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
			N[{root,f/.var->root}]
	]


getrmax[ncof\[Delta]andr_]:=binsearch[Abs@(*Re@*)(ncof\[Delta]andr/."\[Delta]"->#),r,{Log10@Evaluate[("rE")/(2 #)/.Constants`EarthRepl]-1,Log10@Evaluate[("rE")/(2 #)/.Constants`EarthRepl]+1},30][[1]]&


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


(* ::Subsubsection:: *)
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


\[Delta]ofNcap[ncof\[Delta]andr_,NC_]:= Module[{\[Delta]dom={-7,Log10[1.5](*Log10[0.74]*)},N\[Delta]s=10,rmaxintrp,\[Delta]min,NCof\[Delta],\[Delta]minNC,NCof\[Delta]intrp},
(*speed up \[Delta] root finding byinterpolating over rmax first*)

(*Print["test1"];*)
rmaxintrp = Interpolation[Table[{\[Delta]t,Log10@getrmax[ncof\[Delta]andr][10^\[Delta]t]},{\[Delta]t,\[Delta]dom[[1]],\[Delta]dom[[2]],(\[Delta]dom[[2]]-\[Delta]dom[[1]])/N\[Delta]s}],InterpolationOrder->1];

(*Print["test2"];
Print[Plot[rmaxintrp[\[Delta]t],{\[Delta]t,\[Delta]dom[[1]],\[Delta]dom[[2]]}]];
*)
(*\[Delta]min = binsearch[10^rmaxintrp[\[Delta]t],\[Delta]t,{-1,\[Delta]dom[[2]]},30];

Print[\[Delta]min];*)

NCof\[Delta] = 4 \[Pi] NIntegrate[r^2 ncof\[Delta]andr/."\[Delta]"->#,{r,0,Evaluate[10^rmaxintrp[Log10[#]]]}]&;

(*Print["test3"];
*)
(*Print@Table[Log10@NCof\[Delta][10^\[Delta]t],{\[Delta]t,\[Delta]dom[[1]],\[Delta]dom[[2]]}];

Print[((NCof\[Delta][#]-NC)&[\[Delta]t]/.\[Delta]t->10^(\[Delta]dom[[2]]))];

(*Print[ListLinePlot[Table[Log10@NCof\[Delta][\[Delta]t],{\[Delta]t,\[Delta]dom[[1]],\[Delta]dom[[2]]}]]];*)

Print[NC];*)

NCof\[Delta]intrp = Interpolation[Table[{\[Delta]t,NCof\[Delta][10^\[Delta]t]},{\[Delta]t,\[Delta]dom[[1]],\[Delta]dom[[2]],Evaluate[(\[Delta]dom[[2]]-\[Delta]dom[[1]])/N\[Delta]s]}],InterpolationOrder->1];
(*Print["test4"];
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


get\[Delta]fromPDict[pdict_,\[Alpha]Dby\[Alpha]_:1,fD_:0.05]:=Module[{meD,mpD,\[Beta]D,v0,\[Kappa],nD,nCap,NC,\[CapitalGamma]cap,\[CapitalGamma]evap,Ncapof\[Delta]table,\[Delta]t,rmax},

{meD, mpD, \[Beta]D, v0, \[Kappa]}={"meD","mpD","\[Beta]D","v0","\[Kappa]"}/.pdict;

nD =Capture`naDM[fD,meD+mpD];
nCap=getncap[mpD, meD, \[Beta]D, nD, \[Alpha]Dby\[Alpha]];

NC = "Nc"/.pdict/."nD"->nD;
\[CapitalGamma]cap = "dNcMaxdt"/.pdict/."nD"->nD; (*capture rate averaged over Earth*)
\[CapitalGamma]evap = "\[CapitalGamma]total"/.pdict;(*per particle evaporation rate averaged over Earth*)

\[Delta]t=\[Delta]ofNcap[nCap,NC];
rmax=getrmax[nCap][\[Delta]t];
(*\[Delta]t={Print@#[[1]],#[[2]]}[[2]]&@AbsoluteTiming[\[Delta]ofNcap[nCap,NC]];
rmax={Print@#[[1]],#[[2]]}[[2]]&@AbsoluteTiming[getrmax[nCap][\[Delta]t]];*)

<|"\[Delta]"->\[Delta]t,"rmax"->rmax,"NC"->NC,"\[CapitalGamma]cap"->\[CapitalGamma]cap,"\[CapitalGamma]evap"->\[CapitalGamma]evap,"meD"->meD,"mpD"->mpD,"\[Beta]D"->\[Beta]D,"nD"->nD,"v0"->v0,"\[Kappa]"->\[Kappa],"fD"->fD,"\[Alpha]Dby\[Alpha]"->\[Alpha]Dby\[Alpha]|>
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


Clear[Compute\[Delta]DictonScan]
Compute\[Delta]DictonScan[\[Delta]List_,v0ind_:1,\[Alpha]Dby\[Alpha]_:1,fD_:0.05]:= Module[{return\[Delta]List,\[Delta]Listcurrent,PDict,v0,mratio,\[Delta]table},
(*
\[Delta]List - of the form: {<|"\[Delta]"->\[Delta],"file"->"\\path\\to\\file"|>,...}
*)
return\[Delta]List ={};

Monitor[
Do[

\[Delta]Listcurrent = \[Delta]List[[\[Delta]]];

(*PDict= Utilities`ReadIt[\[Delta]Listcurrent["file"]][[;;9,v0ind,;;]];*)
PDict= Utilities`ReadIt[\[Delta]Listcurrent["file"]][[;;,v0ind,;;]];

{v0, mratio} = {"v0","meD"/"mpD"}/.First@First[PDict];

(*\[Delta]table=Flatten[Table[Table[{Log10[("mpD" ("c")^2)/("JpereV")]/.PDict[[i,j]]/.Constants`SIConstRepl,Log10["\[Kappa]"]/.PDict[[i,j]],Log10@"\[Delta]"/.get\[Delta]fromPDict[PDict[[i,j]]]},{i,Dimensions[PDict][[1]]}],{j,Dimensions[PDict][[2]]}],1];*)
(*\[Delta]table=Flatten[Table[Table[get\[Delta]fromPDict[PDict[[i,j]]],{i,Dimensions[PDict][[1]]}],{j,Dimensions[PDict][[2]]}],\[Alpha]Dby\[Alpha],fD];*)
\[Delta]table=Flatten[Table[Table[get\[Delta]fromPDict[PDict[[i,j]],\[Alpha]Dby\[Alpha],fD],{i,Dimensions[PDict][[1]]}],{j,Dimensions[PDict][[2]]}],\[Alpha]Dby\[Alpha],fD];

AppendTo[\[Delta]Listcurrent,<|"\[Delta]table"->\[Delta]table,"v0"->v0,"meDbympD"->mratio|>];
AppendTo[return\[Delta]List,\[Delta]Listcurrent];

,{\[Delta],Length[\[Delta]List]}];
,{\[Delta],j,i}];

Print["Finished running the PDicts for: \n{v0,meD/mpD}=",First@First[{"v0","meD"/"mpD"}/.PDict]];

return\[Delta]List
]


(* ::Subsubsection:: *)
(*Scan over unfixed parameters*)


Clear[ScanOverUnfixedParameters]
ScanOverUnfixedParameters[\[Delta]ScanList_,\[Alpha]Dby\[Alpha]s_:{1,2},fDs_:{0.01,0.05}]:=Module[{\[Delta]Dict},
Do[
Print["Running: ",{v0ind,\[Alpha]Dby\[Alpha],fD}];
\[Delta]Dict=Quiet@Compute\[Delta]DictonScan[\[Delta]ScanList,v0ind,\[Alpha]Dby\[Alpha],fD];
(*Print[\[Delta]Dict];*)
(*SaveIt[NotebookDirectory[]<>StringForm["\[Delta]Dict_``_aDbya_``_fD_``",{v0ind,\[Alpha]Dby\[Alpha],fD}],\[Delta]Dict];*)
Capture`ExportDatFileToDir[\[Delta]Dict,ToString@StringForm["\[Delta]Dict_``_aDbya_``_fD_``",v0ind,\[Alpha]Dby\[Alpha],fD],"\[Delta]Dicts"]
,{v0ind,4},{\[Alpha]Dby\[Alpha],\[Alpha]Dby\[Alpha]s},{fD,fDs}
]
]


(* ::Subsubsection:: *)
(*Get \[Delta] Dict Truth*)


(*Clear[Get\[Delta]DictTruth]
Get\[Delta]DictTruth[\[Delta]Scan_]:=Module[{DEBUG=False,\[Delta]truthtable,\[Delta]outof\[Delta]intable,\[Delta]outof\[Delta]in,\[Delta]truth,\[Delta]Scantruth},

(*take the \[Delta]Scans, and replace \[Delta] by the truth value*)

\[Delta]truthtable ={};

Do[
\[Delta]outof\[Delta]intable =Table[{Log10@"\[Delta]"/.\[Delta]Scan[[j]],Log10@"\[Delta]"/.\[Delta]Scan[[j]]["\[Delta]table"][[i]]},{j,Length[\[Delta]Scan]}];
\[Delta]outof\[Delta]in = Interpolation[\[Delta]outof\[Delta]intable,InterpolationOrder->1];
\[Delta]truth = \[Delta]/.FindRoot[10^\[Delta]outof\[Delta]in[Log10@\[Delta]]-\[Delta],{\[Delta],0.01}];

\[Delta]Scantruth = \[Delta]Scan[[1]]["\[Delta]table"][[i]];
\[Delta]Scantruth["\[Delta]"]=\[Delta]truth;

{\[Delta]Scantruth["rmax"],\[Delta]Scantruth["NC"]}={getrmax[#]["\[Delta]"],Ncapof\[Delta][#,"\[Delta]"]}&@getncap["mpD","meD","\[Beta]D","nD","\[Alpha]Dby\[Alpha]"]/.\[Delta]Scantruth//Quiet;

AppendTo[\[Delta]Scantruth,get\[Phi]DAnalyticEstimates["mpD","meD","\[Beta]D","nD","\[Alpha]Dby\[Alpha]"]/.\[Delta]Scantruth];

If[DEBUG&&(#[[1]]==13&&#[[2]]==-13)&@(Log10@{"mpD" ("c")^2/("JpereV"),"\[Kappa]"}/.\[Delta]Scan[[1]]["\[Delta]table"][[i]]/.Constants`SIConstRepl//N//Round),Print[i];Print[\[Delta]truth];Print[\[Delta]Scantruth]];

AppendTo[\[Delta]truthtable,\[Delta]Scantruth];

,{i,Length[\[Delta]Scan[[1]]["\[Delta]table"]]}];(*assuming all \[Delta]tables have the same ordering size*)

\[Delta]truthtable
]*)


Clear[Get\[Delta]DictTruth]
Get\[Delta]DictTruth[\[Delta]Scan_]:=Module[{DEBUG=False,\[Delta]truthtable,\[Delta]outof\[Delta]intable,\[Delta]outof\[Delta]in,\[Delta]truth,\[Delta]Scantruth},

(*take the \[Delta]Scans, and replace \[Delta] by the truth value*)

\[Delta]truthtable ={};

Monitor[Do[
(*Print[i];*)
\[Delta]outof\[Delta]intable =Table[{Log10@"\[Delta]"/.\[Delta]Scan[[j]],Log10@"\[Delta]"/.Flatten[\[Delta]Scan[[j]]["\[Delta]table"]][[i]]},{j,Length[\[Delta]Scan]}];
(*Print[\[Delta]outof\[Delta]intable];*)
\[Delta]outof\[Delta]in = Interpolation[\[Delta]outof\[Delta]intable,InterpolationOrder->1];
\[Delta]truth = \[Delta]/.FindRoot[10^\[Delta]outof\[Delta]in[Log10@\[Delta]]-\[Delta],{\[Delta],0.01}];

\[Delta]Scantruth = Flatten[\[Delta]Scan[[1]]["\[Delta]table"]][[i]];
\[Delta]Scantruth["\[Delta]"]=\[Delta]truth;

(*Print[\[Delta]truth];
Print[Flatten[\[Delta]Scan[[1]]["\[Delta]table"]][[i]]];
Print["Before"];*)

(*{\[Delta]Scantruth["rmax"],\[Delta]Scantruth["NC"]}={getrmax[#]["\[Delta]"],Ncapof\[Delta][#,"\[Delta]"]}&@getncap["mpD","meD","\[Beta]D","nD","\[Alpha]Dby\[Alpha]"]/.\[Delta]Scantruth(*//Quiet*);*)
{\[Delta]Scantruth["rmax"],\[Delta]Scantruth["NC"]}={getrmax[#][\[Delta]truth],Ncapof\[Delta][#,\[Delta]truth]}&@(getncap["mpD","meD","\[Beta]D","nD","\[Alpha]Dby\[Alpha]"]/.\[Delta]Scantruth(*//Quiet*));

(*Print["After"];*)

AppendTo[\[Delta]Scantruth,get\[Phi]DAnalyticEstimates["mpD","meD","\[Beta]D","nD","\[Alpha]Dby\[Alpha]"]/.\[Delta]Scantruth];

If[DEBUG,If[(#[[1]]==13&&#[[2]]==-13)&@(Log10@{"mpD" ("c")^2/("JpereV"),"\[Kappa]"}/.\[Delta]Scan[[1]]["\[Delta]table"][[i]]/.Constants`SIConstRepl//N//Round),Print[i];Print[\[Delta]truth];Print[\[Delta]Scantruth]]];

AppendTo[\[Delta]truthtable,\[Delta]Scantruth];

,{i,Length[Flatten[\[Delta]Scan[[1]]["\[Delta]table"]]]}],{i}];(*assuming all \[Delta]tables have the same ordering size*)

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
(*get analytic estimate for \[Phi]D in the small charge regime*)


get\[Phi]DAnalyticEstimates[mpD_,meD_,\[Beta]D_,nF_,\[Alpha]Dby\[Alpha]_,\[Phi]g_:Capture`Vgrav] := Module[{niofr,npD,neD,\[Rho]S,\[Lambda]D,nC},
niofr = 1/(mi^2 Sqrt[2 \[Pi]] \[Beta]D^2) nF (mi \[Beta]D)^(3/2) (2 mi vesc \[Beta]D+E^(1/2 mi vesc^2 \[Beta]D) Sqrt[2 \[Pi]] Sqrt[mi \[Beta]D]-E^(1/2 mi vesc^2 \[Beta]D) Sqrt[mi] Sqrt[2 \[Pi]] Sqrt[\[Beta]D] Erf[(Sqrt[mi] vesc Sqrt[\[Beta]D])/Sqrt[2]]);
npD=niofr/.vesc->Sqrt[-((2 Vgp)/mi)]/.mi->mpD;
neD=niofr/.vesc->Sqrt[-((2 Vge)/mi)]/.mi->meD;
\[Rho]S = ("e")/("\[Epsilon]0") Sqrt[\[Alpha]Dby\[Alpha]](npD-neD)/.{Vgp->mpD \[Phi]g[r],Vge->meD \[Phi]g[r]};
(*\[Lambda]D = (("e")^2/("\[Epsilon]0")\[Alpha]Dby\[Alpha](D[npD,Vgp]+D[neD,Vge]))^(-1/2)/.{Vgp->mpD \[Phi]g[r],Vge->mpD \[Phi]g[r]};*)
\[Lambda]D =(-(("e")^2/("\[Epsilon]0"))\[Alpha]Dby\[Alpha](D[npD,Vgp]+D[neD,Vge])/.{Vgp->mpD \[Phi]g[r],Vge->meD \[Phi]g[r]})^(-1/2);
nC = (mpD "G" "ME" )/(("e")^2/("\[Epsilon]0") \[Alpha]Dby\[Alpha]((4 \[Pi])/3 ("rE")^3)) Boole["rE" - r > 0] - (npD-neD)/.Constants`SIConstRepl/.Constants`EarthRepl;

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


(* ::Section::Closed:: *)
(*End*)


End[];


EndPackage[];
