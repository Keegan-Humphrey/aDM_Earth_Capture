(* ::Package:: *)

Needs["Constants`",NotebookDirectory[]<>"../packages/Constants.wl"]
Needs["Utilities`",NotebookDirectory[]<>"../packages/Utilities.wl"]
(*Needs["Dielectrics`",NotebookDirectory[]<>"../packages/Dielectrics.wl"]*)
(*Needs["FormFactors`",NotebookDirectory[]<>"../packages/FormFactors.wl"]*)
(*Needs["EnergyLoss`",NotebookDirectory[]<>"../packages/EnergyLoss.wl"]*)
Needs["Capture`",NotebookDirectory[]<>"../packages/Capture.wl"]
Needs["Screening`",NotebookDirectory[]<>"../packages/Screening.wl"]
Needs["CollisionRate`",NotebookDirectory[]<>"../packages/CollisionRate.wl"]


BeginPackage["ThermTimescales`"];


(* ::Text:: *)
(*The purpose of this package is to compute timescales for plasma interactions and thermalization among the various ambient and captured populations. *)


(* ::Subsubsection:: *)
(*Public Declarations*)


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


(* ::Section:: *)
(*Private*)


Begin["`Private`"];


(* ::Subsection:: *)
(*Define Base Functions*)


(* ::Text:: *)
(*Should change the inputs to relative and target speeds*)


(* ::Subsubsection::Closed:: *)
(*Module for Captured Relaxation Times*)


(* ::Input::Initialization:: *)
Clear[Get\[Tau]C]
Get\[Tau]C[Logm_,v0kms_?NumericQ,\[Delta]t_:10^-4]:=Module[{ncap,ncapofr,Ntot,rmax,rat,ravg,\[Tau]relCintegrand,\[Tau]relC,meDbympD=10^-4,fD=0.05,\[Alpha]Dby\[Alpha]=1,DEBUG=False},
ncap = Screening`getncap["mpD","meD","\[Beta]D","nD",\[Alpha]Dby\[Alpha]]/."nD"->Capture`naDM[fD,"mpD"+"meD"]/."\[Beta]D"->(*3/("mpD"("v0")^2)*)3 /("mpD" ("v0")^2)/.{Private`v->"v0"}/."v0"->v0kms 10^3/."meD"->"mpD"meDbympD/."mpD"->("JpereV")/("c")^2 10^Logm/.Constants`SIConstRepl/.Constants`EarthRepl;

(*Print[ncap];*)

rmax=Screening`getrmax[ncap][\[Delta]t];

ncapofr=Evaluate[ncap/.{"\[Delta]"->\[Delta]t,Private`r->#}]&;



Ntot = Screening`Ncapof\[Delta][ncap,\[Delta]t];

ravg =Re[ (4 \[Pi])/Ntot NIntegrate[r^3 ncap/."\[Delta]"->\[Delta]t/.Private`r->r,{r,0,rmax}]];

\[Tau]relCintegrand=(v /(("m\[Chi]")/2 v^2) CollisionRate`dEkdrcoul )^-1/."nF"->ncap/.{"m\[Chi]"->"mpD","mT"->"mpD"}/.{"\[Alpha]D"->"\[Alpha]","meD"->"mpD"meDbympD}/.{"nD"->Capture`naDM[fD,"mpD"]}/."\[Beta]D"->3/("mpD" ("v0")^2)/.Private`v->v/.v->"vesc" \[Delta]t /."v0"->v0kms 10^3/."mpD"->("JpereV")/("c")^2 10^Logm/.\[Delta]->\[Delta]t/.Constants`SIConstRepl/.Constants`EarthRepl/.{Private`r->#,"\[Delta]"->\[Delta]t}&;
(*Print@Plot[Log10@\[Tau]relCintegrand[10^logr],{logr,0,Log10@rmax},PlotRange->{0,2}];*)
\[Tau]relC  = 1/rmax NIntegrate[\[Tau]relCintegrand[r],{r,0,0.95 rmax}];(*averaged relaxation time through an orbit assuming right through Earth, so this is really an lower bound since some orbits will miss earth*)

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
Get\[Tau]A[Logm_,v0kms_?NumericQ,\[Delta]t_:10^-4,vrat_:10^-4]:=Module[{ncap,Ntot,rmax,rat,ravg,\[Tau]relCintegrand,\[Tau]relA,meDbympD=10^-4,fD=0.05,\[Alpha]Dby\[Alpha]=1},

ncap = Screening`getncap["mpD","meD","\[Beta]D","nD",\[Alpha]Dby\[Alpha]]/."nD"->Capture`naDM[fD,"mpD"]/."\[Beta]D"->3/("mpD" ("v0")^2)/.{Private`v->"v0",v->"v0"}/."v0"->v0kms 10^3/."meD"->"mpD"meDbympD/."mpD"->("JpereV")/("c")^2 10^Logm/.Constants`SIConstRepl/.Constants`EarthRepl;

rmax=Screening`getrmax[ncap][\[Delta]t];

Ntot = Screening`Ncapof\[Delta][ncap,\[Delta]t];

ravg =Re[ (4 \[Pi])/Ntot NIntegrate[r^3 ncap/."\[Delta]"->\[Delta]t/.Private`r->r,{r,0,rmax}]];

\[Tau]relA=(v /(("m\[Chi]")/2 v^2) CollisionRate`dEkdrcoul )^-1/.{"m\[Chi]"->"mpD","mT"->"meD"}/.{"\[Alpha]D"->"\[Alpha]"\[Alpha]Dby\[Alpha],"meD"->"mpD"meDbympD}/."nF"->Capture`naDM[fD,"mpD"]/."mpD"->("JpereV")/("c")^2 10^Logm/.Private`v->v/.v->"vesc"vrat(* \[Delta]t *)/."v0"->v0kms 10^3/.Constants`SIConstRepl/.Constants`EarthRepl;

(*\[Tau]relC  = 1/rmax NIntegrate[\[Tau]relCintegrand[r],{r,0,0.99 rmax}];(*averaged relaxation time through an orbit*)*)

(*Print@Plot[Log10@\[Tau]relCintegrand[10^logr],{logr,0,Log10@rmax},PlotRange->{0,2}];*)

<|"rmax"->rmax,"ravg"->ravg,"Ntot"->Ntot,"\[Tau]relA"->\[Tau]relA|>
]


(* ::Subsection:: *)
(*Timescale Modules*)


(* ::Subsubsection:: *)
(*\[Tau]AA*)


Clear[Get\[Tau]AA]
Get\[Tau]AA[logm_,log\[Kappa]_,pdict_,\[Delta]in_]:=Module[{flatpdict,pdictind,\[Tau]relAdict,\[Tau]relAA,\[Tau]nearEarth},

flatpdict = Flatten[pdict,1];

pdictind=FirstPosition[{Round@Log10@("mpD"//N),Round@Log10@("\[Kappa]"//N)}/.flatpdict,{Round@Log10@(10^logm ("JpereV")/("c")^2/.Constants`SIConstRepl),log\[Kappa]}][[1]];
(*\[Delta]temp="\[Delta]"/.pdict[[\[Delta]dictind]];*)

\[Tau]relAdict=Get\[Tau]A[logm,10^-3 "v0"/.flatpdict[[pdictind]],\[Delta]in,("v0")/("vesc")/.flatpdict[[pdictind]]/.Constants`EarthRepl];
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

\[Tau]relAdict=Quiet@Get\[Tau]C[logm,10^-3 "v0"/.flatpdict[[pdictind]],\[Delta]in];

\[Tau]relAE=Quiet[(v\[Chi](("\[Kappa]")^2/.flatpdict[[pdictind]])Quiet@Max[10^-100,("dEdlofr"/.Capture`GetTotalELFunction[\[Sigma]dicte[[eind]],\[Sigma]dictNuc[[Nucind]]])[1 10^6,v\[Chi]]]/((10^logm ("JpereV")/("c")^2/.Constants`SIConstRepl) v\[Chi]^2/2))^-1/.v\[Chi]->"v0"/.flatpdict[[pdictind]]/.Constants`EarthRepl(*/.\[Tau]relAdict*)];
(*\[Tau]relAE=1;*)
\[Tau]EarthPass = ("rE")/("v0")/.flatpdict[[pdictind]]/.Constants`EarthRepl;

<|"logm"->logm,"log\[Kappa]"->log\[Kappa],"logassrat"->Log10[\[Tau]relAE/\[Tau]EarthPass],"log\[Tau]rel"->Log10@\[Tau]relAE,"log\[Tau]comp"->Log10@\[Tau]EarthPass|>
]


(* ::Subsubsection:: *)
(*\[Tau]CA*)


Clear[Get\[Tau]CA]
Get\[Tau]CA[logm_,log\[Kappa]_,pdict_,\[Delta]in_]:=Module[{flatpdict,pdictind,\[Delta]temp,\[Tau]relCAdict,\[Tau]relCA,\[Tau]evap},

flatpdict = Flatten[pdict,1];

pdictind=FirstPosition[{Round@Log10@("mpD"//N),Round@Log10@("\[Kappa]"//N)}/.flatpdict,{Round@Log10@(10^logm ("JpereV")/("c")^2/.Constants`SIConstRepl),log\[Kappa]}][[1]];
(*\[Delta]temp="\[Delta]"/.pdict[[pdictind]];*)

(*\[Tau]relCdict=Get\[Tau]C[logm,10^-3"v0"/.pdict[[pdictind]],\[Delta]temp];
\[Tau]relCC = "\[Tau]relC"/.\[Tau]relCdict;*)

\[Tau]relCAdict=Get\[Tau]A[logm,10^-3 "v0"/.flatpdict[[pdictind]],\[Delta]in,("v0")/("vesc")/.flatpdict[[pdictind]]/.Constants`EarthRepl];
\[Tau]relCA = "\[Tau]relA"/.\[Tau]relCAdict;

(*\[Tau]evap =("\[Kappa]")^-2 ("\[CapitalGamma]evap")^-1/.flatpdict[[pdictind]];*)
\[Tau]evap = (("\[Kappa]")^-2 /.#)If[StringQ@("\[CapitalGamma]evap"/.#),"\[CapitalGamma]total"/.#,"\[CapitalGamma]evap"/.#]^-1&@flatpdict[[pdictind]];

<|"logm"->logm,"log\[Kappa]"->log\[Kappa],"logassrat"->Log10[\[Tau]relCA/\[Tau]evap],"log\[Tau]rel"->Log10@\[Tau]relCA,"log\[Tau]comp"->Log10@\[Tau]evap|>
]


(* ::Subsection:: *)
(*Pipeline Modules*)


(* ::Subsubsection:: *)
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
Get\[Tau]sFromPDict[pdict_,\[Delta]_,\[Sigma]dicte_,\[Sigma]dictNuc_,\[Alpha]Dby\[Alpha]_:1,fD_:0.05]:= Module[{pdictnew,logm,log\[Kappa],\[Tau]stable,\[Tau]AA,\[Tau]CC,\[Tau]CA,\[Tau]CE,\[Tau]AE,log\[Tau]sdict,rmax,\[CapitalGamma]pC,\[CapitalGamma]pA,NE,NA,\[CapitalGamma]evapE,\[CapitalGamma]capE,\[CapitalGamma]evaptot,\[CapitalGamma]captot,tE,timeenoughforequil,Nctot},

pdictnew = {};

\[Tau]stable=Flatten[Table[

AppendTo[pdictnew,{}];

Table[

{logm ,log\[Kappa]}=Log10@{(("m\[Chi]" ("c")^2)/("JpereV")/.pdict[[mind,\[Kappa]ind]]/.Constants`SIConstRepl),("\[Kappa]"/.pdict[[mind,\[Kappa]ind]])};

Print[{logm ,log\[Kappa]}];

(*can probably do the index search globally and pass it as needed.*)
(*need to update these with the relative velocity rather than the velocity of the probe.*)
\[Tau]AA = Get\[Tau]AA[logm,log\[Kappa],pdict,\[Delta]];
\[Tau]CC = Get\[Tau]CC[logm,log\[Kappa],pdict,\[Delta]];
\[Tau]CA = Get\[Tau]CA[logm,log\[Kappa],pdict,\[Delta]];
\[Tau]CE = Get\[Tau]CE[logm,log\[Kappa],pdict,\[Delta],\[Sigma]dicte,\[Sigma]dictNuc];
\[Tau]AE = Get\[Tau]AE[logm,log\[Kappa],pdict,\[Delta],\[Sigma]dicte,\[Sigma]dictNuc];

(*add these to the pdict as an association*)

log\[Tau]sdict = <|"\[Tau]AA"->"log\[Tau]rel"/.\[Tau]AA,"\[Tau]CC"->"log\[Tau]rel"/.\[Tau]CC,"\[Tau]CA"->"log\[Tau]rel"/.\[Tau]CA,"\[Tau]CE"->"log\[Tau]rel"/.\[Tau]CE,"\[Tau]AE"->"log\[Tau]rel"/.\[Tau]AE,"\[Tau]Epass"->"log\[Tau]comp"/.\[Tau]AE,"\[Tau]evapE"->"log\[Tau]comp"/.\[Tau]CE|>;

Print[log\[Tau]sdict];

rmax = Screening`getrmax[Screening`getncap["mpD","meD","\[Beta]D",Capture`naDM[fD,"mpD"],\[Alpha]Dby\[Alpha]]/.pdict[[mind,\[Kappa]ind]]][\[Delta]];

(*These should really be rates that are integrated over r*)
\[CapitalGamma]pC = Total[10^-#&/@({"\[Tau]CC","\[Tau]CA"}/.log\[Tau]sdict)];
\[CapitalGamma]pA = Total[10^-#&/@({"\[Tau]AA","\[Tau]CA"}/.log\[Tau]sdict)];

NE = Capture`naDM[fD,"mpD"]((4\[Pi])/3 ("rE")^3)/.pdict[[mind,\[Kappa]ind]]/.Constants`EarthRepl;
NA = Capture`naDM[fD,"mpD"]((4\[Pi])/3 rmax^3)/.pdict[[mind,\[Kappa]ind]];
\[CapitalGamma]evapE = NE/\!\(\*SuperscriptBox[\(10\), \("\<\[Tau]evapE\>"\)]\)/.log\[Tau]sdict;

\[CapitalGamma]capE="dNcMaxdt"/.pdict[[mind,\[Kappa]ind]]/."nD"->Capture`naDM[fD,"mpD"]/.pdict[[mind,\[Kappa]ind]];

Nctot = (\[CapitalGamma]capE-\[CapitalGamma]evapE+ NA \[CapitalGamma]pA)/\[CapitalGamma]pC;

Print[Nctot];

tE = 4.543 10^9 (356  24 3600);(*[s] age of the Earth*)

(*these rates may need to be modified slightly to account for the varying rates with r.*)
\[CapitalGamma]evaptot = \[CapitalGamma]pC + \[CapitalGamma]evapE/Nctot;
\[CapitalGamma]captot = NA \[CapitalGamma]pA + \[CapitalGamma]capE;

timeenoughforequil = 1/tE<(\[CapitalGamma]evaptot);

Nctot = If[timeenoughforequil, Nctot, \[CapitalGamma]captot tE];

Print[timeenoughforequil];


AppendTo[pdictnew[[-1]],pdict[[mind,\[Kappa]ind]]];
AssociateTo[pdictnew[[-1,-1]],{"\[CapitalGamma]total"->\[CapitalGamma]evaptot,"log\[Tau]s"->log\[Tau]sdict,"Nc"->Nctot,"\[CapitalGamma]captot"->\[CapitalGamma]captot,"\[CapitalGamma]evaptot"->\[CapitalGamma]evaptot,"\[CapitalGamma]evapE"-> \[CapitalGamma]evapE/Nctot,"\[CapitalGamma]capE"->\[CapitalGamma]capE,"timeenoughforequil"->timeenoughforequil}];

,{\[Kappa]ind,7,7(*Dimensions[pdict][[2]]*)}]
,{mind,4,4(*Dimensions[pdict][[1]]*)}],1];

pdictnew
]


(* ::Section::Closed:: *)
(*End*)


End[];


EndPackage[];
