(* ::Package:: *)

BeginPackage["Constants`"];


(* ::Text:: *)
(*The purpose of this package is to produce replacement tables for commonly used variables in the code. The format is for string variables in any expression*)
(**)
(*input:	 ("e")^2/(4 \[Pi] "\[Epsilon]0""hbar" "c")/.SIConstRepl*)
(*output: 	137^-1*)


(* ::Subsubsection:: *)
(*Public Declarations*)


(* ::Text:: *)
(*Public Functions*)


SIParams::usage = "Computes a replacement table of all relevant constants and parameters from a given electron number density and temperature";
neSIIron::usage = "Computes the number density of conduction electrons in iron from a given mass density";

nIFromDensity::usage = "Computes the number density of Ions from a given mass density";


(* ::Text:: *)
(*Public Replacement Tables*)


SIcoreparams::usage = "Iron parameters in the core (8 conduction electrons)";
SIcrustparams::usage = "Iron parameters in the crust (8 conduction electrons)";

SIConstRepl::usage = "Replacement table for SI constants";

MassDensities::usage = "mass densities in the mantle [kg \!\(\*SuperscriptBox[\(m\), \(-3\)]\)]";


(* ::Text:: *)
(*Misc Constants*)


vescape::usage = "escape velocity at the surface of the Earth";
vescgal::usage = "galactic escape velocity";

REarth::usage = "radius of the Earth";
rcore::usage = "radius of Earth's core";

EarthRepl::usage = "Replacement Table for constants specific to the Earth";

\[Beta]crust::usage = "\[Beta] at Earth's surface";
\[Beta]core::usage = "\[Beta] at Earth's core";


(* ::Section:: *)
(*Private*)


Begin["`Private`"];


(* ::Subsection:: *)
(*Define Functions*)


(* ::Subsubsection:: *)
(*Parameter Replacement Table - SI*)


SIParams[TK_, neSI_,nconduction_]:=Module[{lhbar,lc,lm,l\[Alpha], le,leVperK,l\[Beta],leVm,lne,lqF,lvF,l\[Mu],l\[Omega]p,l\[Chi],l\[Epsilon]0,lkB,lJpereV,lD,lnI,uperkg,lM,lEF},
(*Input:
TK - [K] temperature
neSI - [m^-3] charge number density 
nconduction - [] number of conduction electrons assumed in neSI

Usage:
	Ouputs a replacement rule in natural units for variables used in this notebook. 
*)
	lhbar = 1.055 10^-34; (*[m^2 kg s^-1] ie. [J s]*)
	lc = 3 10^8; (*[m s^-1]*)
	
	lm = 9.11 10^-31; (*[kg]*)
	uperkg = 1.66054 10^(-27); (* [u / kg] atomic mass unit per kg*)
	lM =  55.845 uperkg; (*[kg] nist iron atom mass*)
	
	l\[Alpha] = 1/137;
	(*le = Sqrt[(lhbar lc l\[Alpha])/(4 Pi)];*)
	l\[Epsilon]0 = 8.85 10^-12; (* [F m^-1]*)
	le = Sqrt[4 Pi lhbar lc l\[Alpha] l\[Epsilon]0];
	(*le = 1.602 10^-19; (*[C]*)*)
	lkB = 1.381 10^-23;(*[J K^-1]*)
	lJpereV = 1.602 10^-19;
	
	l\[Beta] = 1/(TK lkB);  (*[J^-1]*)
	
	lne = neSI; (*[m^-3]*)
	lnI = neSI/nconduction; (*[m^-3] ion number density*)
	
	lqF=(3 Pi^2 lne)^(1/3); (*[(m^-1)]*)
	lvF = (lhbar lqF)/lm;(*[m s^-1]*)
	l\[Mu] = (lhbar^2 (3 Pi^2 lne)^(2/3))/(2 lm); (*[J] E_F ~ chemical potential - valid for E_F \[Beta] >> 1 which is true throughout the earth*)
	lEF = l\[Mu];
	l\[Omega]p = Sqrt[( le^2 lne)/lm 1/ l\[Epsilon]0 ]; 
(*[(s^-1)] Plasmon freqency*)
	l\[Chi] = Sqrt[le^2/(Pi lhbar lvF ) 1/(4 \[Pi] l\[Epsilon]0)]; 
	lD = l\[Beta] l\[Mu];
	(*{hbar->lhbar,c->lc,m->lm,e-> le,\[Beta]->l\[Beta],ne->lne,qF->lqF,vF->lvF,\[Mu]->l\[Mu],\[Omega]p->l\[Omega]p,\[Alpha]->l\[Alpha],\[Chi]->l\[Chi],\[Epsilon]0->l\[Epsilon]0,JpereV->lJpereV,D->lD}*)
	{"hbar"->lhbar,"\[HBar]"->lhbar,"c"->lc,"m"->lm,"e"-> le,"\[Beta]"->l\[Beta],"ne"->lne,"qF"->lqF,"vF"->lvF,"\[Mu]"->l\[Mu],"\[Omega]p"->l\[Omega]p,"\[Alpha]"->l\[Alpha],"\[Chi]"->l\[Chi],"\[Epsilon]0"->l\[Epsilon]0,"JpereV"->lJpereV,"kB"->lkB,"D"->lD,"nI"->lnI,"Z"->nconduction,"M"->lM,"EF"->lEF}
	]


(* ::Subsubsection:: *)
(*Iron Conduction Electron Number Density - SI*)


neSIIron[massdensity_,nconduction_]:= Module[ {uperkg,M,densityatoms},

(*Input:
	massdensity - [kg m^-3]
Usage:
	Outputs number density of valence electrons in SI units [m^-3]as: number of valence electrons per atom times number of atoms per unit volume
*)
	uperkg = 1.66054 10^(-27); (* [u / kg] atomic mass unit per kg*)
	M =  55.845 uperkg; (*[kg] nist iron atom mass*)
	densityatoms =massdensity / M; (*density of atoms is the mass density of iron over the mass of an atom*)
	(*nconduction= 8;*)
	(*nconduction= 2;*)
	(*nconduction= 0.5;*)
	nconduction densityatoms 
]


nIFromDensity[massdensity_,M_]:= Module[ {uperkg,densityatoms},

(*Input:
	massdensity - [kg m^-3]
	M - [kg] atom mass
Usage:
	Outputs number density of valence electrons in SI units [m^-3]as: number of valence electrons per atom times number of atoms per unit volume
*)
	uperkg = 1.66054 10^(-27); (* [u / kg] atomic mass unit per kg*)
	densityatoms =massdensity / M; (*density of atoms is the mass density of iron over the mass of an atom*)
	(*nconduction= 8;*)
	(*nconduction= 2;*)
	(*nconduction= 0.5;*)
	densityatoms 
]


(* ::Subsection:: *)
(*Define Public Replacement Tables*)


(* ::Subsection:: *)
(*Misc Constants*)


(* ::Text:: *)
(*Escape velocity from surface of the Earth*)


vescape = 11200; (*[m s^-1]*)


(* ::Text:: *)
(*Galactic escape velocity*)


vescgal = 5.37 10^5;


(* ::Text:: *)
(*Earth and core radii*)


REarth = 6.371 10^6; (*[m]*)
rcore = REarth - 2.885 10^6; (*[m]*)


(* ::Text:: *)
(****** should put all this into a dictionary*)


T0=0;  (*[K] 0 temperature*)
(*Tcrust=290;  (*[K] ambient temperature*)*)
Tcrust=273 + (1000 + 3700)/2;(* average temperature in the Mantle https://education.nationalgeographic.org/resource/mantle/ *)
Tcore=5470;  (*[K] core temperature*)


\[Rho]crust = 7800;(*[kg m^-3] iron density at crust*)
\[Rho]core = 13000;(*[kg m^-3] iron density at earth's core*)


EarthRepl = <|"rE"->REarth,"ME"->5.97 10^24,"vesc"->vescape,"rcore"->rcore,
			"Tcrust"->Tcrust,"Tcore"->Tcore,"\[Beta]crust"->\[Beta]crust,"\[Beta]core"->\[Beta]core,
			"SiO2Frac"->0.447,"MgOFrac"->0.387|>;


MassDensities = <|"\[Rho]M"->4500,"\[Rho]SiO2"->2650,"\[Rho]MgO"->3580|> (*[kg m^-3] mantle and silicate densities*)


(* ::Subsubsection:: *)
(*Compute Outputs*)


(* ::Text:: *)
(*With 8 conduction electrons*)


SIcoreparams = SIParams[Tcore,neSIIron[\[Rho]core,8],8];
SIcrustparams = SIParams[Tcrust,neSIIron[\[Rho]crust,8],8];


SIConsts={"e","m","hbar","c","\[Alpha]","\[Epsilon]0","JpereV","kB"}/.SIcoreparams;
SIConstRepl={"e"->SIConsts[[1]],"m"->SIConsts[[2]],"hbar"->SIConsts[[3]],"\[HBar]"->SIConsts[[3]],"c"->SIConsts[[4]],"\[Alpha]"->SIConsts[[5]],"\[Epsilon]0"->SIConsts[[6]],"JpereV"->SIConsts[[7]],"kB"->SIConsts[[8]],"G"->6.67 10^-11};


\[Beta]crust = 1/(Tcrust "kB")/.SIConstRepl;
\[Beta]core = 1/(Tcore "kB")/.SIConstRepl;


(* ::Section::Closed:: *)
(*End*)


End[];


EndPackage[];


(* ::Text:: *)
(*Below can be used to convert string replacement tables to expression replacement tables*)
(**)
(*- Mathematica whines if you try to do it internally to the package*)


(*SIcoreparams = StrToExpConverter[SIcoreparams];
SIcrustparams = StrToExpConverter[SIcrustparams];
SIConstRepl = StrToExpConverter[SIConstRepl];*)
