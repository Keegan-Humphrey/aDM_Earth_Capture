(* ::Package:: *)

BeginPackage["FormFactors`"];


(* ::Text:: *)
(*The purpose of this package is to define nuclear and atomic form factors for material in the earth, in order to estimate energy loss due to DM scattering off nuclei. *)


(* ::Subsubsection::Closed:: *)
(*Public Declarations*)


(* ::Text:: *)
(*FF definitions*)


AFF::usage = "atomic form factors fit to sum over exponentials";
HFF::usage = "Helm nuclear form factor";
Zeff::usage = "Effective nuclear charge include nuclear and atomic form factors";


(* ::Text:: *)
(*Constant definitions*)


FeNucCoeffs::usage = "Fe form factor parameters";
SiNucCoeffs::usage = "Si form factor parameters";
MgNucCoeffs::usage = "Mg form factor parameters";
ONucCoeffs::usage = "O form factor parameters";

sH::usage = "[m] skin depth used in Helm FF";


(* ::Section:: *)
(*Private*)


Begin["`Private`"];


(* ::Subsubsection:: *)
(*Define AFFs*)


(* ::Text:: *)
(*From : https : // lampx . tugraz . at/~hadley/ss1/crystaldiffraction/atomicformfactors/formfactors . php*)
(*and axion scattering paper*)


(*AFF[q_,coeffs_]:=Module[{},
	(*assumes q is in [\[Angstrom]^-1]]*)
	Sum[coeffs[["a"]][[i]] Exp[-coeffs[["b"]][[i]] (q/(4 \[Pi]))^2],{i,4}]+coeffs[["c"]]
]*)
(*AFF[q_,coeffs_]:=Module[{},
	(*assumes q is in [fm^-1]]*)
	Sum[coeffs[["a"]][[i]] Exp[-coeffs[["b"]][[i]] ((10^5 q)/(4 \[Pi]))^2],{i,4}]+coeffs[["c"]]
]*)

(*testflag= True;*)

AFF[q_?NumericQ,coeffs_]:=Module[{},
	(*assumes q is in [m^-1]] coefficients are given in angstrom*)
	
	(*exparglist = coeffs[["b"]]((10^-10  q)/(4 \[Pi]))^2;
	testflag=False;
	If[!NumericQ@exparglist,testflag=True];
	If[testflag, Print[q];];
	If[testflag, Print[Max[exparglist]];];
	If[testflag, Print[exparglist];testflag=False];
	If[Max[exparglist]<500,Sum[coeffs[["a"]][[i]] Exp[-exparglist[[i]]],{i,4}]+coeffs[["c"]],0]*)
	
	(*Print[Max@Table[coeffs[["b"]][[i]]((10^-10  q)/(4 \[Pi]))^2,{i,4}]];*)
	
	If[((Max@coeffs[["b"]])((10^-10  q)/(4 \[Pi]))^2>500),Return[coeffs[["c"]],Module]];
	
	(*Print[Max@coeffs[["b"]]((10^-10  q)/(4 \[Pi]))^2];*)
	
	Sum[coeffs[["a"]][[i]] Exp[-coeffs[["b"]][[i]]((10^-10  q)/(4 \[Pi]))^2],{i,4}]+coeffs[["c"]]
]


(* ::Subsubsection:: *)
(*Define Helm Nuclear FFs*)


(* ::Text:: *)
(*following analytic handbook*)


(*rn[A_]:= 1.14 A^(1/3) 10^-5 (*[\[Angstrom]] nuclear radius*)
sH = 0.9 10^-5; (*[\[Angstrom]] skin depth*)*)

(*rn[A_]:= 1.14 A^(1/3)(*[fm] nuclear radius*)
sH = 0.9 ; (*[fm] skin depth*)*)

rn[A_]:= 1.14 A^(1/3) 10^-15(*[m] nuclear radius*)
sH = 0.9  10^-15; (*[m] skin depth*)


(* ::Text:: *)
(*define Helm FF in momentum space*)


HFF[q_,coeffs_]:= If[(q sH)^2<10^3, 3 BesselJ[1,q coeffs[["rn"]]]/(q coeffs[["rn"]]) E^(-(q sH)^2/2),0]


(* ::Subsubsection:: *)
(*Define total Subscript[Z, eff](q)*)


Zeff[q_?NumericQ,coeffs_]:= (coeffs[["Z"]]-AFF[q,coeffs])HFF[q,coeffs]
(*Zeff[q_,coeffs_]:= AFF[q,coeffs]*)


(* ::Subsubsection:: *)
(*Define FF parameters*)


(* ::Text:: *)
(*Subscript[m, N] is in kg*)


FeNucCoeffs = <|"a"->{11.7695,7.3573,3.5222,2.3045},"b"->{4.7611,0.3072,15.3535,76.8805},"c"->1.0369,"Z"->26,"A"->56,"rn"->rn[56],"mN"->56 1.66 10^-27,"nI"->1.402 10^29|>;
SiNucCoeffs = <|"a"->{6.2915,3.0353,1.9891,1.541},"b"->{2.4386,32.3337,0.6785,81.6937},"c"->1.1407,"Z"->14,"A"->28,"rn"->rn[28],"mN"->28 1.66 10^-27,"nI"->4.995 10^28|>;
MgNucCoeffs = <|"a"->{5.4204,2.1735,1.2269,2.3073},"b"->{2.8275,79.2611,0.380,7.1937},"c"->0.8584,"Z"->12,"A"->24,"rn"->rn[24],"mN"->24 1.66 10^-27,"nI"->4.362 10^28|>;
ONucCoeffs = <|"a"->{3.0485,2.2868,1.5463,0.867},"b"->{13.2771,5.7011,0.3239,32.9089},"c"->0.2508,"Z"->8,"A"->16,"rn"->rn[16],"mN"->16 1.66 10^-27,"nI"->1.435 10^29|>;


(* ::Text:: *)
(*(*Mg - (1.738 10^-3 10^6)/("mN")/.MgNucCoeffs = 4.362449799196787`*^28*)*)
(*(*Fe - "nI"/.SIcoreparams - 1.402 10^29 *)*)
(*(*Si - https://lampz.tugraz.at/~hadley/memm/materials/silicon/silicon.php*)*)
(*(*Oxygen - 2 ("nI"/.SiNucCoeffs) + ("nI"/.MgNucCoeffs) - 1.435 10^29*)*)
(**)


(* ::Text:: *)
(*Correct Z accounting for fit uncertainty, st. FFs vanish exactly at 0 momentum transfer*)


FeNucCoeffs[["Z"]]=Sum[FeNucCoeffs[["a"]][[i]],{i,4}] + FeNucCoeffs[["c"]];
SiNucCoeffs[["Z"]]=Sum[SiNucCoeffs[["a"]][[i]],{i,4}] + SiNucCoeffs[["c"]];
MgNucCoeffs[["Z"]]=Sum[MgNucCoeffs[["a"]][[i]],{i,4}] + MgNucCoeffs[["c"]];
ONucCoeffs[["Z"]]=Sum[ONucCoeffs[["a"]][[i]],{i,4}] + ONucCoeffs[["c"]];


(* ::Section::Closed:: *)
(*End*)


End[];


EndPackage[];
