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

AFF[q_,coeffs_]:=Module[{},
	(*assumes q is in [m^-1]]*)
	Sum[coeffs[["a"]][[i]] Exp[-coeffs[["b"]][[i]] ((10^-10  q)/(4 \[Pi]))^2],{i,4}]+coeffs[["c"]]
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


HFF[q_,coeffs_]:= 3 BesselJ[1,q coeffs[["rn"]]]/(q coeffs[["rn"]]) E^(-(q sH)^2/2)


(* ::Subsubsection:: *)
(*Define total Subscript[Z, eff](q)*)


Zeff[q_,coeffs_]:= (coeffs[["Z"]]-AFF[q,coeffs])HFF[q,coeffs]


(* ::Subsubsection:: *)
(*Define FF parameters*)


(* ::Text:: *)
(*Subscript[m, N] is in kg*)


FeNucCoeffs = <|"a"->{11.7695,7.3573,3.5222,2.3045},"b"->{4.7611,0.3072,15.3535,76.8805},"c"->1.0369,"Z"->26,"A"->56,"rn"->rn[56],"mN"->56 1.66 10^-27|>;
SiNucCoeffs = <|"a"->{6.2915,3.0353,1.9891,1.541},"b"->{2.4386,32.3337,0.6785,81.6937},"c"->1.1407,"Z"->14,"A"->28,"rn"->rn[28],"mN"->28 1.66 10^-27|>;
MgNucCoeffs = <|"a"->{5.4204,2.1735,1.2269,2.3073},"b"->{2.8275,79.2611,0.380,7.1937},"c"->0.8584,"Z"->12,"A"->24,"rn"->rn[24],"mN"->24 1.66 10^-27|>;
ONucCoeffs = <|"a"->{3.0485,2.2868,1.5463,0.867},"b"->{13.2771,5.7011,0.3239,32.9089},"c"->0.2508,"Z"->8,"A"->16,"rn"->rn[16],"mN"->16 1.66 10^-27|>;


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
