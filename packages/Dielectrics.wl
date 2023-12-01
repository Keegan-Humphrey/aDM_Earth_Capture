(* ::Package:: *)

BeginPackage["Dielectrics`"];


(* ::Text:: *)
(*The purpose of this package is to define the Degenerate Mermin Dielectric at finite energy and momentum transfer.*)


(* ::Subsubsection:: *)
(*Public Declarations*)


(* ::Text:: *)
(*Dimensionless variables*)


uf::usage = "dimensionless variable u in terms of q and \[Omega]";
zf::usage = "dimensionless variable z in terms of q";


(* ::Text:: *)
(*Degenerate*)


uzLinDielectric::usage = "Lindhard Degenerate Dielectric in u z variables";
LinDielectric::usage = "Lindhard Degenerate Dielectric in q \[Omega] variables";

\[Epsilon]M::usage = "Degenerate Mermin Dielectric as a function of dimensionless variables u, z, u\[Nu] (see overleaf)";
\[Epsilon]M\[Omega]q::usage = "Degenerate Mermin Dielectric as a function of dimensionful variables \[Omega], k, \[Nu] (see overleaf)";

\[Epsilon]MNum::usage = "Numerical argument version of \[Epsilon]M (for Numerical Integration)";
\[Epsilon]M\[Omega]qNum::usage = "Numerical argument version of \[Epsilon]M\[Omega]q (for Numerical Integration)";

\[Epsilon]RPA0::usage = "Static limit of RPA dielectric";


(* ::Text:: *)
(*Optical*)


\[Epsilon]Mq0::usage = "Degenerate Mermin Dielectric in the long wavelength limit (matches onto Drude model)";
ELFM0::usage = "Loss function constructed out of \[Epsilon]Mq0";

\[Omega]i::usage = "Plasmon frequency as a function of collision rate and position of Drude oscillator peak";


(* ::Section:: *)
(*Private*)


Begin["`Private`"];


(* ::Subsection:: *)
(*Lindhard Dielectric (Random Phase Approximation) Private Definitions*)


(* ::Text:: *)
(*See the "0 Temperature Limit and canonical Lindhard Function" section of Overleaf aDM notes for details and explanation of the following analysis - including relevant citations. *)


zf[q_]:=q/(2 "qF") 
uf[q_,\[Omega]_]:= \[Omega]/(q "vF")


uzLinDielectric[u_,z_]:=Module[{f1,f2,g,x,\[Chi]},
	(*
		Input:
		z - [] q/(2 qF)
		u - [] \[Omega]/(q qF)
	
	Output:
		\[Epsilon]^Lin - [] Lindhard dielectric function at 0 temperature
	*)
	(*\[Chi] = Sqrt["e"^2/(Pi "hbar" "vF") 1/(4 \[Pi] "\[Epsilon]0")];*)
	g=(1-x^2) Log[Abs[(1+x)/(1-x)]];
	f1= 1/2 + 1/(8 z) ((g/.x->z-u)+(g/.x->z+u));
	f2= Piecewise[{{\[Pi]/2 u,z+u<1},{\[Pi] /(8z)(1-(z-u)^2),Abs[z-u]<1<z+u},{0,Abs[z-u]>1}}];
	
	1 + "\[Chi]"^2/z^2 ( f1 + I f2)
]


LinDielectric[q_,\[Omega]_] := uzLinDielectric[uf[q,\[Omega]],zf[q]] (*Lindhard Dielectric in q, \[Omega] variables*)


\[Epsilon]RPA0[z_]:=Refine[uzLinDielectric[0,z],Assumptions->z>0]


(* ::Subsection::Closed:: *)
(*Mermin Dielectric (Random Time Approximation) Private Definitions*)


(* ::Subsubsection::Closed:: *)
(*Subscript[\[Epsilon], RPA](\[Omega] +I \[Nu])*)


\[CapitalZeta][a_,u\[Nu]_]:= Module[{h1,h2},
	(* Used in the Definition of the Mermin Dielectric in the Degenerate Limit
	There's an implicit {1,I} on output ie. [[1]] corresponds to real part and [[2]] corresponds to imaginary part 
	- easier to do this way and keep track of components than deal with refining arctan and logs*)
	h1=- (ArcTan[(1-a)/u\[Nu]] + ArcTan[(1+a)/u\[Nu]]);
	h2= 1/2 Log[((1+a)^2+u\[Nu]^2)/((1-a)^2+u\[Nu]^2)];
	({a,u\[Nu]} +{{a u\[Nu], 1/2 (1-a^2+u\[Nu]^2)},{1/2 (1-a^2+u\[Nu]^2),-a u\[Nu]}} . {h1,h2})
 ]


(*\[CapitalZeta][a,u\[Nu]]//MatrixForm*) 


\[Epsilon]RPAC[up_,z_,u\[Nu]_]:={1,0} + "\[Chi]"^2/(4 z^3) (\[CapitalZeta][up+z,u\[Nu]]-\[CapitalZeta][up-z,u\[Nu]])


(* ::Subsection::Closed:: *)
(*Mermin Dielectric Public Definitions*)


(* ::Subsubsection::Closed:: *)
(*Degenerate*)


\[Epsilon]M[up_,z_,u\[Nu]_]:= 1 + ((up +I u\[Nu]) ({1,I} . \[Epsilon]RPAC[up,z,u\[Nu]]-1))/(up + I u\[Nu]  ({1,I} . \[Epsilon]RPAC[up,z,u\[Nu]]-1)/(\[Epsilon]RPA0[z]-1))
\[Epsilon]M\[Omega]q[\[Omega]_,k_,\[Nu]_,kF_,vF_]:= \[Epsilon]M[\[Omega]/(vF k),k/(2 kF),\[Nu]/(vF k)] 



\[Epsilon]MNum[up_?NumberQ,z_?NumberQ,u\[Nu]_?NumberQ,params_]:= \[Epsilon]M[up,z,u\[Nu]] /. params
\[Epsilon]M\[Omega]qNum[\[Omega]_?NumberQ,k_?NumberQ,\[Nu]_?NumberQ,params_]:= \[Epsilon]M\[Omega]q[\[Omega],k,\[Nu],"qF","vF"] /. params


(* ::Subsubsection::Closed:: *)
(*Optical Dielectrics*)


(* ::Text:: *)
(*See overleaf for descriptions*)
(**)
(*These are the Long Wavelength Limit of the Mermin Dielectric (which matches onto the Drude Model) and the corresponding loss function. (commented code gives \[Epsilon]Mq0)*)


(*\[Epsilon]M[\[Omega]/(2 qF vF z ),z,\[Nu]/(2 qF vF z )];
Limit[%,z->0]*)


\[Epsilon]Mq0[\[Omega]_,\[Nu]_,\[Omega]p_] := 1 -\[Omega]p^2/(\[Omega] (\[Omega]+I \[Nu]))


ELFM0[\[Omega]_,\[Nu]_,\[Omega]p_]:=Simplify[ComplexExpand[Im[-1/\[Epsilon]Mq0[\[Omega],\[Nu],\[Omega]p]]]]


(* ::Text:: *)
(*See overleaf, the below is the plasmon frequency for a Drude Dielectric with collision rate \[Nu]i and peak position \[Omega]exti*)
(*(find extremum of Drude, solve for value of \[Omega]p, pick positive solution)*)


(*Simplify[D[ELFM0[\[Omega],\[Nu],\[Omega]p],\[Omega]]]
Solve[%==0,\[Omega]]
\[Omega]p^2/.Solve[(\[Omega]/.%[[4]])==\[Omega],\[Omega]p]*)


\[Omega]i[\[Omega]exti_,\[Nu]i_]:=Sqrt[Sqrt[\[Nu]i^2 \[Omega]exti^2 + 4 \[Omega]exti^4] - \[Omega]exti^2]


(* ::Section::Closed:: *)
(*End*)


End[];


EndPackage[];
