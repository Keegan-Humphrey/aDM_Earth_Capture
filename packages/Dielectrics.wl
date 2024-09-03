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

\[CapitalZeta]::usage = "Z used in Mermin Dielectric";

SeriesBoundary::usage = "(u,z) boundary point where we switch to approximate formulae";
\[Epsilon]RPAC::usage = "Complex RPA dielectric";

\[Epsilon]RPA0::usage = "Static limit of RPA dielectric";


(* ::Text:: *)
(*Optical*)


\[Epsilon]Mq0::usage = "Degenerate Mermin Dielectric in the long wavelength limit (matches onto Drude model)";
ELFM0::usage = "Loss function constructed out of \[Epsilon]Mq0";

\[Omega]i::usage = "Plasmon frequency as a function of collision rate and position of Drude oscillator peak";


(* ::Text:: *)
(*Electron Evaporation*)


\[Epsilon]nd::usage = "1st correction to degenerate limit to RPA dielectric, used to compute the evaporation rate from electronic scattering in the nearly degenerate case.";
\[Epsilon]Mnd::usage = "1st correction to degenerate limit to RTA dielectric, used to compute the evaporation rate from electronic scattering in the nearly degenerate case.";


(* ::Text:: *)
(*Structure Function*)


kb::usage = "momentum kinematic boundary";
kinTheta::usage = "Heaviside theta for kinematics";
StructureFunction::usage = "Degenerate electronic structure function via fluctuation dissipation theorem";
StructureFunctionuz::usage = "Degenerate electronic structure function via fluctuation dissipation theorem in u and z variables";
StructureFunctionuz1osc::usage = "Degenerate electronic structure function via fluctuation dissipation theorem in u and z variables and 1 oscillator";


(* ::Text:: *)
(*\[Epsilon] Interpolations*)


InterpolateStruc::usage = "Interpolate the structure function";
InterpolateStrucuz::usage = "Interpolate the structure function in dimensionless u and z variables";
InterpolateStrucuzIndivOsc::usage = "Interpolate the structure function in dimensionless u and z variables (individual oscillators)";
InterpolatezIntegranduzIndivOsc::usage = "Interpolate the integrand of z integral for d\[Sigma]dER in dimensionless u and z variables (individual oscillators)";



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


(* ::Subsection:: *)
(*Mermin Dielectric (Random Time Approximation) Private Definitions*)


(* ::Subsubsection::Closed:: *)
(*Subscript[\[Epsilon], RPA](\[Omega] +I \[Nu])*)


\[CapitalZeta][a_,u\[Nu]_]:= Module[{h1,h2},
	(* Used in the Definition of the Mermin Dielectric in the Degenerate Limit
	There's an implicit {1,I} on output ie. [[1]] corresponds to real part and [[2]] corresponds to imaginary part 
	- easier to do this way and keep track of components than deal with refining arctan and logs*)
	h1= -(ArcTan[(1-a)/u\[Nu]] + ArcTan[(1+a)/u\[Nu]]);
	h2= 1/2 Log[((1+a)^2+u\[Nu]^2)/((1-a)^2+u\[Nu]^2)];
	({a,u\[Nu]} +{{a u\[Nu], 1/2 (1-a^2+u\[Nu]^2)},{1/2 (1-a^2+u\[Nu]^2),-a u\[Nu]}} . {h1,h2})
 ]


(*\[CapitalZeta][a,u\[Nu]]//MatrixForm*) 


(*\[Epsilon]RPAC[up_,z_,u\[Nu]_]:={1,0} + "\[Chi]"^2/(4 z^3) (\[CapitalZeta][up+z,u\[Nu]]-\[CapitalZeta][up-z,u\[Nu]])*)


(*Series[\[CapitalZeta]test[a,u\[Nu]],{a,\[Infinity],2}]
largea=%//Normal
(%/.a->(u+z))-(%/.a->(u-z))//Simplify (*far from z=u*)*)
(*Series[\[CapitalZeta]test[a,u\[Nu]],{a,0,10}];
smalla=%//Normal;
(largea/.a->(u+z))-(smalla/.a->(u-z))//Simplify;
Series[%,{u\[Nu],0,1}]//PowerExpand//Normal//Simplify(*close to z=u*)*)
(*\[Epsilon]RPACnearzequ[u_,z_,u\[Nu]]:={-(1/(315 (u+z)))(-210-10 u^10+80 u^9 z+315 (-2+\[Pi] u\[Nu]) z^2+210 z^4+42 z^6+18 z^8+10 z^10+84 u^5 z (2+3 z^2)-18 u^8 (1+15 z^2)+12 u^7 z (9+40 z^2)-42 u^6 (1+6 z^2+10 z^4)+210 u^4 (-1-z^2+2 z^6)-4 u z^3 (105+42 z^2+27 z^4+20 z^6)-12 u^3 z (-35+21 z^4+40 z^6)+3 u^2 (210-105 \[Pi] u\[Nu]+70 z^4+84 z^6+90 z^8)),-(1/2) \[Pi] (-1+u^2-2 u z+z^2)+1/(315 (u+z)^2)2 u\[Nu] (-105+35 u^12-280 u^11 z-315 z^2+315 z^4+105 z^6+63 z^8+45 z^10+35 z^12-10 u^9 z (27+140 z^2)+5 u^10 (9+182 z^2)+12 u^7 z (-21-30 z^2+140 z^4)+u^8 (63+585 z^2+525 z^4)+42 u^5 z (-5+6 z^2+30 z^4+40 z^6)-21 u^6 (-5-12 z^2+30 z^4+140 z^6)-4 u^3 z^3 (-105-63 z^2+90 z^4+350 z^6)+105 u^4 (3-z^2-6 z^4-6 z^6+5 z^8)-2 u z (315+105 z^4+126 z^6+135 z^8+140 z^10)+u^2 (-315-630 z^2-105 z^4+252 z^6+585 z^8+910 z^10))}(*2nd order expansions of \[CapitalZeta][u+z] in large argument - 10th order in a expansion of \[CapitalZeta][u-z] in small argument*)*)
\[Epsilon]RPACnearzequ[u_,z_,u\[Nu]]:={1/315 (315 (-2+\[Pi] u\[Nu])+210 (u-z)^2+42 (u-z)^4+2 (9+5 (u-z)^2) (u-z)^6) (u-z)+2/(3 (u+z)),-2 u\[Nu]-1/2 \[Pi] (-1+(u-z)^2)+2 u\[Nu] (u-z)^2+2/3 u\[Nu] (u-z)^4+2/5 u\[Nu] (u-z)^6+2/63 u\[Nu] (9+7 (u-z)^2) (u-z)^8-(2 u\[Nu])/(3 (u+z)^2)}(*2nd order expansions of \[CapitalZeta][u+z] in large argument - 10th order in a expansion of \[CapitalZeta][u-z] in small argument*)
\[Epsilon]RPAClargea[u_,z_,u\[Nu]_]:={(4 z)/(3 (-u^2+z^2)),(8 u u\[Nu] z)/(3 (u-z)^2 (u+z)^2)}(*large argument expansion of \[CapitalZeta][u+z]-\[CapitalZeta][u-z]*)
(*\[Epsilon]RPACnearzequ[u_,z_,u\[Nu]]:={-(1/(45045 (u+z)))(-30030-630 u^14+7560 u^13 z+45045 (-2+\[Pi] u\[Nu]) z^2+30030 z^4+6006 z^6+2574 z^8+1430 z^10+910 z^12+630 z^14-910 u^12 (1+45 z^2)+1820 u^11 z (5+72 z^2)+1716 u^7 z (9+40 z^2+70 z^4)+2860 u^9 z (4+35 z^2+126 z^4)-1430 u^10 (1+28 z^2+189 z^4)-858 u^8 (3+45 z^2+175 z^4+315 z^6)-12012 u^5 z (-2-3 z^2+10 z^6+30 z^8)+6006 u^6 (-1-6 z^2-10 z^4+45 z^8)+30030 u^4 (-1-z^2+2 z^6+5 z^8+9 z^10)-4 u z^3 (15015+6006 z^2+3861 z^4+2860 z^6+2275 z^8+1890 z^10)-52 u^3 z (-1155+693 z^4+1320 z^6+1925 z^8+2520 z^10)+13 u^2 (6930-3465 \[Pi] u\[Nu]+2310 z^4+2772 z^6+2970 z^8+3080 z^10+3150 z^12)),-(1/2) \[Pi] (-1+u^2-2 u z+z^2)+1/(45045 (u+z)^2)2 u\[Nu] (-15015+3465 u^16-41580 u^15 z-45045 z^2+45045 z^4+15015 z^6+9009 z^8+6435 z^10+5005 z^12+4095 z^14+3465 z^16+315 u^14 (13+704 z^2)-630 u^13 z (65+1078 z^2)+715 u^10 (9+182 z^2+693 z^4)-1820 u^11 z (22+225 z^2+693 z^4)+455 u^12 (11+387 z^2+2772 z^4)+1430 u^9 z (-27-140 z^2-63 z^4+1386 z^6)-3003 u^6 (-5-12 z^2+30 z^4+140 z^6+225 z^8)+1716 u^7 z (-21-30 z^2+140 z^4+630 z^6+1155 z^8)-429 u^8 (-21-195 z^2-175 z^4+1575 z^6+6930 z^8)-6006 u^5 z (5-6 z^2-30 z^4-40 z^6+15 z^8+210 z^10)-4 u^3 z^3 (-15015-9009 z^2+12870 z^4+50050 z^6+102375 z^8+169785 z^10)+15015 u^4 (3-z^2-6 z^4-6 z^6+5 z^8+33 z^10+84 z^12)-2 u z (45045+15015 z^4+18018 z^6+19305 z^8+20020 z^10+20475 z^12+20790 z^14)+u^2 (-45045-90090 z^2-15015 z^4+36036 z^6+83655 z^8+130130 z^10+176085 z^12+221760 z^14))}(*14th order near z=u*)*)
SeriesBoundary = {10^3.5,10^3}; (*(u,z) boundary point where we switch to approximate formulae*)
NearzequBoundary = 5;(*size of |u-z| where we switch from nearzequ expansion to largea expansion*)
\[Epsilon]RPAC[u_,z_,u\[Nu]_]:={1,0} + "\[Chi]"^2/(4 z^3) Piecewise[{{\[Epsilon]RPACnearzequ[u,z,u\[Nu]],Abs[u-z]<5&&(u>=SeriesBoundary[[1]]||z>=SeriesBoundary[[2]])},{\[Epsilon]RPAClargea[u,z,u\[Nu]],Abs[u-z]>NearzequBoundary&&(u>=SeriesBoundary[[1]]||z>=SeriesBoundary[[2]])}},({{1,0},{0,Sign[u]}} . Abs@(\[CapitalZeta][u+z,u\[Nu]]-\[CapitalZeta][u-z,u\[Nu]]))](*abs is for numerical stability , matrix prefactor is for parity of real and imaginary parts under u->-u*)


(* ::Subsection:: *)
(*Mermin Dielectric Public Definitions*)


(* ::Subsubsection:: *)
(*Degenerate*)


(* ::Text:: *)
(*\[Epsilon]M has a removable branch point at z=1 (k=2 Subscript[k, F]). So we make it piecewise and evaluate it on the limit there (dielectric debug notebook)*)


(*Limit[\[Epsilon]M[u,z,u\[Nu]],z->1]*)
(*\[Epsilon]Matbranchpoint[u_,z_,u\[Nu]_]:= 1+(("\[Chi]")^2 (8+4 (-1+u) u\[Nu] (ArcTan[(2-u)/u\[Nu]]+ArcTan[u/u\[Nu]])+4 (1+u) u\[Nu] (ArcTan[u/u\[Nu]]-ArcTan[(2+u)/u\[Nu]])+(-2 u+u^2-u\[Nu]^2) Log[(u^2+u\[Nu]^2)/(4-4 u+u^2+u\[Nu]^2)]-(2 u+u^2-u\[Nu]^2) Log[((2+u)^2+u\[Nu]^2)/(u^2+u\[Nu]^2)]+2 \[ImaginaryI] (-((-2 u+u^2-u\[Nu]^2) (ArcTan[(2-u)/u\[Nu]]+ArcTan[u/u\[Nu]]))+(-2 u-u^2+u\[Nu]^2) (ArcTan[u/u\[Nu]]-ArcTan[(2+u)/u\[Nu]])+(-1+u) u\[Nu] Log[(u^2+u\[Nu]^2)/(4-4 u+u^2+u\[Nu]^2)]-(1+u) u\[Nu] Log[((2+u)^2+u\[Nu]^2)/(u^2+u\[Nu]^2)])))/(2 (8+2 (-2+u+\[ImaginaryI] u\[Nu]) u\[Nu] ArcTan[(2-u)/u\[Nu]]+4 (u+\[ImaginaryI] u\[Nu]) u\[Nu] ArcTan[u/u\[Nu]]-4 u\[Nu] ArcTan[(2+u)/u\[Nu]]-2 u u\[Nu] ArcTan[(2+u)/u\[Nu]]-2 \[ImaginaryI] u\[Nu]^2 ArcTan[(2+u)/u\[Nu]]-2 \[ImaginaryI] u\[Nu] Log[(u^2+u\[Nu]^2)/(4-4 u+u^2+u\[Nu]^2)]+\[ImaginaryI] u u\[Nu] Log[(u^2+u\[Nu]^2)/(4-4 u+u^2+u\[Nu]^2)]-u\[Nu]^2 Log[(u^2+u\[Nu]^2)/(4-4 u+u^2+u\[Nu]^2)]-2 \[ImaginaryI] u\[Nu] Log[((2+u)^2+u\[Nu]^2)/(u^2+u\[Nu]^2)]-\[ImaginaryI] u u\[Nu] Log[((2+u)^2+u\[Nu]^2)/(u^2+u\[Nu]^2)]+u\[Nu]^2 Log[((2+u)^2+u\[Nu]^2)/(u^2+u\[Nu]^2)]))//FullSimplify*)
\[Epsilon]Matbranchcut[u_,z_,u\[Nu]_]:=1+1/2 ("\[Chi]")^2 (1-I u/u\[Nu] (1+8/(-8+2 (-2+u+I u\[Nu]) u\[Nu] ArcTan[(-2+u)/u\[Nu]]-4 (u+I u\[Nu]) u\[Nu] ArcTan[u/u\[Nu]]+u\[Nu] (2 (2+u+I u\[Nu]) ArcTan[(2+u)/u\[Nu]]+(2 I-I u+u\[Nu]) Log[(u^2+u\[Nu]^2)/((-2+u)^2+u\[Nu]^2)]+I (2+u+I u\[Nu]) Log[1+(4 (1+u))/(u^2+u\[Nu]^2)]))))


(* ::Text:: *)
(*Also want to replace the \[Epsilon] by it's large z limit*)


(*Series[\[Epsilon]M[u,z,u\[Nu]],{z,\[Infinity],6}]//Normal//Simplify;
%/.u\[Nu]->u\[Nu]t/."\[Chi]"->\[Chi]//Simplify;
Refine[%,Assumptions->{Element[{u,z,u\[Nu]t,\[Chi],1+z},Reals]}];
cexpandedlargez=ComplexExpand[%]//Simplify;
cexpandedlargez/.Arg[1+z]->0;
Series[%,{z,\[Infinity],6}]//Simplify*)
(*\[Epsilon]Mlargez[u_,z_,u\[Nu]_]:=1+((1+5 u^2+5 \[ImaginaryI] u u\[Nu]) ("\[Chi]")^2)/(15 z^6)+("\[Chi]")^2/(3 z^4)*)
(*\[Epsilon]Mlargez[u_,z_,u\[Nu]_]:=1+((3+35 u^4+42 \[ImaginaryI] u u\[Nu]+70 \[ImaginaryI] u^3 u\[Nu]-7 u^2 (-6+5 u\[Nu]^2)) ("\[Chi]")^2)/(105 z^8)+((1+5 u^2+5 \[ImaginaryI] u u\[Nu]) ("\[Chi]")^2)/(15 z^6)+("\[Chi]")^2/(3 z^4);*)
(*\[Epsilon]Mlargez[u_,z_,u\[Nu]_]:=1+1/(315 z^10)(("\[Chi]")^2 (5+105 u^6+315 \[ImaginaryI] u^5 u\[Nu]+u^2 (135-483 u\[Nu]^2)-315 u^4 (-1+u\[Nu]^2)-21 \[ImaginaryI] u^3 u\[Nu] (-34+5 u\[Nu]^2)-3 \[ImaginaryI] u u\[Nu] (-45+28 u\[Nu]^2))+3 ("\[Chi]")^2 (3+35 u^4+42 \[ImaginaryI] u u\[Nu]+70 \[ImaginaryI] u^3 u\[Nu]-7 u^2 (-6+5 u\[Nu]^2)) z^2+21 ("\[Chi]")^2 (1+5 u^2+5 \[ImaginaryI] u u\[Nu]) z^4+105 ("\[Chi]")^2 z^6)*)
(*\[Epsilon]Mlargez[u_,z_,u\[Nu]_]:=1+(7 I ("\[Chi]")^2 (u+I u\[Nu]) (1+5 u^2+10 I u u\[Nu]-5 u\[Nu]^2+5 z^2) (-2 z+(-1+z^2) Log[(1+z)/Abs[1-z]]))/(z (-560 I u u\[Nu]^4+140 u\[Nu]^5-210 I u z^6+56 I u u\[Nu]^2 (6+10 u^2+5 z^2)-28 u\[Nu]^3 (6+30 u^2+5 z^2)+4 u\[Nu] (3+35 u^4+7 z^2+35 z^4+7 u^2 (6+5 z^2))+105 I u z^5 (-1+z^2) Log[(1+z)/Abs[1-z]])) (*10th order in 1/z*)*)
(*\[Epsilon]Mlargez[u_,z_,u\[Nu]_]:=1+1/(3465 z^12)("\[Chi]")^2 (35+1155 u^8+4620 \[ImaginaryI] u^7 u\[Nu]-4620 \[ImaginaryI] u^5 u\[Nu] (-5+u\[Nu]^2)-462 u^6 (-14+15 u\[Nu]^2)-132 \[ImaginaryI] u^3 u\[Nu] (-127+154 u\[Nu]^2)+231 u^4 (30-136 u\[Nu]^2+5 u\[Nu]^4)+44 \[ImaginaryI] u u\[Nu] (35-66 u\[Nu]^2+21 u\[Nu]^4)+22 u^2 (70-579 u\[Nu]^2+294 u\[Nu]^4))+(("\[Chi]")^2 (5+105 u^6+315 \[ImaginaryI] u^5 u\[Nu]+u^2 (135-483 u\[Nu]^2)-315 u^4 (-1+u\[Nu]^2)-21 \[ImaginaryI] u^3 u\[Nu] (-34+5 u\[Nu]^2)-3 \[ImaginaryI] u u\[Nu] (-45+28 u\[Nu]^2)))/(315 z^10)+(("\[Chi]")^2 (3+35 u^4+42 \[ImaginaryI] u u\[Nu]+70 \[ImaginaryI] u^3 u\[Nu]-7 u^2 (-6+5 u\[Nu]^2)))/(105 z^8)+(("\[Chi]")^2 (1+5 u^2+5 \[ImaginaryI] u u\[Nu]))/(15 z^6)+("\[Chi]")^2/(3 z^4)(*12th order in 1/z*)*)
(*\[Epsilon]Mlargez[u_,z_,u\[Nu]_]:=1+1/(45045 z^14)("\[Chi]")^2 (315+15015 u^10+75075 \[ImaginaryI] u^9 u\[Nu]-30030 \[ImaginaryI] u^7 u\[Nu] (-22+5 u\[Nu]^2)-15015 u^8 (-9+10 u\[Nu]^2)+3003 \[ImaginaryI] u^5 u\[Nu] (350-504 u\[Nu]^2+5 u\[Nu]^4)+15015 u^6 (18-90 u\[Nu]^2+5 u\[Nu]^4)+286 \[ImaginaryI] u^3 u\[Nu] (1322-4203 u\[Nu]^2+1512 u\[Nu]^4)+429 u^4 (350-3734 u\[Nu]^2+2387 u\[Nu]^4)-13 \[ImaginaryI] u u\[Nu] (-1575+5984 u\[Nu]^2-5412 u\[Nu]^4+924 u\[Nu]^6)-13 u^2 (-1575+23518 u\[Nu]^2-34716 u\[Nu]^4+8316 u\[Nu]^6))+1/(3465 z^12)("\[Chi]")^2 (35+1155 u^8+4620 \[ImaginaryI] u^7 u\[Nu]-4620 \[ImaginaryI] u^5 u\[Nu] (-5+u\[Nu]^2)-462 u^6 (-14+15 u\[Nu]^2)-132 \[ImaginaryI] u^3 u\[Nu] (-127+154 u\[Nu]^2)+231 u^4 (30-136 u\[Nu]^2+5 u\[Nu]^4)+44 \[ImaginaryI] u u\[Nu] (35-66 u\[Nu]^2+21 u\[Nu]^4)+22 u^2 (70-579 u\[Nu]^2+294 u\[Nu]^4))+(("\[Chi]")^2 (5+105 u^6+315 \[ImaginaryI] u^5 u\[Nu]+u^2 (135-483 u\[Nu]^2)-315 u^4 (-1+u\[Nu]^2)-21 \[ImaginaryI] u^3 u\[Nu] (-34+5 u\[Nu]^2)-3 \[ImaginaryI] u u\[Nu] (-45+28 u\[Nu]^2)))/(315 z^10)+(("\[Chi]")^2 (3+35 u^4+42 \[ImaginaryI] u u\[Nu]+70 \[ImaginaryI] u^3 u\[Nu]-7 u^2 (-6+5 u\[Nu]^2)))/(105 z^8)+(("\[Chi]")^2 (1+5 u^2+5 \[ImaginaryI] u u\[Nu]))/(15 z^6)+("\[Chi]")^2/(3 z^4)(*14th order in 1/z*)*)
\[Epsilon]Mlargez[u_,z_,u\[Nu]_]:= 1+(("\[Chi]")^2 (3+42 u^2+35 u^4))/(105 z^8)+(("\[Chi]")^2 (1+5 u^2))/(15 z^6)+(I ("\[Chi]")^2 u u\[Nu])/(3 z^6)+("\[Chi]")^2/(3 z^4)(*8th order with u\[Nu]->u\[Nu]/z*)
(*\[Epsilon]Mlargeu[u_,z_,u\[Nu]_]:=1-("\[Chi]")^2/(3 u^2 z^2)+(I ("\[Chi]")^2 u\[Nu])/(3 u^3 z^2)-(("\[Chi]")^2 (3-5 u\[Nu]^2+5 z^2))/(15 u^4 z^2)+(I ("\[Chi]")^2 u\[Nu] (4 z (17-15 u\[Nu]^2+45 z^2)-3 (-9+5 u\[Nu]^2-15 z^2) (-1+z^2) Log[(-1+z)^2]+3 (-9+5 u\[Nu]^2-15 z^2) (-1+z^2) Log[(1+z)^2]))/(45 u^5 z^2 (4 z+(-1+z^2) Log[(-1+z)^2]-(-1+z^2) Log[(1+z)^2]))*)
\[Epsilon]Mlargeu[u_,z_,u\[Nu]_]:=1-(("\[Chi]")^2 (2+u^2-2 I u u\[Nu]-3 u\[Nu]^2))/(3 u^6)+(8 ("\[Chi]")^2 u\[Nu] (I u+3 u\[Nu]))/(1125 u^6 z^6)+(8 ("\[Chi]")^2 u\[Nu] (I u+3 u\[Nu]))/(525 u^6 z^4)-(("\[Chi]")^2 (15+35 u^4-35 I u^3 u\[Nu]-147 u\[Nu]^2+35 u\[Nu]^4+35 I u u\[Nu] (-2+u\[Nu]^2)-7 u^2 (-3+5 u\[Nu]^2)))/(105 u^6 z^2)-(("\[Chi]")^2 z^2)/(3 u^6)(*6th order in 1/u and 1/z*)
\[Epsilon]Mlargeuandz[u_,z_,u\[Nu]_]:=1+(("\[Chi]")^2 (u+I u\[Nu]))/(3 z^4 (u-(4 I u\[Nu])/(3 z^2 (-2+z Log[z/Abs[1-z]]))));


(*\[Epsilon]M[up_,z_,u\[Nu]_]:= Piecewise[{{\[Epsilon]Matbranchcut[up,z,u\[Nu]],z==1},{\[Epsilon]Mlargez[up,z,u\[Nu]],z>10^6&&up<10^6},{\[Epsilon]Mlargeu[up,z,u\[Nu]],z<10^6&&up>10^6},{\[Epsilon]Mlargeuandz[up,z,u\[Nu]],z>10^6&&up>10^6}},1 + ((up +I u\[Nu]) ({1,I} . \[Epsilon]RPAC[up,z,u\[Nu]]-1))/(up + I u\[Nu]  ({1,I} . \[Epsilon]RPAC[up,z,u\[Nu]]-1)/(\[Epsilon]RPA0[z]-1))]*)
(*\[Epsilon]M[up_,z_,u\[Nu]_]:= Piecewise[{{\[Epsilon]Matbranchcut[up,z,u\[Nu]],z==1},{\[Epsilon]Mlargez[up,z,u\[Nu]],z>10^7},{\[Epsilon]Mlargeu[up,z,u\[Nu]],z<10^7&&up>10^7}},1 + ((up +I u\[Nu]) ({1,I} . \[Epsilon]RPAC[up,z,u\[Nu]]-1))/(up + I u\[Nu]  ({1,I} . \[Epsilon]RPAC[up,z,u\[Nu]]-1)/(\[Epsilon]RPA0[z]-1))]*)
(*\[Epsilon]M[up_,z_,u\[Nu]_]:= Piecewise[{{\[Epsilon]Matbranchcut[up,z,u\[Nu]],z==1},{\[Epsilon]Mlargez[up,z,u\[Nu]],z>10^7||up>10^7}},1 + ((up +I u\[Nu]) ({1,I} . \[Epsilon]RPAC[up,z,u\[Nu]]-1))/(up + I u\[Nu]  ({1,I} . \[Epsilon]RPAC[up,z,u\[Nu]]-1)/(\[Epsilon]RPA0[z]-1))]*)

(*\[Epsilon]M[up_,z_,u\[Nu]_]:=1 + ((up +I u\[Nu]) ({1,I} . \[Epsilon]RPAC[up,z,u\[Nu]]-1))/(up + I u\[Nu]  ({1,I} . \[Epsilon]RPAC[up,z,u\[Nu]]-1)/(\[Epsilon]RPA0[z]-1));*)
\[Epsilon]M[up_,z_,u\[Nu]_]:= Piecewise[{{\[Epsilon]Matbranchcut[up,z,u\[Nu]],z==1}},1 + ((up +I u\[Nu]) ({1,I} . \[Epsilon]RPAC[up,z,u\[Nu]]-1))/(up + I u\[Nu]  ({1,I} . \[Epsilon]RPAC[up,z,u\[Nu]]-1)/(\[Epsilon]RPA0[z]-1))]

\[Epsilon]M\[Omega]q[\[Omega]_,k_,\[Nu]_,kF_,vF_]:= \[Epsilon]M[\[Omega]/(vF k),k/(2 kF),\[Nu]/(vF k)] 


\[Epsilon]MNum[up_?NumberQ,z_?NumberQ,u\[Nu]_?NumberQ,params_]:= Quiet@N[\[Epsilon]M[up,z,u\[Nu]] /. params,50]
\[Epsilon]M\[Omega]qNum[\[Omega]_?NumberQ,k_?NumberQ,\[Nu]_?NumberQ,params_]:= Quiet@N[\[Epsilon]M\[Omega]q[\[Omega],k,\[Nu],"qF","vF"] /. params,50]


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


(* ::Subsubsection::Closed:: *)
(*Structure Function*)


u\[Nu]Fit[z_] := ("\[Nu]i"/(2 "vF" "qF" z)) (*u \[Nu]*)
zb[sign_,\[Omega]_,m\[Chi]_,v\[Chi]_,params_]:= ((m\[Chi] v\[Chi] + <|"+"->1,"-"->(-1)|>[[sign]] Sqrt[2 m\[Chi] (1/2 m\[Chi] v\[Chi]^2 - \[Omega] "\[HBar]")])/(2"qF""\[HBar]")) /.params (*section 10.2 of notes*)
kb[sign_,\[Omega]_,m\[Chi]_,v\[Chi]_,params_:Constants`SIConstRepl]:= ((m\[Chi] v\[Chi] + <|"+"->1,"-"->(-1)|>[[sign]] Sqrt[2 m\[Chi] (1/2 m\[Chi] v\[Chi]^2 - \[Omega] "\[HBar]")])/("\[HBar]")) /.params 
Thetab[\[Omega]_,k_,m\[Chi]_,v\[Chi]_,params_]:=HeavisideTheta[kb["+",\[Omega],m\[Chi],v\[Chi],params]-k]HeavisideTheta[k-kb["-",\[Omega],m\[Chi],v\[Chi],params]]HeavisideTheta[1/2 m\[Chi] v\[Chi]^2 - \[Omega] "\[HBar]"/.params]
TotalThetab[\[Omega]_,k_,m\[Chi]_,v\[Chi]_,params_]:=Sum[Thetab[\[Omega],k,m\[Chi],v\[Chi],params[[i]]],{i,Length[params]}]


kinTheta[k_,\[Omega]_,v\[Chi]_,m\[Chi]_] := HeavisideTheta[N[(k v\[Chi]-("\[HBar]" k^2)/(2 m\[Chi])/.Constants`SIConstRepl)-\[Omega],50]]


(* ::Text:: *)
(*Turns out using an If is faster than Heaviside in this case - THIS IS NOT THE STRUCTURE FUNCTION THIS IS THE ARGUMENT OF THE SCATTERING RATE*)


(*StructureFunction[\[Omega]_,k_?NumberQ,m\[Chi]_,v\[Chi]_,params_,n_:4,truncate_:True,cut_:4]:=Abs@Sum[(2 ("e")^2/(v\[Chi] "\[HBar]" (2 \[Pi])^2 "\[Epsilon]0") (1+HeavisideTheta[cut - "\[Beta]" "\[HBar]" \[Omega]](1/(1-E^(- "\[Beta]" "\[HBar]" \[Omega]))-1))/.params[[i]])(z^n/z /.z->k /(2 "qF")/.params[[i]])If[truncate,Thetab[\[Omega],k,m\[Chi],v\[Chi],params[[i]]],1]Im[-("Ai"/.params[[i]]) / Dielectrics`\[Epsilon]MNum[("\[HBar]"\[Omega])/(2"EF"z)/.z->k /(2 "qF")/.params[[i]], k /(2 "qF")/.params[[i]], u\[Nu]Fit[z]/.z-> k /(2 "qF")/. params[[i]], params[[i]]]],{i,Length[params]}]*)
(*StructureFunction[\[Omega]_,k_?NumberQ,m\[Chi]_,v\[Chi]_,params_,n_:0,truncate_:True,cut_:4]:=Abs@Sum[(2 ("e")^2/(v\[Chi] "\[HBar]" (2 \[Pi])^2 "\[Epsilon]0") (1+If[cut > "\[Beta]" "\[HBar]" \[Omega],(1/(1-E^(- "\[Beta]" "\[HBar]" \[Omega]))-1),0])/.params[[i]])(z^n/z /.z->k /(2 "qF")/.params[[i]])If[truncate,Thetab[\[Omega],k,m\[Chi],v\[Chi],params[[i]]],1]Im[-("Ai"/.params[[i]]) / Dielectrics`\[Epsilon]MNum[("\[HBar]"\[Omega])/(2"EF"z)/.z->k /(2 "qF")/.params[[i]], k /(2 "qF")/.params[[i]], u\[Nu]Fit[z]/.z-> k /(2 "qF")/. params[[i]], params[[i]]]],{i,Length[params]}]*)
(*StructureFunction[\[Omega]_,k_?NumberQ,m\[Chi]_,v\[Chi]_,params_,n_:0,truncate_:True,cut_:4]:=Abs@Sum[Im[-("Ai"/.params[[i]]) / Dielectrics`\[Epsilon]MNum[("\[HBar]"\[Omega])/(2"EF"z)/.z->k /(2 "qF")/.params[[i]], k /(2 "qF")/.params[[i]], u\[Nu]Fit[z]/.z-> k /(2 "qF")/. params[[i]], params[[i]]]],{i,Length[params]}]*)
StructureFunction[\[Omega]_,k_?NumberQ,m\[Chi]_,v\[Chi]_,params_,n_:0,truncate_:True,cut_:4]:=Sum[(2 ("e")^2/(v\[Chi]^2 ("\[HBar]")^2 "ne" (2 \[Pi])^2 "\[Epsilon]0") (If[4 >#,If[#>0.01,1/(1-E^(- #)),1/#-1/2],1]&["\[Beta]" "\[HBar]" \[Omega]])/.params[[i]])(z^n/z /.z->k /(2 "qF")/.params[[i]])If[truncate,kinTheta[\[Omega],k,m\[Chi],v\[Chi]],1]Im[-("Ai"/.params[[i]]) / Dielectrics`\[Epsilon]MNum[("\[HBar]"\[Omega])/(2"EF"z)/.z->k /(2 "qF")/.params[[i]], k /(2 "qF")/.params[[i]], u\[Nu]Fit[z]/.z-> k /(2 "qF")/. params[[i]], params[[i]]]],{i,Length[params]}]


(* ::Text:: *)
(*THIS IS THE STRUCTURE FUNCTION*)


(*StructureFunctionuz[u_,z_?NumberQ,params_,cut_:4]:=(*\[Omega] = (4 kF vF u z) *)Abs@Sum[ (2(1+If[cut > "\[Beta]" "\[HBar]" (4 "qF" "vF" u z),(1/(1-E^(- "\[Beta]" "\[HBar]" (2 "qF" "vF" u z)))-1),0])/.params[[i]])((4"\[Epsilon]0" ("qF")^2 z^2)/(("e")^2)  )Im[-("Ai"/.params[[i]]) / Dielectrics`\[Epsilon]MNum[u, z,("\[Nu]i"/(2 "vF" "qF" z))/.params[[i]],params[[i]]]]/.params[[i]],{i,Length[params]}]*)
(*StructureFunctionuz[u_,z_?NumberQ,params_,cut_:4]:=(*\[Omega] = (4 kF vF u z) *)Abs@Sum[ Im[-("Ai"/.params[[i]]) / Dielectrics`\[Epsilon]MNum[u, z,("\[Nu]i"/(2 "vF" "qF" z))/.params[[i]],params[[i]]]]/.params[[i]],{i,Length[params]}]*)
StructureFunctionuz[u_,z_?NumberQ,params_,cut_:4]:=(*\[Omega] = (4 kF vF u z) *)Abs@Sum[ (2(If[4 >#,If[#>0.01,1/(1-E^(- #)),1/#-1/2],1]&["\[Beta]" "\[HBar]" (4 "qF" "vF" u z)])/.params[[i]])((4"\[Epsilon]0" ("qF")^2 z^2)/(("e")^2)  /.params[[i]])Im[-("Ai"/.params[[i]]) / Dielectrics`\[Epsilon]MNum[u, z,("\[Nu]i"/(2 "vF" "qF" z))/.params[[i]],params[[i]]]],{i,Length[params]}]
StructureFunctionuz1osc[u_,z_?NumberQ,params_,cut_:4]:=(*\[Omega] = (4 kF vF u z) *)Abs@( (2(If[4 >#,If[#>0.01,1/(1-E^(- #)),1/#-1/2],1]&["\[Beta]" "\[HBar]" (4 "qF" "vF" u z)])/.params)((4"\[Epsilon]0" ("qF")^2 z^2)/(("e")^2)  /.params)Im[-("Ai"/.params) / Dielectrics`\[Epsilon]MNum[u, z,("\[Nu]i"/(2 "vF" "qF" z))/.params,params]])
zIntegranduz1osc[u_,z_?NumberQ,params_,cut_:4]:=(*\[Omega] = (4 kF vF u z) *)Abs@( (2(If[4 >#,If[#>0.01,1/(1-E^(- #)),1/#-1/2],1]&["\[Beta]" "\[HBar]" (4 "qF" "vF" u z)])/.params) 1/z Im[-("Ai"/.params) / Dielectrics`\[Epsilon]MNum[u, z,("\[Nu]i"/(2 "vF" "qF" z))/.params,params]]/.params)


(* ::Subsubsection::Closed:: *)
(*Interpolate Structure Function*)


Clear[InterpolateStruc]
InterpolateStruc[params_]:=Module[{mmax,vmax,\[Omega]max,kmax,\[Omega]maxminrat,kmaxminrat,gridsize=200,Interpoints\[Omega],Interpointsk,Interpoints\[Omega]andk,Interpoints\[Epsilon],Inter\[Epsilon]f},
mmax = 10^15 ("JpereV")/("c")^2/.Constants`SIConstRepl;
vmax = 0.2 "c"/.Constants`SIConstRepl;
kmax = (mmax vmax)/("\[HBar]")/.Constants`SIConstRepl;
\[Omega]max= (mmax vmax^2)/(2 "\[HBar]")/.Constants`SIConstRepl;
\[Omega]maxminrat = 10^-30;
kmaxminrat = 10^-20;

Interpoints\[Omega]=10^Subdivide[Log10[\[Omega]maxminrat \[Omega]max],Log10[\[Omega]max],gridsize];
Interpointsk=10^Subdivide[Log10[kmaxminrat kmax],Log10[kmax],gridsize];

Interpoints\[Omega]andk =Table[ Table[{Interpointsk[[i]],Interpoints\[Omega][[j]]},{i,gridsize}],{j,gridsize}];

Interpoints\[Epsilon]=Table[ Table[{{Interpointsk[[i]],Interpoints\[Omega][[j]]},10^-60+StructureFunction[Interpointsk[[i]],Interpoints\[Omega][[j]],mmax,vmax,params,0,False]},{i,gridsize}],{j,gridsize}];

(*\[Epsilon][k,\[Omega]]*)
Inter\[Epsilon]f=Interpolation[Flatten[Log10[Interpoints\[Epsilon]],1],InterpolationOrder->3];

Inter\[Epsilon]f
]


InterpolateStrucuz[params_]:=Module[{mmax,vmax,umax,zmax,umaxminrat,zmaxminrat,gridsize=600,Interpointsu,Interpointsz,Interpoints\[Omega]andk,Interpoints\[Epsilon],Inter\[Epsilon]f},
mmax = 10^15 ("JpereV")/("c")^2/.Constants`SIConstRepl;
vmax = 0.2 "c"/.Constants`SIConstRepl;
umax = Max[vmax/(2 "vF")/.params];
zmax= Max[(mmax vmax)/(2 "qF""\[HBar]")/.params];
umaxminrat = 10^-20;
zmaxminrat = 10^-20;

Interpointsu=10^Subdivide[Log10[umaxminrat umax],Log10[umax],gridsize];
Interpointsz=10^Subdivide[Log10[zmaxminrat zmax],Log10[zmax],gridsize];

(*Interpoints\[Omega]andk =Table[ Table[{Interpoints\[Omega][[i]],Interpointsk[[j]]},{i,gridsize}],{j,gridsize}];*)

Interpoints\[Epsilon]=Table[Table[{{Interpointsu[[i]],Interpointsz[[j]]},StructureFunctionuz[Interpointsu[[i]],Interpointsz[[j]],params]},{i,gridsize}],{j,gridsize}];

(*Print[Flatten[Log10[Interpoints\[Epsilon]],1][[;;5]]];*)

(*\[Epsilon][\[Omega],k]*)
Inter\[Epsilon]f=Interpolation[Flatten[Log10[Interpoints\[Epsilon]],1],InterpolationOrder->3];

Inter\[Epsilon]f
]


InterpolateStrucuzIndivOsc[params_]:=Module[{mmax,vmax,umax,zmax,umaxminrat,zmaxminrat,gridsize=600,Interpointsu,Interpointsz,Interpoints\[Omega]andk,Interpoints\[Epsilon],Inter\[Epsilon]f},
mmax = 10^15 ("JpereV")/("c")^2/.Constants`SIConstRepl;
vmax = 0.2 "c"/.Constants`SIConstRepl;
umax = Max[vmax/(2 "vF")/.params];
zmax= Max[(mmax vmax)/(2 "qF""\[HBar]")/.params];
umaxminrat = 10^-20;
zmaxminrat = 10^-20;

Interpointsu=10^Subdivide[Log10[umaxminrat umax],Log10[umax],gridsize];
Interpointsz=10^Subdivide[Log10[zmaxminrat zmax],Log10[zmax],gridsize];

(*Interpoints\[Omega]andk =Table[ Table[{Interpoints\[Omega][[i]],Interpointsk[[j]]},{i,gridsize}],{j,gridsize}];*)

Interpoints\[Epsilon]=Table[Table[Table[{{Interpointsu[[i]],Interpointsz[[j]]},StructureFunctionuz1osc[Interpointsu[[i]],Interpointsz[[j]],params[[osc]]]},{i,gridsize}],{j,gridsize}],{osc,Length[params]}];

(*Print[Flatten[Log10[Interpoints\[Epsilon]],1][[;;5]]];*)

(*Table[{\[Epsilon][u,z],params}]*)
Inter\[Epsilon]f=Table[{Interpolation[Flatten[Log10[Interpoints\[Epsilon][[osc]]],1],InterpolationOrder->3],params[[osc]]},{osc,Length[params]}];

Inter\[Epsilon]f
]


InterpolatezIntegranduzIndivOsc[params_]:=Module[{mmax,vmax,umax,zmax,umaxminrat,zmaxminrat,gridsize=600,Interpointsu,Interpointsz,Interpoints\[Omega]andk,Interpoints\[Epsilon],Inter\[Epsilon]f},
mmax = 10^15 ("JpereV")/("c")^2/.Constants`SIConstRepl;
vmax = 0.2 "c"/.Constants`SIConstRepl;
umax = Max[vmax/(2 "vF")/.params];
zmax= Max[(mmax vmax)/(2 "qF""\[HBar]")/.params];
umaxminrat = 10^-20;
zmaxminrat = 10^-20;

Interpointsu=10^Subdivide[Log10[umaxminrat umax],Log10[umax],gridsize];
Interpointsz=10^Subdivide[Log10[zmaxminrat zmax],Log10[zmax],gridsize];

(*Interpoints\[Omega]andk =Table[ Table[{Interpoints\[Omega][[i]],Interpointsk[[j]]},{i,gridsize}],{j,gridsize}];*)

Interpoints\[Epsilon]=Table[Table[Table[{{Interpointsu[[i]],Interpointsz[[j]]},zIntegranduz1osc[Interpointsu[[i]],Interpointsz[[j]],params[[osc]]]},{i,gridsize}],{j,gridsize}],{osc,Length[params]}];

(*Print[Flatten[Log10[Interpoints\[Epsilon]],1][[;;5]]];*)

(*Table[{\[Epsilon][u,z],params}]*)
Inter\[Epsilon]f=Table[{Interpolation[Flatten[Log10[Interpoints\[Epsilon][[osc]]],1],InterpolationOrder->3],params[[osc]]},{osc,Length[params]}];

Inter\[Epsilon]f
]


(* ::Subsection:: *)
(*RTA in the nearly degenerate case *)


(* ::Subsubsection:: *)
(*\[Epsilon]nd*)


h[a_,u\[Nu]_]:= 1/2 Log[((1+a)^2+u\[Nu]^2)/((1-a)^2+u\[Nu]^2)] - I (ArcTan[(1-a)/u\[Nu]]+ArcTan[(1+a)/u\[Nu]]) 
\[Epsilon]nd[u_,z_,u\[Nu]_]:= 1 + ("\[Chi]")^2/(4 z^3) (1/("D")-z)(h[u+z,u\[Nu]]-h[u-z,u\[Nu]])


(* ::Subsubsection:: *)
(*\[Epsilon]nd(k,0)*)


\[Epsilon]nd0[z_]:=Module[{h1},
	(*
		Input:
		z - [] q/(2 qF)
		u - [] \[Omega]/(q qF)
	
	Output:
		\[Epsilon]^Lin - [] Lindhard dielectric function at 0 temperature
	*)
	(*\[Chi] = Sqrt["e"^2/(Pi "hbar" "vF") 1/(4 \[Pi] "\[Epsilon]0")];*)
	h1= Log[Abs[(1+z)/(1-z)]];
	
	1 + "\[Chi]"^2/(4 z^3)(1/("D")-z) (2 h1 )
]


(* ::Subsubsection:: *)
(*Mermin \[Epsilon]nd*)


\[Epsilon]Mnd[up_,z_,u\[Nu]_]:= 1 + ((up +I u\[Nu]) ( \[Epsilon]nd[up,z,u\[Nu]]-1))/(up + I u\[Nu]  (\[Epsilon]nd[up,z,u\[Nu]]-1)/(\[Epsilon]nd0[z]-1))
(*Piecewise[{{\[Epsilon]Matbranchcut[up,z,u\[Nu]],z==1}},1 + ((up +I u\[Nu]) ({1,I} . \[Epsilon]RPAC[up,z,u\[Nu]]-1))/(up + I u\[Nu]  ({1,I} . \[Epsilon]RPAC[up,z,u\[Nu]]-1)/(\[Epsilon]RPA0[z]-1))]*)


(* ::Section::Closed:: *)
(*End*)


End[];


EndPackage[];
