(* ::Package:: *)

BeginPackage["Capture`"];


(* ::Subsubsection:: *)
(*Public Declarations*)


Getl::usage = "Get interaction length from the sampling formula [m]";
Vgrav::usage = "Gravitational potential for the earth assuming constant density [SI]";
Updatev::usage = "update velocity according to passed velocity mass and interaction length";


(* ::Chapter:: *)
(*Public*)


Begin["`Private`"];


(* ::Text:: *)
(*Propagation to the earth*)


(*EK\[Infinity][m\[Chi]_,v\[Infinity]_]:= 1/2 m\[Chi] v\[Infinity]^2;
E\[Infinity][m\[Chi]_,v\[Infinity]_,r\[Infinity]_,V_]:=EK\[Infinity][m\[Chi],v\[Infinity]]+V[r\[Infinity]];
(*Vtest[r_]:=G\[Alpha]/r;*)
vperp[v\[Infinity]_,b\[Infinity]_,r_]:=(b\[Infinity] v\[Infinity])/r;
vr[m\[Chi]_,v\[Infinity]_,b\[Infinity]_,r\[Infinity]_,V_,r_]:=Sqrt[2/m\[Chi] (E\[Infinity][m\[Chi],v\[Infinity],r\[Infinity],V]-V[r])-vperp[v\[Infinity],b\[Infinity],r]^2];
br[m\[Chi]_,v\[Infinity]_,b\[Infinity]_,r\[Infinity]_,V_,r_]:=b\[Infinity]/Sqrt[1+(V[r\[Infinity]]-V[r])/EK\[Infinity][m\[Chi],v\[Infinity]]];
Vgrav[r_,m\[Chi]_,] :=- ((m\[Chi] "G""ME")/r);*)


(* ::Text:: *)
(*Update trajectory*)


(*vr[r_,E_,V_,m_,vp_]:=Sqrt[2/m (E-V[r])-vp]
V[r_]:=("\[Alpha]")/r*)


(*Block[{m,r,v={vr,vp1,vp2},l,speed,\[CapitalDelta]t,\[CapitalDelta]r,L,vp,E},
speed = Sqrt[v . v];
\[CapitalDelta]t = l/speed;
\[CapitalDelta]r = \[CapitalDelta]t v[[1]];

E=m/2 speed^2;
vp = v . DiagonalMatrix[{0,1,1}] . v;
(*L = m r vp;*)
v[[1]]=vr[r+\[CapitalDelta]r,E,V,m,vp];
v

]*)


(* ::Text:: *)
(*Get interaction length*)


Getl[m\[Chi]_,v\[Chi]_,\[Kappa]_,f_,Natural_:False] := Module[{\[Lambda],m\[Chi]natconv,v\[Chi]natconv},
(*f is a log log function from kinematic variables to the inverse interaction length
\[Kappa] is kinetic mixing
Natural == True is for natural units for kinematic variables ([eV] for mass and [c] for speed])*)
m\[Chi]natconv=("JpereV")/("c")^2/.Constants`SIConstRepl;
v\[Chi]natconv="c"/.Constants`SIConstRepl;
\[Lambda]=If[Natural,10^-f[Log10[m\[Chi] m\[Chi]natconv],Log10[v\[Chi] v\[Chi]natconv]],10^-f[Log10[m\[Chi]],Log10[v\[Chi]]]];
-(\[Lambda]/\[Kappa]^2)Log[1-Random[]]
]


(* ::Text:: *)
(*Gravitational Potential*)


Vgrav[r_]:=Piecewise[{{-(("G" "ME")/r),r>"rE"/. EarthRepl},{-(("ME" "G")/(2("rE")^3))( 3("rE")^2-r^2),r<"rE"/. EarthRepl}}]/.SIConstRepl/.EarthRepl


Updatev[v\[Chi]_,m\[Chi]_,r_,l_,V_]:=Module[{speed,\[CapitalDelta]t,\[CapitalDelta]r,vr,L,vp,E,reached},
speed = \[Sqrt](v\[Chi] . v\[Chi]);
\[CapitalDelta]t = l/speed;
\[CapitalDelta]r = \[CapitalDelta]t v\[Chi][[1]];
E=(m\[Chi]/2)speed^2+V[r];
vp = \[Sqrt](v\[Chi] . DiagonalMatrix[{0,1,1}] . v\[Chi]); (*Subscript[v, \[Perpendicular]]^2*)
vr=\[Sqrt](Simplify[(2/m\[Chi])E]-(2/m\[Chi])V[Abs[r+\[CapitalDelta]r]]-vp^2);(*updated radial velocity*)
reached=(m\[Chi]/2)speed^2-(m\[Chi]/2)vp^2+V[r] >V[Abs[r+\[CapitalDelta]r]];(*check that Energy is high enough to reach r+\[CapitalDelta]r*)
<|"v\[Chi]"->If[reached,{vr,v\[Chi][[2]],v\[Chi][[3]]},{-v\[Chi][[1]],v\[Chi][[2]],v\[Chi][[3]]}],"r"->Abs[r+\[CapitalDelta]r],"reached?"->reached,"E"->E,"vp"->vp,"\[CapitalDelta]r"->\[CapitalDelta]r,"V[r]"->V[r],"V[|r+\[CapitalDelta]r|]"->V[Abs[r+\[CapitalDelta]r]],"l"->l,"vr0"->v\[Chi][[1]],"m\[Chi]"->N[m\[Chi]]|>
]


(* ::Chapter:: *)
(*End*)


End[];


EndPackage[];
