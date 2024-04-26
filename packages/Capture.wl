(* ::Package:: *)

Needs["Constants`"]


BeginPackage["Capture`"];


(* ::Subsubsection:: *)
(*Public Declarations*)


Getl::usage = "Get interaction length from the sampling formula [m]";
Vgrav::usage = "Gravitational potential for the earth assuming constant density [SI]";

Updatev::usage = "update velocity according to passed velocity mass and interaction length";
Rotatev\[Chi]::usage = "Rotate v\[Chi] about an axis orthogonal to it, and parallel to it";
RotateuAboutvby\[Phi]::usage = "Rotate a vector about another by an angle using Rodrigues formula";

\[Sigma]selector::usage = "Sample to find the interaction process that occurred using the cross-sections of possible processes";

Interpolate\[Epsilon]euztrunced::usage = "Interpolate \[Epsilon] to speed up function evaluation time for the Monte-Carlo";
d\[Sigma]dERd\[CapitalOmega]eofv\[Chi]and\[Epsilon]::usage = "Interpolate over \[Epsilon] for electronic contribution to reduce time for a function evaluation of \[Epsilon]. We do it over the range we will be probing in the Monte-Carlo";
Get\[Xi]\[Omega]andcos\[Theta]::usage = "Sample the distribution of d\[Sigma]dERd\[CapitalOmega] for \[Xi]\[Omega] and cos\[Theta]";
GetNMaxofd\[Sigma]dERd\[CapitalOmega]::usage = "Numerically solve for the maximum of d\[Sigma]dERd\[CapitalOmega] slightly away from the boundaries, and with a small rescaling buffer (for use with the rejection method)";


(* ::Chapter:: *)
(*Public*)


Begin["`Private`"];


(* ::Subsubsection:: *)
(*Propagation to the earth*)


(*EK\[Infinity][m\[Chi]_,v\[Infinity]_]:= 1/2 m\[Chi] v\[Infinity]^2;
E\[Infinity][m\[Chi]_,v\[Infinity]_,r\[Infinity]_,V_]:=EK\[Infinity][m\[Chi],v\[Infinity]]+V[r\[Infinity]];
(*Vtest[r_]:=G\[Alpha]/r;*)
vperp[v\[Infinity]_,b\[Infinity]_,r_]:=(b\[Infinity] v\[Infinity])/r;
vr[m\[Chi]_,v\[Infinity]_,b\[Infinity]_,r\[Infinity]_,V_,r_]:=Sqrt[2/m\[Chi] (E\[Infinity][m\[Chi],v\[Infinity],r\[Infinity],V]-V[r])-vperp[v\[Infinity],b\[Infinity],r]^2];
br[m\[Chi]_,v\[Infinity]_,b\[Infinity]_,r\[Infinity]_,V_,r_]
:=b\[Infinity]/Sqrt[1+(V[r\[Infinity]]-V[r])/EK\[Infinity][m\[Chi],v\[Infinity]]];
Vgrav[r_,m\[Chi]_,] :=- ((m\[Chi] "G""ME")/r);*)


(* ::Subsubsection:: *)
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


(* ::Subsubsection:: *)
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


(* ::Subsubsection:: *)
(*Gravitational Potential*)


Vgrav[r_]:=Piecewise[{{-(("G" "ME")/r),r>"rE"/. Constants`EarthRepl},{-(("ME" "G")/(2("rE")^3))( 3("rE")^2-r^2),r<"rE"/. Constants`EarthRepl}}]/.Constants`SIConstRepl/.Constants`EarthRepl


(* ::Subsubsection:: *)
(*Update velocity according to interaction length and direction of propagation*)


Updatev[v\[Chi]_,m\[Chi]_,r_,l_,V_]:=Module[{speed,\[CapitalDelta]t,\[CapitalDelta]r,vr,L,vp,E,reached},
speed = \[Sqrt](v\[Chi] . v\[Chi]);
\[CapitalDelta]t = l/speed;
\[CapitalDelta]r = \[CapitalDelta]t v\[Chi][[1]];
E=(m\[Chi]/2)speed^2+m\[Chi] V[r];
vp = \[Sqrt](v\[Chi] . DiagonalMatrix[{0,1,1}] . v\[Chi]); (*Subscript[v, \[Perpendicular]]^2*)
(*vr=Sqrt[(Simplify[(2/m\[Chi])E]-2(V[Abs[r+\[CapitalDelta]r]]-V[r])-vp^2)];(*updated radial velocity*)*)
vr=-Sqrt[((2/m\[Chi])(E-m\[Chi] (V[r]+V[Abs[r+\[CapitalDelta]r]]))-vp^2)];(*updated radial velocity*)
reached=(1/2)speed^2-(1/2)vp^2+V[r] >V[Abs[r+\[CapitalDelta]r]];(*check that Energy is high enough to reach r+\[CapitalDelta]r*)
<|"v\[Chi]"->If[reached,{vr,v\[Chi][[2]],v\[Chi][[3]]},{-v\[Chi][[1]],v\[Chi][[2]],v\[Chi][[3]]}],"r"->Abs[r+\[CapitalDelta]r],"vr"->vr,"reached?"->reached,"E"->E,"vp"->vp,"\[CapitalDelta]r"->\[CapitalDelta]r,"V[r]"->V[r],"V[|r+\[CapitalDelta]r|]"->V[Abs[r+\[CapitalDelta]r]],"l"->l,"vr0"->v\[Chi][[1]],"m\[Chi]"->N[m\[Chi]]|>
]


RotateuAboutvby\[Phi][u_,v_,\[Phi]_]:= u Cos[\[Phi]] + v (u . v) (1-Cos[\[Phi]]) - Cross[u,v] Sin[\[Phi]](*Rodrigues formula for rotation of u about v (unit vector) by \[Phi]*)


Rotatev\[Chi][v\[Chi]_,cos\[Theta]_,m\[Chi]_,\[Phi]_,q_]:=Module[{v\[Chi]hat, w,\[CapitalDelta]v\[Chi],v\[Chi]p},
v\[Chi]hat=v\[Chi]/Sqrt[v\[Chi] . v\[Chi]];
w=Cross[v\[Chi]hat,{0,1,0}];
w=w/Sqrt[w . w];
(*Print[w];*)
(*Print[v\[Chi]hat];*)
\[CapitalDelta]v\[Chi]=("\[HBar]" q )/m\[Chi] (w Sqrt[1-cos\[Theta]^2] + v\[Chi]hat cos\[Theta])/.Constants`SIConstRepl;
(*Print[\[CapitalDelta]v\[Chi]];*)

(*Rotate about v\[Chi] by \[Phi]*)
v\[Chi]p = v\[Chi] - RotateuAboutvby\[Phi][\[CapitalDelta]v\[Chi],v\[Chi]hat,\[Phi]]
]


(* ::Subsubsection:: *)
(*Select Interaction Process*)


\[Sigma]selector[\[Sigma]s_]:=Module[{\[Xi]2,Ps,Psums,N\[Sigma],n},

\[Xi]2 = Random[];

N\[Sigma]=Length[\[Sigma]s];
Ps=\[Sigma]s/Total[\[Sigma]s];
Psums=Table[Sum[Ps[[i]],{i,j}],{j,Length[Ps]}];

n=0;
Do[If[ \[Xi]2<=Psums[[i]],n=i;Break[]],{i,1,N\[Sigma]}];
(*{n,\[Xi]2,N[Psums]}*)
n
]


(* ::Subsection:: *)
(*Interpolation*)


Interpolate\[Epsilon]euztrunced[m\[Chi]_,params_,uselog_:True,m_:40,n_:40,\[Delta]_:10^-12]:=Module[{v\[Chi]max,\[Omega]min,\[Omega]max,qmin,qmax,umin,umax,zmin,zmax,InterPu,InterPz,Inter\[Epsilon]Table,Inter\[Epsilon]f},
(*Interpolate over \[Epsilon] for electronic contribution to reduce time for a function evaluation of \[Epsilon]. We do it over the range we will be probing in the Monte-Carlo*)

(*Find the range over which we are interpolating*)
v\[Chi]max=2 10^-2 "c"/.params; (*nonrelativistic limit*)
\[Omega]min= \[Delta] "\[Omega]edgei"/.params;
\[Omega]max = 1/(2"\[HBar]") m\[Chi] v\[Chi]max^2/.params;
qmax=1/("\[HBar]") m\[Chi] v\[Chi]max/.params;
qmin=\[Delta] qmax;(*\[Delta] = Subscript[q, max]/Subscript[q, min] - necessary since the \[Epsilon] has a removable singularity at q=0*)

(*Now convert to u and z variables*)
umin = \[Omega]min/( qmax "vF")/.params;
umax = Min[\[Omega]max/(qmin "vF")/.params,10^3];(*Truncate u at a large value of \[Omega] (where we get numerical instabilities)*)
zmin = qmin/(2 "qF")/.params;
zmax = qmax/(2 "qF")/.params;

(*Log Subdivide the intervals*)
     InterPu= If[uselog,10^Subdivide[Log10[umin],Log10[umax],m],Subdivide[umin,umax,m]];
InterPz= If[uselog,10^Subdivide[Log10[zmin],Log10[zmax],n],Subdivide[zmin,zmax,n]];

Inter\[Epsilon]Table=Table[{{InterPu[[i]],InterPz[[j]]},Abs[Im[-1/Dielectrics`\[Epsilon]M[InterPu[[i]],InterPz[[j]],("\[Nu]i")/(2 "qF""vF"InterPz[[j]] )]]/.params]},{i,m+1},{j,n+1}];
Inter\[Epsilon]f=If[uselog,Interpolation[Flatten[Log10[Inter\[Epsilon]Table],1],InterpolationOrder->4],Interpolation[Flatten[Inter\[Epsilon]Table,1],InterpolationOrder->4]];

<|"\[Epsilon]f"->Inter\[Epsilon]f,"\[Epsilon]Table"->Inter\[Epsilon]Table,"m\[Chi]"->m\[Chi],"v\[Chi]max"->v\[Chi]max,"us"->InterPu,"zs"->InterPz,"\[Omega]range"->{\[Omega]min,\[Omega]max},"qrange"->{qmin,qmax},"urange"->{umin,umax},"zrange"->{zmin,zmax},"params"->params,"Npoints"->{m,n},"\[Delta]"->\[Delta]|>
]


d\[Sigma]dERd\[CapitalOmega]eofv\[Chi]and\[Epsilon][\[Xi]\[Omega]_,cos\[Theta]_,v\[Chi]_,\[Kappa]_,\[Epsilon]Dict_]:=Module[{m\[Chi],\[Omega]max,\[Omega],qpm,\[Xi]\[Theta],prefactor,d\[Sigma]dERd\[CapitalOmega],loguzs},
(* Get d\[Sigma]/(Subscript[dE, R]d\[CapitalOmega]) given an interpolated \[Epsilon] dictionary. Do this as a function of dimensionless variables which will be sampled.
   \[Xi]\[Omega] \[Element] (0,1) - []
  cos\[Theta] \[Epsilon] (0,1) - [] negative disallowed since gas is assumed degenerate: hence no energy can be removed from the Fermi sphere.*)
m\[Chi]=\[Epsilon]Dict["m\[Chi]"];
\[Omega]max =(m\[Chi] v\[Chi]^2)/(2 "\[HBar]")/.\[Epsilon]Dict["params"];
\[Omega]= (\[Epsilon]Dict["\[Omega]range"][[1]]+(cos\[Theta]^2 \[Omega]max-\[Epsilon]Dict["\[Omega]range"][[1]])\[Xi]\[Omega]); (*define \[Omega] in terms of \[Xi]\[Omega]*)
qpm={q->(m\[Chi] v\[Chi])/("\[HBar]") (cos\[Theta]-Sqrt[cos\[Theta]^2-(2\[Omega] "\[HBar]")/(m\[Chi] v\[Chi]^2)]),q->(m\[Chi] v\[Chi])/("\[HBar]") (cos\[Theta]+Sqrt[cos\[Theta]^2-(2\[Omega] "\[HBar]")/(m\[Chi] v\[Chi]^2)])}/.\[Epsilon]Dict["params"];(*solutions for q to \[Omega]-\[Omega](q,cos\[Theta])=0 in the delta function that arises in the q integral. These must be summed over*)
\[Xi]\[Theta]=Sqrt[cos\[Theta]^2-(2 "\[HBar]" \[Omega])/(m\[Chi] v\[Chi]^2)]/.\[Epsilon]Dict["params"];(*the absolute value of the jacobian of the delta function, evaluated on the solution of q, for each of the above solutions*)

loguzs=Table[{Log10[\[Omega]/("vF"q )/.\[Epsilon]Dict["params"]],Log10[q/(2"qF")/.\[Epsilon]Dict["params"]]}/.qpm[[i]]/.\[Epsilon]Dict["params"],{i,2}];
prefactor=( 1/(v\[Chi]^2 "ne") 1/(("\[HBar]")^2 (2 \[Pi])^3) 1/\[Xi]\[Theta] (\[Kappa]^2 ("e")^2)/(4 \[Pi] "\[Epsilon]0") 2/(1-E^(-"\[Beta]" "\[HBar]" \[Omega])))/.\[Epsilon]Dict["params"];

d\[Sigma]dERd\[CapitalOmega] =prefactor Table[10^\[Epsilon]Dict["\[Epsilon]f"][loguzs[[i,1]],loguzs[[i,2]]],{i,2}];

<|"d\[Sigma]dERd\[CapitalOmega]"->d\[Sigma]dERd\[CapitalOmega],"\[Omega]"->\[Omega],"\[Omega]max"->\[Omega]max,"cos\[Theta]"->cos\[Theta],"\[Xi]pm"->\[Xi]\[Theta],"qpm"->qpm,"loguzs"->loguzs,"prefactor"->prefactor|>
]


Get\[Xi]\[Omega]andcos\[Theta][v\[Chi]_,\[Epsilon]Dict_,NMax_]:=Module[{nits=0,\[Xi]s,M\[Xi],d\[Sigma]dERd\[CapitalOmega],\[Omega],Accepted=False,Errors=0},
While[!Accepted,
nits+=1;
\[Xi]s={Random[],Random[]};
M\[Xi]=NMax Random[];
(*\[Sigma]s=Total[d\[Sigma]dERd\[CapitalOmega]eofv\[Chi]and\[Epsilon][\[Xi]s[[1]],\[Xi]s[[2]],v\[Chi],1,\[Epsilon]Dict][["d\[Sigma]dERd\[CapitalOmega]"]]];*)
{d\[Sigma]dERd\[CapitalOmega],\[Omega]}={"d\[Sigma]dERd\[CapitalOmega]","\[Omega]"}/.d\[Sigma]dERd\[CapitalOmega]eofv\[Chi]and\[Epsilon][\[Xi]s[[1]],\[Xi]s[[2]],v\[Chi],1,\[Epsilon]Dict];(*Set \[Kappa] to 1, overall scaling is irrelevant*)
d\[Sigma]dERd\[CapitalOmega]=Total[d\[Sigma]dERd\[CapitalOmega]];
Accepted=M\[Xi]<= d\[Sigma]dERd\[CapitalOmega];
Errors=If[d\[Sigma]dERd\[CapitalOmega]-NMax> 0,Break[];<|"\[Xi]\[Omega]"->\[Xi]s[[1]],"cos\[Theta]"->\[Xi]s[[2]],"M\[Xi]"->M\[Xi],"\[Sigma]s"->d\[Sigma]dERd\[CapitalOmega]|>,0];
];
<|"\[Xi]\[Omega]"->\[Xi]s[[1]],"\[Omega]"->\[Omega],"cos\[Theta]"->\[Xi]s[[2]],"M\[Xi]"->M\[Xi],"d\[Sigma]dERd\[CapitalOmega]"->d\[Sigma]dERd\[CapitalOmega],"Max"->NMax,"Errors"->Errors,"nits"->nits|>
]


(*GetNMaxofd\[Sigma]dERd\[CapitalOmega][\[Epsilon]Dict_,v\[Chi]_,buffer_:1.1 ]:=Module[{},
buffer NMaximize[Total[d\[Sigma]dERd\[CapitalOmega]eofv\[Chi]and\[Epsilon][\[Xi]\[Omega],cos\[Theta],v\[Chi],1,\[Epsilon]Dict][["d\[Sigma]dERd\[CapitalOmega]"]]],{\[Xi]\[Omega],cos\[Theta]}\[Element]Rectangle[{\[Epsilon]Dict["\[Delta]"]^(1/4),\[Epsilon]InterpFeOsc2["\[Delta]"]^(1/4)},{1-\[Epsilon]Dict["\[Delta]"]^(1/4),1-\[Epsilon]Dict["\[Delta]"]^(1/4)}]][[1]] 
]*)
GetNMaxofd\[Sigma]dERd\[CapitalOmega][\[Epsilon]Dict_,buffer_:1.1 ]:=Module[{v\[Chi]max},
v\[Chi]max=10^-3 "c"/.Constants`SIConstRepl;
buffer NMaximize[Total[d\[Sigma]dERd\[CapitalOmega]eofv\[Chi]and\[Epsilon][\[Xi]\[Omega],cos\[Theta],v\[Chi]max,1,\[Epsilon]Dict][["d\[Sigma]dERd\[CapitalOmega]"]]],{\[Xi]\[Omega],cos\[Theta]}\[Element]Rectangle[{\[Epsilon]Dict["\[Delta]"]^(1/4),\[Epsilon]Dict["\[Delta]"]^(1/4)},{1-\[Epsilon]Dict["\[Delta]"]^(1/4),1-\[Epsilon]Dict["\[Delta]"]^(1/4)}]][[1]] 
]


(* ::Chapter:: *)
(*End*)


End[];


EndPackage[];
