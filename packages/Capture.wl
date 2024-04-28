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
d\[Sigma]dERd\[CapitalOmega]eofv\[Chi]and\[Epsilon]::usage = "define d\[Sigma]dERd\[CapitalOmega] for electronic scattering";
Get\[Xi]\[Omega]andcos\[Theta]::usage = "Sample the distribution of d\[Sigma]dERd\[CapitalOmega] for \[Xi]\[Omega] and cos\[Theta]";
GetNMaxofd\[Sigma]dERd\[CapitalOmega]::usage = "Numerically solve for the maximum of d\[Sigma]dERd\[CapitalOmega] slightly away from the boundaries, and with a small rescaling buffer (for use with the rejection method)";

d\[Sigma]d\[CapitalOmega]Nuc0::usage = "define d\[Sigma]d\[CapitalOmega] for 0 temperature nuclear scattering";
GetNMaxofd\[Sigma]d\[CapitalOmega]Nuc0::usage = "Numerically solve for the maximum of d\[Sigma]dERd\[CapitalOmega] slightly away from the boundaries, and with a small rescaling buffer (for use with the rejection method)";
Get\[Omega]andcos\[Theta]Nuc0::usage = "Sample the distribution of d\[Sigma]dERd\[CapitalOmega] for \[Xi]\[Omega] and cos\[Theta]";

ProcessDataAssemblerElectronic::usage = "take data of a electronic material / process and format it for the MC. Also does interpolation of \[Epsilon].";
Interpolated\[Sigma]dERd\[Sigma]Max::usage = "Interpolate the maximum of d\[Sigma]dERd\[CapitalOmega] as a function of v\[Chi]";

Getv\[Chi]\[Infinity]::usage = "Get the initial v\[Chi] at \[Infinity] by sampling a MB and a 2 sphere";
Sample2Sphere::usage = "Point in R3, normalized gives random point on sphere or radius 1";
Getv\[Chi]atEarth::usage = "sample MB and sphere at \[Infinity] to find trajectories that reach earth with b_E less than earth's radius";

testMCprivaterun::usage = "";
RunMC::usage = "Run the Monte-Carlo for one particle given a set of initial conditions";


(* ::Chapter:: *)
(*Public*)


Begin["`Private`"];


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


Vgrav[r_]:=Piecewise[{{-(("G" "ME")/r),r>"rE"/. Constants`EarthRepl},{-(("ME" "G")/(2("rE")^3))( 3("rE")^2-r^2),r<"rE"/. Constants`EarthRepl}},-(("G" "ME")/("rE"/. Constants`EarthRepl))]/.Constants`SIConstRepl/.Constants`EarthRepl


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
(*Electronic Scattering *)


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
GetNMaxofd\[Sigma]dERd\[CapitalOmega][\[Epsilon]Dict_,v\[Chi]_,buffer_:1.1 ]:=Module[{},
(*v\[Chi]max=10^-3 "c"/.Constants`SIConstRepl;*)
buffer NMaximize[Total[d\[Sigma]dERd\[CapitalOmega]eofv\[Chi]and\[Epsilon][\[Xi]\[Omega],cos\[Theta],v\[Chi],1,\[Epsilon]Dict][["d\[Sigma]dERd\[CapitalOmega]"]]],{\[Xi]\[Omega],cos\[Theta]}\[Element]Rectangle[{\[Epsilon]Dict["\[Delta]"]^(1/4),\[Epsilon]Dict["\[Delta]"]^(1/4)},{1-\[Epsilon]Dict["\[Delta]"]^(1/4),1-\[Epsilon]Dict["\[Delta]"]^(1/4)}]][[1]] 
]


(* ::Subsection:: *)
(*Nuclear Scattering*)


d\[Sigma]d\[CapitalOmega]Nuc0[cos\[Theta]_,v\[Chi]_,\[Kappa]_,m\[Chi]_,coeffs_,params_]:=Module[{mN,\[Omega],qpm,\[Xi]\[Theta],Sofq\[Omega],d\[Sigma]d\[CapitalOmega]},
(* Get d\[Sigma]/(d\[CapitalOmega]) for nuclear scattering at 0 temperature. At 0 temperature there is a delta function which gives \[Omega] in terms of q (and hence cos\[Theta] through q_0)
  cos\[Theta] \[Epsilon] (0,1) - [] negative disallowed since gas is assumed 0 T: hence no energy can be removed from the nuclei.*)
mN ="mN"/. coeffs;
\[Omega]=(2 cos\[Theta]^2 mN m\[Chi]^2 v\[Chi]^2)/((mN+m\[Chi])^2 "\[HBar]")/.params; (*define \[Omega] in terms of cos\[Theta]*)
qpm={q->(m\[Chi] v\[Chi])/("\[HBar]") (cos\[Theta]-Sqrt[cos\[Theta]^2-(2\[Omega] "\[HBar]")/(m\[Chi] v\[Chi]^2)]),q->(m\[Chi] v\[Chi])/("\[HBar]") (cos\[Theta]+Sqrt[cos\[Theta]^2-(2\[Omega] "\[HBar]")/(m\[Chi] v\[Chi]^2)])}/.params;(*solutions for q to \[Omega]-\[Omega](q,cos\[Theta])=0 in the delta function that arises in the q integral. These must be summed over*)
\[Xi]\[Theta]=Sqrt[cos\[Theta]^2-(2 "\[HBar]" \[Omega])/(m\[Chi] v\[Chi]^2)]/.params;(*the absolute value of the jacobian of the delta function, evaluated on the roots of q in the argument, for each of the above solutions*)

d\[Sigma]d\[CapitalOmega] =Total@Table[1/(v\[Chi]^2 )  q^2/(("\[HBar]")^2 (2 \[Pi])^3) 1/\[Xi]\[Theta] ((\[Kappa] ("e")^2)/(4 \[Pi] "\[Epsilon]0" q^2))^2 FormFactors`Zeff[q,coeffs]/.qpm[[i]]/.params,{i,2}];

<|"d\[Sigma]d\[CapitalOmega]"->d\[Sigma]d\[CapitalOmega],"\[Omega]"->\[Omega],"cos\[Theta]"->cos\[Theta],"\[Xi]pm"->\[Xi]\[Theta],"qpm"->qpm|>
]


GetNMaxofd\[Sigma]d\[CapitalOmega]Nuc0[v\[Chi]_,m\[Chi]_,coeffs_,params_,buffer_:1.1 ,\[Delta]_:10^-3]:=Module[{},
(*v\[Chi]max=10^-3 "c"/.Constants`SIConstRepl;*)
buffer NMaximize[Total@d\[Sigma]d\[CapitalOmega]Nuc0[cos\[Theta],v\[Chi],1,m\[Chi],coeffs,params][["d\[Sigma]d\[CapitalOmega]"]],cos\[Theta]\[Element]Interval[{\[Delta],1-\[Delta]}]][[1]]
]


Get\[Omega]andcos\[Theta]Nuc0[v\[Chi]_,m\[Chi]_,coeffs_,params_,NMax_]:=Module[{nits=0,cos\[Theta],M\[Xi],d\[Sigma]d\[CapitalOmega],\[Omega],Accepted=False,Errors=0},
While[!Accepted,
nits+=1;
cos\[Theta]=Random[];(*assumes 0 temp so cos\[Theta] > 0*)
M\[Xi]=NMax Random[];
(*\[Sigma]s=Total[d\[Sigma]dERd\[CapitalOmega]eofv\[Chi]and\[Epsilon][\[Xi]s[[1]],\[Xi]s[[2]],v\[Chi],1,\[Epsilon]Dict][["d\[Sigma]dERd\[CapitalOmega]"]]];*)
{d\[Sigma]d\[CapitalOmega],\[Omega]}={"d\[Sigma]d\[CapitalOmega]","\[Omega]"}/.d\[Sigma]d\[CapitalOmega]Nuc0[cos\[Theta],v\[Chi],1,m\[Chi],coeffs,params];(*Set \[Kappa] to 1, overall scaling is irrelevant*)
(*d\[Sigma]d\[CapitalOmega]=Total[d\[Sigma]d\[CapitalOmega]];*)
Accepted=M\[Xi]<= d\[Sigma]d\[CapitalOmega];
Errors=If[d\[Sigma]d\[CapitalOmega]-NMax> 0,Break[];True,0];
];
<|"\[Omega]"->\[Omega],"cos\[Theta]"->cos\[Theta],"M\[Xi]"->M\[Xi],"d\[Sigma]dERd\[CapitalOmega]"->d\[Sigma]d\[CapitalOmega],"Max"->NMax,"Errors"->Errors,"nits"->nits|>
]


(* ::Subsection:: *)
(*Initial Conditions*)


Getv\[Chi]\[Infinity][m\[Chi]_,\[Beta]D_,maxrat_:Sqrt[50]]:=Module[{fMB,NMax,v\[Chi]Max,v\[Chi]MP,nits=0,v\[Chi],P,M\[Xi],Accepted=False},
v\[Chi]MP=Sqrt[2/(m\[Chi] \[Beta]D)];(*most probable velocity in MB*)
v\[Chi]Max = maxrat v\[Chi]MP; (*maximum a simple rescaling of most probable, since higher velocities are exponentially suppressed*)
fMB=(#^2 E^(-\[Beta]D m\[Chi]/2 #^2))&;
NMax=2/(E m\[Chi] \[Beta]D);(*value of fMB at v\[Chi] MP*)
While[!Accepted && nits<30,
nits+=1;
v\[Chi]=Random[] v\[Chi]Max;
(*Print[v\[Chi]];*)
M\[Xi]=Random[]NMax;
P=fMB[v\[Chi]];(*fMB of getting velocity v\[Chi]*)
Accepted=M\[Xi]<= P;
];
<|"v\[Chi]"->v\[Chi],"v\[Chi]Max"->v\[Chi]Max,"M\[Xi]"->M\[Xi],"fMB"->P,"Max"->NMax,"nits"->nits|>
]


Sample2Sphere[]:=Module[{x},
(*Point in R3, normalized gives random point on sphere or radius 1*)
x=Table[1 - 2Random[],{i,3}];
x/Sqrt[x . x]
]


Getv\[Chi]atEarth[m\[Chi]_,\[Beta]D_,r\[Infinity]_:10^1 "rE"/.Constants`EarthRepl,V_:Capture`Vgrav]:=Module[{rE,Accepted=False,nits=0,\[Chi]speed\[Infinity],v\[Chi]\[Infinity],vperp\[Infinity],EK\[Infinity],b\[Infinity],bEarth,vrEarth,succeeded=False,v\[Chi],\[Chi]speedEarth},

rE="rE"/.Constants`EarthRepl;

While[!Accepted &&nits<10000,
nits++;
\[Chi]speed\[Infinity]="v\[Chi]"/.Getv\[Chi]\[Infinity][m\[Chi],\[Beta]D];
v\[Chi]\[Infinity]= \[Chi]speed\[Infinity] Sample2Sphere[];

If[v\[Chi]\[Infinity][[1]]>0,Continue[]];(*neglect anything pointed away*)

vperp\[Infinity]=Sqrt[v\[Chi]\[Infinity] . DiagonalMatrix[{0,1,1}] . v\[Chi]\[Infinity]];
EK\[Infinity]=1/2 m\[Chi] \[Chi]speed\[Infinity]^2;
b\[Infinity]=r\[Infinity] vperp\[Infinity]/\[Chi]speed\[Infinity];
bEarth = b\[Infinity]/Sqrt[1 + (m\[Chi](V[r\[Infinity]]-V[rE]))/EK\[Infinity]];
If[bEarth<rE,succeeded=True;Break[]]
];
vrEarth = Sqrt[v\[Chi]\[Infinity][[1]]^2+2(V[r\[Infinity]]-V[rE])];
Print[v\[Chi]\[Infinity][[1]]^2];
Print[2(V[r\[Infinity]]-V[rE])];
v\[Chi]={-vrEarth,v\[Chi]\[Infinity][[2]],v\[Chi]\[Infinity][[3]]};
\[Chi]speedEarth = Sqrt[v\[Chi] . v\[Chi]];
<|"v\[Chi]\[Infinity]"->v\[Chi]\[Infinity],"v\[Chi]"->v\[Chi],"\[Chi]speedEarth"->\[Chi]speedEarth,"m\[Chi]"->m\[Chi],"\[Beta]D"->\[Beta]D,"\[Chi]speed\[Infinity]"->\[Chi]speed\[Infinity],"bEarth"->bEarth,"nits"->nits,"succeeded"->succeeded|>
]


(* ::Subsection:: *)
(*MC PreProcessing*)


(*ProcessDataAssemblerElectronic[m\[Chi]_,params_,d\[Sigma]dERd\[CapitalOmega]sampler_,\[Lambda]Dict_,ProcNum_:1,uselog_:True,m_:200,n_:200,\[Delta]_:10^-7]:= Module[{\[Epsilon]Dict,d\[Sigma]dERd\[CapitalOmega]Max,d\[Sigma]dERd\[CapitalOmega]sample,\[Sigma]},
\[Epsilon]Dict= Interpolate\[Epsilon]euztrunced[m\[Chi],params,uselog,m,n,\[Delta]];
d\[Sigma]dERd\[CapitalOmega]Max = GetNMaxofd\[Sigma]dERd\[CapitalOmega][\[Epsilon]Dict,#]&;
d\[Sigma]dERd\[CapitalOmega]sample = Get\[Xi]\[Omega]andcos\[Theta][#,\[Epsilon]Dict,d\[Sigma]dERd\[CapitalOmega]Max[#]]&;
\[Sigma]=1/("ne"/.\[Lambda]Dict["fitparams"][[1]]) 10^("f"/.\[Lambda]Dict)[Log10[m\[Chi]],Log10[#]]&;
<|"ProcNum"->ProcNum,"\[Epsilon]Dict"->\[Epsilon]Dict,"d\[Sigma]dERd\[CapitalOmega]Max"->d\[Sigma]dERd\[CapitalOmega]Max,"d\[Sigma]dERd\[CapitalOmega]sample"->d\[Sigma]dERd\[CapitalOmega]sample,"\[Sigma]"->\[Sigma],"\[Lambda]"->\[Lambda]Dict["f"]|>
]*)
ProcessDataAssemblerElectronic[m\[Chi]_,params_,d\[Sigma]dERd\[CapitalOmega]sampler_,\[Lambda]Dict_,Save_:False,fname_:"procdata",ProcNum_:1,uselog_:True,m_:200,n_:200,\[Delta]_:10^-7]:= Module[{\[Epsilon]Dict,d\[Sigma]dERd\[CapitalOmega]Maxloglog,d\[Sigma]dERd\[CapitalOmega]Max,d\[Sigma]dERd\[CapitalOmega]sample,\[Sigma],outdict},
\[Epsilon]Dict= Interpolate\[Epsilon]euztrunced[m\[Chi],params,uselog,m,n,\[Delta]];
d\[Sigma]dERd\[CapitalOmega]Maxloglog = "Maxf"/.Interpolated\[Sigma]dERd\[Sigma]Max[GetNMaxofd\[Sigma]dERd\[CapitalOmega][\[Epsilon]Dict,#]&];
d\[Sigma]dERd\[CapitalOmega]Max=(10^d\[Sigma]dERd\[CapitalOmega]Maxloglog[Log10[#]])&;
d\[Sigma]dERd\[CapitalOmega]sample = Get\[Xi]\[Omega]andcos\[Theta][#,\[Epsilon]Dict,d\[Sigma]dERd\[CapitalOmega]Max[#]]&;
\[Sigma]=1/("ne"/.\[Lambda]Dict["fitparams"][[1]]) 10^("f"/.\[Lambda]Dict)[Log10[m\[Chi]],Log10[#]]&;
outdict=<|"ProcNum"->ProcNum,"\[Epsilon]Dict"->\[Epsilon]Dict,"d\[Sigma]dERd\[CapitalOmega]Max"->d\[Sigma]dERd\[CapitalOmega]Max,"d\[Sigma]dERd\[CapitalOmega]sample"->d\[Sigma]dERd\[CapitalOmega]sample,"\[Sigma]"->\[Sigma],"\[Lambda]"->\[Lambda]Dict["f"]|>;
If[Save,Utilities`SaveIt[NotebookDirectory[]<>fname,outdict]];(*Note that we have to save the object internally in the package to use it in the package. It is saved as private so can only be accessed from in here.*)
outdict 
]


Interpolated\[Sigma]dERd\[Sigma]Max[d\[Sigma]dERd\[CapitalOmega]Maxf_,m_:20,uselog_:True]:=Module[{v\[Chi]Max,v\[Chi]Min,v\[Chi]s,MaxInterTable,MaxInterf},
v\[Chi]Max = 2 10^-1 "c" /.Constants`SIConstRepl;
v\[Chi]Min = 1/10 Constants`vescape;
v\[Chi]s = N@If[uselog,10^Subdivide[Log10[v\[Chi]Min],Log10[v\[Chi]Max],m],Subdivide[v\[Chi]Min,v\[Chi]Max,m]];
MaxInterTable = Table[{v\[Chi]s[[i]],d\[Sigma]dERd\[CapitalOmega]Maxf[v\[Chi]s[[i]]]},{i,m}];
MaxInterf=If[uselog,Interpolation[Log10[MaxInterTable],InterpolationOrder->4],Interpolation[MaxInterTable,InterpolationOrder->4]];
<|"Maxf"->MaxInterf,"MaxTable"->MaxInterTable,"v\[Chi]s"->v\[Chi]s,"v\[Chi]range"->{v\[Chi]Min,v\[Chi]Max}|>
]


(* ::Subsection:: *)
(*Run MC*)


testMCprivaterun[params_,\[Lambda]_,fileloc_]:=Module[{FeDataet},
FeDataet=ProcessDataAssemblerElectronic[10 "m"/.Constants`SIConstRepl,params,Get\[Xi]\[Omega]andcos\[Theta],\[Lambda]];
Utilities`SaveIt[fileloc,FeDataet];
(*FeDataet = Utilities`ReadIt[fileloc];*)
("d\[Sigma]dERd\[CapitalOmega]sample"/.FeDataet)[10^5];
(*This is getting convoluted, we might want to just dump all this into a notebook so we don't have to deal with all the scope crap.*)
RunMC[10 "m"/.Constants`SIConstRepl,10^-6.5,("v\[Chi]"/.Getv\[Chi]atEarth[ 10"m"/.Constants`SIConstRepl,(0.001"JpereV"/.Constants`SIConstRepl)^-1]),{"\[Lambda]"["f"]/.FeDataet},{"\[Sigma]"/.FeDataet},{"d\[Sigma]dERd\[CapitalOmega]Max"/.FeDataet},{"d\[Sigma]dERd\[CapitalOmega]sample"/.FeDataet}]
]


RunMC[m\[Chi]_,\[Kappa]_,v\[Chi]0_,\[Lambda]s_,\[Sigma]s_,d\[Sigma]dERd\[CapitalOmega]Maxes_,d\[Sigma]dERd\[CapitalOmega]sample_,V_:Vgrav]:=Module[{v\[Chi],r,V,nits,ncaptured,nescaped,ProcNum,l,reached,ER,cos\[Theta],\[Phi],q},
(*
Psuedo Code
initialize cross-sections, interaction lengths, and cross-sections for each process
Get initial conditions: v\[Chi] (with perp and radial components) at earth's radius
While in the Earth
	Get interaction length in the Earth
	Update radius according to interaction length and velocity
	Update velocity according to potential change 
	select process - omit for now
	sample to find \[Omega] and cos\[Theta] 
	update velocity
*)

{v\[Chi],r} = {v\[Chi]0,"rE"/.EarthRepl};

nits=0;
ncaptured=0;
nescaped=0;
While[nits<2,
(*nits<30,*)
nits++;
Print["N iterations: ",nits];
ProcNum=\[Sigma]selector[Table[\[Sigma]s[[i]][Sqrt[v\[Chi] . v\[Chi]]],{i,Length[\[Sigma]s]}]];
Print["Process number: ",ProcNum];
(*l = Getl[m\[Chi],Sqrt[v\[Chi].v\[Chi]],\[Kappa],\[Lambda]s[[ProcNum]]["f"]];*)
l = Getl[m\[Chi],Sqrt[v\[Chi] . v\[Chi]],\[Kappa],\[Lambda]s[[ProcNum]]];
Print["l is ",l];
{v\[Chi],r,reached}={"v\[Chi]","r","reached?"}/.Updatev[v\[Chi],m\[Chi],r,l,V];
Print["v\[Chi] is ",v\[Chi]];
Print["|v\[Chi]| is ",Sqrt[v\[Chi] . v\[Chi]]];
Print["r is ",r];
Print[reached];
If[r>"rE"/.Constants`EarthRepl  ,nescaped++;Break[]];
{ER,cos\[Theta]}={"\[HBar]""\[Omega]"/.Constants`SIConstRepl,"cos\[Theta]"}/.d\[Sigma]dERd\[CapitalOmega]sample[[ProcNum]][Sqrt[v\[Chi] . v\[Chi]]];
\[Phi]=2 \[Pi] Random[];
q=Sqrt[2 m\[Chi] ER ]/("\[HBar]")/.Constants`SIConstRepl;
v\[Chi] = Rotatev\[Chi][v\[Chi],cos\[Theta],m\[Chi],\[Phi],q];
If[Sqrt[v\[Chi] . v\[Chi]]<"vesc"/.Constants`EarthRepl,ncaptured++;Print["Captured with: \!\(\*SubscriptBox[\(v\), \(\[Chi]\)]\)=",Sqrt[v\[Chi] . v\[Chi]]];Break[]];
(*Print["v\[Chi] after scatter is ",v\[Chi]];*)
];
<|"ER"->ER,"v\[Chi]"->v\[Chi],"|v\[Chi]|"->Sqrt[v\[Chi] . v\[Chi]],"\[Phi]"->\[Phi],"l"->l,"cos\[Theta]"->cos\[Theta],"q"->q,"r"->r,"reached?"->reached,"\[Kappa]"->N@\[Kappa],"nits"->nits,"ncaptured"->ncaptured,"nescaped"->nescaped|>
]


(* ::Chapter:: *)
(*End*)


End[];


EndPackage[];
