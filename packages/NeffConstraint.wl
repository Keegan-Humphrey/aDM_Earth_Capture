(* ::Package:: *)

BeginPackage["NeffConstraint`"];


(* ::Text:: *)
(*The purpose of this package is to define modules and constants used in the computation of the Subscript[\[CapitalDelta]N, eff] constraint on \[Kappa]*)


(* ::Subsubsection:: *)
(*Public Declarations*)


CosmoParams::usage = "Produce a replacement table (string to N) for parameters used in the computation of the \!\(\*SubscriptBox[\(\[CapitalDelta]N\), \(eff\)]\) constraint on \[Kappa]";
\[CapitalLambda]CDMrepl::usage = "Density parameters for \[CapitalLambda]CDM cosmology";

H::usage = "Hubble as a function of temperature";
aH::usage = "\!\(\*OverscriptBox[\(a\), \(.\)]\) as a function of temperature";

ne::usage = "electron number density as a function of temperature";
Z1::usage = "classical single particle partition function";
\[CapitalLambda]::usage = "Debeye energy of SM electrons as a function of temperature";

nD0::usage = "DM number density today as a function of temperature";
nDR::usage = "DM number density at recombination as a function of temperature";

NumIntegrand::usage = "integrand used in ESM Integral";
NonNumIntegrand::usage = "integrand used in ESM Integral in symbolic form";
ESMIntegral::usage = "integral to compute energy transfer / interaction rates w electrons as a function of temperature";

Ettotheld\[CapitalGamma]dEt::usage = "\!\(\*SuperscriptBox[\(E\), \(l\)]\)\!\(\*FractionBox[\(d\[CapitalGamma]\), SubscriptBox[\(dE\), \(t\)]]\)";

InterpolateEtransRate::usage = "Interpolate over temperature to find energy transfer rate as a funtion of mass";
InterpolateTotEtrans::usage = "Interpolate over mass to find energy transfered to DS";
Interpolatey::usage = "Find the CMB distortion comptonization parameter";

Compute\[CapitalDelta]NeffConstraint::usage = "Compute \[CapitalDelta]Neff constraint and produce contour plot of \[Kappa] bound in \!\(\*SubscriptBox[\(m\), \(D\)]\),\!\(\*SubscriptBox[\(f\), \(D\)]\) plane";


(* ::Section:: *)
(*Private*)


Begin["`Private`"];


CosmoParams[]:= Module[{Td,ad,Tdt,T0,aR,TR,aB,TB,h,eVm,c,\[Rho]critbaumanncm,\[Rho]critbaumannm,\[Rho]crit,H0tinv,H0linv,H0inv,H0,n0,nd,nr,mSMtot,me,\[Alpha],g,gstar,pth,kF,EF,EFd,Degen,\[Rho]DM0},

(*Cosmo Constants - temperatures and scale factors*)
Td = 5 10^5; (*[eV] annihilation temperature*)
ad =(2 10^9)^-1; (*[] scale factor (~(1/z)) at annihilation*)
Tdt = 2.5 10^-4; (*[eV] Temperature of the Universe today*)
T0 = 2.5 10^-4; (*[eV] Temperature of the Universe today*)
aR = (1100)^-1;(*[] scale factor ~ (1/z) at recombination*)
TR = (Tdt/aR); (*[eV] Temperature at recombination*)
aB = 10^-7; (*[] scale factor when brehmstrahlung becomes inefficient*)
TB = T0/aB; (*[eV] temperature when brehmstrahlung becomes inefficient*)

(*Conversion*)
eVm = N[200 10^-15 10^6]; (*[eV m] = 1 follows from 1 = 200 MeV fm*)
c = 3 10^8; (*[m s^-1]*)

(*Critical Density*)
h= 0.67; (*baumann 1.2.64*)
\[Rho]critbaumanncm = 1.1 10^-5 h^2; (*[protons cm^-3]*) 
\[Rho]critbaumannm = \[Rho]critbaumanncm  (10^6); (*[protons m^-3]*)
\[Rho]crit = \[Rho]critbaumannm eVm^3; (*[protons eV^3]*)

(*Hubble*)
H0tinv = 4.4 10^17; (*s*)
H0linv = c H0tinv; (*m*)
H0inv = eVm^-1 H0linv; (*eV^-1*)
H0= (H0inv^-1); (*eV*) 

(*Number density*)
n0 = "\[CapitalOmega]b" \[Rho]crit /.\[CapitalLambda]CDMrepl;  (*[eV^3] proton / electron number density today*)
nd = n0 (Td/Tdt)^3; (*[eV^3] number density at decoupling*)
nr = n0 (TR/Tdt)^3; (*[eV^3] number density at recombination*)

(*SM parameters*)
mSMtot= me + 2000 me;
me = 5 10^5; (*[eV] electron mass*)
\[Alpha]=1/137;
g = 4; (*[] number of dofs in SM electron*)
gstar = 3.38;(*[] gstar after decoupling*)

(*Thermal Quantities*)
pth =( (me Tdt)/(2 \[Pi]))^(1/2); (*[eV] thermal de-broglie momentum today *)
kF = N[(3 \[Pi]^2 n0)^(1/3)]; (*[eV] fermi momentum today*)
EF = kF^2/(2 me);
EFd = EF Td/Tdt; (*[eV]fermi energy at decoupling*)
Degen = EF/Tdt; (*Degeneracy of gas today*)

(*DM energy density*)
\[Rho]DM0= "\[CapitalOmega]c"/"\[CapitalOmega]b" n0 mSMtot/.\[CapitalLambda]CDMrepl;

{"Td"->Td,"ad"->ad,"Tdt"->Tdt,"T0"->T0,"aR"->aR,"TR"->TR,"aB"->aB,"TB"->TB,"eVm"->eVm,"c"->c,"h"->h,"\[Rho]crit"->\[Rho]crit,"H0"->H0,"n0"->n0,"nd"->nd,"nr"->nr,"me"->me,"\[Alpha]"->\[Alpha],"g"->g,"gstar"->gstar,"pth"->pth,"kF"->kF,"EF"->EF,"EFd"->EFd,"Degen"->Degen,"\[Rho]DM0"->\[Rho]DM0}
]


(* ::Text:: *)
(*Define hubble and scale speed*)


(*Density Parameters for Hubble*)
(*\[CapitalLambda]CDMrepl= {"\[CapitalOmega]k" ->0.01,"\[CapitalOmega]r"->9.4 10^-5, "\[CapitalOmega]m"->0.32,"\[CapitalOmega]\[CapitalLambda]"->0.68,"\[CapitalOmega]b"->0.05,"\[CapitalOmega]c"->0.27}; (*[] \[CapitalLambda]CDM density parameters*)*)
\[CapitalLambda]CDMrepl= {"\[CapitalOmega]k" ->0.00,"\[CapitalOmega]r"->9.4 10^-5, "\[CapitalOmega]m"->0.32,"\[CapitalOmega]\[CapitalLambda]"->0.68,"\[CapitalOmega]b"->0.05,"\[CapitalOmega]c"->0.27}; (*[] \[CapitalLambda]CDM density parameters*)
H[\[Beta]_] := "H0" ("\[CapitalOmega]r" "a"^-4+"\[CapitalOmega]m" "a"^-3+"\[CapitalOmega]k" "a"^-2 + "\[CapitalOmega]\[CapitalLambda]")^(1/2)/.\[CapitalLambda]CDMrepl/."a"->("Tdt" \[Beta])
aH[\[Beta]_] :="a" "H0" ("\[CapitalOmega]r" "a"^-4+"\[CapitalOmega]m" "a"^-3+"\[CapitalOmega]k" "a"^-2 + "\[CapitalOmega]\[CapitalLambda]")^(1/2)/.\[CapitalLambda]CDMrepl/."a"->("Tdt" \[Beta])


(* ::Text:: *)
(*Number density today*)


ne[\[Beta]_]:= "n0" (1/(\[Beta] "Tdt"))^3 (*[eV^3] number density as a function of temperature*)


(* ::Text:: *)
(*Thermal Quantities*)


Z1[\[Beta]_] := ("g" "pth"^3)/"n0" (1/("Tdt" \[Beta]))^(-3/2)(*[] single particle classical partition function - normalization of distribution Subscript[n, e]*)


(* ::Text:: *)
(*more number density*)


(*mSMtot= me + 2000 me*)
(*\[Rho]DM0= "\[CapitalOmega]c"/"\[CapitalOmega]b" "n0" "mSMtot"/.\[CapitalLambda]CDMrepl*)
nD0[fD_,mDMtot_] :=  "\[Rho]DM0"/mDMtot fD
nDR[fD_,mDMtot_] := 1/"aR"^3 nD0[fD,mDMtot]


(* ::Text:: *)
(*Debeye energy - Assuming electron mass*)


\[CapitalLambda][\[Beta]_,m_]:= 4 \[Pi] "\[Alpha]" (2\[Pi] "n0")/(m "Tdt" ) (1/(\[Beta] "Tdt"))^2


(* ::Text:: *)
(*Determine interaction rate as a function of temperature *)


(*NumIntegrand[mDt_,mt_,\[Kappa]_,\[Xi]\[Beta]_?NumberQ,\[Xi]E_?NumberQ,params_,l_:1]:= Module[{Integrand,E1ms,E1},
(*l - == 1 => energy transfer rate integrand, == 0 => interaction rate (over Subscript[n, D]) integrand*)
 E1ms=Integrate[Ett^l/(Ett+"\[CapitalLambda]")^2 (("ms" + "ESMt")^2+("ms"+(Ett-"ESMt"))^2),Ett];
	E1=Simplify[E1ms/."ms"->("m"^2+"mD"^2)/(2 "mD")];
	Integrand =Simplify[Simplify[(E1/.Ett->(2 "mD")/(("ESMt"+"mD")^2/("ESMt"^2-"m"^2)-1))-(E1/.Ett->0)]];
2/\[Pi] ("\[Alpha]"^2 \[Kappa]^2)/mDt (E^(-"\[Beta]" "ESMt" + "\[Beta]" "\[Mu]")/aH["\[Beta]"]^l Integrand /."\[Mu]"->"m" - 1/"\[Beta]" Log["Qt"]/.{"\[CapitalLambda]"->\[CapitalLambda]["\[Beta]",mt],"Qt"->Z1["\[Beta]"]}/."ESMt"->1/"\[Beta]" \[Xi]E + "m"/."\[Beta]"->\[Xi]\[Beta]/"m"/.{"m"->mt,"mD"->mDt})/.params
]*)
NumIntegrand[mDt_,mt_,\[Kappa]_,\[Xi]\[Beta]_?NumberQ,\[Xi]E_?NumberQ,params_,l_:1]:= Module[{P\[CapitalMu],prefactor,IntegrandEt,boundsE,boundsEt,thermalrepl,\[CapitalLambda]repl,msrepl,\[Xi]Erepl,\[Xi]\[Beta]repl,IntegrandESM},
(*l == 1 => energy transfer rate integrand (per particle per efold)
l == 0 => interaction rate (over Subscript[n, D]) integrand*)
P\[CapitalMu]=("ms"+"ESM")^2+("ms"+(Et-"ESM"))^2;
prefactor = ("n0")/(("T0" "\[Beta]")^3 aH["\[Beta]"]^l) 1/("Z1") ((4 \[Pi] "\[Alpha]")^2 \[Kappa]^2)/("mD" (2 \[Pi])^3) ("g")/4 ;(* prefactor of eqn (629) of notes over nD*)
IntegrandEt=E^(-"\[Beta]"("ESM"-"m")) P\[CapitalMu]/(Et + "\[CapitalLambda]")^2 Et^l;
boundsEt={0,(2 "mD")/(("ESM"+"mD")^2/(("ESM")^2-("m")^2)-1) };
thermalrepl={"\[CapitalLambda]"->"\[CapitalLambda]0" (1/("T0" "\[Beta]"))^2,"Z1"->"g" ((2 \[Pi] "\[Beta]")/( "m"))^(-3/2)};
\[CapitalLambda]repl={"\[CapitalLambda]0"->4\[Pi] "\[Alpha]" (2 \[Pi] "n0")/("m" "T0")};
msrepl = {"ms"->(("m")^2+("mD")^2)/(2 "mD")};
\[Xi]Erepl={"ESM"->\[Xi]E /"\[Beta]" + "m"};
\[Xi]\[Beta]repl="\[Beta]"->\[Xi]\[Beta]/"m";

IntegrandESM=Integrate[IntegrandEt,Et];
prefactor((IntegrandESM/.Et->boundsEt[[2]])-(IntegrandESM/.Et->boundsEt[[1]]))/.msrepl/.thermalrepl/.\[CapitalLambda]repl/.\[Xi]Erepl/.\[Xi]\[Beta]repl/.{"m"->mt,"mD"->mDt}/.params
]


(*NonNumIntegrand[mDt_,mt_,\[Kappa]_,\[Xi]\[Beta]_,\[Xi]E_,l_:1]:= Module[{Integrand,E1ms,E1},
(*l - == 1 => energy transfer rate integrand, == 0 => interaction rate (over Subscript[n, D]) integrand*)
 E1ms=Integrate[Ett^l/(Ett+"\[CapitalLambda]")^2 (("ms" + "ESMt")^2+("ms"+(Ett-"ESMt"))^2),Ett];
	E1=Simplify[E1ms/."ms"->("m"^2+"mD"^2)/(2 "mD")];
	Integrand =Simplify[Simplify[(E1/.Ett->(2 "mD")/(("ESMt"+"mD")^2/("ESMt"^2-"m"^2)-1))-(E1/.Ett->0)]];
2/\[Pi] ("\[Alpha]"^2 \[Kappa]^2)/mDt (E^(-"\[Beta]" "ESMt" + "\[Beta]" "\[Mu]")/aH["\[Beta]"]^l Integrand /."\[Mu]"->"m" - 1/"\[Beta]" Log["Qt"]/.{"\[CapitalLambda]"->\[CapitalLambda]["\[Beta]",mt],"Qt"->Z1["\[Beta]"]}/."ESMt"->1/"\[Beta]" \[Xi]E + "m"/."\[Beta]"->\[Xi]\[Beta]/"m"/.{"m"->mt,"mD"->mDt})
]*)
NonNumIntegrand[mDt_,mt_,\[Kappa]_,\[Xi]\[Beta]_,\[Xi]E_,l_:1]:= Module[{P\[CapitalMu],prefactor,IntegrandEt,boundsE,boundsEt,thermalrepl,\[CapitalLambda]repl,msrepl,\[Xi]Erepl,\[Xi]\[Beta]repl,IntegrandESM},
(*l == 1 => energy transfer rate integrand (per particle per efold)
l == 0 => interaction rate (over Subscript[n, D]) integrand*)
P\[CapitalMu]=("ms"+"ESM")^2+("ms"+(Et-"ESM"))^2;
prefactor = ("n0")/(("T0" "\[Beta]")^3 aH["\[Beta]"]^l) 1/("Z1") ((4 \[Pi] "\[Alpha]")^2 \[Kappa]^2)/("mD" (2 \[Pi])^3) ("g")/4 ;(* prefactor of eqn (629) of notes over nD*)
IntegrandEt=E^(-"\[Beta]"("ESM"-"m")) P\[CapitalMu]/(Et + "\[CapitalLambda]")^2 Et^l;
boundsEt={0,(2 "mD")/(("ESM"+"mD")^2/(("ESM")^2-("m")^2)-1) };
thermalrepl={"\[CapitalLambda]"->"\[CapitalLambda]0" (1/("T0" "\[Beta]"))^2,"Z1"->"g" ((2 \[Pi] "\[Beta]")/( "m"))^(-3/2)};
\[CapitalLambda]repl={"\[CapitalLambda]0"->4\[Pi] "\[Alpha]" (2 \[Pi] "n0")/("m" "T0")};
msrepl = {"ms"->(("m")^2+("mD")^2)/(2 "mD")};
\[Xi]Erepl={"ESM"->\[Xi]E /"\[Beta]" + "m"};
\[Xi]\[Beta]repl="\[Beta]"->\[Xi]\[Beta]/"m";

IntegrandESM=Integrate[IntegrandEt,Et];
prefactor((IntegrandESM/.Et->boundsEt[[2]])-(IntegrandESM/.Et->boundsEt[[1]]))/.msrepl/.thermalrepl/.\[CapitalLambda]repl/.\[Xi]Erepl/.\[Xi]\[Beta]repl/.{"m"->mt,"mD"->mDt}
]


Ettotheld\[CapitalGamma]dEt[mD_,m_,\[Kappa]_,\[Beta]_,Et_,l_:1]:=Module[{P\[CapitalMu],prefactor,EintsIntegrand,boundsE,boundsEt,thermalrepl,\[CapitalLambda]repl,msrepl,\[Xi]Erepl,\[Xi]\[Beta]repl,IntegrandESM,pthav,Etint,},
P\[CapitalMu]=("ms"+ESM)^2+("ms"+(Et-ESM))^2;
prefactor = ("n0")/("T0" \[Beta])^3 1/("Z1") ((4 \[Pi] "\[Alpha]")^2 \[Kappa]^2)/(mD (2 \[Pi])^3) ("g")/4 ;(* prefactor of eqn (629) of notes over nD*)
EintsIntegrand=E^(-\[Beta](ESM-m)) P\[CapitalMu]/(Et + "\[CapitalLambda]")^2 Et^l;
boundsE= {Et/2+1/2 Sqrt[((Et + 2 mD)(2 m^2+Et mD))/mD],\[Infinity]};
(*boundsEt={0,\[Infinity] };*)
(*shorthandrepl={"ms"->(m^2+mD^2)/(2 mD),"\[CapitalLambda]"->4\[Pi] "\[Alpha]" (2 \[Pi] "n0")/(m "T0") (1/("T0" \[Beta]))^2,"Z1"->"g" ((2 \[Pi] \[Beta])/ m)^(-3/2)};*)
thermalrepl={"\[CapitalLambda]"->"\[CapitalLambda]0" (1/("T0" \[Beta]))^2,"Z1"->"g" ((2 \[Pi] \[Beta])/ m)^(-3/2)};
\[CapitalLambda]repl={"\[CapitalLambda]0"->4\[Pi] "\[Alpha]" (2 \[Pi] "n0")/(m "T0")};
msrepl = {"ms"->(m^2+mD^2)/(2 mD)};
pthav=3 Sqrt[2] Sqrt[(2 m)/(\[Beta] 2 \[Pi])]; (*average thermal momentum for boltzman distributed species (ie for SM electrons)*)
Etint=0-(Integrate[prefactor EintsIntegrand,ESM]/.ESM->boundsE[[1]])//Simplify;
<|"Ettotheld\[CapitalGamma]dEt"->Etint,"msrepl"->msrepl,"\[CapitalLambda]repl"->\[CapitalLambda]repl,"thermalrepl"->thermalrepl|>
]


(*P\[CapitalMu]=("ms"+ESM)^2+("ms"+(Et-ESM))^2;
prefactor = ("n0")/("T0" \[Beta])^31/("Z1")((4 \[Pi]"\[Alpha]")^2 \[Kappa]^2)/(mD (2 \[Pi])^3) ("g")/4 ;(* prefactor of eqn (629) of notes over nD*)
integrand=E^(-\[Beta](ESM-m))P\[CapitalMu]/(Et + "\[CapitalLambda]")^2Et^l
boundsE= {m,\[Infinity]};
boundsEt={0,(2 mD)/((ESM+mD)^2/(ESM^2-m^2)-1) };
(*shorthandrepl={"ms"->(m^2+mD^2)/(2 mD),"\[CapitalLambda]"->4\[Pi] "\[Alpha]" (2 \[Pi] "n0")/(m "T0") (1/("T0" \[Beta]))^2,"Z1"->"g" ((2 \[Pi] \[Beta])/ m)^(-3/2)};*)
thermalrepl={"\[CapitalLambda]"->"\[CapitalLambda]0" (1/("T0" \[Beta]))^2,"Z1"->"g" ((2 \[Pi] \[Beta])/ m)^(-3/2)};
\[CapitalLambda]repl={"\[CapitalLambda]0"->4\[Pi] "\[Alpha]" (2 \[Pi] "n0")/(m "T0")}
msrepl = {"ms"->(m^2+mD^2)/(2 mD)}
pthav=3 Sqrt[2] Sqrt[(2 m)/(\[Beta] 2 \[Pi])] (*average thermal momentum for boltzman distributed species (ie for SM electrons)*)
Integrate[integrand/. l->0,Et];
(*Etintegral=((%/.Et->boundsEt[[2]])-(%/.Et->boundsEt[[1]]))*)
Etintegral= Et  integrand/.l->0/.Et->pthav^2/(2 m)*)


ESMIntegral[mD_,m_,\[Kappa]_,\[Xi]\[Beta]_?NumberQ,params_,l_:1]:=(m/.params)/\[Xi]\[Beta] NIntegrate[ NumIntegrand[mD,m,\[Kappa],\[Xi]\[Beta],\[Xi]E,params,l],{\[Xi]E,0,\[Infinity]}](*extra factor comes from jacobian*)


(* ::Text:: *)
(*Interpolate over temperature, find total energy transferred to dark sector at recombination*)


InterpolateEtransRate[mD_,m_,\[Kappa]_,params_,l_:1]:=Module[{n=13,Inter\[Xi]\[Beta]s,InterPoints,RateLogFnof\[Beta],\[CapitalDelta]Etransfered,Td,TR,Tdt},

Td = "Td"/.params; (*this is Subscript[T, 0]*)
TR = "TR" /.params;
Tdt = "Tdt"/.params;

(*Interpolate to find the energy transfer rate between e e decoupling and recombination*)
Inter\[Xi]\[Beta]s= Subdivide[Log10[m/Td],Log10[m/TR],n]; (*points to be samples, we have made \[Beta] dimensionless (\[Xi]\[Beta]=m \[Beta])*)
InterPoints=Table[{Inter\[Xi]\[Beta]s[[i]],Log10[Abs[ESMIntegral[mD,m,\[Kappa],10^Inter\[Xi]\[Beta]s[[i]],params,l]]]},{i,Length[Inter\[Xi]\[Beta]s]}] ;(*Abs since numerical instability at early times*)
RateLogFnof\[Beta]=Interpolation[InterPoints];
\[CapitalDelta]Etransfered = NIntegrate[Tdt 10^RateLogFnof\[Beta][Log10[\[Beta] m]],{\[Beta],1/Td,1/TR}];
(*RateFnof\[Beta][\[Beta]_] := 10^RateLogFnof\[Beta][Log10[\[Beta] m]];*)
(*<|"f"->RateFnof\[Beta],"LogLogf"->RateLogFnof\[Beta],"\[Xi]\[Beta]"->Inter\[Xi]\[Beta]s,"Points"->InterPoints|>*)
<|"f"->RateLogFnof\[Beta],"\[CapitalDelta]Etot"->\[CapitalDelta]Etransfered,"\[Xi]\[Beta]"->Inter\[Xi]\[Beta]s,"Points"->InterPoints,"mD"->mD,"m"->m,"\[Kappa]"->\[Kappa],"params"->params|>
]


(*Interpolatey[mD_,m_,mpD_,fD_,\[Kappa]_,params_,l_:1]:=Module[{n=13,Inter\[Xi]\[Beta]s,InterPoints,RateLogFnof\[Beta],\[CapitalDelta]Etransfered,Td,TR,Tdt},

Td = "Td"/.params; (*this is Subscript[T, 0]*)
TR = "TR" /.params;
Tdt = "Tdt"/.params;

(*Interpolate to find the energy transfer rate between e e decoupling and recombination*)
Inter\[Xi]\[Beta]s= Subdivide[Log10[m/Td],Log10[m/TR],n]; (*points to be samples, we have made \[Beta] dimensionless (\[Xi]\[Beta]=m \[Beta])*)
InterPoints=Table[{Inter\[Xi]\[Beta]s[[i]],Log10[Abs[ESMIntegral[mD,m,\[Kappa],10^Inter\[Xi]\[Beta]s[[i]],params,l]]]},{i,Length[Inter\[Xi]\[Beta]s]}] ;(*Abs since numerical instability at early times*)
RateLogFnof\[Beta]=Interpolation[InterPoints];
\[CapitalDelta]Etransfered = NIntegrate[Tdt 10^RateLogFnof\[Beta][Log10[\[Beta] m]],{\[Beta],1/Td,1/TR}];
(*RateFnof\[Beta][\[Beta]_] := 10^RateLogFnof\[Beta][Log10[\[Beta] m]];*)
(*<|"f"->RateFnof\[Beta],"LogLogf"->RateLogFnof\[Beta],"\[Xi]\[Beta]"->Inter\[Xi]\[Beta]s,"Points"->InterPoints|>*)
<|"f"->RateLogFnof\[Beta],"\[CapitalDelta]Etot"->\[CapitalDelta]Etransfered,"\[Xi]\[Beta]"->Inter\[Xi]\[Beta]s,"Points"->InterPoints,"mD"->mD,"m"->m,"\[Kappa]"->\[Kappa],"params"->params|>
]*)
Interpolatey[mD_,m_,mpD_,fD_,\[Kappa]_,params_,l_:1]:=Module[{n=13,Inter\[Xi]\[Beta]s,InterPoints,yLogFnof\[Beta],yprefactor,y,TB,TR,T0},

TB = "TB"/.params; (*brehmstrahlung inefficiency*)
TR = "TR" /.params; (*recombination*)
T0 = "T0"/.params; (*today*)

(*Interpolate to find the energy transfer rate between e e decoupling and recombination*)
Inter\[Xi]\[Beta]s= Subdivide[Log10[m/TB],Log10[m/TR],n]; (*points to be samples, we have made \[Beta] dimensionless (\[Xi]\[Beta]=m \[Beta])*)
(*Print[Inter\[Xi]\[Beta]s];*)
yprefactor= T0 /4  (nD0[fD,mpD]/(T0^3 \[Beta]^3))(\[Pi]^2/30 2)^-1 \[Beta]^4/.params;
(*Print[Table[(yprefactor /.\[Beta]->10^Inter\[Xi]\[Beta]s[[i]]/m),{i,Length[Inter\[Xi]\[Beta]s]}]];*)
InterPoints=Table[{Inter\[Xi]\[Beta]s[[i]],Log10[Abs[(yprefactor /.\[Beta]->10^Inter\[Xi]\[Beta]s[[i]]/m)ESMIntegral[mD,m,\[Kappa],10^Inter\[Xi]\[Beta]s[[i]],params,l]]]},{i,Length[Inter\[Xi]\[Beta]s]}] ;
Print[InterPoints];(*Abs since numerical instability at early times*)
(*computed in the notes*)
yLogFnof\[Beta]=Interpolation[InterPoints];
y = NIntegrate[10^yLogFnof\[Beta][Log10[\[Beta] m]],{\[Beta],1/TB,1/TR}];
(*RateFnof\[Beta][\[Beta]_] := 10^RateLogFnof\[Beta][Log10[\[Beta] m]];*)
(*<|"f"->RateFnof\[Beta],"LogLogf"->RateLogFnof\[Beta],"\[Xi]\[Beta]"->Inter\[Xi]\[Beta]s,"Points"->InterPoints|>*)
<|"f"->yLogFnof\[Beta],"y"->y,"\[Xi]\[Beta]"->Inter\[Xi]\[Beta]s,"Points"->InterPoints,"mD"->mD,"m"->m,"\[Kappa]"->\[Kappa],"params"->params|>
]


(* ::Text:: *)
(*Interpolate over a given mass range*)


InterpolateTotEtrans[mD_:{10^-3 ("me"/.CosmoParams[]),10^4 ("me"/.CosmoParams[])},m_:("me"/.CosmoParams[]),\[Kappa]_:1,params_,l_:1]:=Module[{n=Round[Log10[mD[[2]]]-Log10[mD[[1]]]]+1,Interms,TableofInters,\[CapitalDelta]EPoints,TotLogFnofm},
(*Interpolate to find the energy transfer rate between e e decoupling and recombination*)
Interms= N[Subdivide[Log10[mD[[2]]],Log10[mD[[1]]],n]]; (*points to be samples, we have made \[Beta] dimensionless (\[Xi]\[Beta]=m \[Beta])*)
TableofInters=Table[InterpolateEtransRate[10^Interms[[i]],m,\[Kappa],params,l],{i,Length[Interms]}];(*table of interpolations*)
\[CapitalDelta]EPoints=Table[{Interms[[i]],Log10[TableofInters[[i]][["\[CapitalDelta]Etot"]]]},{i,Length[Interms]}] ;(*Abs since numerical instability at early times*)
TotLogFnofm=Interpolation[\[CapitalDelta]EPoints];
(*RateFnof\[Beta][\[Beta]_] := 10^RateLogFnof\[Beta][Log10[\[Beta] m]];*)
(*<|"f"->RateFnof\[Beta],"LogLogf"->RateLogFnof\[Beta],"\[Xi]\[Beta]"->Inter\[Xi]\[Beta]s,"Points"->InterPoints|>*)
<|"f"->TotLogFnofm,"Inters"->TableofInters,"Points"->\[CapitalDelta]EPoints,"Interms"->Interms,"mD"->mD,"m"->m,"\[Kappa]"->\[Kappa],"params"->params|>
]


(* ::Text:: *)
(*Neff constraint and dark temperature*)


(*Compute\[CapitalDelta]NeffConstraint[\[CapitalDelta]EtotDict_]:=Module[{gstar,\[CapitalDelta]Neff,constraint,TDRconst},
gstar=2 (TDRt/"TR")^4/.\[CapitalDelta]EtotDict[["params"]];(* []effective gstar for dark photon at recombination*)
\[CapitalDelta]Neff = 8/7 (11/4)^(4/3) gstar; (*[]*)
constraint = 0.1; (*[]constraint on \[CapitalDelta]Neff*)
TDRconst = TDRt/.Solve[\[CapitalDelta]Neff==constraint,TDRt][[4]] ;(*[eV]upper bound on TD at recombination*)
Print[Plot[\[CapitalDelta]EtotDict[["f"]][mDM+6],{mDM,Log10[\[CapitalDelta]EtotDict[["mD"]][[1]]]-6+6,Log10[\[CapitalDelta]EtotDict[["mD"]][[2]]]-6},PlotLabel->"\!\(\*SubscriptBox[\(Log\), \(10\)]\)[Energy Transfered to DS] [eV]",AxesLabel->{"\!\(\*SubscriptBox[\(Log\), \(10\)]\)[\!\(\*SubscriptBox[\(m\), \(D\)]\)] [MeV]"}]];
Print[ContourPlot[Log10[\[Pi]^2/45 TDRconst^3/nDR[10^fDt,10^(mDM+6)]/.\[CapitalDelta]EtotDict[["params"]]],{fDt,-5,0},{mDM,Log10[\[CapitalDelta]EtotDict[["mD"]][[1]]]-6,Log10[\[CapitalDelta]EtotDict[["mD"]][[2]]]-6},ContourLabels->True,FrameLabel->{"\!\(\*SubscriptBox[\(Log\), \(10\)]\)[\!\(\*SubscriptBox[\(f\), \(D\)]\)] []","\!\(\*SubscriptBox[\(Log\), \(10\)]\)[\!\(\*SubscriptBox[\(m\), \(D\)]\)] [MeV]"},PlotLabel->"\!\(\*SubscriptBox[\(Log\), \(10\)]\)[\!\(\*FractionBox[SubscriptBox[\(\[Rho]\), SubscriptBox[\(\[Gamma]\), \(D\)]], SubscriptBox[\(\[Rho]\), \(aDM\)]]\)]"]];
Print[ContourPlot[1/2 Log10[(3TDRconst)/10^\[CapitalDelta]EtotDict[["f"]][6+mDM] (1+\[Pi]^2/45 TDRconst^3/nDR[10^fDt,10^(6+mDM)])/.\[CapitalDelta]EtotDict[["params"]]],{fDt,-5,0},{mDM,Log10[\[CapitalDelta]EtotDict[["mD"]][[1]]]-6,Log10[\[CapitalDelta]EtotDict[["mD"]][[2]]]-6},ContourLabels->True,PlotLabel->"\!\(\*SubscriptBox[\(Log\), \(10\)]\)[\[Kappa]] [] - \!\(\*SubscriptBox[\(\[CapitalDelta]N\), \(eff\)]\) bound",FrameLabel->{"\!\(\*SubscriptBox[\(Log\), \(10\)]\)[\!\(\*SubscriptBox[\(f\), \(D\)]\)] []","\!\(\*SubscriptBox[\(Log\), \(10\)]\)[\!\(\*SubscriptBox[\(m\), \(D\)]\)] [MeV]"}]]
(*Print[ContourPlot[1/2Log10[(3TDRconst)/10^\[CapitalDelta]EtotDict[["f"]][6+mDM](1+\[Pi]^2/45TDRconst^3/nDR[10^fDt,10^(6+mDM)])],{fDt,-5,0},{mDM,Log10[10^-3 me]-6,Log10[10^4 me]-6},ContourLabels->True,PlotLabel->"Subscript[Log, 10][\[Kappa]] [] - Subscript[\[CapitalDelta]N, eff] bound",FrameLabel->{"Subscript[Log, 10][Subscript[f, D]] []","Subscript[Log, 10][Subscript[m, D]] [MeV]"}]]*)
]*)
Compute\[CapitalDelta]NeffConstraint[\[CapitalDelta]EtotDict_,mpD_]:=Module[{gstar,\[CapitalDelta]Neff,constraint,TDRconst},
gstar=2 (TDRt/"TR")^4/.\[CapitalDelta]EtotDict[["params"]];(* []effective gstar for dark photon at recombination*)
\[CapitalDelta]Neff = 8/7 (11/4)^(4/3) gstar; (*[]*)
constraint = 0.1; (*[]constraint on \[CapitalDelta]Neff*)
TDRconst = TDRt/.Solve[\[CapitalDelta]Neff==constraint,TDRt][[4]] ;(*[eV]upper bound on TD at recombination*)
Print[Plot[\[CapitalDelta]EtotDict[["f"]][mDM+6],{mDM,Log10[\[CapitalDelta]EtotDict[["mD"]][[1]]]-6,Log10[\[CapitalDelta]EtotDict[["mD"]][[2]]]-6},PlotLabel->"\!\(\*SubscriptBox[\(Log\), \(10\)]\)[Energy Transfered to DS] [eV]",AxesLabel->{"\!\(\*SubscriptBox[\(Log\), \(10\)]\)[\!\(\*SubscriptBox[\(m\), \(D\)]\)] [MeV]"}]];
Print[ContourPlot[Log10[\[Pi]^2/45 TDRconst^3/nDR[10^fDt,10^(6+mDM)]/.\[CapitalDelta]EtotDict[["params"]]],{fDt,-5,0},{mDM,Log10[\[CapitalDelta]EtotDict[["mD"]][[1]]]-6,Log10[\[CapitalDelta]EtotDict[["mD"]][[2]]]-6},ContourLabels->True,FrameLabel->{"\!\(\*SubscriptBox[\(Log\), \(10\)]\)[\!\(\*SubscriptBox[\(f\), \(D\)]\)] []","\!\(\*SubscriptBox[\(Log\), \(10\)]\)[\!\(\*SubscriptBox[\(m\), \(pD\)]\)] [MeV]"},PlotLabel->"\!\(\*SubscriptBox[\(Log\), \(10\)]\)[\!\(\*FractionBox[SubscriptBox[\(\[Rho]\), SubscriptBox[\(\[Gamma]\), \(D\)]], SubscriptBox[\(\[Rho]\), \(aDM\)]]\)]"]];
Print[ContourPlot[1/2 Log10[(3TDRconst)/10^\[CapitalDelta]EtotDict[["f"]][6+mDM] (1+\[Pi]^2/45 TDRconst^3/nDR[10^fDt,mpD])/.\[CapitalDelta]EtotDict[["params"]]],{fDt,-5,0},{mDM,Log10[\[CapitalDelta]EtotDict[["mD"]][[1]]]-6,Log10[\[CapitalDelta]EtotDict[["mD"]][[2]]]-6},ContourLabels->True,PlotLabel->"\!\(\*SubscriptBox[\(Log\), \(10\)]\)[\[Kappa]] [] - \!\(\*SubscriptBox[\(\[CapitalDelta]N\), \(eff\)]\) bound",FrameLabel->{"\!\(\*SubscriptBox[\(Log\), \(10\)]\)[\!\(\*SubscriptBox[\(f\), \(D\)]\)] []","\!\(\*SubscriptBox[\(Log\), \(10\)]\)[\!\(\*SubscriptBox[\(m\), \(D\)]\)] [MeV]"}]]
(*Print[ContourPlot[1/2Log10[(3TDRconst)/10^\[CapitalDelta]EtotDict[["f"]][6+mDM](1+\[Pi]^2/45TDRconst^3/nDR[10^fDt,10^(6+mDM)])],{fDt,-5,0},{mDM,Log10[10^-3 me]-6,Log10[10^4 me]-6},ContourLabels->True,PlotLabel->"Subscript[Log, 10][\[Kappa]] [] - Subscript[\[CapitalDelta]N, eff] bound",FrameLabel->{"Subscript[Log, 10][Subscript[f, D]] []","Subscript[Log, 10][Subscript[m, D]] [MeV]"}]]*)
]


(* ::Section::Closed:: *)
(*End*)


End[];


EndPackage[];
