(* ::Package:: *)

BeginPackage["NeffConstraint`"];


(* ::Text:: *)
(*The purpose of this package is to define modules and constants used in the computation of the Subscript[\[CapitalDelta]N, eff] constraint on \[Kappa]*)


(* ::Subsubsection:: *)
(*Public Declarations*)


(* ::Section:: *)
(*Private*)


Begin["`Private`"];


(* ::Text:: *)
(*Cosmo Constants - temperatures and scale factors*)


Td = 5 10^5; (*[eV] annihilation temperature*)
ad =(2 10^9)^-1; (*[] scale factor (~(1/z)) at annihilation*)
Tdt = 2.5 10^-4 (*[eV] Temperature of the Universe today*)
T0 = 2.5 10^-4 (*[eV] Temperature of the Universe today*)
aR = (1100)^-1;(*[] scale factor ~ (1/z) at recombination*)
TR = (Tdt/aR) (*[eV] Temperature at recombination*)


(* ::Text:: *)
(*Thermal Quantities*)


me = 5 10^5; (*[eV] electron mass*)
\[Alpha]=1/137;
pth =( (me Tdt)/(2 \[Pi]))^(1/2); (*[eV] thermal de-broglie momentum*)
g = 4; (*[] dofs in SM electron*)
Z1[\[Beta]_] := (g pth^3)/n0 (1/(Tdt \[Beta]))^(-3/2)(*[] single particle classical partition function*)
kF = N[(3 \[Pi]^2 n0)^(1/3)] (*[eV] fermi momentum today*)
EF = kF^2/(2 me)
EFd = EF Td/Tdt (*[eV]fermi energy at decoupling*)
Degen = EF/Tdt (*Degeneracy of gas today - degeneracy is independent of time*)



(* ::Text:: *)
(*Conversion*)


eVm = N[200 10^-15 10^6] (*[eV m] = 1 follows from 1 = 200 MeV fm*)


(* ::Text:: *)
(*Hubble*)


H0tinv = 4.4 10^17; (*s*)
H0linv = 3 10^8 H0tinv; (*m*)
H0inv = (200 10^-15 10^6)^-1 H0linv; (*eV^-1*)
H0= (H0inv^-1) (*eV*) 


(* ::Text:: *)
(*Critical density*)


h= 0.67 (*baumann 1.2.64*)
\[Rho]critbaumanncm = 1.1 10^-5 h^2 (*[protons cm^-3]*) 
\[Rho]critbaumannm = \[Rho]critbaumanncm  (10^6) (*[protons m^-3]*)


(* ::Text:: *)
(*Define hubble and scale speed*)


\[CapitalLambda]CDMrepl= {\[CapitalOmega]k ->0.01,\[CapitalOmega]r->9.4 10^-5, \[CapitalOmega]m->0.32,\[CapitalOmega]\[CapitalLambda]->0.68,\[CapitalOmega]b->0.05,\[CapitalOmega]c->0.27}; (*[] \[CapitalLambda]CDM density parameters*)
H[\[Beta]_] :=H0 (\[CapitalOmega]r a^-4+\[CapitalOmega]m a^-3+\[CapitalOmega]k a^-2 + \[CapitalOmega]\[CapitalLambda])^(1/2)/.\[CapitalLambda]CDMrepl/.a->(Tdt \[Beta])
aH[\[Beta]_] :=a H0 (\[CapitalOmega]r a^-4+\[CapitalOmega]m a^-3+\[CapitalOmega]k a^-2 + \[CapitalOmega]\[CapitalLambda])^(1/2)/.\[CapitalLambda]CDMrepl/.a->(Tdt \[Beta])


(* ::Text:: *)
(*Number density today*)


n0 = \[CapitalOmega]b \[Rho]critbaumannm eVm^3 /.\[CapitalLambda]CDMrepl  (*[eV^3] proton / electron number density today*)
ne[\[Beta]_]:= n0 (1/(\[Beta] Tdt))^3 (*[eV^3] number density as a function of temperature*)
nd = n0 (Td/Tdt)^3 (*[eV^3] number density at decoupling*)
nr = n0 ((Td/Tdt)^3) (*[eV^3] number density at decoupling*)


(* ::Text:: *)
(*Determine interaction rate as a function of temperature*)



Integrate[Ett/(Ett+\[CapitalLambda])^2 ((ms + ESMt)^2+(ms+(Ett-ESMt))^2),Ett];
Simplify[%/.ms->(m^2+mD^2)/(2 mD)];
tempintegrand =Simplify[Simplify[(%/.Ett->(2 mD)/((ESMt+mD)^2/(ESMt^2-m^2)-1))-(%/.Ett->0)]];
NumIntegrand[mDt_,mt_,\[Kappa]_,\[Xi]\[Beta]_?NumberQ,\[Xi]E_?NumberQ]:=2/\[Pi] (\[Alpha]^2 \[Kappa]^2)/mDt (E^(-\[Beta] ESMt + \[Beta] \[Mu])/aH[\[Beta]] tempintegrand /.\[Mu]->m - 1/\[Beta] Log[Qt]/.{\[CapitalLambda]->\[CapitalLambda][\[Beta]],Qt->Z1[\[Beta]]}/.ESMt->1/\[Beta] \[Xi]E + m/.\[Beta]->\[Xi]\[Beta]/m/.{m->mt,mD->mDt})
ESMIntegral[mD_,m_,\[Kappa]_,\[Xi]\[Beta]_?NumberQ]:=NIntegrate[m/\[Xi]\[Beta] NumIntegrand[mD,m,\[Kappa],\[Xi]\[Beta],\[Xi]E],{\[Xi]E,0,\[Infinity]}](*extra factor comes from jacobian*)


(* ::Text:: *)
(*Interpolate over temperature, find total energy transfered to dark sector at recombination*)


InterpolateEtransRate[mD_,m_,\[Kappa]_]:=Module[{n=13,Inter\[Xi]\[Beta]s,InterPoints,RateLogFnof\[Beta],\[CapitalDelta]Etransfered},
(*Interpolate to find the energy transfer rate between e e decoupling and recombination*)
Inter\[Xi]\[Beta]s= Subdivide[Log10[m/Td],Log10[m/TR],n]; (*points to be samples, we have made \[Beta] dimensionless (\[Xi]\[Beta]=m \[Beta])*)
InterPoints=Table[{Inter\[Xi]\[Beta]s[[i]],Log10[Abs[ESMIntegral[mD,m,\[Kappa],10^Inter\[Xi]\[Beta]s[[i]]]]]},{i,Length[Inter\[Xi]\[Beta]s]}] ;(*Abs since numerical instability at early times*)
RateLogFnof\[Beta]=Interpolation[InterPoints];
\[CapitalDelta]Etransfered = NIntegrate[Tdt 10^RateLogFnof\[Beta][Log10[\[Beta] me]],{\[Beta],1/Td,1/TR}];
(*RateFnof\[Beta][\[Beta]_] := 10^RateLogFnof\[Beta][Log10[\[Beta] m]];*)
(*<|"f"->RateFnof\[Beta],"LogLogf"->RateLogFnof\[Beta],"\[Xi]\[Beta]"->Inter\[Xi]\[Beta]s,"Points"->InterPoints|>*)
<|"f"->RateLogFnof\[Beta],"\[CapitalDelta]Etot"->\[CapitalDelta]Etransfered,"\[Xi]\[Beta]"->Inter\[Xi]\[Beta]s,"Points"->InterPoints,"mD"->mD,"m"->m,"\[Kappa]"->\[Kappa]|>
]


(* ::Text:: *)
(*Interpolate over a given mass range*)


InterpolateTotEtrans[mD_:{10^-3 me,10^4 me},m_:me,\[Kappa]_:1]:=Module[{n=Round[Log10[mD[[2]]]-Log10[mD[[1]]]]+1,Interms,TableofInters,\[CapitalDelta]EPoints,TotLogFnofm},
(*Interpolate to find the energy transfer rate between e e decoupling and recombination*)
Interms= N[Subdivide[Log10[mD[[2]]],Log10[mD[[1]]],n]]; (*points to be samples, we have made \[Beta] dimensionless (\[Xi]\[Beta]=m \[Beta])*)
TableofInters=Table[InterpolateEtransRate[10^Interms[[i]],m,\[Kappa]],{i,Length[Interms]}];(*table of interpolations*)
\[CapitalDelta]EPoints=Table[{Interms[[i]],Log10[TableofInters[[i]][["\[CapitalDelta]Etot"]]]},{i,Length[Interms]}] ;(*Abs since numerical instability at early times*)
TotLogFnofm=Interpolation[\[CapitalDelta]EPoints];
(*RateFnof\[Beta][\[Beta]_] := 10^RateLogFnof\[Beta][Log10[\[Beta] m]];*)
(*<|"f"->RateFnof\[Beta],"LogLogf"->RateLogFnof\[Beta],"\[Xi]\[Beta]"->Inter\[Xi]\[Beta]s,"Points"->InterPoints|>*)
<|"f"->TotLogFnofm,"Inters"->TableofInters,"Points"->\[CapitalDelta]EPoints,"Interms"->Interms,"mD"->mD,"m"->m,"\[Kappa]"->\[Kappa]|>
]


(* ::Section::Closed:: *)
(*End*)


End[];


EndPackage[];
