(* ::Package:: *)

Needs["FormFactors`"]
Needs["Constants`"]
Needs["Utilities`"]
Needs["Dielectrics`"]


(*FrontEndTokenExecute@"Save"*)


BeginPackage["EnergyLoss`"];


(* ::Subsubsection:: *)
(*Public Declarations*)


(* ::Text:: *)
(*Kinematics*)


Kinematics::usage = "Make a mesh table of DM masses and velocity defining kinematics for EL tables";


(* ::Text:: *)
(*Electronic*)


LossIntegrand::usage = "Loss integrand from optical fits";

EnergyLossMSI::usage = "Compute the energy loss for 1 oscillator using theory collision rate";
(*EnergyLossMSIFit::usage = "Compute the energy loss for 1 oscillator using fit parameters";*)
EnergyLossMSIFitEnhanced::usage = "Compute the energy loss for 1 oscillator using fit parameters";

(*EnergyLossMSIFitEnhanced::usage = "Compute the energy loss for 1 oscillator using fit parameters with Temperature enhancement";
uMIntegralFitEnhanced::usage = "test";
*)
(*BandGapTruncation::usage = "";*)

(*EnergyLossMSIFitSumOscillators::usage = "Compute EL from optical fits of a material";*)
EnergyLossMSIFitSumOscillatorsEnhanced::usage = "Compute EL from optical fits of a material";

(*uMIntegralFitEnhanced::usage="";
zMIntegralFitEnhanced::usage="";*)
(*\[Epsilon]test::usage="";*)

EnergyLossTableAndInter::usage = "Create an EL table and interpolation function from theory collision rate";
(*EnergyLossTableAndInterFIT::usage = "Create an EL table and interpolation function from optical fits of a material";
EnergyLossTableAndInterFITTotalParams::usage = "Create an EL table and interpolation function from optical fits of a material using 1 repl table";*)
EnergyLossTableAndInterFITTotalParamsEnhanced::usage = "Create an EL table and interpolation function from optical fits of a material using 1 repl table including the temperature enhancement from FD Theorem";
(*EnergyLossTableAndInterFITTotalParamsEnhancedOscillators::usage = "Create an EL table and interpolation function from optical fits of a material using 1 repl table including the temperature enhancement from FD Theorem, and giving information for each oscillator";*)


(* ::Text:: *)
(*Nuclear*)


SNuc::usage = "Nuclear Structure function, using degenerate limit, and billiard ball atom approximation with nuclear and atomic form factors";
VCoul::usage = "Coulomb potential, excluding effective kinetic mixing and species charge";

\[Omega]IntegralNuc::usage = "\[Omega] integral appearing in nuclear Energy Loss integral at Finite Temperature";
ELIntegralNuc::usage = "Integral for the nuclear Energy Loss at Finite Temperature";
ELIntegralNuc0::usage = "Integral for the nuclear Energy Loss at 0 Temperature";

EnergyLossNuc::usage = "Compute energy loss at finite temperature from nuclear scattering";
EnergyLossNuc0::usage = "Compute energy loss at 0 temperature from nuclear scattering";

EnergyLossTableAndInterNuc0::usage = "Create an EL table and interpolation function for 0 Temperature Nuclear Scattering";
EnergyLossTableAndInterNuc::usage = "Create an EL table and interpolation function for 0 Temperature Nuclear Scattering";


(* ::Text:: *)
(*d\[Sigma]/(d Subscript[E, R]) - Nuclear*)


dPd\[Omega]Nuc0::usage = "Derivative of Transition rate wrt \[Omega] for Nuclear Interactions(proportional to \!\(\*FormBox[FractionBox[\(d\[Sigma]\), \(d\\\ \*SubscriptBox[\(E\), \(R\)]\)],
TraditionalForm]\)) at 0 Temperature";
dPd\[Omega]Nuc0Num::usage = "Derivative of Transition rate wrt \[Omega] for Nuclear Interactions(proportional to \!\(\*FormBox[FractionBox[\(d\[Sigma]\), \(d\\\ \*SubscriptBox[\(E\), \(R\)]\)],
TraditionalForm]\)) at 0 Temperature";

d\[Sigma]dERNuc0::usage = "\!\(\*FormBox[FractionBox[\(d\[Sigma]\), \(d\\\ \*SubscriptBox[\(E\), \(R\)]\)],
TraditionalForm]\) for Nuclear Interactions at 0 Temperature";


(* ::Text:: *)
(*d\[Sigma]/(d Subscript[E, R]) - Electronic*)


dPd\[Omega]e::usage = "Calculate dPd\[Omega] from electronic contributions (using optical data fits)";
dPd\[Omega]eNum::usage = "Numerically Calculate dPd\[Omega] from electronic contributions (using optical data fits)";

(*zIntegrald\[Sigma]dER::usage = "Calculate the integral over z appearing in dPd\[Omega]";*)

d\[Sigma]dERe::usage = "Calculate the distribution of cross-section with recoil energy for electrons";
d\[Sigma]dEReNum::usage = "Numerically Calculate the distribution of cross-section with recoil energy for electrons";


(* ::Text:: *)
(*Interpolate over d\[Sigma]dERe*)


Interpolated\[Sigma]dERe::usage = "Interpolate d\[Sigma]dER over \[Omega] for electronic contribution ";

\[Omega]of\[Xi]::usage = "Get \[Omega](\[Xi]) (where \[Xi] is the dimensionless parameter used in the interpolation)";

\[Omega]maxofv\[Chi]::usage = "Get the maximum energy transfer allowable by DM kinematics (units of [s^-1])";
\[Omega]of\[Xi]andv\[Chi]::usage = "\[Omega]of\[Xi] as a function of DM kinematics";
\[Xi]of\[Omega]andv\[Chi]::usage = "inverse of \[Omega]of\[Xi]andv\[Chi]";

\[Omega]maxNuc::usage = "Get the maximum energy transfer allowable by DM kinematics (units of [s^-1]) for Nuclear";
\[Omega]of\[Xi]Nuc::usage = "\[Omega]of\[Xi] as a function of DM kinematics";
\[Xi]of\[Omega]Nuc::usage = "inverse of \[Omega]of\[Xi]andv\[Chi]";

Getv\[Chi]and\[Xi]List::usage = "Get points to interpolate d\[Sigma]dER over at fixed v\[Chi]";

Interpolatev\[Chi]and\[Omega]\[Sigma]dERe\[Xi]::usage = "Interpolate d\[Sigma]dER over v\[Chi] and \[Omega](\[Xi]) for electronic";
Interpolatev\[Chi]and\[Omega]\[Sigma]dERNuc\[Xi]::usage = "Interpolate d\[Sigma]dER over v\[Chi] and \[Omega](\[Xi]) for nuclear";


(* ::Text:: *)
(*Plot Energy Loss Interpolation Function*)


PlotEL::usage = "Plot the EL interpolation function of m\[Chi] and v\[Chi] (in Natural units)";


(* ::Chapter:: *)
(*Private - Electronic*)


Begin["`Private`"];


(* ::Subsubsection::Closed:: *)
(*Kinematics - Public*)


Kinematics[m\[Chi]_,v\[Chi]_,meshdims_]:=Module[{lm\[Chi]list,lv\[Chi]list},
	lm\[Chi]list = E^Subdivide[Log[m\[Chi][[1]]],Log[m\[Chi][[2]]],meshdims[["m\[Chi]"]]-1];
	lv\[Chi]list = E^Subdivide[Log[v\[Chi][[1]]],Log[v\[Chi][[2]]],meshdims[["v\[Chi]"]]-1];
	(*{m\[Chi]list,v\[Chi]list}*)
	<|"m\[Chi]"->lm\[Chi]list,"v\[Chi]"->lv\[Chi]list|>
]


(* ::Section:: *)
(*Electonic*)


(* ::Subsection:: *)
(*EL Integrals - Electronic*)


(* ::Subsubsection::Closed:: *)
(*Integration Bounds*)


upp[z_?NumberQ,m\[Chi]_?NumberQ,v\[Chi]_?NumberQ,params_]:=v\[Chi]/"vF" + z "m"/m\[Chi]/.params (*u_+' - upper boundary of EL Integrals*)

upm[z_?NumberQ,m\[Chi]_?NumberQ,v\[Chi]_?NumberQ,params_]:=v\[Chi]/"vF" - z "m"/m\[Chi]/.params (* u_-' - lower boundary of EL Integrals*)


zt[m\[Chi]_?NumberQ,v\[Chi]_?NumberQ,params_]:=((m\[Chi] v\[Chi])/("m" "vF"))/.params (* 0 of upm*)


(* ::Subsubsection::Closed:: *)
(*Using Theoretical Collision Rate*)


uMIntegral[z_?NumberQ,m\[Chi]_?NumberQ,v\[Chi]_?NumberQ,params_] := NIntegrate[u Im[-1/Dielectrics`\[Epsilon]MNum[u,z,Dielectrics`u\[Nu][z]/.params,params]] ,{u,-upp[z,m\[Chi],v\[Chi],params],upm[z,m\[Chi],v\[Chi],params]}]

zMIntegral[m\[Chi]_?NumberQ,v\[Chi]_?NumberQ,params_] := NIntegrate[z uMIntegral[z,m\[Chi],v\[Chi],params],{z,0,\[Infinity]}]


(* ::Subsubsection::Closed:: *)
(*Using Optical Fits Parameters*)


u\[Nu]Fit[z_] := ("\[Nu]i"/(2 "vF" "qF" z)) (*u \[Nu]*)


uMIntegralFit[z_?NumberQ, m\[Chi]_?NumberQ, v\[Chi]_?NumberQ, params_] :=
    NIntegrate[u Im[-1 / Dielectrics`\[Epsilon]MNum[u, z, u\[Nu]Fit[z] /. params, params]], {u,
         -upp[z, m\[Chi], v\[Chi], params], upm[z, m\[Chi], v\[Chi], params]}]

zMIntegralFit[m\[Chi]_?NumberQ, v\[Chi]_?NumberQ, params_] :=
    NIntegrate[z uMIntegralFit[z, m\[Chi], v\[Chi], params], {z, 0, \[Infinity]}]


(* ::Text:: *)
(*Loss integrand from Fits*)


LossIntegrand[strtotalparams_,q_,\[Omega]_]:=Log10[Sum[( ("Ai")/("JpereV") (("qF" "vF")/("c"))^2 8 /Pi  ("e")^2/(4 \[Pi] "\[Epsilon]0")/.strtotalparams[[i]])Im[(-Dielectrics`uf[q,\[Omega]]Dielectrics`zf[q]/.strtotalparams)/Dielectrics`\[Epsilon]MNum[Dielectrics`uf[q,\[Omega]]/.strtotalparams[[i]],Dielectrics`zf[q]/.strtotalparams[[i]],("\[Nu]i")/(2 "vF" "qF" Dielectrics`zf[q])/.strtotalparams[[i]],strtotalparams[[i]]]],{i,Length[strtotalparams]}]]


(* ::Subsubsection::Closed:: *)
(*Using Optical Fits Parameters - With Temperature Enhancement*)


EnhancementFactorPW[u_,z_,IRCut_]:=\!\(\*
TagBox[GridBox[{
{"\[Piecewise]", GridBox[{
{"1", 
RowBox[{"1", "<=", 
RowBox[{"\"\<EF\>\"", " ", "\"\<\[Beta]\>\"", " ", "u", " ", "z"}]}]},
{
FractionBox["1", 
RowBox[{"1", "-", 
SuperscriptBox["E", 
RowBox[{
RowBox[{"-", "4"}], "\"\<\[Beta]\>\"", "\"\<EF\>\"", " ", "u", " ", "z"}]]}]], 
RowBox[{"IRCut", " ", "<=", 
RowBox[{"\"\<EF\>\"", " ", "\"\<\[Beta]\>\"", " ", "u", " ", "z"}], "<=", "1"}]},
{
RowBox[{
FractionBox["1", 
RowBox[{"4", " ", "\"\<EF\>\"", " ", "\"\<\[Beta]\>\"", " ", "z", " ", "u"}]], "+", 
FractionBox["1", "2"]}], 
RowBox[{
RowBox[{"\"\<EF\>\"", " ", "\"\<\[Beta]\>\"", " ", "u", " ", "z"}], " ", "<=", " ", "IRCut"}]}
},
AllowedDimensions->{2, Automatic},
Editable->True,
GridBoxAlignment->{"Columns" -> {{Left}}, "ColumnsIndexed" -> {}, "Rows" -> {{Baseline}}, "RowsIndexed" -> {}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "ColumnsIndexed" -> {}, "Rows" -> {{1.}}, "RowsIndexed" -> {}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.84]}, Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}, "RowsIndexed" -> {}},
Selectable->True]}
},
GridBoxAlignment->{"Columns" -> {{Left}}, "ColumnsIndexed" -> {}, "Rows" -> {{Baseline}}, "RowsIndexed" -> {}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "ColumnsIndexed" -> {}, "Rows" -> {{1.}}, "RowsIndexed" -> {}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.35]}, Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}, "RowsIndexed" -> {}}],
"Piecewise",
DeleteWithContents->True,
Editable->False,
SelectWithContents->True,
Selectable->False,
StripWrapperBoxes->True]\)(*Cut at low and high energies for stability (use limits)*)


BandGapTruncation[\[Omega]_,l_]:= If[l<0,HeavisideTheta["\[Omega]BG"-\[Omega]],1]


(*uMIntegralFitEnhanced[z_?NumberQ, m\[Chi]_?NumberQ, v\[Chi]_?NumberQ, params_] :=
    NIntegrate[u (EnhancementFactorPW[u,z,10^-3]/.params) Im[-1 / Dielectrics`\[Epsilon]MNum[u, z, u\[Nu]Fit[z] /. params, params]], {u,
         -upp[z, m\[Chi], v\[Chi], params], upm[z, m\[Chi], v\[Chi], params]}] (* only include below silicate band gap \[Omega]bg*)
    
zMIntegralFitEnhanced[m\[Chi]_?NumberQ, v\[Chi]_?NumberQ, params_] :=
    NIntegrate[z uMIntegralFitEnhanced[z, m\[Chi], v\[Chi], params], {z, 0, \[Infinity]}]*)
    
(*uMIntegralFitEnhanced[z_?NumberQ, m\[Chi]_?NumberQ, v\[Chi]_?NumberQ, params_,l_:1] :=
    NIntegrate[u^Abs[l] (BandGapTruncation[u z 2 "qF" "vF"/.params,l]/.params)(EnhancementFactorPW[u,z,10^-3]/.params) Im[-1 / Dielectrics`\[Epsilon]MNum[u, z, u\[Nu]Fit[z] /. params, params]], {u,
         -upp[z, m\[Chi], v\[Chi], params], upm[z, m\[Chi], v\[Chi], params]}] (* only include below silicate band gap \[Omega]bg*)
(* l - [] energy moment of the interaction per unit length as a function of energy distribution *)
    
zMIntegralFitEnhanced[m\[Chi]_?NumberQ, v\[Chi]_?NumberQ, params_,l_:1] :=
    NIntegrate[z^Abs[l] uMIntegralFitEnhanced[z, m\[Chi], v\[Chi], params,l], {z, 0, \[Infinity]}]*)
    
    (*uMIntegralFitEnhanced[z_?NumberQ, m\[Chi]_?NumberQ, v\[Chi]_?NumberQ, params_,l_:1] :=
    NIntegrate[u^l (BandGapTruncation[u z 2 "qF" "vF"/.params,l]/.params)(EnhancementFactorPW[u,z,10^-3]/.params) Im[-1 / Dielectrics`\[Epsilon]MNum[u, z, u\[Nu]Fit[z] /. params, params]], {u,
         -upp[z, m\[Chi], v\[Chi], params], upm[z, m\[Chi], v\[Chi], params]}] (* l=-1 => only include below silicate band gap \[Omega]bg*)
(* l - [] energy moment of the interaction per unit length as a function of energy distribution *)
    
zMIntegralFitEnhanced[m\[Chi]_?NumberQ, v\[Chi]_?NumberQ, params_,l_:1] :=
    NIntegrate[z^l uMIntegralFitEnhanced[z, m\[Chi], v\[Chi], params,l], {z, 0, \[Infinity]}]*)
    
    
 (*uMIntegralFitEnhanced[z_?NumberQ, m\[Chi]_?NumberQ, v\[Chi]_?NumberQ, params_,l_:1] :=
    NIntegrate[u^l (BandGapTruncation[u z 2 "qF" "vF"/.params,l]/.params)(EnhancementFactorPW[u,z,10^-3]/.params) Im[-1 / Dielectrics`\[Epsilon]MNum[u, z, u\[Nu]Fit[z] /. params, params]], {u,
         0, upm[z, m\[Chi], v\[Chi], params]}] (* l=-1 => only include below silicate band gap \[Omega]bg*)
(* l - [] energy moment of the interaction per unit length as a function of energy distribution *)
    
zMIntegralFitEnhanced[m\[Chi]_?NumberQ, v\[Chi]_?NumberQ, params_,l_:1] :=
    NIntegrate[z^l uMIntegralFitEnhanced[z, m\[Chi], v\[Chi], params,l], {z, 0, zt[m\[Chi],v\[Chi],params]}]
 *)   
    (*\[Epsilon]test[u_,z_,params_]:= Im[-1 / Dielectrics`\[Epsilon]MNum[u, z, u\[Nu]Fit[z] /. params, params]]*)
    
    uMIntegralFitEnhanced[z_?NumberQ, m\[Chi]_?NumberQ, v\[Chi]_?NumberQ, params_,l_:1] :=
    NIntegrate[u^Abs[l] (BandGapTruncation[u z 2 "qF" "vF"/.params,l]/.params)(EnhancementFactorPW[u,z,10^-3]/.params) Im[-1 / Dielectrics`\[Epsilon]MNum[u, z, u\[Nu]Fit[z] /. params, params]], {u,
         0, upm[z, m\[Chi], v\[Chi], params]}] (* l=-1 => only include below silicate band gap \[Omega]bg*)
(* l - [] energy moment of the interaction per unit length as a function of energy distribution *)
    
zMIntegralFitEnhanced[m\[Chi]_?NumberQ, v\[Chi]_?NumberQ, params_,l_:1] :=
    NIntegrate[z^Abs[l] uMIntegralFitEnhanced[z, m\[Chi], v\[Chi], params,l], {z, 0, zt[m\[Chi],v\[Chi],params]}]
    


(* ::Subsection:: *)
(*EL Single Oscillator - Public*)


(* ::Subsubsection::Closed:: *)
(*Using Theoretical Collision Rate*)


EnergyLossMSI[m\[Chi]_, v\[Chi]_, params_] :=
    Module[{lL, lPrefactor, ldEdr},
        lL = zMIntegral[m\[Chi], v\[Chi], params];
        lPrefactor = 1 / "JpereV" (("qF" "vF") / v\[Chi]) ^ 2 8 / Pi "e"^2 / (4 \[Pi] "\[Epsilon]0") /. params;
        ldEdr = lPrefactor Total[lL];
        Print["For m\[Chi] = ", m\[Chi], " , v\[Chi] = ", v\[Chi], "\nL:\t", Total[lL] "\nPrefactor:\t",
             lPrefactor, "\ndEdr:\t", ldEdr];
        <|"m\[Chi]" -> m\[Chi], "v\[Chi]" -> v\[Chi], "params" -> params, "L" -> lL,
             "Prefactor" -> lPrefactor, "dEdr" -> ldEdr|>
    ]


(* ::Subsubsection::Closed:: *)
(*Using Optical Fits Parameters*)


(*EnergyLossMSIFit[m\[Chi]_, v\[Chi]_, fitrepl_, consts_, T_] :=
    Module[{lL, lPrefactor, ldEdr, fitparams},
        fitparams = Join[Constants`SIParams[T, (("m" "\[Epsilon]0" "\[Omega]i"^2) / ("e"^2) /. fitrepl 
            /. consts)], fitrepl]; (*this extracts the number density from the fitted
             plasmon frequency*)
        lL = zMIntegralFit[m\[Chi], v\[Chi], fitparams];
        lPrefactor = "Ai" / "JpereV" (("qF" "vF") / v\[Chi])^2 8 / Pi "e"^2 / (4 \[Pi] "\[Epsilon]0") /. fitparams /. fitrepl;
        ldEdr = lPrefactor Total[lL];
        Print["For m\[Chi] = ", m\[Chi], " , v\[Chi] = ", v\[Chi], "\nL:\t", Total[lL] "\nPrefactor:\t",
             lPrefactor, "\ndEdr:\t", ldEdr];
        <|"m\[Chi]" -> m\[Chi], "v\[Chi]" -> v\[Chi], "fitparams" -> fitparams, "fitrepl"
             -> fitrepl, "L" -> lL, "Prefactor" -> lPrefactor, "dEdr" -> ldEdr|>
    ]

EnergyLossMSIFit[m\[Chi]_, v\[Chi]_, fitparams_] :=
    Module[
        {lL, lPrefactor, ldEdr},
        lL = zMIntegralFit[m\[Chi], v\[Chi], fitparams];
        lPrefactor = "Ai" / "JpereV" (("qF" "vF") / v\[Chi])^2 8 / Pi "e"^2 / (4 \[Pi] "\[Epsilon]0") /. fitparams;
        ldEdr = lPrefactor Total[lL];
        Print["For m\[Chi] = ", m\[Chi], " , v\[Chi] = ", v\[Chi], "\nL:\t", Total[lL] "\nPrefactor:\t",
             lPrefactor, "\ndEdr:\t", ldEdr];
        <|"m\[Chi]" -> m\[Chi], "v\[Chi]" -> v\[Chi], "fitparams" -> fitparams, "L" -> lL,
             "Prefactor" -> lPrefactor, "dEdr" -> ldEdr|>
    ]*)


(* ::Subsubsection::Closed:: *)
(*Using Optical Fits Parameters - With Temperature Enhancement*)


(*EnergyLossMSIFitEnhanced[m\[Chi]_, v\[Chi]_, fitrepl_, consts_, T_] :=
    Module[{lL, lPrefactor, ldEdr, fitparams},
        fitparams = Join[Constants`SIParams[T, (("m" "\[Epsilon]0" "\[Omega]i"^2) / ("e"^2) /. fitrepl 
            /. consts)], fitrepl]; (*this extracts the number density from the fitted
             plasmon frequency*)
        lL = zMIntegralFitEnhanced[m\[Chi], v\[Chi], fitparams];
        lPrefactor = "Ai" / "JpereV" (("qF" "vF") / v\[Chi])^2 8 / Pi "e"^2 / (4 \[Pi] "\[Epsilon]0") /. fitparams /. fitrepl;
        ldEdr = lPrefactor Total[lL];
        Print["For m\[Chi] = ", m\[Chi], " , v\[Chi] = ", v\[Chi], "\nL:\t", Total[lL] "\nPrefactor:\t",
             lPrefactor, "\ndEdr:\t", ldEdr];
        <|"m\[Chi]" -> m\[Chi], "v\[Chi]" -> v\[Chi], "fitparams" -> fitparams, "fitrepl"
             -> fitrepl, "L" -> lL, "Prefactor" -> lPrefactor, "dEdr" -> ldEdr|>
    ]*)

EnergyLossMSIFitEnhanced[m\[Chi]_, v\[Chi]_, fitparams_,l_:1] :=
    Module[
        {lL, lPrefactor, lthEmoment},
        lL = zMIntegralFitEnhanced[m\[Chi], v\[Chi], fitparams,l];
        (*lPrefactor = "Ai" / "JpereV" (("qF" "vF") / v\[Chi])^2 8 / Pi "e"^2 / (4 \[Pi] "\[Epsilon]0") /. fitparams;*)
        lPrefactor = ("\[HBar]")^(Abs[l]-1) ("Ai" / ("JpereV")^Abs[l]) (("qF" "vF")^(Abs[l]+1)/v\[Chi]^2 (2^(Abs[l]+2) / Pi) ("e"^2 / (4 \[Pi] "\[Epsilon]0"))) /. fitparams;
        lthEmoment = lPrefactor Total[lL];
        Print["For m\[Chi] = ", m\[Chi], " , v\[Chi] = ", v\[Chi], "\nL:\t", Total[lL] "\nPrefactor:\t",
             lPrefactor, StringForm["\n \!\(\*SuperscriptBox[\(``\), \(th\)]\) E moment:\t",l], lthEmoment];
        <|"m\[Chi]" -> m\[Chi], "v\[Chi]" -> v\[Chi], "fitparams" -> fitparams, "L" -> lL,
             "Prefactor" -> lPrefactor, "dEdr" -> lthEmoment|>
    ]
    
   (* EnergyLossMSIFitEnhanced[m\[Chi]_, v\[Chi]_, fitparams_,l_:1] :=
    Module[
        {lL, lPrefactor, ldEdr},
        lL = zMIntegralFitEnhanced[m\[Chi], v\[Chi], fitparams,l];
        (*lPrefactor = "Ai" / "JpereV" (("qF" "vF") / v\[Chi])^2 8 / Pi "e"^2 / (4 \[Pi] "\[Epsilon]0") /. fitparams;*)
        lPrefactor = ("\[HBar]")^(l-1) ("Ai" / ("JpereV")^l) (("qF" "vF")^(l+1)/v\[Chi]^2 (2^(l+2) / Pi) ("e"^2 / (4 \[Pi] "\[Epsilon]0"))) /. fitparams;
        ldEdr = lPrefactor Total[lL];
        Print["For m\[Chi] = ", m\[Chi], " , v\[Chi] = ", v\[Chi], "\nL:\t", Total[lL] "\nPrefactor:\t",
             lPrefactor, "\ndEdr:\t", ldEdr];
        <|"m\[Chi]" -> m\[Chi], "v\[Chi]" -> v\[Chi], "fitparams" -> fitparams, "L" -> lL,
             "Prefactor" -> lPrefactor, "dEdr" -> ldEdr|>
    ]*)


(* ::Subsection:: *)
(*EL Sum over Oscillators - Public*)


(* ::Subsubsection::Closed:: *)
(*Without Temperature Enhancement*)


(*EnergyLossMSIFitSumOscillators[m\[Chi]_, v\[Chi]_, fitreplNested_, consts_, T_] :=
    Module[{ELlist},
        ELlist = Table[EnergyLossMSIFit[m\[Chi], v\[Chi], fitreplNested[[i]], consts,T]["dEdr"], {i, Length[fitreplNested]}];
        Print["dEdr total: ", Total[ELlist]];
        <|"m\[Chi]" -> m\[Chi], "v\[Chi]" -> v\[Chi], "fitparams" -> fitreplNested, "consts"
             -> consts, "T" -> T, "ELlist" -> ELlist, "dEdr" -> Total[ELlist]|>
    ]

EnergyLossMSIFitSumOscillators[m\[Chi]_, v\[Chi]_, fitparams_] :=
    Module[{ELlist},
        ELlist = Table[EnergyLossMSIFit[m\[Chi], v\[Chi], fitparams[[i]]]["dEdr"], {i, Length[fitparams]}];
        Print["dEdr total: ", Total[ELlist]];
        <|"m\[Chi]" -> m\[Chi], "v\[Chi]" -> v\[Chi], "fitparams" -> fitparams, "ELlist" 
            -> ELlist, "dEdr" -> Total[ELlist]|>
    ]*)


(* ::Subsubsection::Closed:: *)
(*With Temperature Enhancement*)


(*EnergyLossMSIFitSumOscillatorsEnhanced[m\[Chi]_, v\[Chi]_, fitreplNested_, consts_, T_] :=
    Module[{ELlist},
        ELlist = Table[EnergyLossMSIFitEnhanced[m\[Chi], v\[Chi], fitreplNested[[i]], consts,T]["dEdr"], {i, Length[fitreplNested]}];
        Print["dEdr total: ", Total[ELlist]];
        <|"m\[Chi]" -> m\[Chi], "v\[Chi]" -> v\[Chi], "fitparams" -> fitreplNested, "consts"
             -> consts, "T" -> T, "ELlist" -> ELlist, "dEdr" -> Total[ELlist]|>
    ]*)

EnergyLossMSIFitSumOscillatorsEnhanced[m\[Chi]_, v\[Chi]_, fitparams_,l_:1] :=
    Module[{ELlist},
        ELlist = Table[EnergyLossMSIFitEnhanced[m\[Chi], v\[Chi], fitparams[[i]],l]["dEdr"], {i, Length[fitparams]}];
        Print["dEdr total: ", Total[ELlist]];
        <|"m\[Chi]" -> m\[Chi], "v\[Chi]" -> v\[Chi], "fitparams" -> fitparams, "ELlist" 
            -> ELlist, "dEdr" -> Total[ELlist]|>
    ]


(* ::Subsection:: *)
(*EL Interpolation and Table - Public*)


(* ::Subsubsection::Closed:: *)
(*Using Theoretical Collision Rate*)


EnergyLossTableAndInter[m\[Chi]_ : {5 10^5, 10^9}, v\[Chi]_ : {10 ^ -4, 10 ^ -2
    }, meshdims_ : <|"m\[Chi]" -> 5, "v\[Chi]" -> 5|>, \[Epsilon]_, params_, save_:True, fname_
    :"outputdict"] :=
    Module[
        {kinlist, EnergyLossTable, m\[Chi]Mesh, v\[Chi]Mesh, EnergyLossMesh, InterpolationTable,
             LogInterpolation, outputdict}
        ,
(*
	  m\[Chi] - [eV] max and min of masses in grid
	v\[Chi] - [c] max and min of velocities in grid
	meshdims - size of equally (in log space) spaced grids of m\[Chi] and v\[Chi] 
    to compute energy loss over
	
	Output energy loss is in [eV m^-1]
	*)
        kinlist = Kinematics[m\[Chi], v\[Chi], meshdims];
        EnergyLossTable = {};
        Do[AppendTo[EnergyLossTable, EnergyLossMSI[(kinlist[["m\[Chi]"]])[[
            i]], (kinlist[["v\[Chi]"]])[[j]], \[Epsilon], params]], {i, meshdims[["m\[Chi]"]]}, {j, 
            meshdims[["v\[Chi]"]]}];
        (*SaveIt[NotebookDirectory[]<>fname,EnergyLossTable];*)
        (*EnergyLossTable = ReadIt[NotebookDirectory[]<>fname];*)
        m\[Chi]Mesh = ArrayReshape[Table[EnergyLossTable[[i]][["m\[Chi]"]], {i,
             Length[EnergyLossTable]}], {meshdims[["m\[Chi]"]], meshdims[["v\[Chi]"]]}];
        v\[Chi]Mesh = ArrayReshape[Table[EnergyLossTable[[i]][["v\[Chi]"]], {i,
             Length[EnergyLossTable]}], {meshdims[["m\[Chi]"]], meshdims[["v\[Chi]"]]}];
        EnergyLossMesh = ArrayReshape[Table[EnergyLossTable[[i]][["dEdr"
            ]], {i, Length[EnergyLossTable]}], {meshdims[["m\[Chi]"]], meshdims[["v\[Chi]"]]
            }];
        InterpolationTable = Table[{{m\[Chi]Mesh[[i, j]], v\[Chi]Mesh[[i, j]]},
             EnergyLossMesh[[i, j]]}, {i, meshdims[["m\[Chi]"]]}, {j, meshdims[["v\[Chi]"]]
            }];
        LogInterpolation = Interpolation[Flatten[Log[10, InterpolationTable
            ], 1]];
        outputdict = <|"f" -> LogInterpolation, "m\[Chi]Mesh" -> m\[Chi]Mesh, "v\[Chi]Mesh"
             -> v\[Chi]Mesh, "EnergyLossMesh" -> EnergyLossMesh, "params" -> params, "fname"
             -> fname, "kinlist" -> kinlist, "meshdims" -> meshdims, "\[Epsilon]" -> ToString[
            \[Epsilon]]|>;
        If[save,
            Utilities`SaveIt[NotebookDirectory[] <> fname, outputdict]
        ];
        outputdict
    ]


(* ::Subsubsection::Closed:: *)
(*Using Optical Fits Parameters*)


(*EnergyLossTableAndInterFIT[m\[Chi]_ : {5 10^5, 10^9}, v\[Chi]_ : {10 ^ -4, 10 ^
     -2}, meshdims_ : <|"m\[Chi]" -> 5, "v\[Chi]" -> 5|>, fitreplNested_, consts_, 
    T_, save_:True, fname_:"outputdict"] :=
    Module[
        {kinlist, EnergyLossTable, m\[Chi]Mesh, v\[Chi]Mesh, EnergyLossMesh, InterpolationTable,
             LogInterpolation, outputdict, count, times, timesMesh}
        ,
(*
	  m\[Chi] - [kg] max and min of masses in grid
	v\[Chi] - [m s^-1] max and min of velocities in grid
	meshdims - size of equally (in log space) spaced grids of m\[Chi] and v\[Chi] 
    
    to compute energy loss over
	
	Output energy loss is in [eV m^-1]
	*)
        kinlist = Kinematics[m\[Chi], v\[Chi], meshdims];
        count = 0;
        EnergyLossTable = {};
        times = {};
        Do[
            count++;
            Print[count, " of ", meshdims[["m\[Chi]"]] meshdims[["v\[Chi]"]]];
            AppendTo[times, AbsoluteTiming[AppendTo[EnergyLossTable, 
                EnergyLossMSIFitSumOscillators[(kinlist[["m\[Chi]"]])[[i]], (kinlist[["v\[Chi]"
                ]])[[j]], fitreplNested, consts, T]]][[1]]]
            ,
            {i, meshdims[["m\[Chi]"]]}
            ,
            {j, meshdims[["v\[Chi]"]]}
        ];
        m\[Chi]Mesh = ArrayReshape[Table[EnergyLossTable[[i]][["m\[Chi]"]], {i,
             Length[EnergyLossTable]}], {meshdims[["m\[Chi]"]], meshdims[["v\[Chi]"]]}];
        v\[Chi]Mesh = ArrayReshape[Table[EnergyLossTable[[i]][["v\[Chi]"]], {i,
             Length[EnergyLossTable]}], {meshdims[["m\[Chi]"]], meshdims[["v\[Chi]"]]}];
        EnergyLossMesh = ArrayReshape[Table[EnergyLossTable[[i]][["dEdr"
            ]], {i, Length[EnergyLossTable]}], {meshdims[["m\[Chi]"]], meshdims[["v\[Chi]"]]
            }];
        timesMesh = ArrayReshape[times, {meshdims[["m\[Chi]"]], meshdims[[
            "v\[Chi]"]]}];
        InterpolationTable = Table[{{m\[Chi]Mesh[[i, j]], v\[Chi]Mesh[[i, j]]},
             EnergyLossMesh[[i, j]]}, {i, meshdims[["m\[Chi]"]]}, {j, meshdims[["v\[Chi]"]]
            }];
        LogInterpolation = Interpolation[Flatten[Log[10, InterpolationTable
            ], 1]];
        outputdict = <|"f" -> LogInterpolation, "m\[Chi]Mesh" -> m\[Chi]Mesh, "v\[Chi]Mesh"
             -> v\[Chi]Mesh, "EnergyLossMesh" -> EnergyLossMesh, "consts" -> consts, "fitparams"
             -> fitreplNested, "T" -> T, "fname" -> fname, "kinlist" -> kinlist, 
            "meshdims" -> meshdims, "timesMesh" -> timesMesh, "totalTime" -> Total[
            timesMesh, 2]|>;
        If[save,
            Utilities`SaveIt[NotebookDirectory[] <> fname, outputdict]
        ];
        outputdict
    ]
*)


(*EnergyLossTableAndInterFITTotalParams[m\[Chi]_ : {5 10^5, 10^9}, v\[Chi]_ : {10 ^ -4, 10 ^
     -2}, meshdims_ : <|"m\[Chi]" -> 5, "v\[Chi]" -> 5|>, fitreplNested_, save_:True,
     fname_:"outputdict"] :=
    Module[
        {kinlist, EnergyLossTable, m\[Chi]Mesh, v\[Chi]Mesh, EnergyLossMesh, InterpolationTable,
             LogInterpolation, outputdict, count, times, timesMesh}
        ,
(*
	  m\[Chi] - [kg] max and min of masses in grid
	v\[Chi] - [m s^-1] max and min of velocities in grid
	meshdims - size of equally (in log space) spaced grids of m\[Chi] and v\[Chi] 
    
    to compute energy loss over
	
	Output energy loss is in [eV m^-1]
	*)
        kinlist = Kinematics[m\[Chi], v\[Chi], meshdims];
        count = 0;
        EnergyLossTable = {};
        times = {};
        Do[
            count++;
            Print[count, " of ", meshdims[["m\[Chi]"]] meshdims[["v\[Chi]"]]];
            AppendTo[times, AbsoluteTiming[AppendTo[EnergyLossTable, 
                EnergyLossMSIFitSumOscillators[(kinlist[["m\[Chi]"]])[[i]], (kinlist[["v\[Chi]"
                ]])[[j]], fitreplNested]]][[1]]]
            ,
            {i, meshdims[["m\[Chi]"]]}
            ,
            {j, meshdims[["v\[Chi]"]]}
        ];
        m\[Chi]Mesh = ArrayReshape[Table[EnergyLossTable[[i]][["m\[Chi]"]], {i,
             Length[EnergyLossTable]}], {meshdims[["m\[Chi]"]], meshdims[["v\[Chi]"]]}];
        v\[Chi]Mesh = ArrayReshape[Table[EnergyLossTable[[i]][["v\[Chi]"]], {i,
             Length[EnergyLossTable]}], {meshdims[["m\[Chi]"]], meshdims[["v\[Chi]"]]}];
        EnergyLossMesh = ArrayReshape[Table[EnergyLossTable[[i]][["dEdr"
            ]], {i, Length[EnergyLossTable]}], {meshdims[["m\[Chi]"]], meshdims[["v\[Chi]"]]
            }];
        timesMesh = ArrayReshape[times, {meshdims[["m\[Chi]"]], meshdims[[
            "v\[Chi]"]]}];
        InterpolationTable = Table[{{m\[Chi]Mesh[[i, j]], v\[Chi]Mesh[[i, j]]},
             EnergyLossMesh[[i, j]]}, {i, meshdims[["m\[Chi]"]]}, {j, meshdims[["v\[Chi]"]]
            }];
        LogInterpolation = Interpolation[Flatten[Log[10, InterpolationTable
            ], 1]];
        outputdict = <|"f" -> LogInterpolation, "m\[Chi]Mesh" -> m\[Chi]Mesh, "v\[Chi]Mesh"
             -> v\[Chi]Mesh, "EnergyLossMesh" -> EnergyLossMesh, "fitparams" -> fitreplNested,
             "fname" -> fname, "kinlist" -> kinlist, "meshdims" -> meshdims, "timesMesh"
             -> timesMesh, "totalTime" -> Total[timesMesh, 2]|>;
        If[save,
            Utilities`SaveIt[NotebookDirectory[] <> fname, outputdict]
        ];
        outputdict
    ]*)


(* ::Subsubsection::Closed:: *)
(*Using Optical Fits Parameters - With Temperature Enhancement*)


(*CORRECT ONE WHICH WORKS - DOESN'T GIVE INDIVIDUAL OSCILLATOR INFO*)
EnergyLossTableAndInterFITTotalParamsEnhanced[m\[Chi]_ : {5 10^5, 10^9}, v\[Chi]_ : {10 ^ -4, 10 ^
     -2}, meshdims_ : <|"m\[Chi]" -> 5, "v\[Chi]" -> 5|>, fitreplNested_, save_:True,
     fname_:"outputdict",l_:1] :=
    Module[
        {kinlist, EnergyLossTable, m\[Chi]Mesh, v\[Chi]Mesh, EnergyLossMesh, InterpolationTable,
             LogInterpolation, outputdict, count, times, timesMesh}
        ,
(*
	  m\[Chi] - [kg] max and min of masses in grid
	v\[Chi] - [m s^-1] max and min of velocities in grid
	meshdims - size of equally (in log space) spaced grids of m\[Chi] and v\[Chi] 
    l - [] which moment of interaction per unit length (1/vdP/Subscript[dE, R])to calculate 
           (l=1 gives Energy Loss per unit length [eV m^-1]], 
            l=0 gives inverse interaction length [m^-1])
    
    to compute energy loss over
	
	Output energy loss is in [eV m^-1]
	*)
        kinlist = Kinematics[m\[Chi], v\[Chi], meshdims];
        count = 0;
        EnergyLossTable = {};
        times = {};
        Do[
            count++;
            Print[count, " of ", meshdims[["m\[Chi]"]] meshdims[["v\[Chi]"]]];
            AppendTo[times, AbsoluteTiming[AppendTo[EnergyLossTable, 
                EnergyLossMSIFitSumOscillatorsEnhanced[(kinlist[["m\[Chi]"]])[[i]], (kinlist[["v\[Chi]"
                ]])[[j]], fitreplNested,l]]][[1]]]
            ,
            {i, meshdims[["m\[Chi]"]]}
            ,
            {j, meshdims[["v\[Chi]"]]}
        ];
        m\[Chi]Mesh = ArrayReshape[Table[EnergyLossTable[[i]][["m\[Chi]"]], {i,
             Length[EnergyLossTable]}], {meshdims[["m\[Chi]"]], meshdims[["v\[Chi]"]]}];
        v\[Chi]Mesh = ArrayReshape[Table[EnergyLossTable[[i]][["v\[Chi]"]], {i,
             Length[EnergyLossTable]}], {meshdims[["m\[Chi]"]], meshdims[["v\[Chi]"]]}];
        EnergyLossMesh = ArrayReshape[Table[EnergyLossTable[[i]][["dEdr"
            ]], {i, Length[EnergyLossTable]}], {meshdims[["m\[Chi]"]], meshdims[["v\[Chi]"]]
            }];
        timesMesh = ArrayReshape[times, {meshdims[["m\[Chi]"]], meshdims[[
            "v\[Chi]"]]}];
        InterpolationTable = Table[{{m\[Chi]Mesh[[i, j]], v\[Chi]Mesh[[i, j]]},
             EnergyLossMesh[[i, j]]}, {i, meshdims[["m\[Chi]"]]}, {j, meshdims[["v\[Chi]"]]
            }];
        LogInterpolation = Interpolation[Flatten[Log[10, InterpolationTable
            ], 1]];
        outputdict = <|"f" -> LogInterpolation, "m\[Chi]Mesh" -> m\[Chi]Mesh, "v\[Chi]Mesh"
             -> v\[Chi]Mesh, "EnergyLossMesh" -> EnergyLossMesh, "fitparams" -> fitreplNested,
             "fname" -> fname, "kinlist" -> kinlist, "meshdims" -> meshdims, "timesMesh"
             -> timesMesh, "totalTime" -> Total[timesMesh, 2]|>;
        If[save,
            Utilities`SaveIt[NotebookDirectory[] <> fname, outputdict]
        ];
        outputdict
    ]
    
(*MODIFYING TO GET INDIVIDUAL OSCILLATOR FITS*)
(*EnergyLossTableAndInterFITTotalParamsEnhancedOscillators[m\[Chi]_ : {5 10^5, 10^9}, v\[Chi]_ : {10 ^ -4, 10 ^
     -2}, meshdims_ : <|"m\[Chi]" -> 5, "v\[Chi]" -> 5|>, fitreplNested_, save_:True,
     fname_:"outputdict",l_:1] :=
    Module[
        {kinlist, EnergyLossTable, m\[Chi]Mesh, v\[Chi]Mesh, EnergyLossMesh, InterpolationTable,
             LogInterpolation, EnergyLossMeshELlist, LogInterpolationELlist, outputdict,
              count, times, timesMesh}
        ,
(*
	  m\[Chi] - [kg] max and min of masses in grid
	v\[Chi] - [m s^-1] max and min of velocities in grid
	meshdims - size of equally (in log space) spaced grids of m\[Chi] and v\[Chi] to compute energy loss / interaction length over
    l - [] which moment of interaction per unit length (1/vdP/Subscript[dE, R])to calculate 
           (l=1 gives Energy Loss per unit length [eV m^-1]], 
            l=0 gives inverse interaction length [m^-1])
	
	Output is in [eV m^-1] or [m^-1], Oscillator suffix is for individual oscillators in the fit
		EnergyLossMesh(Oscillators) - table of energy loss / interaction length over the kinematic mesh
		f(Oscillators) - Log Log interpolation of EnergyLossMesh
	*)
        kinlist = Kinematics[m\[Chi], v\[Chi], meshdims];
        count = 0;
        EnergyLossTable = {};
        times = {};
        Do[
            count++;
            Print[count, " of ", meshdims[["m\[Chi]"]] meshdims[["v\[Chi]"]]];
            AppendTo[times, AbsoluteTiming[AppendTo[EnergyLossTable, 
                EnergyLossMSIFitSumOscillatorsEnhanced[(kinlist[["m\[Chi]"]])[[i]], (kinlist[["v\[Chi]"
                ]])[[j]], fitreplNested,l]]][[1]]]
            ,
            {i, meshdims[["m\[Chi]"]]}
            ,
            {j, meshdims[["v\[Chi]"]]}
        ];
        
        
        (*Save Times*)
        timesMesh = ArrayReshape[times, {meshdims[["m\[Chi]"]], meshdims[[
            "v\[Chi]"]]}];
        
        (*Compute kinematic meshes*)
        m\[Chi]Mesh = ArrayReshape[Table[EnergyLossTable[[i]][["m\[Chi]"]], {i,
             Length[EnergyLossTable]}], {meshdims[["m\[Chi]"]], meshdims[["v\[Chi]"]]}];
        v\[Chi]Mesh = ArrayReshape[Table[EnergyLossTable[[i]][["v\[Chi]"]], {i,
             Length[EnergyLossTable]}], {meshdims[["m\[Chi]"]], meshdims[["v\[Chi]"]]}];
        
        (*Compute total *)
        EnergyLossMesh = ArrayReshape[Table[EnergyLossTable[[i]][["dEdr"
            ]], {i, Length[EnergyLossTable]}], {meshdims[["m\[Chi]"]], meshdims[["v\[Chi]"]]
            }];
        LogInterpolation=Utilities`InterpolateTable[m\[Chi]Mesh,v\[Chi]Mesh,meshdims,EnergyLossMesh];
        
        (*Compute Individual Oscillators*)
        EnergyLossMeshELlist = Table[ArrayReshape[Table[EnergyLossTable[[i]][["ELlist"
            ]][[j]], {i, Length[EnergyLossTable]}], {meshdims[["m\[Chi]"]], meshdims[["v\[Chi]"]]
            }],{j,Length[fitreplNested]}];    
        LogInterpolationELlist=Table[Utilities`InterpolateTable[m\[Chi]Mesh,v\[Chi]Mesh,meshdims,EnergyLossMeshELlist[[j]]],{j,Length[fitreplNested]}];
        
        
        outputdict = <|"f" -> LogInterpolation, "fOscillators"->LogInterpolationELlist, "m\[Chi]Mesh" -> m\[Chi]Mesh, "v\[Chi]Mesh"
             -> v\[Chi]Mesh, "EnergyLossMesh" -> EnergyLossMesh, "EnergyLossMeshOscillators"->EnergyLossMeshELlist, 
             "fitparams" -> fitreplNested, "fname" -> fname, "kinlist" -> kinlist, "meshdims" -> meshdims, "timesMesh"
             -> timesMesh, "totalTime" -> Total[timesMesh, 2]|>;
        If[save,
            Utilities`SaveIt[NotebookDirectory[] <> fname, outputdict]
        ];
        outputdict
    ]*)


(*EnergyLossTableAndInterFITTotalParamsEnhanced[m\[Chi]_ : {5 10^5, 10^9}, v\[Chi]_ : {10 ^ -4, 10 ^
     -2}, meshdims_ : <|"m\[Chi]" -> 5, "v\[Chi]" -> 5|>, fitreplNested_, save_:True,
     fname_:"outputdict",\[Beta]_] :=
    Module[
        {kinlist, EnergyLossTable, m\[Chi]Mesh, v\[Chi]Mesh, EnergyLossMesh, InterpolationTable,
             LogInterpolation, outputdict, count, times, timesMesh}
        ,
(*
	  m\[Chi] - [kg] max and min of masses in grid
	v\[Chi] - [m s^-1] max and min of velocities in grid
	meshdims - size of equally (in log space) spaced grids of m\[Chi] and v\[Chi] 
    
    to compute energy loss over
	
	Output energy loss is in [eV m^-1]
	*)
        kinlist = Kinematics[m\[Chi], v\[Chi], meshdims];
        count = 0;
        EnergyLossTable = {};
        times = {};
        (*Do[fitreplNested[[i]]=Utilities`ReplaceParams[fitreplNested[[i]],\[Beta],"\[Beta]"],{i,Length[fitreplNested]}];*)
        Do[
            count++;
            Print[count, " of ", meshdims[["m\[Chi]"]] meshdims[["v\[Chi]"]]];
            AppendTo[times, AbsoluteTiming[AppendTo[EnergyLossTable, 
                EnergyLossMSIFitSumOscillatorsEnhanced[(kinlist[["m\[Chi]"]])[[i]], (kinlist[["v\[Chi]"
                ]])[[j]], fitreplNested]]][[1]]]
            ,
            {i, meshdims[["m\[Chi]"]]}
            ,
            {j, meshdims[["v\[Chi]"]]}
        ];
        m\[Chi]Mesh = ArrayReshape[Table[EnergyLossTable[[i]][["m\[Chi]"]], {i,
             Length[EnergyLossTable]}], {meshdims[["m\[Chi]"]], meshdims[["v\[Chi]"]]}];
        v\[Chi]Mesh = ArrayReshape[Table[EnergyLossTable[[i]][["v\[Chi]"]], {i,
             Length[EnergyLossTable]}], {meshdims[["m\[Chi]"]], meshdims[["v\[Chi]"]]}];
        EnergyLossMesh = ArrayReshape[Table[EnergyLossTable[[i]][["dEdr"
            ]], {i, Length[EnergyLossTable]}], {meshdims[["m\[Chi]"]], meshdims[["v\[Chi]"]]
            }];
        timesMesh = ArrayReshape[times, {meshdims[["m\[Chi]"]], meshdims[[
            "v\[Chi]"]]}];
        InterpolationTable = Table[{{m\[Chi]Mesh[[i, j]], v\[Chi]Mesh[[i, j]]},
             EnergyLossMesh[[i, j]]}, {i, meshdims[["m\[Chi]"]]}, {j, meshdims[["v\[Chi]"]]
            }];
        LogInterpolation = Interpolation[Flatten[Log[10, InterpolationTable
            ], 1]];
        outputdict = <|"f" -> LogInterpolation, "m\[Chi]Mesh" -> m\[Chi]Mesh, "v\[Chi]Mesh"
             -> v\[Chi]Mesh, "EnergyLossMesh" -> EnergyLossMesh, "fitparams" -> fitreplNested,
             "fname" -> fname, "kinlist" -> kinlist, "meshdims" -> meshdims, "timesMesh"
             -> timesMesh, "totalTime" -> Total[timesMesh, 2]|>;
        If[save,
            Utilities`SaveIt[NotebookDirectory[] <> fname, outputdict]
        ];
        outputdict
    ]*)


(* ::Section::Closed:: *)
(*Nuclear*)


(* ::Subsection::Closed:: *)
(*EL Integrals - Nuclear - Finite Temperature*)


(* ::Subsubsection::Closed:: *)
(*Coulomb potential*)


VCoul[q_,params_]:= ((("e")^2)/(q^2 "\[Epsilon]0"))/.params (*no kappa included here*)


(* ::Subsubsection::Closed:: *)
(*Nuclear Structure Function - Finite Temperature*)


(* ::Text:: *)
(*Below coeffs corresponds to the replacement table of coefficients of the nuclear and atomic form factors, as well as constants describing the atomic species*)


pb[sign_,q_,\[Omega]_]:= ("mN" \[Omega])/(q "\[HBar]") + <|"+"->1,"-"->(-1)|>[[sign]] q/2 (*centre of mass momentum integral bounds*)
BoltzmannExp[p_,\[Beta]_]:= Exp[-\[Beta] (("\[HBar]")^2 p^2)/(2 "mN")]


SNuc[q_,\[Omega]_,\[Beta]_,coeffs_]:= -(FormFactors`Zeff[q,coeffs]^2/(2 \[Pi] q "\[HBar]")) (2 \[Pi])^(3/2) Sqrt["mN" \[Beta]] (BoltzmannExp[-pb["+",q,\[Omega]],\[Beta]]-BoltzmannExp[-pb["-",q,\[Omega]],\[Beta]])/. coeffs


(* ::Subsubsection::Closed:: *)
(*Nuclear Average Energy Loss Integral - Finite Temperature*)


(* ::Text:: *)
(*Below params correspond to the physical constants replacement table*)


\[Omega]b[sign_,m\[Chi]_,v\[Chi]_,q_]:= q v\[Chi] + <|"+"->1,"-"->(-1)|>[[sign]] (("\[HBar]" q^2)/(2 m\[Chi]))(*\[Omega] integral bounds*)


\[Omega]IntegralNuc[q_?NumberQ,m\[Chi]_?NumberQ,v\[Chi]_?NumberQ,\[Beta]_,coeffs_,params_]:= NIntegrate[\[Omega] SNuc[q,\[Omega],\[Beta],coeffs]/.params,{\[Omega],-\[Omega]b["+",m\[Chi],v\[Chi],q]/.params,\[Omega]b["-",m\[Chi],v\[Chi],q]/.params}]


ELIntegralNuc[m\[Chi]_?NumberQ,v\[Chi]_?NumberQ,coeffs_,params_,statpars_]:= ("nI"/.statpars) NIntegrate[q/v\[Chi]^2 VCoul[q,params]^2/(2 \[Pi])^2 \[Omega]IntegralNuc[q,m\[Chi],v\[Chi],("\[Beta]"/.statpars),coeffs,params],{q,0,\[Infinity]}] (*doesn't include the \[Kappa]^2*)


(* ::Subsubsection::Closed:: *)
(*EL Integrals - Nuclear 0 Temperature*)


ELIntegralNuc0[m\[Chi]_?NumberQ,v\[Chi]_?NumberQ,coeffs_,params_,statpars_:<|"nI"->10^23|>]:= ("nI" /.statpars)/(2 \[Pi] v\[Chi]^2) NIntegrate[q^3/(2("mN"/. coeffs)) VCoul[q,params]^2 FormFactors`Zeff[q,coeffs]^2,{q,0,(2 (("mN")^-1+m\[Chi]^-1)^-1 v\[Chi])/("\[HBar]")/. params /.coeffs}]


(* ::Subsection::Closed:: *)
(*EL Module - Nuclear*)


EnergyLossNuc0[m\[Chi]_, v\[Chi]_,coeffs_,params_,statpars_:<|"nI"->10^23|>,print_:True] :=
    Module[{ldEdr},
    (*output in [eV m^-1]*)
        ldEdr = (1/("JpereV")/.params) ELIntegralNuc0[m\[Chi], v\[Chi],coeffs,params,statpars];
        If[print,
           Print["For m\[Chi] = ", N[m\[Chi]], " , v\[Chi] = ", v\[Chi], "\ndEdr:\t", ldEdr];
        ];
        <|"m\[Chi]" -> N[m\[Chi]], "v\[Chi]" -> v\[Chi], "params" -> params,"coeffs"->coeffs,"statpars"->statpars,"dEdr" -> ldEdr|>
    ]


EnergyLossNuc[m\[Chi]_, v\[Chi]_,coeffs_,params_,statpars_,print_:True] :=
    Module[{ldEdr},
    (*output in [eV m^-1]*)
        ldEdr = (1/("JpereV")/.params) ELIntegralNuc[m\[Chi],v\[Chi],coeffs,params,statpars];
        (*Print["For m\[Chi] = ", N[m\[Chi]], " , v\[Chi] = ", v\[Chi], "\ndEdr:\t", ldEdr];*)
        If[print,
           Print["For m\[Chi] = ", N[m\[Chi]], " , v\[Chi] = ", v\[Chi], "\ndEdr:\t", ldEdr];
        ];
        <|"m\[Chi]" -> N[m\[Chi]], "v\[Chi]" -> v\[Chi], "params" -> params,"coeffs"->coeffs,"statpars"->statpars,"dEdr" -> ldEdr|>
    ]


(* ::Subsection::Closed:: *)
(*EL Interpolation and Table - Nuclear*)


(* ::Subsubsection::Closed:: *)
(*0 Temperature*)


EnergyLossTableAndInterNuc0[m\[Chi]_ : {5 10^5, 10^9}, v\[Chi]_ : {10^-4, 10^-2
    }, meshdims_ : <|"m\[Chi]" -> 5, "v\[Chi]" -> 5|>, coeffs_,params_,statpars_, save_:True, fname_
    :"outputdict",label_,print_] :=
    Module[
        {kinlist, EnergyLossTable, m\[Chi]Mesh, v\[Chi]Mesh, EnergyLossMesh, InterpolationTable,
             LogInterpolation, outputdict}
        ,
(*
	  m\[Chi] - [eV] max and min of masses in grid
	v\[Chi] - [c] max and min of velocities in grid
	meshdims - size of equally (in log space) spaced grids of m\[Chi] and v\[Chi] 
    to compute energy loss over
	
	Output energy loss is in [eV m^-1]
	*)
        kinlist = Kinematics[m\[Chi], v\[Chi], meshdims];
        EnergyLossTable = {};
        Do[AppendTo[EnergyLossTable, EnergyLossNuc0[(kinlist[["m\[Chi]"]])[[
            i]], (kinlist[["v\[Chi]"]])[[j]],coeffs,params,statpars,print]], {i, meshdims[["m\[Chi]"]]}, {j, 
            meshdims[["v\[Chi]"]]}];
        (*SaveIt[NotebookDirectory[]<>fname,EnergyLossTable];*)
        (*EnergyLossTable = ReadIt[NotebookDirectory[]<>fname];*)
        m\[Chi]Mesh = ArrayReshape[Table[EnergyLossTable[[i]][["m\[Chi]"]], {i,
             Length[EnergyLossTable]}], {meshdims[["m\[Chi]"]], meshdims[["v\[Chi]"]]}];
        v\[Chi]Mesh = ArrayReshape[Table[EnergyLossTable[[i]][["v\[Chi]"]], {i,
             Length[EnergyLossTable]}], {meshdims[["m\[Chi]"]], meshdims[["v\[Chi]"]]}];
        EnergyLossMesh = ArrayReshape[Table[EnergyLossTable[[i]][["dEdr"
            ]], {i, Length[EnergyLossTable]}], {meshdims[["m\[Chi]"]], meshdims[["v\[Chi]"]]
            }];
        InterpolationTable = Table[{{m\[Chi]Mesh[[i, j]], v\[Chi]Mesh[[i, j]]},
             EnergyLossMesh[[i, j]]}, {i, meshdims[["m\[Chi]"]]}, {j, meshdims[["v\[Chi]"]]
            }];
        LogInterpolation = Interpolation[Flatten[Log[10, InterpolationTable
            ], 1]];
        outputdict = <|"f" -> LogInterpolation, "m\[Chi]Mesh" -> m\[Chi]Mesh, "v\[Chi]Mesh"
             -> v\[Chi]Mesh, "EnergyLossMesh" -> EnergyLossMesh, "params" -> params, "fname"
             -> fname, "kinlist" -> kinlist, "meshdims" -> meshdims, "\[Epsilon]" -> ToString[
            label]|>;
        If[save,
            Utilities`SaveIt[NotebookDirectory[] <> fname, outputdict]
        ];
        outputdict
    ]


(* ::Subsubsection::Closed:: *)
(*Finite Temperature*)


EnergyLossTableAndInterNuc[m\[Chi]_ : {5 10^5, 10^9}, v\[Chi]_ : {10^-4, 10^-2
    }, meshdims_ : <|"m\[Chi]" -> 5, "v\[Chi]" -> 5|>, coeffs_,params_,statpars_, save_:True, fname_
    :"outputdict",label_,print_] :=
    Module[
        {kinlist, EnergyLossTable, m\[Chi]Mesh, v\[Chi]Mesh, EnergyLossMesh, InterpolationTable,
             LogInterpolation, outputdict}
        ,
(*
	  m\[Chi] - [eV] max and min of masses in grid
	v\[Chi] - [c] max and min of velocities in grid
	meshdims - size of equally (in log space) spaced grids of m\[Chi] and v\[Chi] 
    to compute energy loss over
	
	Output energy loss is in [eV m^-1]
	*)
        kinlist = Kinematics[m\[Chi], v\[Chi], meshdims];
        EnergyLossTable = {};
        Do[AppendTo[EnergyLossTable, EnergyLossNuc[(kinlist[["m\[Chi]"]])[[
            i]], (kinlist[["v\[Chi]"]])[[j]],coeffs,params,statpars,print]], {i, meshdims[["m\[Chi]"]]}, {j, 
            meshdims[["v\[Chi]"]]}];
        (*SaveIt[NotebookDirectory[]<>fname,EnergyLossTable];*)
        (*EnergyLossTable = ReadIt[NotebookDirectory[]<>fname];*)
        m\[Chi]Mesh = ArrayReshape[Table[EnergyLossTable[[i]][["m\[Chi]"]], {i,
             Length[EnergyLossTable]}], {meshdims[["m\[Chi]"]], meshdims[["v\[Chi]"]]}];
        v\[Chi]Mesh = ArrayReshape[Table[EnergyLossTable[[i]][["v\[Chi]"]], {i,
             Length[EnergyLossTable]}], {meshdims[["m\[Chi]"]], meshdims[["v\[Chi]"]]}];
        EnergyLossMesh = ArrayReshape[Table[EnergyLossTable[[i]][["dEdr"
            ]], {i, Length[EnergyLossTable]}], {meshdims[["m\[Chi]"]], meshdims[["v\[Chi]"]]
            }];
        InterpolationTable = Table[{{m\[Chi]Mesh[[i, j]], v\[Chi]Mesh[[i, j]]},
             EnergyLossMesh[[i, j]]}, {i, meshdims[["m\[Chi]"]]}, {j, meshdims[["v\[Chi]"]]
            }];
        LogInterpolation = Interpolation[Flatten[Log[10, InterpolationTable
            ], 1]];
        outputdict = <|"f" -> LogInterpolation, "m\[Chi]Mesh" -> m\[Chi]Mesh, "v\[Chi]Mesh"
             -> v\[Chi]Mesh, "EnergyLossMesh" -> EnergyLossMesh, "params" -> params, "fname"
             -> fname, "kinlist" -> kinlist, "meshdims" -> meshdims, "\[Epsilon]" -> ToString[
            label]|>;
        If[save,
            Utilities`SaveIt[NotebookDirectory[] <> fname, outputdict]
        ];
        outputdict
    ]


(* ::Section:: *)
(*d\[Sigma]/(d Subscript[E, R]) *)


(* ::Subsection:: *)
(*0 Temperature - Nuclear*)


qb[sign_,m\[Chi]_,v\[Chi]_,\[Omega]_]:= ((m\[Chi] v\[Chi] + <|"+"->1,"-"->(-1)|>[[sign]] Sqrt[2 m\[Chi] (1/2 m\[Chi] v\[Chi]^2 - \[Omega] "\[HBar]")])/("\[HBar]")) (*centre of mass momentum integral bounds*)


dPd\[Omega]Nuc0[\[Omega]_,m\[Chi]_,v\[Chi]_,coeffs_,params_,statpars_]:= ("nI" /.statpars)/(2 \[Omega] "\[HBar]") q^2/("\[HBar]" (2 \[Pi])^2 v\[Chi]) VCoul[q,params]^2 FormFactors`Zeff[q,coeffs]^2 HeavisideTheta[qb["+",m\[Chi],v\[Chi],\[Omega]]-q]HeavisideTheta[q-qb["-",m\[Chi],v\[Chi],\[Omega]]]/. q -> Sqrt[((2 coeffs[["mN"]] \[Omega])/("\[HBar]")) ]/.params

dPd\[Omega]Nuc0Num[\[Omega]_,m\[Chi]_,v\[Chi]_,coeffs_,params_,statpars_]:= ("nI" /.statpars)/(2 \[Omega] "\[HBar]") q^2/("\[HBar]" (2 \[Pi])^2 v\[Chi]) VCoul[q,params]^2 FormFactors`Zeff[q,coeffs]^2 HeavisideTheta[qb["+",m\[Chi],v\[Chi],\[Omega]]-q]HeavisideTheta[q-qb["-",m\[Chi],v\[Chi],\[Omega]]]/. q -> Sqrt[((2 coeffs[["mN"]] \[Omega])/("\[HBar]")) ]/.params


d\[Sigma]dERNuc0[\[Omega]_,m\[Chi]_,v\[Chi]_,coeffs_,params_]:= 1/(2 \[Omega] ("\[HBar]")^3) q^2/((2 \[Pi])^2 v\[Chi]^2) VCoul[q,params]^2 FormFactors`Zeff[q,coeffs]^2 HeavisideTheta[qb["+",m\[Chi],v\[Chi],\[Omega]]-q]HeavisideTheta[q-qb["-",m\[Chi],v\[Chi],\[Omega]]]/. q -> Sqrt[((2 coeffs[["mN"]] \[Omega])/("\[HBar]")) ]/.params


(* ::Subsection:: *)
(*Electronic*)


zb[sign_,\[Omega]_,m\[Chi]_,v\[Chi]_,params_]:= ((m\[Chi] v\[Chi] + <|"+"->1,"-"->(-1)|>[[sign]] Sqrt[2 m\[Chi] (1/2 m\[Chi] v\[Chi]^2 - \[Omega] "\[HBar]")])/(2"qF""\[HBar]")) /.params


zIntegrald\[Sigma]dER[\[Omega]_,m\[Chi]_,v\[Chi]_,params_]:=Sum[NIntegrate[1/z Im[-("Ai"/.params[[i]]) / Dielectrics`\[Epsilon]MNum[("\[HBar]"\[Omega])/(2"EF"z)/.params[[i]], z, u\[Nu]Fit[z] /. params[[i]], params[[i]]]],{z,zb["-",\[Omega],m\[Chi],v\[Chi],params[[i]]],zb["+",\[Omega],m\[Chi],v\[Chi],params[[i]]]}],{i,Length[params]}] (*taking \[Kappa] = 1*)
zIntegrald\[Sigma]dER1oscillator[\[Omega]_,m\[Chi]_,v\[Chi]_,paramsi_]:=NIntegrate[1/z Im[-("Ai"/.paramsi) / Dielectrics`\[Epsilon]MNum[("\[HBar]"\[Omega])/(2"EF"z)/.paramsi, z, u\[Nu]Fit[z] /. paramsi, paramsi]],{z,zb["-",\[Omega],m\[Chi],v\[Chi],paramsi],zb["+",\[Omega],m\[Chi],v\[Chi],paramsi]}](*taking \[Kappa] = 1*)


(*dPd\[Omega]e[\[Omega]_,m\[Chi]_,v\[Chi]_,params_]:=( ("e")^2/(v\[Chi] "\[HBar]" (2 \[Pi])^2 "\[Epsilon]0")2/(1 - E^(-"\[Beta]" "\[HBar]" \[Omega]))/.params)NIntegrate[1/zSum[Im[-("Ai"/.params[[i]]) / Dielectrics`\[Epsilon]MNum[("\[HBar]"\[Omega])/(2"EF"z)/.params[[i]], z, u\[Nu]Fit[z] /. params[[i]], params[[i]]]],{i,Length[params]}],{z,zb["-",\[Omega],m\[Chi],v\[Chi],params],zb["+",\[Omega],m\[Chi],v\[Chi],params]}] (*taking \[Kappa] = 1*)

dPd\[Omega]e[\[Omega]_,m\[Chi]_,v\[Chi]_,params_,\[Beta]_]:=( ("e")^2/(v\[Chi] "\[HBar]" (2 \[Pi])^2 "\[Epsilon]0")2/(1 - E^(-\[Beta] "\[HBar]" \[Omega]))/.params)NIntegrate[1/zSum[Im[-("Ai"/.params[[i]]) / Dielectrics`\[Epsilon]MNum[("\[HBar]"\[Omega])/(2"EF"z)/.params[[i]], z, u\[Nu]Fit[z] /. params[[i]], params[[i]]]],{i,Length[params]}],{z,zb["-",\[Omega],m\[Chi],v\[Chi],params],zb["+",\[Omega],m\[Chi],v\[Chi],params]}] (*taking \[Kappa] = 1*)
*)
(*dPd\[Omega]e[\[Omega]_,m\[Chi]_,v\[Chi]_,params_]:= Sum[( ("e")^2/(v\[Chi] "\[HBar]" (2 \[Pi])^2 "\[Epsilon]0")2/(1 - E^(-"\[Beta]" "\[HBar]" \[Omega]))/.params[[i]])zIntegrald\[Sigma]dER1oscillator[\[Omega],m\[Chi],v\[Chi],params[[i]]],{i,Length[params]}]

dPd\[Omega]e[\[Omega]_,m\[Chi]_,v\[Chi]_,params_,\[Beta]_]:=Sum[( ("e")^2/(v\[Chi] "\[HBar]" (2 \[Pi])^2 "\[Epsilon]0")2/(1 - E^(-\[Beta] "\[HBar]" \[Omega]))/.params[[i]])zIntegrald\[Sigma]dER1oscillator[\[Omega],m\[Chi],v\[Chi],params[[i]]],{i,Length[params]}]
*)

dPd\[Omega]e[\[Omega]_,m\[Chi]_,v\[Chi]_,params_,cut_:4]:= Sum[( ("e")^2/(v\[Chi] "\[HBar]" (2 \[Pi])^2 "\[Epsilon]0") (1+HeavisideTheta[cut - "\[Beta]" "\[HBar]" \[Omega]](1/(1-E^(- "\[Beta]" "\[HBar]" \[Omega]))-1))/.params[[i]])zIntegrald\[Sigma]dER1oscillator[\[Omega],m\[Chi],v\[Chi],params[[i]]],{i,Length[params]}]
dPd\[Omega]e[\[Omega]_,m\[Chi]_,v\[Chi]_,params_,\[Beta]_,cut_:4]:=Sum[( ("e")^2/(v\[Chi] "\[HBar]" (2 \[Pi])^2 "\[Epsilon]0") (1+HeavisideTheta[cut - \[Beta] "\[HBar]" \[Omega]](1/(1-E^(-\[Beta] "\[HBar]" \[Omega]))-1))/.params[[i]])zIntegrald\[Sigma]dER1oscillator[\[Omega],m\[Chi],v\[Chi],params[[i]]],{i,Length[params]}]

dPd\[Omega]eNum[\[Omega]_?NumberQ,m\[Chi]_,v\[Chi]_,params_,cut_:4]:= Sum[( ("e")^2/(v\[Chi] "\[HBar]" (2 \[Pi])^2 "\[Epsilon]0") (1+HeavisideTheta[cut - "\[Beta]" "\[HBar]" \[Omega]](1/(1-E^(- "\[Beta]" "\[HBar]" \[Omega]))-1))/.params[[i]])zIntegrald\[Sigma]dER1oscillator[\[Omega],m\[Chi],v\[Chi],params[[i]]],{i,Length[params]}]


d\[Sigma]dERe[\[Omega]_,m\[Chi]_,v\[Chi]_,params_,cut_:4]:= Sum[( ("e")^2/(v\[Chi]^2 ("\[HBar]")^2 "ne" (2 \[Pi])^2 "\[Epsilon]0") (1+HeavisideTheta[cut - "\[Beta]" "\[HBar]" \[Omega]](1/(1-E^(- "\[Beta]" "\[HBar]" \[Omega]))-1))/.params[[i]])zIntegrald\[Sigma]dER1oscillator[\[Omega],m\[Chi],v\[Chi],params[[i]]],{i,Length[params]}]
d\[Sigma]dERe[\[Omega]_,m\[Chi]_,v\[Chi]_,params_,\[Beta]_,cut_:4]:=Sum[( ("e")^2/(v\[Chi]^2 ("\[HBar]")^2 "ne" (2 \[Pi])^2 "\[Epsilon]0") (1+HeavisideTheta[cut- \[Beta] "\[HBar]" \[Omega]](1/(1-E^(-\[Beta] "\[HBar]" \[Omega]))-1))/.params[[i]])zIntegrald\[Sigma]dER1oscillator[\[Omega],m\[Chi],v\[Chi],params[[i]]],{i,Length[params]}]

d\[Sigma]dEReNum[\[Omega]_?NumberQ,m\[Chi]_,v\[Chi]_,params_,\[Beta]_]:=Sum[( ("e")^2/(v\[Chi]^2 ("\[HBar]")^2 "ne" (2 \[Pi])^2 "\[Epsilon]0") (1+HeavisideTheta[4 - \[Beta] "\[HBar]" \[Omega]](1/(1-E^(-\[Beta] "\[HBar]" \[Omega]))-1))/.params[[i]])zIntegrald\[Sigma]dER1oscillator[\[Omega],m\[Chi],v\[Chi],params[[i]]],{i,Length[params]}]


(* ::Subsubsection:: *)
(*Interpolate d\[Sigma]dER - Nuclear and Electronic*)


(* ::Text:: *)
(*Interpolate over \[Omega]*)


Interpolated\[Sigma]dERe[m\[Chi]_,v\[Chi]_,params_]:=Module[{ERmin,ERmax,InterP\[Omega],\[Epsilon],n=20,Interd\[Sigma]Table,Interd\[Sigma]f},
(*Interpolate over d\[Sigma]dER for electronic contribution (to avoid doing the integral over momentum transfer at every evaluation and to speed up the numerical integration over it to find the normalization)*)
(*local parameters*)
ERmin= "\[Omega]edgei"/.params;
ERmax = 1/(2"\[HBar]") m\[Chi] v\[Chi]^2/.params;
InterP\[Omega]=  10^Subdivide[Log10[ERmin],Log10[ERmax],n];
\[Epsilon]=10^-10;(*Subdivide takes us slightly over ERmax, which throws errors. So subtract off a small fraction \[Epsilon] of ER max*)

Interd\[Sigma]Table=Table[{InterP\[Omega][[i]],d\[Sigma]dERe[InterP\[Omega][[i]](1-\[Epsilon]),m\[Chi],v\[Chi],{params},"\[Beta]"/.params]HeavisideTheta[InterP\[Omega][[i]]-ERmin]},{i,n+1}];
Interd\[Sigma]f=Interpolation[Interd\[Sigma]Table];
<|"d\[Sigma]f"->Interd\[Sigma]f,"d\[Sigma]Table"->Interd\[Sigma]Table,"Interpoints"->InterP\[Omega]|>
]


(* ::Text:: *)
(*Kinematics and Changes of Variable*)


\[Omega]of\[Xi][\[Xi]_,\[Omega]min_,\[Omega]max_]:=\[Omega]min + \[Xi] (\[Omega]max - \[Omega]min)


(* ::Text:: *)
(*For electronic*)


\[Omega]maxofv\[Chi][v\[Chi]_,m\[Chi]_,params_]:= 1/(2"\[HBar]") m\[Chi] v\[Chi]^2/.params
\[Omega]of\[Xi]andv\[Chi][\[Xi]_,\[Omega]min_,v\[Chi]_,m\[Chi]_,params_]:=\[Omega]min + \[Xi] (\[Omega]maxofv\[Chi][v\[Chi],m\[Chi],params] - \[Omega]min)
\[Xi]of\[Omega]andv\[Chi][\[Omega]_,\[Omega]min_,v\[Chi]_,m\[Chi]_,params_]:=(\[Omega]-\[Omega]min )/(\[Omega]maxofv\[Chi][v\[Chi],m\[Chi],params] - \[Omega]min)


(* ::Text:: *)
(*Interpolate over the v\[Chi] as well - Electronic*)


Getv\[Chi]and\[Xi]List[v\[Chi]_,n_:20,\[Epsilon]_:10^-4,uselog_:False]:=Module[{Interpoints\[Xi]},
Interpoints\[Xi]=If[uselog,10^Subdivide[Log10[\[Epsilon]],Log10[1-\[Epsilon]],n],Subdivide[\[Epsilon],1-\[Epsilon],n]];(*move slightly away from the boundaries so we can do a log-log interpolation (d\[Sigma]/Subscript[dE, R]->0 at Subscript[\[Omega], max], and \[Xi]->0 at Subscript[\[Omega], min])*)
(*Interpoints\[Xi]=10^Subdivide[Log10[\[Epsilon]],Log10[1-\[Epsilon]],n];*)
(*Interpoints\[Xi]=Subdivide[\[Epsilon],1-\[Epsilon],n];*)(*move slightly away from the boundaries so we can do a log-log interpolation (d\[Sigma]/Subscript[dE, R]->0 at Subscript[\[Omega], max], and \[Xi]->0 at Subscript[\[Omega], min])*)
Table[{v\[Chi],Interpoints\[Xi][[i]]},{i,n+1}]
]


(*Interpolatev\[Chi]and\[Omega]\[Sigma]dERe\[Xi][m\[Chi]_:0.5 10^6 ("JpereV")/("c")^2/.Constants`SIConstRepl,v\[Chi]_:{4 10^-4,10^-2}"c"/.Constants`SIConstRepl,params_,m_:5,n_:60,\[Epsilon]_:10^-4,cut_:4,uselog_:False]:=Module[{\[Omega]min,\[Omega]max,Interpointsv\[Chi],Interpointsv\[Chi]and\[Xi],Interd\[Sigma]Table,Interd\[Sigma]f,Inter\[Sigma]Table,Inter\[Sigma]f},
(*Interpolate over d\[Sigma]dER for electronic contribution (to avoid doing the integral over momentum transfer at every evaluation and to speed up the numerical integration over it to find the normalization)*)

\[Omega]min= "\[Omega]edgei"/.params;


Interpointsv\[Chi]=10^Subdivide[Log10[v\[Chi][[1]]],Log10[v\[Chi][[2]]],m];
Interpointsv\[Chi]and\[Xi]=SortBy[Table[Getv\[Chi]and\[Xi]List[Interpointsv\[Chi][[i]],n,\[Epsilon],uselog],{i,m+1}],Last];

(*Interpolate to get d\[Sigma]/dER*)
Interd\[Sigma]Table=Table[{Interpointsv\[Chi]and\[Xi][[i,j]],d\[Sigma]dERe[\[Omega]of\[Xi][Interpointsv\[Chi]and\[Xi][[i,j]][[2]],\[Omega]min,\[Omega]maxofv\[Chi][Interpointsv\[Chi]and\[Xi][[i,j]][[1]],m\[Chi],params]],m\[Chi],Interpointsv\[Chi]and\[Xi][[i,j]][[1]],{params},"\[Beta]"/.params,cut]},{i,m+1},{j,n+1}];
Interd\[Sigma]f=Interpolation[Flatten[Log10[Interd\[Sigma]Table],1],InterpolationOrder->4];

(*Interpolate to get \[Sigma]*)
Inter\[Sigma]Table =Table[{Interpointsv\[Chi][[i]],("\[HBar]"/.params)(\[Omega]maxofv\[Chi][Interpointsv\[Chi][[i]],m\[Chi],params]-\[Omega]min)NIntegrate[10^Interd\[Sigma]f[Log10[Interpointsv\[Chi][[i]]],Log10[\[Xi]]],{\[Xi],\[Epsilon],1-\[Epsilon]}]},{i,m+1}]; (*The prefactor on the Integral is the Jacobian of the transformation E_R -> \[Xi]*)
Inter\[Sigma]f = Interpolation[Log10[Inter\[Sigma]Table],InterpolationOrder->4];
(*Inter\[Sigma]f = Interpolation[Inter\[Sigma]Table,InterpolationOrder->4];*)

<|"d\[Sigma]f"->Interd\[Sigma]f,"d\[Sigma]Table"->Interd\[Sigma]Table,"Interpoints"->Interpointsv\[Chi]and\[Xi],"Interpointsv\[Chi]"->Interpointsv\[Chi],"\[Sigma]f"->Inter\[Sigma]f,"\[Sigma]Table"->Inter\[Sigma]Table,"m\[Chi]"->m\[Chi],"v\[Chi]"->Log10[v\[Chi]],"\[Xi]"->Log10[{\[Epsilon],1-\[Epsilon]}],"params"->params,"\[Omega]min"->\[Omega]min|>
]*)


Interpolatev\[Chi]and\[Omega]\[Sigma]dERe\[Xi][m\[Chi]_:0.5 10^6 ("JpereV")/("c")^2/.Constants`SIConstRepl,params_,m_:5,n_:60,\[Epsilon]_:10^-4,cut_:4,uselog_:False,get\[Sigma]too_,regionindex_:1,materialname_:"material",Truncatev\[Chi]_:False]:=Module[{\[Omega]min,v\[Chi],Interpointsv\[Chi],Interpointsv\[Chi]and\[Xi],Interd\[Sigma]Table,Interd\[Sigma]f,Inter\[Sigma]Table,Inter\[Sigma]f,InterdEdlTable,InterdEdlf,Inter\[Sigma]of\[Xi]Table,Inter\[Sigma]of\[Xi]f,v\[Chi]maxindexfor\[Sigma]int},
(*Interpolate over d\[Sigma]dER for electronic contribution (to avoid doing the integral over momentum transfer at every evaluation and to speed up the numerical integration over it to find the normalization)*)

\[Omega]min= "\[Omega]edgei"/.params; (*the cutoff used in the data fitting*)

v\[Chi] = {3Sqrt[("\[HBar]" \[Omega]min)/(2 m\[Chi])], 2 10^-1 "c"}/.Constants`SIConstRepl;(*interpolate from the minimum velocity (up to an order 1 fudge factor (3) for numerical stability) to the boundary of the non-relativistic regime*)

Interpointsv\[Chi]=10^Subdivide[Log10[v\[Chi][[1]]],Log10[v\[Chi][[2]]],m];
Interpointsv\[Chi]and\[Xi]=SortBy[Table[Getv\[Chi]and\[Xi]List[Interpointsv\[Chi][[i]],n,\[Epsilon],uselog],{i,m+1}],Last];

(*Interpolate to get d\[Sigma]/dER*)
Interd\[Sigma]Table=Table[{Interpointsv\[Chi]and\[Xi][[i,j]],d\[Sigma]dERe[\[Omega]of\[Xi][Interpointsv\[Chi]and\[Xi][[i,j]][[2]],\[Omega]min,\[Omega]maxofv\[Chi][Interpointsv\[Chi]and\[Xi][[i,j]][[1]],m\[Chi],params]],m\[Chi],Interpointsv\[Chi]and\[Xi][[i,j]][[1]],{params},"\[Beta]"/.params,cut]},{i,m+1},{j,n+1}];
Interd\[Sigma]f=Interpolation[Flatten[Log10[Interd\[Sigma]Table],1],InterpolationOrder->4];
Print["d\[Sigma] interpolation done"];

(*Interpolate to get \[Sigma]*)
(*If[get\[Sigma]too,
Inter\[Sigma]Table =Table[{Interpointsv\[Chi][[i]],("\[HBar]"/.params)(\[Omega]maxofv\[Chi][Interpointsv\[Chi][[i]],m\[Chi],params]-\[Omega]min)NIntegrate[10^Interd\[Sigma]f[Log10[Interpointsv\[Chi][[i]]],Log10[\[Xi]]],{\[Xi],\[Epsilon],1-\[Epsilon]}]},{i,m+1}]; (*The prefactor on the Integral is the Jacobian of the transformation E_R -> \[Xi]*)
Inter\[Sigma]f = Interpolation[Log10[Inter\[Sigma]Table],InterpolationOrder->4];
Print["\[Sigma] interpolation done"];
Return[<|"d\[Sigma]f"->Interd\[Sigma]f,"d\[Sigma]Table"->Interd\[Sigma]Table,"Interpoints"->Interpointsv\[Chi]and\[Xi],"Interpointsv\[Chi]"->Interpointsv\[Chi],"\[Sigma]f"->Inter\[Sigma]f,"\[Sigma]Table"->Inter\[Sigma]Table,"m\[Chi]"->m\[Chi],"v\[Chi]"->Log10[v\[Chi]],"\[Xi]"->Log10[{\[Epsilon],1-\[Epsilon]}],"params"->params,"\[Omega]min"->\[Omega]min|>],
Return[<|"d\[Sigma]f"->Interd\[Sigma]f,"d\[Sigma]Table"->Interd\[Sigma]Table,"Interpoints"->Interpointsv\[Chi]and\[Xi],"Interpointsv\[Chi]"->Interpointsv\[Chi],"m\[Chi]"->m\[Chi],"v\[Chi]"->Log10[v\[Chi]],"\[Xi]"->Log10[{\[Epsilon],1-\[Epsilon]}],"params"->params,"\[Omega]min"->\[Omega]min|>]
];*)
Inter\[Sigma]Table =Table[{Interpointsv\[Chi][[i]],("\[HBar]"/.params)(\[Omega]maxofv\[Chi][Interpointsv\[Chi][[i]],m\[Chi],params]-\[Omega]min)NIntegrate[10^Interd\[Sigma]f[Log10[Interpointsv\[Chi][[i]]],Log10[\[Xi]]],{\[Xi],\[Epsilon],1-\[Epsilon]}]},{i,m+1}]; (*The prefactor on the Integral is the Jacobian of the transformation E_R -> \[Xi]*)
Inter\[Sigma]f = Interpolation[Log10[Inter\[Sigma]Table],InterpolationOrder->4];
Print["\[Sigma] interpolation done"];

InterdEdlTable =Table[{Interpointsv\[Chi][[i]],("ne"/.params)("\[HBar]"/.params)^2 (\[Omega]maxofv\[Chi][Interpointsv\[Chi][[i]],m\[Chi],params]-\[Omega]min)NIntegrate[\[Omega]of\[Xi]andv\[Chi][\[Xi],\[Omega]min,Interpointsv\[Chi][[i]],m\[Chi],params] 10^Interd\[Sigma]f[Log10[Interpointsv\[Chi][[i]]],Log10[\[Xi]]],{\[Xi],\[Epsilon],1-\[Epsilon]}]},{i,m+1}]; (*The prefactor on the Integral is the Jacobian of the transformation E_R -> \[Xi]*)
InterdEdlf= Interpolation[Log10[InterdEdlTable],InterpolationOrder->4];
Print["dEdlf interpolation done"];

(*Print[N@Interpointsv\[Chi]];
Print@Position[Interpointsv\[Chi],_?(#<10^6&)];
Print[Position[Interpointsv\[Chi],_?(#<10^6&)][[-1,1]]];
Print[Interpointsv\[Chi][[Position[Interpointsv\[Chi],_?(#<10^6&)][[-1,1]]]]];*)

v\[Chi]maxindexfor\[Sigma]int=If[Truncatev\[Chi],Position[Interpointsv\[Chi],_?(#<10^6&)][[-1,1]],m+1];

(*Inter\[Sigma]of\[Xi]Table=Table[{Interpointsv\[Chi]and\[Xi][[i,j]],("\[HBar]"/.params)(\[Omega]maxofv\[Chi][Interpointsv\[Chi]and\[Xi][[i,j]][[1]],m\[Chi],params]-\[Omega]min)Abs[Re[NIntegrate[10^Interd\[Sigma]f[Log10[Interpointsv\[Chi]and\[Xi][[i,j]][[1]]],Log10[\[Xi]]],{\[Xi],Interpointsv\[Chi]and\[Xi][[i,j]][[2]],1-\[Epsilon]}]]]},{i,m+1},{j,n}]; *)
Inter\[Sigma]of\[Xi]Table=Table[{Interpointsv\[Chi]and\[Xi][[i,j]],("\[HBar]"/.params)(\[Omega]maxofv\[Chi][Interpointsv\[Chi]and\[Xi][[i,j]][[1]],m\[Chi],params]-\[Omega]min)Abs[Re[NIntegrate[10^Interd\[Sigma]f[Log10[Interpointsv\[Chi]and\[Xi][[i,j]][[1]]],Log10[\[Xi]]],{\[Xi],Interpointsv\[Chi]and\[Xi][[i,j]][[2]],1-\[Epsilon]}]]]},{i,v\[Chi]maxindexfor\[Sigma]int},{j,n}]; 
Inter\[Sigma]of\[Xi]f=Interpolation[Flatten[Log10[Inter\[Sigma]of\[Xi]Table],1],InterpolationOrder->4];
(*Inter\[Sigma]of\[Xi]f=Interpolation[Flatten[Inter\[Sigma]of\[Xi]Table,1],InterpolationOrder->4];*)
Print["\[Sigma] of \[Xi] interpolation done"];

(* I THINK I AM DOUBLE COUNTING FACTORS OF \[HBar] AND THE JACOBIAN. CHECK THIS CAREFULLYX*)

(*Inter\[Sigma]Table =Table[{Interpointsv\[Chi][[i]],("\[HBar]"/.params)(\[Omega]maxofv\[Chi][Interpointsv\[Chi][[i]],m\[Chi],params]-\[Omega]min)NIntegrate[10^Interd\[Sigma]f[Log10[Interpointsv\[Chi][[i]]],Log10[\[Xi]]],{\[Xi],\[Epsilon],1-\[Epsilon]}]},{i,m+1}]; (*The prefactor on the Integral is the Jacobian of the transformation E_R -> \[Xi]*)
Inter\[Sigma]f = Interpolation[Log10[Inter\[Sigma]Table],InterpolationOrder->4];
Print["\[Sigma] interpolation done"];
*)

(*<|"d\[Sigma]f"->Interd\[Sigma]f,"d\[Sigma]Table"->Interd\[Sigma]Table,"Interpoints"->Interpointsv\[Chi]and\[Xi],"Interpointsv\[Chi]"->Interpointsv\[Chi],"\[Sigma]f"->Inter\[Sigma]f,"\[Sigma]Table"->Inter\[Sigma]Table,"dEdlf"->InterdEdlf,"dEdlTable"->InterdEdlTable,"\[Sigma]of\[Xi]Table"->Inter\[Sigma]of\[Xi]Table,"\[Sigma]of\[Xi]f"->Inter\[Sigma]of\[Xi]f,"m\[Chi]"->m\[Chi],"v\[Chi]"->Log10[v\[Chi]],"\[Xi]"->Log10[{\[Epsilon],1-\[Epsilon]}],"params"->params,"\[Omega]min"->\[Omega]min,"regionindex"->regionindex,"materialname"->materialname|>*)
<|"d\[Sigma]f"->Interd\[Sigma]f,"Interpoints"->Interpointsv\[Chi]and\[Xi],"Interpointsv\[Chi]"->Interpointsv\[Chi],"\[Sigma]f"->Inter\[Sigma]f,"dEdlf"->InterdEdlf,"\[Sigma]of\[Xi]f"->Inter\[Sigma]of\[Xi]f,"m\[Chi]"->m\[Chi],"v\[Chi]"->Log10[v\[Chi]],"\[Xi]"->Log10[{\[Epsilon],1-\[Epsilon]}],"params"->params,"\[Omega]min"->\[Omega]min,"regionindex"->regionindex,"materialname"->materialname|>
]


(* ::Text:: *)
(*Interpolate over the v\[Chi] as well - Nuclear*)


\[Omega]maxNuc[v\[Chi]_,m\[Chi]_,mN_,params_]:= (2 v\[Chi]^2)/("\[HBar]") (m\[Chi]^2 mN)/(m\[Chi] + mN)^2/.params
\[Omega]of\[Xi]Nuc[\[Xi]_,\[Omega]min_,v\[Chi]_,m\[Chi]_,mN_,params_]:=\[Omega]min + \[Xi] (\[Omega]maxNuc[v\[Chi],m\[Chi],mN,params] - \[Omega]min)
\[Xi]of\[Omega]Nuc[\[Omega]_,\[Omega]min_,v\[Chi]_,m\[Chi]_,mN_,params_]:=(\[Omega]-\[Omega]min )/(\[Omega]maxNuc[v\[Chi],m\[Chi],mN,params] - \[Omega]min)


(* ::Text:: *)
(***** should turn these into pure functions and put them in the module*)


(*Interpolatev\[Chi]and\[Omega]\[Sigma]dERNuc\[Xi][m\[Chi]_:0.5 10^6 ("JpereV")/("c")^2/.Constants`SIConstRepl,v\[Chi]_:{4 10^-4,10^-2}"c"/.Constants`SIConstRepl,FFcoeffs_,params_,m_:5,n_:60,\[Epsilon]_:10^-4]:=Module[{\[Omega]min,mN,Interpointsv\[Chi],Interpointsv\[Chi]and\[Xi],Interd\[Sigma]Table,Interd\[Sigma]f},
(*Interpolate over d\[Sigma]dER for electronic contribution (to avoid doing the integral over momentum transfer at every evaluation and to speed up the numerical integration over it to find the normalization)*)

\[Omega]min=0;
mN = "mN"/.FFcoeffs;

Interpointsv\[Chi]=10^Subdivide[Log10[v\[Chi][[1]]],Log10[v\[Chi][[2]]],m];
Interpointsv\[Chi]and\[Xi]=SortBy[Table[Getv\[Chi]and\[Xi]List[Interpointsv\[Chi][[i]],n,\[Epsilon],True],{i,m+1}],Last];

Interd\[Sigma]Table=Table[{Interpointsv\[Chi]and\[Xi][[i,j]],d\[Sigma]dERNuc0[\[Omega]of\[Xi][Interpointsv\[Chi]and\[Xi][[i,j]][[2]],\[Omega]min,\[Omega]maxNuc[Interpointsv\[Chi]and\[Xi][[i,j]][[1]],m\[Chi],mN,params]],m\[Chi],Interpointsv\[Chi]and\[Xi][[i,j]][[1]],FFcoeffs,params]},{i,m+1},{j,n+1}];
Interd\[Sigma]f=Interpolation[Flatten[Log10[Interd\[Sigma]Table],1],InterpolationOrder->4];
<|"d\[Sigma]f"->Interd\[Sigma]f,"d\[Sigma]Table"->Interd\[Sigma]Table,"Interpoints"->Interpointsv\[Chi]and\[Xi],"Interpointsv\[Chi]"->Interpointsv\[Chi],"m\[Chi]"->m\[Chi],"v\[Chi]"->Log10[v\[Chi]],"\[Xi]"->Log10[{\[Epsilon],1-\[Epsilon]}],"NucleusParams"->FFcoeffs,"params"->params,"\[Omega]min"->\[Omega]min|>
]*)


Interpolatev\[Chi]and\[Omega]\[Sigma]dERNuc\[Xi][m\[Chi]_:0.5 10^6 ("JpereV")/("c")^2/.Constants`SIConstRepl,FFcoeffs_,regionindex_:1,materialname_:"material",Truncatev\[Chi]_:False,m_:20,n_:60,v\[Chi]_:{10^-4 "vesc",2 10^-1"c"}/.Constants`SIConstRepl/.Constants`EarthRepl,params_:Constants`SIConstRepl,\[Epsilon]_:10^-4,get\[Sigma]too_:True]:=Module[{\[Omega]min,mN,Interpointsv\[Chi],Interpointsv\[Chi]and\[Xi],Interd\[Sigma]Table,Interd\[Sigma]f,Inter\[Sigma]Table,Inter\[Sigma]f,InterdEdlTable,InterdEdlf,Inter\[Sigma]of\[Xi]Table,Inter\[Sigma]of\[Xi]f,v\[Chi]maxindexfor\[Sigma]int},
(*Interpolate over d\[Sigma]dER and \[Sigma] for Nuclear contribution (to avoid doing the integral over momentum transfer at every evaluation and to speed up the numerical integration over it to find the normalization)*)

\[Omega]min=0;
mN = "mN"/.FFcoeffs;

Interpointsv\[Chi]=10^Subdivide[Log10[v\[Chi][[1]]],Log10[v\[Chi][[2]]],m];
Interpointsv\[Chi]and\[Xi]=SortBy[Table[Getv\[Chi]and\[Xi]List[Interpointsv\[Chi][[i]],n,\[Epsilon],True],{i,m+1}],Last];

(*Print[Interpointsv\[Chi]];
Print[Interpointsv\[Chi]and\[Xi]];*)

(*Interpolate to get d\[Sigma]/dER*)
Interd\[Sigma]Table=Table[{Interpointsv\[Chi]and\[Xi][[i,j]],Max[d\[Sigma]dERNuc0[\[Omega]of\[Xi][Interpointsv\[Chi]and\[Xi][[i,j]][[2]],\[Omega]min,\[Omega]maxNuc[Interpointsv\[Chi]and\[Xi][[i,j]][[1]],m\[Chi],mN,params]],m\[Chi],Interpointsv\[Chi]and\[Xi][[i,j]][[1]],FFcoeffs,params],10^-50]},{i,m+1},{j,n+1}];
Interd\[Sigma]f=Interpolation[Flatten[Log10[Interd\[Sigma]Table],1],InterpolationOrder->4];
Print["d\[Sigma] interpolation done"];

(*Print[Interd\[Sigma]f[5,-3]];
Print[Interd\[Sigma]f["Domain"][[1,2]]];
Print[DensityPlot[Interd\[Sigma]f[v\[Chi],\[Xi]],{v\[Chi],Interd\[Sigma]f["Domain"][[1,1]],Interd\[Sigma]f["Domain"][[1,2]]},{\[Xi],Interd\[Sigma]f["Domain"][[2,1]],Interd\[Sigma]f["Domain"][[2,2]]},PlotLegends->Automatic]];
Print[Interd\[Sigma]Table];*)
(*If[get\[Sigma]too,
(*Interpolate to get \[Sigma]*)
Inter\[Sigma]Table =Table[{Interpointsv\[Chi][[i]],("\[HBar]"/.params)(\[Omega]maxNuc[Interpointsv\[Chi][[i]],m\[Chi],mN,params]-\[Omega]min)NIntegrate[10^Interd\[Sigma]f[Log10[Interpointsv\[Chi][[i]]],Log10[\[Xi]]],{\[Xi],\[Epsilon],1-\[Epsilon]}]},{i,m+1}]; (*The prefactor on the Integral is the Jacobian of the transformation E_R -> \[Xi]*)
Inter\[Sigma]f = Interpolation[Log10[Inter\[Sigma]Table],InterpolationOrder->4];
];
*)
Inter\[Sigma]Table =Table[{Interpointsv\[Chi][[i]],("\[HBar]"/.params)(\[Omega]maxNuc[Interpointsv\[Chi][[i]],m\[Chi],mN,params]-\[Omega]min)NIntegrate[10^Interd\[Sigma]f[Log10[Interpointsv\[Chi][[i]]],Log10[\[Xi]]],{\[Xi],\[Epsilon],1-\[Epsilon]}]},{i,m+1}]; (*The prefactor on the Integral is the Jacobian of the transformation E_R -> \[Xi]*)
Inter\[Sigma]f = Interpolation[Log10[Inter\[Sigma]Table],InterpolationOrder->4];
Print["\[Sigma] interpolation done"];

InterdEdlTable =Table[{Interpointsv\[Chi][[i]],("nI"/.FFcoeffs)("\[HBar]"/.params)^2 (\[Omega]maxNuc[Interpointsv\[Chi][[i]],m\[Chi],mN,params]-\[Omega]min)NIntegrate[\[Omega]of\[Xi]Nuc[\[Xi],\[Omega]min,Interpointsv\[Chi][[i]],m\[Chi],mN,params] 10^Interd\[Sigma]f[Log10[Interpointsv\[Chi][[i]]],Log10[\[Xi]]],{\[Xi],\[Epsilon],1-\[Epsilon]}]},{i,m+1}]; (*The prefactor on the Integral is the Jacobian of the transformation E_R -> \[Xi]*)
InterdEdlf= Interpolation[Log10[InterdEdlTable],InterpolationOrder->4];
Print["dEdlf interpolation done"];

v\[Chi]maxindexfor\[Sigma]int=If[Truncatev\[Chi],Position[Interpointsv\[Chi],_?(#<10^6&)][[-1,1]],m+1];

Inter\[Sigma]of\[Xi]Table=Table[{Interpointsv\[Chi]and\[Xi][[i,j]],("\[HBar]"/.params)(\[Omega]maxNuc[Interpointsv\[Chi]and\[Xi][[i,j]][[1]],m\[Chi],mN,params]-\[Omega]min)Abs[Re[NIntegrate[10^Interd\[Sigma]f[Log10[Interpointsv\[Chi]and\[Xi][[i,j]][[1]]],Log10[\[Xi]]],{\[Xi],Interpointsv\[Chi]and\[Xi][[i,j]][[2]],1-\[Epsilon]}]]]},{i,v\[Chi]maxindexfor\[Sigma]int},{j,n}]; 
Inter\[Sigma]of\[Xi]f=Interpolation[Flatten[Log10[Inter\[Sigma]of\[Xi]Table],1],InterpolationOrder->4];
(*Inter\[Sigma]of\[Xi]f=Interpolation[Flatten[Inter\[Sigma]of\[Xi]Table,1],InterpolationOrder->4];*)
Print["\[Sigma] of \[Xi] interpolation done"];

(*<|"d\[Sigma]f"->Interd\[Sigma]f,"d\[Sigma]Table"->Interd\[Sigma]Table,"Interpoints"->Interpointsv\[Chi]and\[Xi],"\[Sigma]f"->Inter\[Sigma]f,"\[Sigma]Table"->Inter\[Sigma]Table,"dEdlf"->InterdEdlf,"dEdlTable"->InterdEdlTable,"\[Sigma]of\[Xi]Table"->Inter\[Sigma]of\[Xi]Table,"\[Sigma]of\[Xi]f"->Inter\[Sigma]of\[Xi]f,"Interpointsv\[Chi]"->Interpointsv\[Chi],"m\[Chi]"->m\[Chi],"mN"->mN,"v\[Chi]"->Log10[v\[Chi]],"\[Xi]"->Log10[{\[Epsilon],1-\[Epsilon]}],"NucleusParams"->FFcoeffs,"params"->params,"\[Omega]min"->\[Omega]min,"regionindex"->regionindex,"materialname"->materialname|>*)
<|"d\[Sigma]f"->Interd\[Sigma]f,"Interpoints"->Interpointsv\[Chi]and\[Xi],"\[Sigma]f"->Inter\[Sigma]f,"dEdlf"->InterdEdlf,"\[Sigma]of\[Xi]f"->Inter\[Sigma]of\[Xi]f,"Interpointsv\[Chi]"->Interpointsv\[Chi],"m\[Chi]"->m\[Chi],"mN"->mN,"v\[Chi]"->Log10[v\[Chi]],"\[Xi]"->Log10[{\[Epsilon],1-\[Epsilon]}],"NucleusParams"->FFcoeffs,"params"->params,"\[Omega]min"->\[Omega]min,"regionindex"->regionindex,"materialname"->materialname|>
]


(* ::Section:: *)
(*Plotting*)


(* ::Subsection:: *)
(*Plot Energy Loss Interpolation Function - Public*)


(*PlotEL[ELdict_,material_:"Earth"]:=Module[{},
  (* ELdict - association returned by EnergyLossTableAndInterFIT (assumes SI) *)
  Print[Show[{DensityPlot[ELdict[["f"]][m+6 - Log10[("c"^2)/("JpereV")/.Constants`SIConstRepl],v +Log10["c"/.Constants`SIConstRepl]],
  		{m,Log[10,First[ELdict[["kinlist"]][["m\[Chi]"]]]]-6 + Log10[("c")^2/("JpereV")/.Constants`SIConstRepl],Log[10,Last[ELdict[["kinlist"]][["m\[Chi]"]]]]-6 + Log10[("c")^2/("JpereV")/.Constants`SIConstRepl]},
  		{v,Log[10,First[ELdict[["kinlist"]][["v\[Chi]"]]]]-Log10["c"/.Constants`SIConstRepl],Log[10,Last[ELdict[["kinlist"]][["v\[Chi]"]]]]-Log10["c"/.Constants`SIConstRepl]},
  		PlotRange->All,FrameLabel->{"\!\(\*SubscriptBox[\(Log\), \(10\)]\)[m\[Chi]] [MeV]","\!\(\*SubscriptBox[\(Log\), \(10\)]\)[v\[Chi]] [c]"},PlotLabel->StringForm["\!\(\*SubscriptBox[\(Log\), \(10\)]\)[Energy Loss] [eV \!\(\*SuperscriptBox[\(m\), \(-1\)]\)] - ``",material]],
  ContourPlot[ELdict[["f"]][m+6- Log10[("c")^2/("JpereV")/.Constants`SIConstRepl],v+Log10["c"/.Constants`SIConstRepl]],
  		{m,Log[10,First[ELdict[["kinlist"]][["m\[Chi]"]]]]-6 + Log10[("c")^2/("JpereV")/.Constants`SIConstRepl],Log[10,Last[ELdict[["kinlist"]][["m\[Chi]"]]]]-6 + Log10[("c")^2/("JpereV")/.Constants`SIConstRepl]},
  		{v,Log[10,First[ELdict[["kinlist"]][["v\[Chi]"]]]]-Log10["c"/.Constants`SIConstRepl],Log[10,Last[ELdict[["kinlist"]][["v\[Chi]"]]]]-Log10["c"/.Constants`SIConstRepl]},
  		PlotRange->All,ContourLabels->True,ContourShading->False]}]]
  ]*)
  
  PlotEL[ELdict_,material_:"Earth",ChangeTitle_:False,Title_]:=Module[{plottitle},
  (* ELdict - association returned by EnergyLossTableAndInterFIT (assumes SI) *)
  plottitle=If[ChangeTitle,Title,StringForm["\!\(\*SubscriptBox[\(Log\), \(10\)]\)[Energy Loss] [eV \!\(\*SuperscriptBox[\(m\), \(-1\)]\)] - ``",material]];
  Print[Show[{DensityPlot[Re[ELdict[["f"]][m+6 - Log10[("c"^2)/("JpereV")/.Constants`SIConstRepl],v +Log10["c"/.Constants`SIConstRepl]]],
  		{m,Log[10,First[ELdict[["kinlist"]][["m\[Chi]"]]]]-6 + Log10[("c")^2/("JpereV")/.Constants`SIConstRepl],Log[10,Last[ELdict[["kinlist"]][["m\[Chi]"]]]]-6 + Log10[("c")^2/("JpereV")/.Constants`SIConstRepl]},
  		{v,Log[10,First[ELdict[["kinlist"]][["v\[Chi]"]]]]-Log10["c"/.Constants`SIConstRepl],Log[10,Last[ELdict[["kinlist"]][["v\[Chi]"]]]]-Log10["c"/.Constants`SIConstRepl]},
  		PlotRange->All,FrameLabel->{"\!\(\*SubscriptBox[\(Log\), \(10\)]\)[m\[Chi]] [MeV]","\!\(\*SubscriptBox[\(Log\), \(10\)]\)[v\[Chi]] [c]"},PlotLabel->plottitle],
  ContourPlot[Re[ELdict[["f"]][m+6- Log10[("c")^2/("JpereV")/.Constants`SIConstRepl],v+Log10["c"/.Constants`SIConstRepl]]],
  		{m,Log[10,First[ELdict[["kinlist"]][["m\[Chi]"]]]]-6 + Log10[("c")^2/("JpereV")/.Constants`SIConstRepl],Log[10,Last[ELdict[["kinlist"]][["m\[Chi]"]]]]-6 + Log10[("c")^2/("JpereV")/.Constants`SIConstRepl]},
  		{v,Log[10,First[ELdict[["kinlist"]][["v\[Chi]"]]]]-Log10["c"/.Constants`SIConstRepl],Log[10,Last[ELdict[["kinlist"]][["v\[Chi]"]]]]-Log10["c"/.Constants`SIConstRepl]},
  		PlotRange->All,ContourLabels->True,ContourShading->False]}]]
  ]


(* ::Chapter::Closed:: *)
(*End*)


End[];


EndPackage[];
