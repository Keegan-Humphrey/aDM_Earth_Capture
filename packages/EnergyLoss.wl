(* ::Package:: *)

Needs["FormFactors`"]
Needs["Constants`"]
Needs["Utilities`"]
Needs["Dielectrics`"]


FrontEndTokenExecute@"Save"


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

EnergyLossMSIFitEnhanced::usage = "Compute the energy loss for 1 oscillator using fit parameters";

EnergyLossMSIFitSumOscillatorsEnhanced::usage = "Compute EL from optical fits of a material";

EnergyLossTableAndInter::usage = "Create an EL table and interpolation function from theory collision rate";
EnergyLossTableAndInterFITTotalParamsEnhanced::usage = "Create an EL table and interpolation function from optical fits of a material using 1 repl table including the temperature enhancement from FD Theorem";


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
d\[Sigma]dEReInterp::usage = "d\[Sigma]dERe with interpolated structure function";
d\[Sigma]dEReNum::usage = "Numerically Calculate the distribution of cross-section with recoil energy for electrons";

d\[Sigma]dERend::usage = "Numerically Calculate the distribution of cross-section with recoil energy for ups scattering from electrons";



(* ::Text:: *)
(*Interpolate over d\[Sigma]dERe*)


Interpolated\[Sigma]dERe::usage = "Interpolate d\[Sigma]dER over \[Omega] for electronic contribution ";

\[Omega]of\[Xi]::usage = "Get \[Omega](\[Xi]) (where \[Xi] is the dimensionless parameter used in the interpolation)";

\[Omega]maxofv\[Chi]::usage = "Get the maximum energy transfer allowable by DM kinematics (units of [s^-1])";
\[Omega]of\[Xi]andv\[Chi]::usage = "\[Omega]of\[Xi] as a function of DM kinematics";
\[Xi]of\[Omega]andv\[Chi]::usage = "inverse of \[Omega]of\[Xi]andv\[Chi]";

\[Omega]maxofv\[Chi]nd::usage = "Get the maximum energy transfer allowable by DM kinematics (units of [s^-1])";
\[Omega]of\[Xi]andv\[Chi]nd::usage = "\[Omega]of\[Xi] as a function of DM kinematics";
\[Xi]of\[Omega]andv\[Chi]nd::usage = "inverse of \[Omega]of\[Xi]andv\[Chi]";


\[Omega]maxNuc::usage = "Get the maximum energy transfer allowable by DM kinematics (units of [s^-1]) for Nuclear";
\[Omega]of\[Xi]Nuc::usage = "\[Omega]of\[Xi] as a function of DM kinematics";
\[Xi]of\[Omega]Nuc::usage = "inverse of \[Omega]of\[Xi]andv\[Chi]";

Getv\[Chi]and\[Xi]List::usage = "Get points to interpolate d\[Sigma]dER over at fixed v\[Chi]";

Interpolatev\[Chi]and\[Omega]\[Sigma]dERe\[Xi]::usage = "Interpolate d\[Sigma]dER over v\[Chi] and \[Omega](\[Xi]) for electronic";
Interpolatev\[Chi]and\[Omega]\[Sigma]dERe\[Xi]nd::usage = "Interpolate d\[Sigma]dERnd over v\[Chi] and \[Omega](\[Xi]) for electronic";
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


(* ::Subsubsection:: *)
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


(* ::Subsubsection:: *)
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


(* ::Section:: *)
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


(* ::Subsection:: *)
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


(* ::Subsubsection:: *)
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


zIntegrald\[Sigma]dER1oscillatorfnofbounds[\[Omega]_,m\[Chi]_,v\[Chi]_,paramsi_,lb_,ub_]:=NIntegrate[1/z Im[-("Ai"/.paramsi) / Dielectrics`\[Epsilon]MNum[(("\[HBar]"\[Omega])/(2"EF"z)/.paramsi), z, (u\[Nu]Fit[z] /. paramsi), paramsi]],{z,lb,ub},Method->{"LocalAdaptive"},MaxRecursion->12,PrecisionGoal->3,AccuracyGoal->3]

zIntegrald\[Sigma]dER1oscillator[\[Omega]_,m\[Chi]_,v\[Chi]_,paramsi_]:=zIntegrald\[Sigma]dER1oscillatorfnofbounds[\[Omega],m\[Chi],v\[Chi],paramsi,zb["-",\[Omega],m\[Chi],v\[Chi],paramsi],zb["+",\[Omega],m\[Chi],v\[Chi],paramsi]](*taking \[Kappa] = 1*)


dPd\[Omega]e[\[Omega]_,m\[Chi]_,v\[Chi]_,params_,cut_:4]:= Sum[( ("e")^2/(v\[Chi] "\[HBar]" (2 \[Pi])^2 "\[Epsilon]0") If[4 >#,If[#>0.01,1/(1-E^(- #)),1/#-1/2],1]&[ "\[Beta]" "\[HBar]" \[Omega]]/.params[[i]])zIntegrald\[Sigma]dER1oscillator[\[Omega],m\[Chi],v\[Chi],params[[i]]],{i,Length[params]}]
dPd\[Omega]e[\[Omega]_,m\[Chi]_,v\[Chi]_,params_,\[Beta]_,cut_:4]:=Sum[( ("e")^2/(v\[Chi] "\[HBar]" (2 \[Pi])^2 "\[Epsilon]0") If[4 >#,If[#>0.01,1/(1-E^(- #)),1/#-1/2],1]&["\[Beta]" "\[HBar]" \[Omega]]/.params[[i]])zIntegrald\[Sigma]dER1oscillator[\[Omega],m\[Chi],v\[Chi],params[[i]]],{i,Length[params]}]

dPd\[Omega]eNum[\[Omega]_?NumberQ,m\[Chi]_,v\[Chi]_,params_,cut_:4]:= Sum[( ("e")^2/(v\[Chi] "\[HBar]" (2 \[Pi])^2 "\[Epsilon]0") If[4 >#,If[#>0.01,1/(1-E^(- #)),1/#-1/2],1]&[ "\[Beta]" "\[HBar]" \[Omega]]/.params[[i]])zIntegrald\[Sigma]dER1oscillator[\[Omega],m\[Chi],v\[Chi],params[[i]]],{i,Length[params]}]


(*d\[Sigma]dERe[\[Omega]_,m\[Chi]_,v\[Chi]_,params_,cut_:4]:= Sum[( ("e")^2/(v\[Chi]^2 ("\[HBar]")^2 "ne" (2 \[Pi])^2 "\[Epsilon]0") (1+HeavisideTheta[cut - "\[Beta]" "\[HBar]" \[Omega]](1/(1-E^(- "\[Beta]" "\[HBar]" \[Omega]))-1))/.params[[i]])zIntegrald\[Sigma]dER1oscillator[\[Omega],m\[Chi],v\[Chi],params[[i]]],{i,Length[params]}]*)
d\[Sigma]dERe[\[Omega]_,m\[Chi]_,v\[Chi]_,params_,\[Beta]_,cut_:4]:=Sum[( ("e")^2/(v\[Chi]^2 ("\[HBar]")^2 "ne" (2 \[Pi])^2 "\[Epsilon]0") If[4 >#,If[#>0.01,1/(1-E^(- #)),1/#-1/2],1]&[ "\[Beta]" "\[HBar]" \[Omega]]/.params[[i]])zIntegrald\[Sigma]dER1oscillator[\[Omega],m\[Chi],v\[Chi],params[[i]]],{i,Length[params]}]

d\[Sigma]dEReInterp[\[Omega]_,m\[Chi]_,v\[Chi]_,StrucTable_]:= Module[{params,func},

Sum[(("e")^2/(v\[Chi]^2 ("\[HBar]")^2 "ne" (2 \[Pi])^2 "\[Epsilon]0")/.StrucTable[[i,2]])NIntegrate[10^StrucTable[[i,1]][Log10@(\[Omega]/(2 "qF""vF"z)/.StrucTable[[i,2]]),Log10@z],{z,zb["-",\[Omega],m\[Chi],v\[Chi],StrucTable[[i,2]]],zb["+",\[Omega],m\[Chi],v\[Chi],StrucTable[[i,2]]]}],{i,Length[StrucTable]}]
]

d\[Sigma]dEReNum[\[Omega]_?NumberQ,m\[Chi]_,v\[Chi]_,params_,\[Beta]_]:=Sum[( ("e")^2/(v\[Chi]^2 ("\[HBar]")^2 "ne" (2 \[Pi])^2 "\[Epsilon]0") (1+HeavisideTheta[4 - \[Beta] "\[HBar]" \[Omega]](1/(1-E^(-\[Beta] "\[HBar]" \[Omega]))-1))/.params[[i]])zIntegrald\[Sigma]dER1oscillator[\[Omega],m\[Chi],v\[Chi],params[[i]]],{i,Length[params]}]


(* ::Subsubsection:: *)
(*electronic - evaporation*)


(*zIntegrald\[Sigma]dER1oscillatorfnofboundsnd[\[Omega]_,m\[Chi]_,v\[Chi]_,paramsi_,lb_,ub_]:=NIntegrate[1/z Im[("Ai"/.paramsi) Dielectrics`\[Epsilon]MndNum[(("\[HBar]"\[Omega])/(2"EF"z)/.paramsi), z, (u\[Nu]Fit[z] /. paramsi), paramsi] / Abs[Dielectrics`\[Epsilon]MNum[(("\[HBar]"\[Omega])/(2"EF"z)/.paramsi), z, (u\[Nu]Fit[z] /. paramsi), paramsi]]^2],{z,lb,ub},Method->{"LocalAdaptive"},MaxRecursion->12,PrecisionGoal->3,AccuracyGoal->3]*)
zIntegrald\[Sigma]dER1oscillatorfnofboundsnd[\[Omega]_,m\[Chi]_,v\[Chi]_,paramsi_,lb_,ub_]:=NIntegrate[1/z Im[-("Ai"/.paramsi) / Dielectrics`\[Epsilon]MndNum[(("\[HBar]"\[Omega])/(2"EF"z)/.paramsi), z, (u\[Nu]Fit[z] /. paramsi), paramsi]],{z,lb,ub},Method->{"LocalAdaptive"},MaxRecursion->12,PrecisionGoal->3,AccuracyGoal->3]

Clear[zIntegrald\[Sigma]dER1oscillatornd]
zIntegrald\[Sigma]dER1oscillatornd[\[Omega]_,m\[Chi]_,v\[Chi]_,paramsi_]:=Module[{zm,zp,ndlim},
(*compute the d\[Sigma]dERnd z integral 
integration region is the intersection of (0,1/D) and (Subscript[z, -],Subscript[z, +])*)
zm=-zb["-",-\[Omega],m\[Chi],v\[Chi],paramsi];
zp=zb["+",-\[Omega],m\[Chi],v\[Chi],paramsi];
ndlim=1/("D")/.paramsi;
(*Print[zm," ",zp," ",ndlim];*)
If[ndlim<=zm,Return[10^-100,Module](*a negligible value w no log troubles when interpolating*),
	If[ndlim<=zp,Return[zIntegrald\[Sigma]dER1oscillatorfnofboundsnd[\[Omega],m\[Chi],v\[Chi],paramsi,zm,ndlim],Module],
		Return[zIntegrald\[Sigma]dER1oscillatorfnofboundsnd[\[Omega],m\[Chi],v\[Chi],paramsi,zm,zp],Module]
	];
];
]


d\[Sigma]dERend[\[Omega]_,m\[Chi]_,v\[Chi]_,params_,\[Beta]_,cut_:4]:=Sum[( ("e")^2/(v\[Chi]^2 ("\[HBar]")^2 "ne" (2 \[Pi])^2 "\[Epsilon]0") If[cut >#,If[#>0.01,1/(1-E^(- #)),1/#-1/2],1]&[ "\[Beta]" "\[HBar]" \[Omega]]/.params[[o]])zIntegrald\[Sigma]dER1oscillatornd[\[Omega],m\[Chi],v\[Chi],params[[o]]],{o,Length[params]}]


\[Omega]maxofv\[Chi]nd[v\[Chi]_?NumericQ,m\[Chi]_?NumericQ,params_]:= ((2 "qF" v\[Chi])/("D") + (2 (("qF")^2) "\[HBar]" )/(("D")^2 m\[Chi]))/.params
\[Omega]of\[Xi]andv\[Chi]nd[\[Xi]_?NumericQ,\[Omega]min_,v\[Chi]_,m\[Chi]_,params_]:=\[Omega]min + \[Xi] (\[Omega]maxofv\[Chi]nd[v\[Chi],m\[Chi],params] - \[Omega]min)
\[Xi]of\[Omega]andv\[Chi]nd[\[Omega]_,\[Omega]min_,v\[Chi]_,m\[Chi]_,params_]:=(\[Omega]-\[Omega]min )/(\[Omega]maxofv\[Chi]nd[v\[Chi],m\[Chi],params] - \[Omega]min)


(*\[Omega]maxofv\[Chi][v\[Chi]_,m\[Chi]_,params_]:= 1/(2"\[HBar]") m\[Chi] v\[Chi]^2/.params*)


(*Getv\[Chi]and\[Xi]Listnd[v\[Chi]_,n_:20,\[Epsilon]_:10^-4,uselog_:True]:=Module[{Interpoints\[Xi]},
Interpoints\[Xi]=If[uselog,10^Subdivide[Log10[\[Epsilon]],Log10[1-\[Epsilon]],n],Subdivide[\[Epsilon],1-\[Epsilon],n]];(*move slightly away from the boundaries so we can do a log-log interpolation (d\[Sigma]/Subscript[dE, R]->0 at Subscript[\[Omega], max], and \[Xi]->0 at Subscript[\[Omega], min])*)
(*Interpoints\[Xi]=10^Subdivide[Log10[\[Epsilon]],Log10[1-\[Epsilon]],n];*)
(*Interpoints\[Xi]=Subdivide[\[Epsilon],1-\[Epsilon],n];*)(*move slightly away from the boundaries so we can do a log-log interpolation (d\[Sigma]/Subscript[dE, R]->0 at Subscript[\[Omega], max], and \[Xi]->0 at Subscript[\[Omega], min])*)
Table[{v\[Chi],Interpoints\[Xi][[i]]},{i,n+1}]
]*)


Clear[Interpolatev\[Chi]and\[Omega]\[Sigma]dERe\[Xi]nd]
Interpolatev\[Chi]and\[Omega]\[Sigma]dERe\[Xi]nd[m\[Chi]_:0.5 10^6 ("JpereV")/("c")^2/.Constants`SIConstRepl,params_,m_:5,n_:60,(*\[Epsilon]_:10^-6*)\[Epsilon]_:10^-6,cut_:4,uselog_:True,get\[Sigma]too_,regionindex_:1,materialname_:"material",Truncatev\[Chi]_:False]:=Module[{\[Omega]min,v\[Chi],Interpointsv\[Chi],Interpointsv\[Chi]and\[Xi],Interd\[Sigma]Table,Interd\[Sigma]f,Inter\[Sigma]Table,Inter\[Sigma]f,InterdEdlTable,InterdEdlf,Inter\[Sigma]of\[Xi]Table,Inter\[Sigma]of\[Xi]f,v\[Chi]maxindexfor\[Sigma]int},
(*Interpolate over d\[Sigma]dERnd FOR UP SCATTERING for electronic contribution (to avoid doing the integral over momentum transfer at every evaluation and to speed up the numerical integration over it to find the normalization)*)

(*
zIntegrald\[Sigma]dER1oscillatornd returns non-zero value iff 1/D>=(m\[Chi] v\[Chi]-Sqrt[2] Sqrt[m\[Chi] ((m\[Chi] v\[Chi]^2)/2-\[Omega] \[HBar])])/(2 kF \[HBar]) = Subscript[z, -].
RHS is decreasing with Subscript[v, \[Chi]], and increasing with \[Omega]. 
The smallest value of the RHS is given by maximizing \[Omega] and minimizing v\[Chi]

The upper bound on |\[Omega]| st this is satisfied is (2 Subscript[k, F] v\[Chi])/D +(2 (Subscript[k, F]^2) \[HBar] )/(D^2 m\[Chi]) -> so we pass this to Getv\[Chi]and\[Xi]Listnd to find the maximum value of \[Xi]

The lower bound on v\[Chi] st this energy transfer is enough to evaporate is Subscript[v, esc] - (2 \[HBar] kF)/(D m\[Chi])

using the correct upper bound for \[Omega] (\[Omega]maxofv\[Chi]nd) it is going to much smaller \[Omega] as a result of the much lower \[Omega]max. Tantamount to going to lower \[Epsilon]. NIntegrate stuggles so keeps iterating. Just stick to the original \[Omega]max, and set disallowed \[Omega] to negligible values (10^-100).
*)

\[Omega]min= "\[Omega]edgei"/.params; (*the cutoff used in the data fitting*)

(*v\[Chi] = {3Sqrt[("\[HBar]" \[Omega]min)/(2 m\[Chi])], 2 10^-1 "c"}/.Constants`SIConstRepl;(*interpolate from the minimum velocity (up to an order 1 fudge factor (3) for numerical stability) to the boundary of the non-relativistic regime*)*)
(*v\[Chi] = {("D" \[Omega]min)/(2 "qF")+("qF" "\[HBar]")/("D" m\[Chi]), 2 10^-1 "c"}/.Constants`SIConstRepl/.params;(*interpolate from the minimum velocity to the boundary of the non-relativistic regime*)*)
v\[Chi] = {"vesc"-(2 "\[HBar]" "qF")/("D" m\[Chi]), 2 10^-1 "c"}/.Constants`SIConstRepl/.Constants`EarthRepl/.params;(*interpolate from the minimum velocity to the boundary of the non-relativistic regime*)
(*Print[N@v\[Chi]];*)

If[v\[Chi][[1]]>v\[Chi][[2]],Print["minimum energy requires relativistic probes, we won't consider these"];(*Return["NA",Module]*)Return[<|"NA"->Null|>,Module]];

Interpointsv\[Chi]=10^Subdivide[Log10[v\[Chi][[1]]],Log10[v\[Chi][[2]]],m];
Interpointsv\[Chi]and\[Xi]=SortBy[Table[Getv\[Chi]and\[Xi]List[Interpointsv\[Chi][[i]],n,\[Epsilon],uselog],{i,m+1}],Last];

(*Print[N@Interpointsv\[Chi]];
Print[N@Interpointsv\[Chi]and\[Xi]];*)

(*Interpolate to get d\[Sigma]/dER*)
(*Interd\[Sigma]Table=Table[{Interpointsv\[Chi]and\[Xi][[i,j]],d\[Sigma]dERend[\[Omega]of\[Xi][Interpointsv\[Chi]and\[Xi][[i,j]][[2]],\[Omega]min,\[Omega]maxofv\[Chi][Interpointsv\[Chi]and\[Xi][[i,j]][[1]],m\[Chi],params]],m\[Chi],Interpointsv\[Chi]and\[Xi][[i,j]][[1]],{params},"\[Beta]"/.params,cut]},{i,m+1},{j,n+1}];*)
Interd\[Sigma]Table=Table[{Interpointsv\[Chi]and\[Xi][[i,j]],d\[Sigma]dERend[\[Omega]of\[Xi][Interpointsv\[Chi]and\[Xi][[i,j]][[2]],\[Omega]min,\[Omega]maxofv\[Chi]nd[Interpointsv\[Chi]and\[Xi][[i,j]][[1]],m\[Chi],params]],m\[Chi],Interpointsv\[Chi]and\[Xi][[i,j]][[1]],{params},"\[Beta]"/.params,cut]},{i,m+1},{j,n+1}];
Interd\[Sigma]f=Interpolation[Flatten[Log10[Interd\[Sigma]Table],1],InterpolationOrder->3];
Print["d\[Sigma] interpolation done"];

(*Inter\[Sigma]Table =Table[{Interpointsv\[Chi][[i]],("\[HBar]"/.params)(\[Omega]maxofv\[Chi][Interpointsv\[Chi][[i]],m\[Chi],params]-\[Omega]min)NIntegrate[10^Interd\[Sigma]f[Log10[Interpointsv\[Chi][[i]]],Log10[\[Xi]]],{\[Xi],\[Epsilon],1-\[Epsilon]},Method->{"LocalAdaptive"},PrecisionGoal->3,MaxRecursion->8,AccuracyGoal->3]},{i,m+1}]; (*The prefactor on the Integral is the Jacobian of the transformation E_R -> \[Xi]*)*)
Inter\[Sigma]Table =Table[{Interpointsv\[Chi][[i]],("\[HBar]"/.params)(\[Omega]maxofv\[Chi]nd[Interpointsv\[Chi][[i]],m\[Chi],params]-\[Omega]min)NIntegrate[10^Interd\[Sigma]f[Log10[Interpointsv\[Chi][[i]]],Log10[\[Xi]]],{\[Xi],\[Epsilon],1-\[Epsilon]},Method->{"LocalAdaptive"},PrecisionGoal->3,MaxRecursion->8,AccuracyGoal->3]},{i,m+1}]; (*The prefactor on the Integral is the Jacobian of the transformation E_R -> \[Xi]*)
Inter\[Sigma]f = Interpolation[Log10[Inter\[Sigma]Table],InterpolationOrder->4];
Print["\[Sigma] interpolation done"];

(*\[Omega]of\[Xi][\[Xi]_,\[Omega]min_,\[Omega]max_]:=\[Omega]min + \[Xi] (\[Omega]max - \[Omega]min)*)

(*InterdEdlTable =Table[{Interpointsv\[Chi][[i]],("ne"/.params)("\[HBar]"/.params)^2 (\[Omega]maxofv\[Chi][Interpointsv\[Chi][[i]],m\[Chi],params]-\[Omega]min)NIntegrate[\[Omega]of\[Xi]andv\[Chi][\[Xi],\[Omega]min,Interpointsv\[Chi][[i]],m\[Chi],params] 10^Interd\[Sigma]f[Log10[Interpointsv\[Chi][[i]]],Log10[\[Xi]]],{\[Xi],\[Epsilon],1-\[Epsilon]},Method->{"LocalAdaptive"},PrecisionGoal->3,MaxRecursion->8,AccuracyGoal->3]},{i,m+1}]; (*The prefactor on the Integral is the Jacobian of the transformation E_R -> \[Xi]*)*)
InterdEdlTable =Table[{Interpointsv\[Chi][[i]],("ne"/.params)("\[HBar]"/.params)^2 (\[Omega]maxofv\[Chi]nd[Interpointsv\[Chi][[i]],m\[Chi],params]-\[Omega]min)NIntegrate[\[Omega]of\[Xi]andv\[Chi]nd[\[Xi],\[Omega]min,Interpointsv\[Chi][[i]],m\[Chi],params] 10^Interd\[Sigma]f[Log10[Interpointsv\[Chi][[i]]],Log10[\[Xi]]],{\[Xi],\[Epsilon],1-\[Epsilon]},Method->{"LocalAdaptive"},PrecisionGoal->3,MaxRecursion->8,AccuracyGoal->3]},{i,m+1}]; (*The prefactor on the Integral is the Jacobian of the transformation E_R -> \[Xi]*)
InterdEdlf= Interpolation[Log10[InterdEdlTable],InterpolationOrder->4];
Print["dEdlf interpolation done"];

v\[Chi]maxindexfor\[Sigma]int=If[Truncatev\[Chi],Position[Interpointsv\[Chi],_?(#<10^6&)][[-1,1]],m+1];

(*Inter\[Sigma]of\[Xi]Table=Table[{Interpointsv\[Chi]and\[Xi][[i,j]],10^-49+("\[HBar]"/.params)(\[Omega]maxofv\[Chi][Interpointsv\[Chi]and\[Xi][[i,j]][[1]],m\[Chi],params]-\[Omega]min)Abs[Re[NIntegrate[10^Interd\[Sigma]f[Log10[Interpointsv\[Chi]and\[Xi][[i,j]][[1]]],Log10[\[Xi]]],{\[Xi],Interpointsv\[Chi]and\[Xi][[i,j]][[2]],1-\[Epsilon]},Method->{"LocalAdaptive"},PrecisionGoal->3,MaxRecursion->8,AccuracyGoal->3]]]},{i,v\[Chi]maxindexfor\[Sigma]int},{j,n+1}]; *)
Inter\[Sigma]of\[Xi]Table=Table[{Interpointsv\[Chi]and\[Xi][[i,j]],10^-100+("\[HBar]"/.params)(\[Omega]maxofv\[Chi]nd[Interpointsv\[Chi]and\[Xi][[i,j]][[1]],m\[Chi],params]-\[Omega]min)Abs[Re[NIntegrate[10^Interd\[Sigma]f[Log10[Interpointsv\[Chi]and\[Xi][[i,j]][[1]]],Log10[\[Xi]]],{\[Xi],Interpointsv\[Chi]and\[Xi][[i,j]][[2]],1-\[Epsilon]},Method->{"LocalAdaptive"},PrecisionGoal->3,MaxRecursion->8,AccuracyGoal->3]]]},{i,v\[Chi]maxindexfor\[Sigma]int},{j,n+1}]; 
Inter\[Sigma]of\[Xi]f=Interpolation[Flatten[Log10[Inter\[Sigma]of\[Xi]Table],1],InterpolationOrder->4];

Print["\[Sigma] of \[Xi] interpolation done"];

<|"d\[Sigma]f"->Interd\[Sigma]f,"Interpoints"->Interpointsv\[Chi]and\[Xi],"Interpointsv\[Chi]"->Interpointsv\[Chi],"\[Sigma]f"->Inter\[Sigma]f,"dEdlf"->InterdEdlf,"\[Sigma]of\[Xi]f"->Inter\[Sigma]of\[Xi]f,"\[Sigma]of\[Xi]table"->Inter\[Sigma]of\[Xi]Table,"m\[Chi]"->m\[Chi],"v\[Chi]"->Log10[v\[Chi]],"\[Xi]"->Log10[{\[Epsilon],1-\[Epsilon]}],"params"->params,"\[Omega]min"->\[Omega]min,"regionindex"->regionindex,"materialname"->materialname|>
]


(*Clear[Interpolatev\[Chi]and\[Omega]\[Sigma]dERe\[Xi]nd]
Interpolatev\[Chi]and\[Omega]\[Sigma]dERe\[Xi]nd[m\[Chi]_:0.5 10^6 ("JpereV")/("c")^2/.Constants`SIConstRepl,params_,m_:5,n_:60,\[Epsilon]_:10^-6,cut_:4,uselog_:True,get\[Sigma]too_,regionindex_:1,materialname_:"material",Truncatev\[Chi]_:False]:=Module[{\[Omega]min,v\[Chi],Interpointsv\[Chi],Interpointsv\[Chi]and\[Xi],Interd\[Sigma]Table,Interd\[Sigma]f,Inter\[Sigma]Table,Inter\[Sigma]f,InterdEdlTable,InterdEdlf,Inter\[Sigma]of\[Xi]Table,Inter\[Sigma]of\[Xi]f,v\[Chi]maxindexfor\[Sigma]int},
(*Interpolate over d\[Sigma]dERnd FOR UP SCATTERING for electronic contribution (to avoid doing the integral over momentum transfer at every evaluation and to speed up the numerical integration over it to find the normalization)*)

(*
zIntegrald\[Sigma]dER1oscillatornd returns non-zero value iff 1/D>=(m\[Chi] v\[Chi]-Sqrt[2] Sqrt[m\[Chi] ((m\[Chi] v\[Chi]^2)/2-\[Omega] \[HBar])])/(2 kF \[HBar]) = Subscript[z, -].
RHS is decreasing with Subscript[v, \[Chi]], and increasing with \[Omega]. 
The smallest value of the RHS is given by maximizing \[Omega] and minimizing v\[Chi]

The lower bound on v\[Chi] st this is satisfied is (D \[Omega])/(2 Subscript[k, F])+(Subscript[k, F] \[HBar])/(D m\[Chi])
The upper bound on \[Omega] st this is satisfied is (2 Subscript[k, F] v\[Chi])/D - (2 (Subscript[k, F]^2) \[HBar] )/(D^2 m\[Chi]) -> so we pass this to Getv\[Chi]and\[Xi]Listnd to find the maximum value of \[Xi]

Significantly slower than above since it is going to much smaller \[Omega] as a result of the much lower \[Omega]max. Tantamount to going to lower \[Epsilon]. NIntegrate stuggles so keeps iterating. Just stick to the original \[Omega]max, and set disallowed \[Omega] to negligible values.
*)

\[Omega]min= "\[Omega]edgei"/.params; (*the cutoff used in the data fitting*)

(*v\[Chi] = {3Sqrt[("\[HBar]" \[Omega]min)/(2 m\[Chi])], 2 10^-1 "c"}/.Constants`SIConstRepl;(*interpolate from the minimum velocity (up to an order 1 fudge factor (3) for numerical stability) to the boundary of the non-relativistic regime*)*)
v\[Chi] = {("D" \[Omega]min)/(2 "qF")+("qF" "\[HBar]")/("D" m\[Chi]), 2 10^-1 "c"}/.Constants`SIConstRepl/.params;(*interpolate from the minimum velocity to the boundary of the non-relativistic regime*)

(*Print[\[Omega]maxofv\[Chi]nd[v\[Chi],m\[Chi],params]];*)

(*Print[N@v\[Chi]];*)

If[v\[Chi][[1]]>v\[Chi][[2]],Print["minimum energy requires relativistic probes, we won't consider these"];Return["NA",Module]];

Interpointsv\[Chi]=10^Subdivide[Log10[v\[Chi][[1]]],Log10[v\[Chi][[2]]],m];
Interpointsv\[Chi]and\[Xi]=SortBy[Table[Getv\[Chi]and\[Xi]List[Interpointsv\[Chi][[i]],n,\[Epsilon],uselog],{i,m+1}],Last];

(*Print[N@Interpointsv\[Chi]];
Print[N@Interpointsv\[Chi]and\[Xi]];*)

(*Interpolate to get d\[Sigma]/dER*)
(*Print[AbsoluteTiming@With[{i=5,j=5},10^-49+d\[Sigma]dERend[\[Omega]of\[Xi][Interpointsv\[Chi]and\[Xi][[i,j]][[2]],\[Omega]min,\[Omega]maxofv\[Chi]nd[Interpointsv\[Chi]and\[Xi][[i,j]][[1]],m\[Chi],params]],m\[Chi],Interpointsv\[Chi]and\[Xi][[i,j]][[1]],{params},"\[Beta]"/.params,cut]]];
Print[AbsoluteTiming@With[{i=5,j=5},10^-49+d\[Sigma]dERend[\[Omega]of\[Xi][Interpointsv\[Chi]and\[Xi][[i,j]][[2]],\[Omega]min,\[Omega]maxofv\[Chi][Interpointsv\[Chi]and\[Xi][[i,j]][[1]],m\[Chi],params]],m\[Chi],Interpointsv\[Chi]and\[Xi][[i,j]][[1]],{params},"\[Beta]"/.params,cut]]];*)
(*Do[Print@AbsoluteTiming[d\[Sigma]dERend[\[Omega]of\[Xi][Interpointsv\[Chi]and\[Xi][[i,j]][[2]],\[Omega]min,\[Omega]maxofv\[Chi][Interpointsv\[Chi]and\[Xi][[i,j]][[1]],m\[Chi],params]],m\[Chi],Interpointsv\[Chi]and\[Xi][[i,j]][[1]],{params},"\[Beta]"/.params,cut]][[1]],{i,m+1},{j,n+1}];
Do[Print@AbsoluteTiming[d\[Sigma]dERend[\[Omega]of\[Xi][Interpointsv\[Chi]and\[Xi][[i,j]][[2]],\[Omega]min,\[Omega]maxofv\[Chi]nd[Interpointsv\[Chi]and\[Xi][[i,j]][[1]],m\[Chi],params]],m\[Chi],Interpointsv\[Chi]and\[Xi][[i,j]][[1]],{params},"\[Beta]"/.params,cut]][[1]],{i,m+1},{j,n+1}];*)
(*Do[Print[AbsoluteTiming[d\[Sigma]dERend[\[Omega]of\[Xi][Interpointsv\[Chi]and\[Xi][[i,j]][[2]],\[Omega]min,\[Omega]maxofv\[Chi][Interpointsv\[Chi]and\[Xi][[i,j]][[1]],m\[Chi],params]],m\[Chi],Interpointsv\[Chi]and\[Xi][[i,j]][[1]],{params},"\[Beta]"/.params,cut]][[1]],
" ",AbsoluteTiming[d\[Sigma]dERend[\[Omega]of\[Xi][Interpointsv\[Chi]and\[Xi][[i,j]][[2]],\[Omega]min,\[Omega]maxofv\[Chi]nd[Interpointsv\[Chi]and\[Xi][[i,j]][[1]],m\[Chi],params]],m\[Chi],Interpointsv\[Chi]and\[Xi][[i,j]][[1]],{params},"\[Beta]"/.params,cut]][[1]]],{i,m+1},{j,n+1}];*)
(*Print["donzo"]*)
(*Interd\[Sigma]Table=Table[{Interpointsv\[Chi]and\[Xi][[i,j]],10^-49+d\[Sigma]dERend[\[Omega]of\[Xi][Interpointsv\[Chi]and\[Xi][[i,j]][[2]],\[Omega]min,\[Omega]maxofv\[Chi]nd[Interpointsv\[Chi]and\[Xi][[i,j]][[1]],m\[Chi],params]],m\[Chi],Interpointsv\[Chi]and\[Xi][[i,j]][[1]],{params},"\[Beta]"/.params,cut]},{i,m+1},{j,n+1}];*)
Interd\[Sigma]Table=Monitor[Table[{Interpointsv\[Chi]and\[Xi][[i,j]],10^-100+d\[Sigma]dERend[\[Omega]of\[Xi][Interpointsv\[Chi]and\[Xi][[i,j]][[2]],\[Omega]min,\[Omega]maxofv\[Chi]nd[Interpointsv\[Chi]and\[Xi][[i,j]][[1]],m\[Chi],params]],m\[Chi],Interpointsv\[Chi]and\[Xi][[i,j]][[1]],{params},"\[Beta]"/.params,cut]},{i,m+1},{j,n+1}],{i,j}];
(*Print[Interd\[Sigma]Table];*)
Print["donzo"];
Interd\[Sigma]f=Interpolation[Flatten[Log10[Interd\[Sigma]Table],1],InterpolationOrder->3];
Print["d\[Sigma] interpolation done"];

Inter\[Sigma]Table =Table[{Interpointsv\[Chi][[i]],("\[HBar]"/.params)(\[Omega]maxofv\[Chi]nd[Interpointsv\[Chi][[i]],m\[Chi],params]-\[Omega]min)NIntegrate[10^Interd\[Sigma]f[Log10[Interpointsv\[Chi][[i]]],Log10[\[Xi]]],{\[Xi],\[Epsilon],1-\[Epsilon]},Method->{"LocalAdaptive"}]},{i,m+1}]; (*The prefactor on the Integral is the Jacobian of the transformation E_R -> \[Xi]*)
Inter\[Sigma]f = Interpolation[Log10[Inter\[Sigma]Table],InterpolationOrder->4];
Print["\[Sigma] interpolation done"];

(*\[Omega]of\[Xi][\[Xi]_,\[Omega]min_,\[Omega]max_]:=\[Omega]min + \[Xi] (\[Omega]max - \[Omega]min)*)

InterdEdlTable =Table[{Interpointsv\[Chi][[i]],("ne"/.params)("\[HBar]"/.params)^2 (\[Omega]maxofv\[Chi]nd[Interpointsv\[Chi][[i]],m\[Chi],params]-\[Omega]min)NIntegrate[\[Omega]of\[Xi]andv\[Chi]nd[\[Xi],\[Omega]min,Interpointsv\[Chi][[i]],m\[Chi],params] 10^Interd\[Sigma]f[Log10[Interpointsv\[Chi][[i]]],Log10[\[Xi]]],{\[Xi],\[Epsilon],1-\[Epsilon]},Method->{"LocalAdaptive"}]},{i,m+1}]; (*The prefactor on the Integral is the Jacobian of the transformation E_R -> \[Xi]*)
InterdEdlf= Interpolation[Log10[InterdEdlTable],InterpolationOrder->4];
Print["dEdlf interpolation done"];

v\[Chi]maxindexfor\[Sigma]int=If[Truncatev\[Chi],Position[Interpointsv\[Chi],_?(#<10^6&)][[-1,1]],m+1];

Inter\[Sigma]of\[Xi]Table=Table[{Interpointsv\[Chi]and\[Xi][[i,j]],10^-49+("\[HBar]"/.params)(\[Omega]maxofv\[Chi]nd[Interpointsv\[Chi]and\[Xi][[i,j]][[1]],m\[Chi],params]-\[Omega]min)Abs[Re[NIntegrate[10^Interd\[Sigma]f[Log10[Interpointsv\[Chi]and\[Xi][[i,j]][[1]]],Log10[\[Xi]]],{\[Xi],Interpointsv\[Chi]and\[Xi][[i,j]][[2]],1-\[Epsilon]},Method->{"LocalAdaptive"},PrecisionGoal->3,MaxRecursion->8,AccuracyGoal->3]]]},{i,v\[Chi]maxindexfor\[Sigma]int},{j,n+1}]; 
Inter\[Sigma]of\[Xi]f=Interpolation[Flatten[Log10[Inter\[Sigma]of\[Xi]Table],1],InterpolationOrder->4];
Print["\[Sigma] of \[Xi] interpolation done"];

<|"d\[Sigma]f"->Interd\[Sigma]f,"Interpoints"->Interpointsv\[Chi]and\[Xi],"Interpointsv\[Chi]"->Interpointsv\[Chi],"\[Sigma]f"->Inter\[Sigma]f,"dEdlf"->InterdEdlf,"\[Sigma]of\[Xi]f"->Inter\[Sigma]of\[Xi]f,"\[Sigma]of\[Xi]table"->Inter\[Sigma]of\[Xi]Table,"m\[Chi]"->m\[Chi],"v\[Chi]"->Log10[v\[Chi]],"\[Xi]"->Log10[{\[Epsilon],1-\[Epsilon]}],"params"->params,"\[Omega]min"->\[Omega]min,"regionindex"->regionindex,"materialname"->materialname|>
]*)


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
\[Omega]of\[Xi]andv\[Chi][\[Xi]_?NumericQ,\[Omega]min_,v\[Chi]_,m\[Chi]_,params_]:=\[Omega]min + \[Xi] (\[Omega]maxofv\[Chi][v\[Chi],m\[Chi],params] - \[Omega]min)
\[Xi]of\[Omega]andv\[Chi][\[Omega]_,\[Omega]min_,v\[Chi]_,m\[Chi]_,params_]:=(\[Omega]-\[Omega]min )/(\[Omega]maxofv\[Chi][v\[Chi],m\[Chi],params] - \[Omega]min)


(* ::Text:: *)
(*Interpolate over the v\[Chi] as well - Electronic*)


Getv\[Chi]and\[Xi]List[v\[Chi]_,n_:20,\[Epsilon]_:10^-4,uselog_:False]:=Module[{Interpoints\[Xi]},
Interpoints\[Xi]=If[uselog,10^Subdivide[Log10[\[Epsilon]],Log10[1-\[Epsilon]],n],Subdivide[\[Epsilon],1-\[Epsilon],n]];(*move slightly away from the boundaries so we can do a log-log interpolation (d\[Sigma]/Subscript[dE, R]->0 at Subscript[\[Omega], max], and \[Xi]->0 at Subscript[\[Omega], min])*)
(*Interpoints\[Xi]=10^Subdivide[Log10[\[Epsilon]],Log10[1-\[Epsilon]],n];*)
(*Interpoints\[Xi]=Subdivide[\[Epsilon],1-\[Epsilon],n];*)(*move slightly away from the boundaries so we can do a log-log interpolation (d\[Sigma]/Subscript[dE, R]->0 at Subscript[\[Omega], max], and \[Xi]->0 at Subscript[\[Omega], min])*)
Table[{v\[Chi],Interpoints\[Xi][[i]]},{i,n+1}]
]


Interpolatev\[Chi]and\[Omega]\[Sigma]dERe\[Xi][m\[Chi]_:0.5 10^6 ("JpereV")/("c")^2/.Constants`SIConstRepl,params_,m_:5,n_:60,\[Epsilon]_:10^-6,cut_:4,uselog_:False,get\[Sigma]too_,regionindex_:1,materialname_:"material",Truncatev\[Chi]_:False]:=Module[{\[Omega]min,v\[Chi],Interpointsv\[Chi],Interpointsv\[Chi]and\[Xi],Interd\[Sigma]Table,Interd\[Sigma]f,Inter\[Sigma]Table,Inter\[Sigma]f,InterdEdlTable,InterdEdlf,Inter\[Sigma]of\[Xi]Table,Inter\[Sigma]of\[Xi]f,v\[Chi]maxindexfor\[Sigma]int},
(*Interpolate over d\[Sigma]dER for electronic contribution (to avoid doing the integral over momentum transfer at every evaluation and to speed up the numerical integration over it to find the normalization)*)

\[Omega]min= "\[Omega]edgei"/.params; (*the cutoff used in the data fitting*)

v\[Chi] = {3Sqrt[("\[HBar]" \[Omega]min)/(2 m\[Chi])], 2 10^-1 "c"}/.Constants`SIConstRepl;(*interpolate from the minimum velocity (up to an order 1 fudge factor (3) for numerical stability) to the boundary of the non-relativistic regime*)

If[v\[Chi][[1]]>v\[Chi][[2]],Print["minimum energy requires relativistic probes, we won't consider these"];Return["NA",Module]];

Interpointsv\[Chi]=10^Subdivide[Log10[v\[Chi][[1]]],Log10[v\[Chi][[2]]],m];
Interpointsv\[Chi]and\[Xi]=SortBy[Table[Getv\[Chi]and\[Xi]List[Interpointsv\[Chi][[i]],n,\[Epsilon],uselog],{i,m+1}],Last];

(*Interpolate to get d\[Sigma]/dER*)
Interd\[Sigma]Table=Table[{Interpointsv\[Chi]and\[Xi][[i,j]],d\[Sigma]dERe[\[Omega]of\[Xi][Interpointsv\[Chi]and\[Xi][[i,j]][[2]],\[Omega]min,\[Omega]maxofv\[Chi][Interpointsv\[Chi]and\[Xi][[i,j]][[1]],m\[Chi],params]],m\[Chi],Interpointsv\[Chi]and\[Xi][[i,j]][[1]],{params},"\[Beta]"/.params,cut]},{i,m+1},{j,n+1}];
Interd\[Sigma]f=Interpolation[Flatten[Log10[Interd\[Sigma]Table],1],InterpolationOrder->4];
Print["d\[Sigma] interpolation done"];

Inter\[Sigma]Table =Table[{Interpointsv\[Chi][[i]],("\[HBar]"/.params)(\[Omega]maxofv\[Chi][Interpointsv\[Chi][[i]],m\[Chi],params]-\[Omega]min)NIntegrate[10^Interd\[Sigma]f[Log10[Interpointsv\[Chi][[i]]],Log10[\[Xi]]],{\[Xi],\[Epsilon],1-\[Epsilon]},Method->{"LocalAdaptive"}]},{i,m+1}]; (*The prefactor on the Integral is the Jacobian of the transformation E_R -> \[Xi]*)
Inter\[Sigma]f = Interpolation[Log10[Inter\[Sigma]Table],InterpolationOrder->4];
Print["\[Sigma] interpolation done"];

InterdEdlTable =Table[{Interpointsv\[Chi][[i]],("ne"/.params)("\[HBar]"/.params)^2 (\[Omega]maxofv\[Chi][Interpointsv\[Chi][[i]],m\[Chi],params]-\[Omega]min)NIntegrate[\[Omega]of\[Xi]andv\[Chi][\[Xi],\[Omega]min,Interpointsv\[Chi][[i]],m\[Chi],params] 10^Interd\[Sigma]f[Log10[Interpointsv\[Chi][[i]]],Log10[\[Xi]]],{\[Xi],\[Epsilon],1-\[Epsilon]},Method->{"LocalAdaptive"}]},{i,m+1}]; (*The prefactor on the Integral is the Jacobian of the transformation E_R -> \[Xi]*)
InterdEdlf= Interpolation[Log10[InterdEdlTable],InterpolationOrder->4];
Print["dEdlf interpolation done"];

v\[Chi]maxindexfor\[Sigma]int=If[Truncatev\[Chi],Position[Interpointsv\[Chi],_?(#<10^6&)][[-1,1]],m+1];

Inter\[Sigma]of\[Xi]Table=Table[{Interpointsv\[Chi]and\[Xi][[i,j]],10^-49+("\[HBar]"/.params)(\[Omega]maxofv\[Chi][Interpointsv\[Chi]and\[Xi][[i,j]][[1]],m\[Chi],params]-\[Omega]min)Abs[Re[NIntegrate[10^Interd\[Sigma]f[Log10[Interpointsv\[Chi]and\[Xi][[i,j]][[1]]],Log10[\[Xi]]],{\[Xi],Interpointsv\[Chi]and\[Xi][[i,j]][[2]],1-\[Epsilon]},Method->{"LocalAdaptive"},PrecisionGoal->3,MaxRecursion->8,AccuracyGoal->3]]]},{i,v\[Chi]maxindexfor\[Sigma]int},{j,n+1}]; 
Inter\[Sigma]of\[Xi]f=Interpolation[Flatten[Log10[Inter\[Sigma]of\[Xi]Table],1],InterpolationOrder->4];

Print["\[Sigma] of \[Xi] interpolation done"];

<|"d\[Sigma]f"->Interd\[Sigma]f,"Interpoints"->Interpointsv\[Chi]and\[Xi],"Interpointsv\[Chi]"->Interpointsv\[Chi],"\[Sigma]f"->Inter\[Sigma]f,"dEdlf"->InterdEdlf,"\[Sigma]of\[Xi]f"->Inter\[Sigma]of\[Xi]f,"\[Sigma]of\[Xi]table"->Inter\[Sigma]of\[Xi]Table,"m\[Chi]"->m\[Chi],"v\[Chi]"->Log10[v\[Chi]],"\[Xi]"->Log10[{\[Epsilon],1-\[Epsilon]}],"params"->params,"\[Omega]min"->\[Omega]min,"regionindex"->regionindex,"materialname"->materialname|>
]


(* ::Text:: *)
(*Interpolate over the v\[Chi] as well - Nuclear*)


\[Omega]maxNuc[v\[Chi]_,m\[Chi]_,mN_,params_]:= (2 v\[Chi]^2)/("\[HBar]") (m\[Chi]^2 mN)/(m\[Chi] + mN)^2/.params
\[Omega]of\[Xi]Nuc[\[Xi]_?NumericQ,\[Omega]min_,v\[Chi]_,m\[Chi]_,mN_,params_]:=\[Omega]min + \[Xi] (\[Omega]maxNuc[v\[Chi],m\[Chi],mN,params] - \[Omega]min)
\[Xi]of\[Omega]Nuc[\[Omega]_,\[Omega]min_,v\[Chi]_,m\[Chi]_,mN_,params_]:=(\[Omega]-\[Omega]min )/(\[Omega]maxNuc[v\[Chi],m\[Chi],mN,params] - \[Omega]min)


(* ::Text:: *)
(***** should turn these into pure functions and put them in the module*)


Interpolatev\[Chi]and\[Omega]\[Sigma]dERNuc\[Xi][m\[Chi]_:0.5 10^6 ("JpereV")/("c")^2/.Constants`SIConstRepl,FFcoeffs_,regionindex_:1,materialname_:"material",\[Epsilon]_:10^-6,Truncatev\[Chi]_:False,m_:20,n_:60,v\[Chi]_:{10^-4 "vesc",2 10^-1"c"}/.Constants`SIConstRepl/.Constants`EarthRepl,params_:Constants`SIConstRepl]:=Module[{\[Omega]min,mN,Interpointsv\[Chi],Interpointsv\[Chi]and\[Xi],Interd\[Sigma]Table,Interd\[Sigma]f,Inter\[Sigma]Table,Inter\[Sigma]f,InterdEdlTable,InterdEdlf,Inter\[Sigma]of\[Xi]Table,Inter\[Sigma]of\[Xi]f,v\[Chi]maxindexfor\[Sigma]int},
(*Interpolate over d\[Sigma]dER and \[Sigma] for Nuclear contribution (to avoid doing the integral over momentum transfer at every evaluation and to speed up the numerical integration over it to find the normalization)*)

\[Omega]min=0;
mN = "mN"/.FFcoeffs;

Interpointsv\[Chi]=10^Subdivide[Log10[v\[Chi][[1]]],Log10[v\[Chi][[2]]],m];
Interpointsv\[Chi]and\[Xi]=SortBy[Table[Getv\[Chi]and\[Xi]List[Interpointsv\[Chi][[i]],n,\[Epsilon],True],{i,m+1}],Last];

(*Interpolate to get d\[Sigma]/dER*)
Interd\[Sigma]Table=Table[{Interpointsv\[Chi]and\[Xi][[i,j]],Max[d\[Sigma]dERNuc0[\[Omega]of\[Xi][Interpointsv\[Chi]and\[Xi][[i,j]][[2]],\[Omega]min,\[Omega]maxNuc[Interpointsv\[Chi]and\[Xi][[i,j]][[1]],m\[Chi],mN,params]],m\[Chi],Interpointsv\[Chi]and\[Xi][[i,j]][[1]],FFcoeffs,params],10^-50]},{i,m+1},{j,n+1}];
Interd\[Sigma]f=Interpolation[Flatten[Log10[Interd\[Sigma]Table],1],InterpolationOrder->4];
Print["d\[Sigma] interpolation done"];

Inter\[Sigma]Table =Table[{Interpointsv\[Chi][[i]],("\[HBar]"/.params)(\[Omega]maxNuc[Interpointsv\[Chi][[i]],m\[Chi],mN,params]-\[Omega]min)NIntegrate[10^Interd\[Sigma]f[Log10[Interpointsv\[Chi][[i]]],Log10[\[Xi]]],{\[Xi],\[Epsilon],1-\[Epsilon]},Method->"LocalAdaptive"]},{i,m+1}]; (*The prefactor on the Integral is the Jacobian of the transformation E_R -> \[Xi]*)
Inter\[Sigma]f = Interpolation[Log10[Inter\[Sigma]Table],InterpolationOrder->4];
Print["\[Sigma] interpolation done"];

InterdEdlTable =Table[{Interpointsv\[Chi][[i]],("nI"/.FFcoeffs)("\[HBar]"/.params)^2 (\[Omega]maxNuc[Interpointsv\[Chi][[i]],m\[Chi],mN,params]-\[Omega]min)NIntegrate[\[Omega]of\[Xi]Nuc[\[Xi],\[Omega]min,Interpointsv\[Chi][[i]],m\[Chi],mN,params] 10^Interd\[Sigma]f[Log10[Interpointsv\[Chi][[i]]],Log10[\[Xi]]],{\[Xi],\[Epsilon],1-\[Epsilon]},Method->"LocalAdaptive"]},{i,m+1}]; (*The prefactor on the Integral is the Jacobian of the transformation E_R -> \[Xi]*)
InterdEdlf= Interpolation[Log10[InterdEdlTable],InterpolationOrder->4];
Print["dEdlf interpolation done"];

v\[Chi]maxindexfor\[Sigma]int=If[Truncatev\[Chi],Position[Interpointsv\[Chi],_?(#<10^6&)][[-1,1]],m+1];

Inter\[Sigma]of\[Xi]Table=Table[{Interpointsv\[Chi]and\[Xi][[i,j]],10^-49+("\[HBar]"/.params)(\[Omega]maxNuc[Interpointsv\[Chi]and\[Xi][[i,j]][[1]],m\[Chi],mN,params]-\[Omega]min)Abs[Re[NIntegrate[10^Interd\[Sigma]f[Log10[Interpointsv\[Chi]and\[Xi][[i,j]][[1]]],Log10[\[Xi]]],{\[Xi],Interpointsv\[Chi]and\[Xi][[i,j]][[2]],1-\[Epsilon]},Method->"LocalAdaptive"]]]},{i,v\[Chi]maxindexfor\[Sigma]int},{j,n+1}]; 
Inter\[Sigma]of\[Xi]f=Interpolation[Flatten[Log10[Inter\[Sigma]of\[Xi]Table],1],InterpolationOrder->4];
(*Inter\[Sigma]of\[Xi]f=Interpolation[Flatten[Inter\[Sigma]of\[Xi]Table,1],InterpolationOrder->4];*)
Print["\[Sigma] of \[Xi] interpolation done"];

<|"d\[Sigma]f"->Interd\[Sigma]f,"Interpoints"->Interpointsv\[Chi]and\[Xi],"\[Sigma]f"->Inter\[Sigma]f,"dEdlf"->InterdEdlf,"\[Sigma]of\[Xi]f"->Inter\[Sigma]of\[Xi]f,"Interpointsv\[Chi]"->Interpointsv\[Chi],"m\[Chi]"->m\[Chi],"mN"->mN,"v\[Chi]"->Log10[v\[Chi]],"\[Xi]"->Log10[{\[Epsilon],1-\[Epsilon]}],"NucleusParams"->FFcoeffs,"params"->params,"\[Omega]min"->\[Omega]min,"regionindex"->regionindex,"materialname"->materialname|>
]


(* ::Section:: *)
(*Plotting*)


(* ::Subsection:: *)
(*Plot Energy Loss Interpolation Function - Public*)


 
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
