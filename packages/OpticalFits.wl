(* ::Package:: *)

BeginPackage["OpticalFits`"];


(* ::Subsubsection:: *)
(*Public Declarations*)


Aofq::usage = "q dependent Ai coefficients used with truncation fitting";



(* ::Section:: *)
(*Private*)


Begin["`Private`"];


(* ::Subsection:: *)
(*Truncation coefficients*)


(* ::Text:: *)
(***** These need to be changed to "" variables*)


(* ::Subsubsection::Closed:: *)
(*Optical Integral*)


(* ::Text:: *)
(*Comments give derivation of the below expression*)


(*Integrate[ELFM0[\[Omega],\[Nu],\[Omega]p],\[Omega]]*)


(*Integrate[\[Omega]edgei ELFM0[\[Omega]edgei,\[Nu]i,\[Omega]i],\[Omega]edgei]
Limit[%,\[Omega]edgei->\[Infinity]]*)


(* ::Text:: *)
(*Optical integral used in finite q Ai coefficients*)


AqOpticalIntegral[params_]:= (("\[Omega]i")^2 (-(\[Pi]/(2 Sqrt[1/(("\[Nu]i")^2-2 ("\[Omega]i")^2-"\[Nu]i" Sqrt[("\[Nu]i")^2-4 ("\[Omega]i")^2])]))+\[Pi]/(2 Sqrt[1/(("\[Nu]i")^2-2 ("\[Omega]i")^2+"\[Nu]i" Sqrt[("\[Nu]i")^2-4 ("\[Omega]i")^2])])))/(Sqrt[2] Sqrt[("\[Nu]i")^2-4 ("\[Omega]i")^2]) -1/(Sqrt[2] Sqrt[("\[Nu]i")^2-4 ("\[Omega]i")^2]) ("\[Omega]i")^2 (-Sqrt[("\[Nu]i")^2-2 ("\[Omega]i")^2-"\[Nu]i" Sqrt[("\[Nu]i")^2-4 ("\[Omega]i")^2]] ArcTan[(Sqrt[2] "\[Omega]edgei")/Sqrt[("\[Nu]i")^2-2 ("\[Omega]i")^2-"\[Nu]i" Sqrt[("\[Nu]i")^2-4 ("\[Omega]i")^2]]]  +Sqrt[("\[Nu]i")^2-2 ("\[Omega]i")^2+"\[Nu]i" Sqrt[("\[Nu]i")^2-4 ("\[Omega]i")^2]] ArcTan[(Sqrt[2] "\[Omega]edgei")/Sqrt[("\[Nu]i")^2-2 ("\[Omega]i")^2+"\[Nu]i" Sqrt[("\[Nu]i")^2-4 ("\[Omega]i")^2]]])/. params


(* ::Subsubsection::Closed:: *)
(*Finite q Integral*)


(*u\[Omega]edge[z_,params_] := ( "\[Omega]edgei")/(2 "vF" "qF" z)/.params
u\[Nu]fit[z_,params_] := ( "\[Nu]i")/(2 "vF" "qF" z)/.params*)


(* ::Text:: *)
(*Finite q integral cutoff at high u for numerical stability (see overleaf)*)


Truncated\[Omega]Integral[z_,params_,cutoff_:1000]:=Module[{u\[Omega]edge,u\[Nu]Fit},
	u\[Omega]edge = ( "\[Omega]edgei")/(2 "vF" "qF" z)/.params;
	u\[Nu]Fit = ( "\[Nu]i")/(2 "vF" "qF" z)/.params;
	((2 z "qF" "vF")^2/.params)(NIntegrate[u Im[-1/Dielectrics`\[Epsilon]MNum[u,z,(u\[Nu]Fit),params]] 
	,{u,(u\[Omega]edge),cutoff}]+(("\[Nu]i" ("\[Chi]")^2)/(6 cutoff "qF" "vF" z^3) /.params)) 
	]


(* ::Subsubsection::Closed:: *)
(*Finite q Ai coefficients*)


Aofq[z_,params_] :=  ((("Ai"/.params) AqOpticalIntegral[params])/Truncated\[Omega]Integral[z,params,1000]) 


(* ::Section:: *)
(*End*)


End[];


EndPackage[];
