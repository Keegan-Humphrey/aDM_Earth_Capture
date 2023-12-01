(* ::Package:: *)

Needs["Dielectrics`"]


BeginPackage["CollisionRate`"];


(* ::Text:: *)
(*The purpose of this package is to define modules for the calculation of collision rates of electrons in solids analytically*)


(* ::Subsubsection::Closed:: *)
(*Public Declarations*)


(* ::Text:: *)
(*Collision Rates from White Dwarf Paper*)


\[Nu]::usage = "Theoretically calculated collision rate (see overleaf for references)";
u\[Nu]::usage = "Dimensionless version of \[Nu] (for use with \[Epsilon]M)";


(* ::Text:: *)
(*Collision Rate a la Ziman and Marder*)


\[Nu]ZM::usage = "Collision Rate from elastic scattering off nuclei and inelastic single phonon scattering";


(* ::Section::Closed:: *)
(*Private*)


Begin["Private`"];


(* ::Subsubsection::Closed:: *)
(*\[Nu] (Collision Rate)*)


\[Nu][]:= Module[{x,J,\[Omega]px,y,\[Nu]0,\[Nu]},
			(*Courtesy of: https://doi.org/10.1103/PhysRevE.91.023102 , 
							https://doi.org/10.1017/S0263034606060733
				
				- See Overleaf for details*)
			x = "vF"/"c";
			\[Omega]px =Sqrt[(4 Pi "e"^2 "ne")/("m" (1 + x^2)^(1/2)) 1/(4 \[Pi] "\[Epsilon]0")];
			y =Sqrt[3]"\[Beta]" "hbar" \[Omega]px;
		J= (y^3/(3(1+0.07414 y)^3) Log[2.810/y - (0.810 x^2)/(y (1 + x^2))+ 1]+Pi^5/6 y^4/(13.91 + y)^4 )(1 + 6/(5 x^2) + 2/(5 x^4));
			\[Nu]0 = (J  3 )/(2 "hbar" "m" "c"^2  "\[Beta]"^2) Sqrt[("\[Alpha]" x^3)/(Pi^3  (1 + x^2)^(5/2))];
			 \[Nu]=\[Nu]0/Sqrt[1 + 0.2/("\[Beta]" "\[Mu]")];
			\[Nu]
]


u\[Nu][z_]:= Module[{\[Nu]t,u\[Nu]},
			\[Nu]t = \[Nu][];
			u\[Nu] = \[Nu]t/(2 "vF" "qF" z);
			u\[Nu]
]


(* ::Subsubsection::Closed:: *)
(*Ziman Marder Collision Rates*)


(* ::Text:: *)
(*Collision Rate due to elastic scattering off Ions*)


\[Nu]ZimanFaberLin[params_] :=( (4 "m" ("e")^4 ("Z")^2)/(\[Pi] ("\[HBar]")^3 ("\[Epsilon]0")^2)/.params)("nI""qF"/.params)NIntegrate[z^3 (1/((2"qF")^2  z^2 Dielectrics`\[Epsilon]RPA0[z])/.params)^2  ,{z,0,1} ]


(* ::Text:: *)
(*Collision Rate due to single phonon inelastic scattering (temperature dependent)*)


\[Nu]PhononMarderLinT[params_,c_,\[CapitalTheta]_,T_] :=((4\[Pi] ("e")^2 "ne")/("\[Epsilon]0" "m")/.params)( (3 \[Pi] "\[Epsilon]0")/(4\[Pi] (("e")^2) "\[HBar]" (("vF")^2) )/.params)("nI"/.params)(1/(4 ("qF")^4)/.params)(("\[HBar]")/(2 "M" c)/.params)(("kB" T)/("\[HBar]" c )/.params)^5 NIntegrate[z^4/(E^z-1) ((("e")^2 "Z")/("\[Epsilon]0")/(((z "kB" T)/("\[HBar]" c))^2 Dielectrics`\[Epsilon]RPA0[1/(2"qF") ((z "kB" T)/("\[HBar]" c))])/.params)^2  ,{z,0,2 \[CapitalTheta]/T/.params} ]


\[Nu]ZM[params_,c_,\[CapitalTheta]_,T_]:=\[Nu]ZimanFaberLin[params]+\[Nu]PhononMarderLinT[params,c,\[CapitalTheta],T];


(* ::Section::Closed:: *)
(*End*)


End[];


EndPackage[];
