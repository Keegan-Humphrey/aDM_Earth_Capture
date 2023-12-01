(* ::Package:: *)

BeginPackage["Utilities`"];


(* ::Text:: *)
(*Utility Functions used in the code*)


(* ::Subsubsection:: *)
(*Public Declarations*)


(* ::Text:: *)
(*Save and Read Objects to file*)


SaveIt::usage = "Saves an arbitrary object in Mathematica to File ";
ReadIt::usage = "Reads in a saved file";


(* ::Text:: *)
(*Convert string replacement tables to expressions*)


StrToExpConverter::usage = "Convert a string->expr replacement table to expr->expr";
ExpToStrConverter::usage = "Convert a expr->expr replacement table to string->expr";


(* ::Text:: *)
(*Change coefficients in a replacement table*)


RescaleAi::usage = "Rescale Ai coefficient in a replacement table";
RescaleParams::usage = "Rescale coefficients in a replacement table";
AddToParams::usage = "Add to coefficients in a replacement table";
ReplaceParams::usage = "Replace parameter in a repalacement table";


(* ::Text:: *)
(*Contour Plot with fit replacement table*)


ContourTotalParams::usage = "Contour Plot of integrand of EL integrals using fit replacement table";
ContourTotalParams\[Omega]q::usage = "Contour Plot of integrand of EL integrals using fit replacement table in \[Omega] q variables";


(* ::Section:: *)
(*Private*)


Begin["Private`"];


(* ::Subsubsection::Closed:: *)
(*Parallel / best practices tips*)


(* ::Text:: *)
(*Monitor - keeps track*)
(*use Evaluation > Kernel ... to add, and choose kernels to use, this lets use run processes in parallel*)
(*change sleep time to let mathematica processes run for longer*)
(*Note - doesn't work when the keyboard is closed*)


(* ::Subsubsection::Closed:: *)
(*Save and Read Functions*)


(* ::Text:: *)
(*Many thanks to David Curtin for these super useful modules! :)*)


SaveIt[filename_,expr_]:=Module[{output},output=Export[filename<>".dat",ToString[expr//InputForm],"String"];
output];
SaveIt[varnamestring_]:=Module[{output},output=Export[varnamestring<>".dat",ToString[ToExpression[varnamestring]//InputForm],"String"];
output];


ReadIt[filename_]:=Module[{output},output=ToExpression[Import[StringReplace[filename,".dat"->""]<>".dat","String"]];
output
];


(* ::Subsubsection::Closed:: *)
(*String vs Expression Replacement Tables*)


StrToExpConverter[table_]:=Module[{temptable},
		temptable = table /. Rule -> List;
		Do[temptable[[i,1]] = ToExpression[temptable[[i,1]]];
			temptable[[i]]=temptable[[i]]/.List->Rule,
			{i,Length[temptable]}
		];
		temptable
]


ExpToStrConverter[table_]:=Module[{temptable},
		temptable = table /. Rule -> List;
		Do[temptable[[i,1]] = ToString[temptable[[i,1]]];
			temptable[[i]]=temptable[[i]]/.List->Rule,
			{i,Length[temptable]}
		];
		temptable
]


(* ::Subsubsection:: *)
(*Rescale Fit Coefficients*)


RescaleAi[params_,Cr_]:=Module[{list,Ailocation},
	(*
	params - parameter replacement table with string keys
	Cr - rescaling of coordinates (x)
	rescalepar - entry to be rescaled
	*)
	list = params/. Rule->List;
	Ailocation=Position[list,"Ai"][[1,1]];
	list[[Ailocation,2]]*=Cr;
	Table[list[[i]]/.List->Rule,{i,Length[list]}]
]

RescaleParams[params_,Cr_,par_]:=Module[{list,location},
	(*
	params - parameter replacement table with string keys
	Cr - rescaling of coordinates (x)
	par - entry to be rescaled
	*)
	list = params/. Rule->List;
	location=Position[list,par][[1,1]];
	list[[location,2]]*=Cr;
	Table[list[[i]]/.List->Rule,{i,Length[list]}]
]

AddToParams[params_,Cr_,par_]:=Module[{list,location},
	(*
	params - parameter replacement table with string keys
	Cr - rescaling of coordinates (x)
	par - entry to be added to
	*)
	list = params/. Rule->List;
	location=Position[list,par][[1,1]];
	list[[location,2]]+=Cr;
	Table[list[[i]]/.List->Rule,{i,Length[list]}]
]

ReplaceParams[params_,ent_,par_]:=Module[{list,location},
	(*
	params - parameter replacement table with string keys
	ent - new entry for the par 
	par - entry to be replaced
	*)
	list = params/. Rule->List;
	location=Position[list,par][[1,1]];
	list[[location,2]]=ent;
	Table[list[[i]]/.List->Rule,{i,Length[list]}]
]


(* ::Subsection:: *)
(*Plotting*)


(* ::Subsubsection:: *)
(*Plot \[Epsilon] contour via total params*)


ContourTotalParams[strtotalparams_]:=Module[{},
(*expects argument to be string replacement table*)
Print[ContourPlot[Log10[Sum[( ("Ai")/("JpereV") (("qF" "vF")/("c"))^2 8 /Pi  ("e")^2/(4 \[Pi] "\[Epsilon]0")/.strtotalparams[[i]])Im[(-u z)/Dielectrics`\[Epsilon]MNum[u,z,("\[Nu]i")/(2 "vF" "qF" z)/.strtotalparams[[i]],strtotalparams[[i]]]],{i,Length[strtotalparams]}]],{z,0.1,10},{u,0.1,10},FrameLabel->{{"u",""},{"z","\!\(\*SubscriptBox[\(Log\), \(10\)]\)[(Loss Integrand)\!\(\*SuperscriptBox[SubscriptBox[\(\[Beta]\), \(\[Chi]\)], \(2\)]\)]"}},ContourLabels->True,Axes->True]]
]


ContourTotalParams\[Omega]q[strtotalparams_,material_:""]:=Module[{},
(*expects argument to be string replacement table*)
Print[ContourPlot[EnergyLoss`LossIntegrand[strtotalparams,q,\[Omega] (("\[HBar]")/("JpereV"))^-1/.Constants`SIConstRepl],{q,0.1 10^10,10  10^10},{\[Omega],0.1,50},Contours->{9,8,7,6,5,4,3,2,1},FrameLabel->{{"\[Omega] [eV]",""},{"q [\!\(\*SuperscriptBox[\(m\), \(-1\)]\)]",StringForm["\!\(\*SubscriptBox[\(Log\), \(10\)]\)[(Loss Integrand)\!\(\*SuperscriptBox[SubscriptBox[\(\[Beta]\), \(\[Chi]\)], \(2\)]\)\!\(\*SubscriptBox[\(]\), \(``\)]\)",material]}},ContourLabels->True,Axes->True]]
]


(* ::Section::Closed:: *)
(*End*)


End[];


EndPackage[];
