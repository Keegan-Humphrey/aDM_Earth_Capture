(* ::Package:: *)

BeginPackage["Kappamin`"];


(* ::Subsubsection:: *)
(*Public Declarations*)


(* ::Text:: *)
(*Interpolate and compute \[Kappa]min*)


\[Kappa]minInterpolationfnofr::usage = "compute \[Kappa]min using output of EnergyLoss`EnergyLossTableAndInter*";


(* ::Section:: *)
(*Private*)


Begin["`Private`"];


(* ::Subsubsection::Closed:: *)
(*Earth Measurements*)


(*REarth = 6.371 10^6;
rcore = REarth - 2.885 10^6;*)


(* ::Subsubsection::Closed:: *)
(*Evolve Energy Through the Earth*)


(* ::Text:: *)
(*Where f is not a function of r*)


Eat2R[m\[Chi]_,v\[Chi]_,f_,\[Kappa]_,plot_:False,rmax_:2 Constants`REarth]:=Module[{},
		(*{rmax=2 REarth},*)
		(* rmax - [m] default is earth's radius (times 2)
		Function of \[Kappa], solves an ODE for E(r), returning the value at the Earth's diameter. Used this to find \[Kappa] min (given that E increase when it undershoots)
		*)
		sol=NDSolve[{Ek'[r]==-\[Kappa]^2 10^f[m\[Chi],Log[10,Sqrt[Max[2 10^-m\[Chi] Ek[r],(10^((f@Domain[])[[2,1]]))^2]]]],Ek[0]==1/2 10^m\[Chi](10^v\[Chi])^2 },Ek,{r,0,2rmax}];
		If[plot, 
			Print[Plot[(Ek/.First[sol])[r],{r,0,2rmax},AxesLabel->{"r [m]","E\[Chi](r) [eV]"}]]
		];
		Ek[r]/.First[sol]/.r->2rmax
] (*WITHOUT RMAX*)


(* ::Text:: *)
(*Where f is a function of r*)


Eat2Rfnofr[m\[Chi]_,v\[Chi]_,f_,\[Kappa]_,domain_,plot_:False,rmax_:Constants`REarth]:=Module[{},
		(* rmax - [m] default is earth's radius
		Function of \[Kappa], solves an ODE for E(r), returning the value at the Earth's diameter. Used this to find \[Kappa] min (given that E increase when it undershoots)
		*)
		sol=NDSolve[{Ek'[r]==-\[Kappa]^2 ("JpereV"/.Constants`SIConstRepl) 10^f[m\[Chi],Log[10,Sqrt[Max[2/10^m\[Chi] Ek[r],(10^((domain)[[2,1]]))^2]]],r],Ek[0]==1/2 10^m\[Chi]  (10^v\[Chi])^2 },Ek,{r,0,2 rmax}];
		If[plot, 
			Print[Plot[(Ek/.First[sol])[r]/("JpereV"/.Constants`SIConstRepl),{r,0,2rmax},AxesLabel->{"r [m]","E\[Chi](r) [eV]"}]];
			Print[Plot[f[m\[Chi],Log[10,Sqrt[Max[2/10^m\[Chi](Ek/.First[sol])[r],(10^((domain)[[2,1]]))^2]]],r],{r,0,2rmax},AxesLabel->{"r [m]","\!\(\*FractionBox[\(-1\), SuperscriptBox[\(\[Kappa]\), \(2\)]]\)\!\(\*FractionBox[\(dE\), \(dr\)]\)(r) [eV \!\(\*SuperscriptBox[\(m\), \(-1\)]\)]"}]]
		];
		Ek[r]/.First[sol]/.r->2 rmax
]


(* ::Subsubsection::Closed:: *)
(*Binary Search to find \[Kappa]min*)


(* ::Text:: *)
(*Where  f is not a function of r*)


\[Kappa]min[m\[Chi]_,v\[Chi]_,f_,domain_,vmin_:0,n_:20,ofmmin_:-20,ofmmax_:-1]:=Module[{\[Kappa],interval,mid},
			(*This function finds the root of E(2R,\[Kappa]) via binary search with a fixed number of iterations
			Input:
				ofmmin - the order of magnitude of \[Kappa] to initialize the search (\[Kappa] >= 10^offmin)
			  ofmmax - the maximum order of magnitude of \[Kappa] 
			  n - number of iterations

			Output:
					\[Kappa]min - minimum value of \[Kappa] for capture of DM with given kinematics in the earth
					E(2R,\[Kappa]) - the value of the energy at the given precision of \[Kappa] 
			*)

		interval={ofmmin,ofmmax};
		mid=0;
		 Do[mid = (interval[[2]]+interval[[1]])/2;
			interval =If[Eat2R[m\[Chi],v\[Chi],f,10^mid,domain]- 1/2 m\[Chi] vmin^2>0,{mid,interval[[2]]},{interval[[1]],mid}],
			{i,n}
		];
		\[Kappa] = N[10^mid];
			{N[\[Kappa]],Eat2R[m\[Chi],v\[Chi],f,\[Kappa],domain]}
	]


(* ::Text:: *)
(*Where f is a function of r*)


\[Kappa]minfnofr[m\[Chi]_,v\[Chi]_,f_,domain_,vmin_:0,n_:20,ofmmin_:-20,ofmmax_:-1]:=Module[{\[Kappa],interval,mid},
			
			(*This function finds the root of E(2R,\[Kappa]) via binary search with a fixed number of iterations
			Input:
				ofmmin - the order of magnitude of \[Kappa] to initialize the search (\[Kappa] >= 10^offmin)
			  ofmmax - the maximum order of magnitude of \[Kappa] 
			  n - number of iterations
			  vmin - [m s^-1] velocity being searched for at 2 R (eg. escape velocity on surface of earth)

			Output:
					\[Kappa]min - [] minimum value of \[Kappa] for capture of DM with given kinematics in the earth
					E(2R,\[Kappa]) - [] fraction of energy different from target
			*)

		interval={ofmmin,ofmmax};
		mid=0;
		 Do[mid = (interval[[2]]+interval[[1]])/2;
			interval =If[Eat2Rfnofr[m\[Chi],v\[Chi],f,10^mid,domain]- 1/2 10^m\[Chi] vmin^2>0,{mid,interval[[2]]},{interval[[1]],mid}],
			{i,n}
		];
		\[Kappa] = N[10^mid];
		{N[\[Kappa]],(Eat2Rfnofr[m\[Chi],v\[Chi],f,\[Kappa],domain]-1/2 10^m\[Chi] vmin^2)/(1/2 10^m\[Chi] (10^v\[Chi])^2)}
	]


(* ::Subsubsection:: *)
(*Interpolation table for \[Kappa]min - Public*)


\[Kappa]minInterpolationfnofr[kindict_,f_,vescape_,domain_,plot_:False,title_:"Earth"] := Module[{table,interptable,interpolation},
		(* f - fn of r to run over returning Log10[EL] (taking non-log arguments)
			kindict - *)
table = Table[Table[First[\[Kappa]minfnofr[Log[10,kindict[["kinlist"]][["m\[Chi]"]][[i]]],Log[10,kindict[["kinlist"]][["v\[Chi]"]][[j]]],f,domain,vescape]],{i,kindict[["meshdims"]][["m\[Chi]"]]}],{j,kindict[["meshdims"]][["v\[Chi]"]]}];
	
	interptable = Table[{{Log[10,kindict[["m\[Chi]Mesh"]][[i,j]]],Log[10,kindict[["v\[Chi]Mesh"]][[i,j]]]},Log[10,table[[j,i]]]},{i,kindict[["meshdims"]][["m\[Chi]"]]
},{j,kindict[["meshdims"]][["v\[Chi]"]]}];	
	interpolation =  Interpolation[Flatten[interptable,1]];
	If[plot,
Print[Show[{DensityPlot[interpolation[m+6-Log10[("c"^2)/("JpereV")/.Constants`SIConstRepl],v+Log10["c"/.Constants`SIConstRepl]],
	{m,Log[10,First[kindict[["kinlist"]][["m\[Chi]"]]]]-6+Log10[("c"^2)/("JpereV")/.Constants`SIConstRepl],Log[10,Last[kindict[["kinlist"]][["m\[Chi]"]]]]-6+Log10[("c"^2)/("JpereV")/.Constants`SIConstRepl]},
		{v,Log[10,First[kindict[["kinlist"]][["v\[Chi]"]]]]-Log10["c"/.Constants`SIConstRepl],Log[10,Last[kindict[["kinlist"]][["v\[Chi]"]]]]-Log10["c"/.Constants`SIConstRepl]},PlotRange->All,FrameLabel->{"\!\(\*SubscriptBox[\(Log\), \(10\)]\)[m\[Chi]] [MeV]","\!\(\*SubscriptBox[\(Log\), \(10\)]\)[v\[Chi]] [c]"},PlotLabel->StringForm["\!\(\*SubscriptBox[\(Log\), \(10\)]\)[\!\(\*SubscriptBox[\(\[Kappa]\), \(min\)]\)] [] - ``",title]],
ContourPlot[interpolation[m+6- Log10[("c"^2)/("JpereV")/.Constants`SIConstRepl],v+Log10["c"/.Constants`SIConstRepl]],{m,Log[10,First[kindict[["kinlist"]][["m\[Chi]"]]]]-6+Log10[("c"^2)/("JpereV")/.Constants`SIConstRepl],Log[10,Last[kindict[["kinlist"]][["m\[Chi]"]]]]-6+Log10[("c"^2)/("JpereV")/.Constants`SIConstRepl]},
	{v,Log[10,First[kindict[["kinlist"]][["v\[Chi]"]]]]-Log10["c"/.Constants`SIConstRepl],Log[10,Last[kindict[["kinlist"]][["v\[Chi]"]]]]-Log10["c"/.Constants`SIConstRepl]},PlotRange->All,ContourLabels->True,ContourShading->False]}]]
	];
	
	<|"table"->table,"f\[Kappa]min"->interpolation,"ELdict"->kindict,"fEL"->f,"vescape"->vescape|>
]


(* ::Section::Closed:: *)
(*End*)


End[];


EndPackage[];
