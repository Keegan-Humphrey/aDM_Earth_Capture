(* ::Package:: *)

Needs["Constants`"]


BeginPackage["Capture`"];


(* ::Subsection:: *)
(*Public Declarations*)


(* ::Subsubsection:: *)
(*Utilities*)


Profiler::usage = "";


v0to\[Beta]D::usage = "";


CrustOrCore::usage = "";
SortProcessesByRegion::usage = "";
Region\[CapitalTheta]::usage = "";
ltor::usage = "";


ExportToDir::usage = "";


PlotPDist::usage = "";
PlotPDistwPf::usage = "";


PlotKEchange::usage = "";
PlotKEchange\[Kappa]Contour::usage = "";


PlotNc\[Kappa]Contour::usage = "";


PlotCaptureRate\[Kappa]Contour::usage = "";


PlotEvaporationRate\[Kappa]Contour::usage = "";


PlotELfrac::usage = "";
PlotELNucoverE::usage = "";


PlotCaptureRateComparison::usage = "";


(* ::Subsubsection:: *)
(*MC Processes*)


Vgrav::usage = "";


dc::usage = "";
dE::usage = "";
dM::usage = "";
GetTotalELFunction::usage = "";
GeteELFunction::usage = "";
GetNucELFunction::usage = "";
EK::usage = "";
\[CapitalDelta]rTotby\[Kappa]squared::usage = "";


GetWeakCaptureTrajectoryinv::usage = "";
GetWeakCaptureTrajectoryinvFE::usage = "";


(* ::Subsubsection:: *)
(*Weak Capture Probability*)


GetWeakCaptureProbabilityfromv::usage = "";


GetWeakCaptureProbabilityNucfromv::usage = "";


GetWeakCaptureInteractionLength::usage = "";
GetTotaleWeakCaptureInteractionLength::usage = "";


GetWeakCaptureInteractionLengthNuc::usage = "";
GetTotalNucWeakCaptureInteractionLength::usage = "";


(* ::Subsubsection:: *)
(*Initial Conditions*)


Getv\[Chi]Max::usage = "";


Getv\[Chi]atEarth::usage = "";


(* ::Subsubsection:: *)
(*Get Capture Probability*)


GetCaptureProbability::usage = "";
IterateGetCaptureProbability::usage = "";


GetCaptureRate::usage = "";
IterateGetCaptureRate::usage = "";


InterpolatePc::usage = "";


(* ::Subsubsection:: *)
(*Get Total Capture Rate*)


nD0\[CapitalLambda]CDM::usage = "";
naDM::usage = "";


GetNcMax::usage = "";


GetEvaporationRate::usage = "";
IterateGetEvaporationRate::usage = "";


GetEvaporationRateend::usage = "";
IterateGetEvaporationRateend::usage = "";


GetPenetrationDepthfnofv::usage = "";
GetPenetrationDepth::usage = "";
EvaporationwMFP::usage = "";
EvaporationVolume::usage = "";


GetinvMFT::usage = "";
GetTotalinvMFT::usage = "";
GetPejection::usage = "";
GetAppendixEvaporationRate::usage = "";
GetTotalAppendixEvaporationRate::usage = "";


GetNc::usage = "";


ScanGetCaptureRateandGetNc::usage = "";


ScanGetCapture::usage = "";


(* ::Subsubsection:: *)
(*Get Significant Capture Regimes*)


VCoulombE::usage = "";
KEchange::usage = "";
IterateKEchange::usage = "";


(* ::Chapter:: *)
(*Public*)


Begin["`Private`"];


(* ::Subsection::Closed:: *)
(*Utilities*)


(* ::Subsubsection::Closed:: *)
(*Profiler*)


(* ::Input:: *)
(*Clear[Profiler]*)
(*Profiler[fun_,inputs_List]:=TableForm[#&/@Join@@@({AbsoluteTiming[fun@@#],#}&/@inputs)]*)


(* ::Subsubsection::Closed:: *)
(*Convert from v0 to \[Beta]D*)


(* ::Input::Initialization:: *)
Clear[v0to\[Beta]D]
v0to\[Beta]D[v0_,meD_,mpD_]:=Module[{\[Beta]D,vpDmean,veDmean},
(*
v0 is the velocity corresponding to the temperature of the aDM particles. By equipartition: 3 Subscript[T, D] = 1/2(Subscript[m, Subscript[e, D]]Subscript[v^2, Subscript[e, D]]+Subscript[m, Subscript[p, D]]Subscript[v^2, Subscript[p, D]]) so we define 3 Subscript[T, D] = Overscript[m, _] Subscript[v^2, 0] 
*)
\[Beta]D = 6/((meD + mpD)v0^2);
vpDmean = Sqrt[3/(mpD \[Beta]D)];
veDmean = Sqrt[3/(meD \[Beta]D)]; 

<|"\[Beta]D"->\[Beta]D,"v0"->v0,"vpDmean"->vpDmean,"veDmean"->veDmean|>
]


(* ::Subsubsection::Closed:: *)
(*Region handling*)


(* ::Input::Initialization:: *)
CrustOrCore[r_]:=Piecewise[{{1,("rcore"/.Constants`EarthRepl)<r<=("rE"/.Constants`EarthRepl)},{2,r<("rcore"/.Constants`EarthRepl)}},3](*1 => in the crust, 2 => in the core, 3 => not in the Earth*)


(* ::Input::Initialization:: *)
SortProcessesByRegion[associationlist_,key_]:=Module[{sortedlist},
sortedlist={{},{},{}}; (*1 corresponds to crust, 2 to core, 3 to outside earth (empty)*)
Do[AppendTo[sortedlist[[("region"/.associationlist[[i]])]],key/.associationlist[[i]]],{i,Length[associationlist]}];
sortedlist
]


(* ::Input::Initialization:: *)
Clear[Region\[CapitalTheta]]
Region\[CapitalTheta][r_,regionindex_]:=Piecewise[{{1,("rcore"/.Constants`EarthRepl)<r<=("rE"/.Constants`EarthRepl)},{1,r<("rcore"/.Constants`EarthRepl)},{1,r>("rE"/.Constants`EarthRepl)}}[[regionindex;;regionindex]],0]


(* ::Input::Initialization:: *)
ltor[l_,bE_]:=Sqrt[bE^2+l^2]


(* ::Subsubsection::Closed:: *)
(*Export to Directory*)


(* ::Input::Initialization:: *)
Clear[ExportToDir]
ExportToDir[obj_,name_,subsubdirectory_:"subsubdirectory",subdirectory_:DateString[{"Day","_","Month","_","Year"}],ndupsmax_:500]:=Module[{directory,rootname,extension,path,ndups},
directory=NotebookDirectory[]<>"Plots/"<>subdirectory<>"/"<>subsubdirectory;
{rootname,extension}=StringSplit[ToString@name,"."];

If[!DirectoryQ[directory],CreateDirectory[directory]];

ndups=0;
path=directory<>"/"<>rootname<>"_"<>ToString[#]<>"."<>extension&;
While[FileExistsQ[path[ndups]]&&ndups<=ndupsmax,
ndups++;
];

If[(ndups<ndupsmax),Export[path[ndups],obj];Return[<|"path"->path[ndups],"directory"->directory|>,Module],Print["Too many files with this name in the specified directory."]];

]


(* ::Input:: *)
(**)


(* ::Subsubsection::Closed:: *)
(*Plot Subscript[P, cap]*)


(* ::Input::Initialization:: *)
(*Clear[PlotPDist]
PlotPDist[PDict_]:=Module[{Pfunc,Pfs,\[Kappa],mpD,v0},
mpD=Union["mpD"/.PDict][[1]]("c")^2/("JpereV")/.Constants`SIConstRepl;
v0=Union["v0"/.PDict][[1]];
Print[mpD];
Print[v0];
Pfs=InterpolatePc[PDict];
(*Total*)
Do[
Pfunc="Pf"/.Pfs[[p]];
\[Kappa]="\[Kappa]"/.Pfs[[p]];
Print[Show[{DensityPlot[Pfunc[v\[Chi],Log10@bE],{v\[Chi],Pfunc["Domain"][[1,1]],Pfunc["Domain"][[1,2]]},{bE,10^Pfunc["Domain"][[2,1]],10^Pfunc["Domain"][[2,2]]},PlotLegends->Automatic,AxesLabel->Automatic,PlotRange->All,FrameLabel->{"Subscript[Log, 10]v\[Chi] [ms^-1]","Subscript[b, E] [m]"},PlotLabel->"Subscript[Log, 10]Subscript[P, cap] - Subscript[Log, 10]"<>ToString[StringForm["\[Kappa] = ``",Log10@\[Kappa]]]<>", Subscript[Log, 10]Subscript[m, Subscript[p, D]]"<>ToString[StringForm[" = `` [eV]",Log10@mpD]]<>ToString[StringForm[", v0 = `` [km/s]",v0/1000]]],ContourPlot[Pfunc[v\[Chi],Log10@bE],{v\[Chi],Pfunc["Domain"][[1,1]],Pfunc["Domain"][[1,2]]},{bE,10^Pfunc["Domain"][[2,1]],10^Pfunc["Domain"][[2,2]]},ContourShading->False,ContourLabels->True]}]];

(*Electronic*)
Pfunc="PfE"/.Pfs[[p]];
Print[Show[{DensityPlot[Pfunc[v\[Chi],Log10@bE],{v\[Chi],Pfunc["Domain"][[1,1]],Pfunc["Domain"][[1,2]]},{bE,10^Pfunc["Domain"][[2,1]],10^Pfunc["Domain"][[2,2]]},PlotLegends->Automatic,AxesLabel->Automatic,PlotRange->All,FrameLabel->{"Subscript[Log, 10]v\[Chi] [ms^-1]","Subscript[b, E] [m]"},PlotLabel->"Subscript[Log, 10]Subscript[P, cap,e] - Subscript[Log, 10]"<>ToString[StringForm["\[Kappa] = ``",Log10@\[Kappa]]]<>", Subscript[Log, 10]Subscript[m, Subscript[p, D]]"<>ToString[StringForm[" = `` [eV]",Log10@mpD]]<>ToString[StringForm[", v0 = `` [km/s]",v0/1000]]],ContourPlot[Pfunc[v\[Chi],Log10@bE],{v\[Chi],Pfunc["Domain"][[1,1]],Pfunc["Domain"][[1,2]]},{bE,10^Pfunc["Domain"][[2,1]],10^Pfunc["Domain"][[2,2]]},ContourShading->False,ContourLabels->True]}]];

(*Nuclear*)
Pfunc="PfNuc"/.Pfs[[p]];
Print[Show[{DensityPlot[Pfunc[v\[Chi],Log10@bE],{v\[Chi],Pfunc["Domain"][[1,1]],Pfunc["Domain"][[1,2]]},{bE,10^Pfunc["Domain"][[2,1]],10^Pfunc["Domain"][[2,2]]},PlotLegends->Automatic,AxesLabel->Automatic,PlotRange->All,FrameLabel->{"Subscript[Log, 10]v\[Chi] [ms^-1]","Subscript[b, E] [m]"},PlotLabel->"Subscript[Log, 10]Subscript[P, cap,Nuc] - Subscript[Log, 10]"<>ToString[StringForm["\[Kappa] = ``",Log10@\[Kappa]]]<>", Subscript[Log, 10]Subscript[m, Subscript[p, D]]"<>ToString[StringForm[" = `` [eV]",Log10@mpD]]<>ToString[StringForm[", v0 = `` [km/s]",v0/1000]]],ContourPlot[Pfunc[v\[Chi],Log10@bE],{v\[Chi],Pfunc["Domain"][[1,1]],Pfunc["Domain"][[1,2]]},{bE,10^Pfunc["Domain"][[2,1]],10^Pfunc["Domain"][[2,2]]},ContourShading->False,ContourLabels->True]}]];
,{p,Length[Pfs]}]
]*)


(* ::Input::Initialization:: *)
Clear[PlotPDist]
PlotPDist[PDict_]:=Module[{Pfunc,Pfs,\[Kappa],mpD,v0,plot},
mpD=Union["mpD"/.PDict][[1]] ("c")^2/("JpereV")/.Constants`SIConstRepl;
v0=Union["v0"/.PDict][[1]];

Pfs=InterpolatePc[PDict];
(*Total*)
Do[
Pfunc="Pf"/.Pfs[[p]];
\[Kappa]="\[Kappa]"/.Pfs[[p]];
plot=Show[{DensityPlot[Pfunc[v\[Chi],Log10@bE],{v\[Chi],Pfunc["Domain"][[1,1]],Pfunc["Domain"][[1,2]]},{bE,10^Pfunc["Domain"][[2,1]],10^Pfunc["Domain"][[2,2]]},PlotLegends->Automatic,AxesLabel->Automatic,PlotRange->All,FrameLabel->{"\!\(\*SubscriptBox[\(Log\), \(10\)]\)v\[Chi] [\!\(\*SuperscriptBox[\(ms\), \(-1\)]\)]","\!\(\*SubscriptBox[\(b\), \(E\)]\) [m]"},PlotLabel->"\!\(\*SubscriptBox[\(Log\), \(10\)]\)\!\(\*SubscriptBox[\(P\), \(cap\)]\) - \!\(\*SubscriptBox[\(Log\), \(10\)]\)"<>ToString[StringForm["\[Kappa] = ``",Log10@\[Kappa]]]<>", \!\(\*SubscriptBox[\(Log\), \(10\)]\)\!\(\*SubscriptBox[\(m\), SubscriptBox[\(p\), \(D\)]]\)"<>ToString[StringForm[" = `` [eV]",Log10@mpD]]<>ToString[StringForm[", v0 = `` [km/s]",v0/1000]]],ContourPlot[Pfunc[v\[Chi],Log10@bE],{v\[Chi],Pfunc["Domain"][[1,1]],Pfunc["Domain"][[1,2]]},{bE,10^Pfunc["Domain"][[2,1]],10^Pfunc["Domain"][[2,2]]},ContourShading->False,ContourLabels->True]}];
(*ExportToDir[plot,"Pcap_total.png","Pcap_total"];*)
ExportToDir[plot,ToString@StringForm["Pcap_total_Log10mpD_``_Log10\[Kappa]_``.png",Round@Log10@mpD,#]&@Round@Log10@\[Kappa],ToString@StringForm["Pcap_total_Log10mpD_``",Round@Log10@mpD]];

(*Print[StringQ@ToString@StringForm["Pcap_total_Log10mpD_``_Log10\[Kappa]_``.png",Round@Log10@mpD,Round@Log10@\[Kappa]]];
Print[ToString@StringForm["Pcap_total_Log10mpD_``_Log10\[Kappa]_``.png",Round@Log10@mpD,Round@Log10@\[Kappa]]];
Print[StringQ@ToString@StringForm["Pcap_total_Log10mpD_``",Round@Log10@mpD]];
Print[ToString@StringForm["Pcap_total_Log10mpD_``",Round@Log10@mpD]];*)

(*Electronic*)
Pfunc="PfE"/.Pfs[[p]];
plot=Show[{DensityPlot[Pfunc[v\[Chi],Log10@bE],{v\[Chi],Pfunc["Domain"][[1,1]],Pfunc["Domain"][[1,2]]},{bE,10^Pfunc["Domain"][[2,1]],10^Pfunc["Domain"][[2,2]]},PlotLegends->Automatic,AxesLabel->Automatic,PlotRange->All,FrameLabel->{"\!\(\*SubscriptBox[\(Log\), \(10\)]\)v\[Chi] [\!\(\*SuperscriptBox[\(ms\), \(-1\)]\)]","\!\(\*SubscriptBox[\(b\), \(E\)]\) [m]"},PlotLabel->"\!\(\*SubscriptBox[\(Log\), \(10\)]\)\!\(\*SubscriptBox[\(P\), \(cap, e\)]\) - \!\(\*SubscriptBox[\(Log\), \(10\)]\)"<>ToString[StringForm["\[Kappa] = ``",Log10@\[Kappa]]]<>", \!\(\*SubscriptBox[\(Log\), \(10\)]\)\!\(\*SubscriptBox[\(m\), SubscriptBox[\(p\), \(D\)]]\)"<>ToString[StringForm[" = `` [eV]",Log10@mpD]]<>ToString[StringForm[", v0 = `` [km/s]",v0/1000]]],ContourPlot[Pfunc[v\[Chi],Log10@bE],{v\[Chi],Pfunc["Domain"][[1,1]],Pfunc["Domain"][[1,2]]},{bE,10^Pfunc["Domain"][[2,1]],10^Pfunc["Domain"][[2,2]]},ContourShading->False,ContourLabels->True]}];
(*ExportToDir[plot,"Pcap_electronic.png","Pcap_electronic"];*)
ExportToDir[plot,StringForm["Pcap_electronic_Log10mpD_``_Log10\[Kappa]_``.png",Round[Log10@mpD],#]&@Round@Log10@\[Kappa],ToString@StringForm["Pcap_electronic_Log10mpD_``",Round[Log10@mpD]]];

(*Nuclear*)
Pfunc="PfNuc"/.Pfs[[p]];
plot=Show[{DensityPlot[Pfunc[v\[Chi],Log10@bE],{v\[Chi],Pfunc["Domain"][[1,1]],Pfunc["Domain"][[1,2]]},{bE,10^Pfunc["Domain"][[2,1]],10^Pfunc["Domain"][[2,2]]},PlotLegends->Automatic,AxesLabel->Automatic,PlotRange->All,FrameLabel->{"\!\(\*SubscriptBox[\(Log\), \(10\)]\)v\[Chi] [\!\(\*SuperscriptBox[\(ms\), \(-1\)]\)]","\!\(\*SubscriptBox[\(b\), \(E\)]\) [m]"},PlotLabel->"\!\(\*SubscriptBox[\(Log\), \(10\)]\)\!\(\*SubscriptBox[\(P\), \(cap, Nuc\)]\) - \!\(\*SubscriptBox[\(Log\), \(10\)]\)"<>ToString[StringForm["\[Kappa] = ``",Log10@\[Kappa]]]<>", \!\(\*SubscriptBox[\(Log\), \(10\)]\)\!\(\*SubscriptBox[\(m\), SubscriptBox[\(p\), \(D\)]]\)"<>ToString[StringForm[" = `` [eV]",Log10@mpD]]<>ToString[StringForm[", v0 = `` [km/s]",v0/1000]]],ContourPlot[Pfunc[v\[Chi],Log10@bE],{v\[Chi],Pfunc["Domain"][[1,1]],Pfunc["Domain"][[1,2]]},{bE,10^Pfunc["Domain"][[2,1]],10^Pfunc["Domain"][[2,2]]},ContourShading->False,ContourLabels->True]}];
(*ExportToDir[plot,"Pcap_nuclear.png","Pcap_nuclear"];*)
ExportToDir[plot,ToString@StringForm["Pcap_nuclear_Log10mpD_``_Log10\[Kappa]_``.png",Round@Log10@mpD,#]&@Round@Log10@\[Kappa],ToString@StringForm["Pcap_nuclear_Log10mpD_``",Round@Log10@mpD]];


ClearAll[Pfunc,plot];
,{p,Length[Pfs]}];
Clear[Pfs];
]


(* ::Input::Initialization:: *)
Clear[PlotPDistwPf]
PlotPDistwPf[PDict_]:=Module[{Pfunc,Pfs,\[Kappa],mpD,v0,plot},
mpD=Union["mpD"/.PDict][[1]] ("c")^2/("JpereV")/.Constants`SIConstRepl;
v0=Union["v0"/.PDict][[1]];

(*Pfs=InterpolatePc[PDict];*)
Pfs = PDict;
(*Total*)
Do[
Pfunc="Pf"/.Pfs[[p]];
\[Kappa]="\[Kappa]"/.Pfs[[p]];
plot=Show[{DensityPlot[Pfunc[v\[Chi],Log10@bE],{v\[Chi],Pfunc["Domain"][[1,1]],Pfunc["Domain"][[1,2]]},{bE,10^Pfunc["Domain"][[2,1]],10^Pfunc["Domain"][[2,2]]},PlotLegends->Automatic,AxesLabel->Automatic,PlotRange->All,FrameLabel->{"\!\(\*SubscriptBox[\(Log\), \(10\)]\)v\[Chi] [\!\(\*SuperscriptBox[\(ms\), \(-1\)]\)]","\!\(\*SubscriptBox[\(b\), \(E\)]\) [m]"},PlotLabel->"\!\(\*SubscriptBox[\(Log\), \(10\)]\)\!\(\*SubscriptBox[\(P\), \(cap\)]\) - \!\(\*SubscriptBox[\(Log\), \(10\)]\)"<>ToString[StringForm["\[Kappa] = ``",Log10@\[Kappa]]]<>", \!\(\*SubscriptBox[\(Log\), \(10\)]\)\!\(\*SubscriptBox[\(m\), SubscriptBox[\(p\), \(D\)]]\)"<>ToString[StringForm[" = `` [eV]",Log10@mpD]]<>ToString[StringForm[", v0 = `` [km/s]",v0/1000]]],ContourPlot[Pfunc[v\[Chi],Log10@bE],{v\[Chi],Pfunc["Domain"][[1,1]],Pfunc["Domain"][[1,2]]},{bE,10^Pfunc["Domain"][[2,1]],10^Pfunc["Domain"][[2,2]]},ContourShading->False,ContourLabels->True]}];
(*ExportToDir[plot,"Pcap_total.png","Pcap_total"];*)
ExportToDir[plot,ToString@StringForm["Pcap_total_Log10mpD_``_Log10\[Kappa]_``.png",Round@Log10@mpD,#]&@Round@Log10@\[Kappa],ToString@StringForm["Pcap_total_Log10mpD_``",Round@Log10@mpD]];

(*Print[StringQ@ToString@StringForm["Pcap_total_Log10mpD_``_Log10\[Kappa]_``.png",Round@Log10@mpD,Round@Log10@\[Kappa]]];
Print[ToString@StringForm["Pcap_total_Log10mpD_``_Log10\[Kappa]_``.png",Round@Log10@mpD,Round@Log10@\[Kappa]]];
Print[StringQ@ToString@StringForm["Pcap_total_Log10mpD_``",Round@Log10@mpD]];
Print[ToString@StringForm["Pcap_total_Log10mpD_``",Round@Log10@mpD]];*)

(*Electronic*)
Pfunc="PfE"/.Pfs[[p]];
plot=Show[{DensityPlot[Pfunc[v\[Chi],Log10@bE],{v\[Chi],Pfunc["Domain"][[1,1]],Pfunc["Domain"][[1,2]]},{bE,10^Pfunc["Domain"][[2,1]],10^Pfunc["Domain"][[2,2]]},PlotLegends->Automatic,AxesLabel->Automatic,PlotRange->All,FrameLabel->{"\!\(\*SubscriptBox[\(Log\), \(10\)]\)v\[Chi] [\!\(\*SuperscriptBox[\(ms\), \(-1\)]\)]","\!\(\*SubscriptBox[\(b\), \(E\)]\) [m]"},PlotLabel->"\!\(\*SubscriptBox[\(Log\), \(10\)]\)\!\(\*SubscriptBox[\(P\), \(cap, e\)]\) - \!\(\*SubscriptBox[\(Log\), \(10\)]\)"<>ToString[StringForm["\[Kappa] = ``",Log10@\[Kappa]]]<>", \!\(\*SubscriptBox[\(Log\), \(10\)]\)\!\(\*SubscriptBox[\(m\), SubscriptBox[\(p\), \(D\)]]\)"<>ToString[StringForm[" = `` [eV]",Log10@mpD]]<>ToString[StringForm[", v0 = `` [km/s]",v0/1000]]],ContourPlot[Pfunc[v\[Chi],Log10@bE],{v\[Chi],Pfunc["Domain"][[1,1]],Pfunc["Domain"][[1,2]]},{bE,10^Pfunc["Domain"][[2,1]],10^Pfunc["Domain"][[2,2]]},ContourShading->False,ContourLabels->True]}];
(*ExportToDir[plot,"Pcap_electronic.png","Pcap_electronic"];*)
ExportToDir[plot,StringForm["Pcap_electronic_Log10mpD_``_Log10\[Kappa]_``.png",Round[Log10@mpD],#]&@Round@Log10@\[Kappa],ToString@StringForm["Pcap_electronic_Log10mpD_``",Round[Log10@mpD]]];

(*Nuclear*)
Pfunc="PfNuc"/.Pfs[[p]];
plot=Show[{DensityPlot[Pfunc[v\[Chi],Log10@bE],{v\[Chi],Pfunc["Domain"][[1,1]],Pfunc["Domain"][[1,2]]},{bE,10^Pfunc["Domain"][[2,1]],10^Pfunc["Domain"][[2,2]]},PlotLegends->Automatic,AxesLabel->Automatic,PlotRange->All,FrameLabel->{"\!\(\*SubscriptBox[\(Log\), \(10\)]\)v\[Chi] [\!\(\*SuperscriptBox[\(ms\), \(-1\)]\)]","\!\(\*SubscriptBox[\(b\), \(E\)]\) [m]"},PlotLabel->"\!\(\*SubscriptBox[\(Log\), \(10\)]\)\!\(\*SubscriptBox[\(P\), \(cap, Nuc\)]\) - \!\(\*SubscriptBox[\(Log\), \(10\)]\)"<>ToString[StringForm["\[Kappa] = ``",Log10@\[Kappa]]]<>", \!\(\*SubscriptBox[\(Log\), \(10\)]\)\!\(\*SubscriptBox[\(m\), SubscriptBox[\(p\), \(D\)]]\)"<>ToString[StringForm[" = `` [eV]",Log10@mpD]]<>ToString[StringForm[", v0 = `` [km/s]",v0/1000]]],ContourPlot[Pfunc[v\[Chi],Log10@bE],{v\[Chi],Pfunc["Domain"][[1,1]],Pfunc["Domain"][[1,2]]},{bE,10^Pfunc["Domain"][[2,1]],10^Pfunc["Domain"][[2,2]]},ContourShading->False,ContourLabels->True]}];
(*ExportToDir[plot,"Pcap_nuclear.png","Pcap_nuclear"];*)
ExportToDir[plot,ToString@StringForm["Pcap_nuclear_Log10mpD_``_Log10\[Kappa]_``.png",Round@Log10@mpD,#]&@Round@Log10@\[Kappa],ToString@StringForm["Pcap_nuclear_Log10mpD_``",Round@Log10@mpD]];


ClearAll[Pfunc,plot];
,{p,Length[Pfs]}];
Clear[Pfs];
]


(* ::Input::Initialization:: *)
(*Clear[PlotPDistEvsNuc]
PlotPDistEvsNuc[PDict_]:=Module[{Pfunc,\[Kappa],mpD,v0},
mpD=Union["mpD"/.PDict][[1]]("c")^2/("JpereV")/.Constants`SIConstRepl;
v0=Union["v0"/.PDict][[1]];
(*Print[mpD];
Print[v0];
PfEs=InterpolatePc[PDict];*)
Do[
Pfunc="PfE"/.PDict[[p]];
\[Kappa]="\[Kappa]"/.PDict[[p]];
Print[Show[{DensityPlot[Pfunc[v\[Chi],bE],{v\[Chi],Pfunc["Domain"][[1,1]],Pfunc["Domain"][[1,2]]},{bE,Pfunc["Domain"][[2,1]],Pfunc["Domain"][[2,2]]},PlotLegends->Automatic,AxesLabel->Automatic,PlotRange->All,FrameLabel->{"Subscript[Log, 10]v\[Chi]","Subscript[Log, 10]Subscript[b, E]"},PlotLabel->"Subscript[P, cap,e] - Subscript[Log, 10]"<>ToString[StringForm["\[Kappa] = ``",Log10@\[Kappa]]]<>", Subscript[Log, 10]Subscript[m, Subscript[p, D]]"<>ToString[StringForm[" = `` [eV]",Log10@mpD]]<>ToString[StringForm[", v0 = `` [m/s]",v0]]],ContourPlot[Pfunc[v\[Chi],bE],{v\[Chi],Pfunc["Domain"][[1,1]],Pfunc["Domain"][[1,2]]},{bE,Pfunc["Domain"][[2,1]],Pfunc["Domain"][[2,2]]},ContourShading->False,ContourLabels->True]}]];
,{p,Length[PDict]}];

Do[
Pfunc="PfNuc"/.PDict[[p]];
\[Kappa]="\[Kappa]"/.PDict[[p]];
Print[Show[{DensityPlot[Pfunc[v\[Chi],bE],{v\[Chi],Pfunc["Domain"][[1,1]],Pfunc["Domain"][[1,2]]},{bE,Pfunc["Domain"][[2,1]],Pfunc["Domain"][[2,2]]},PlotLegends->Automatic,AxesLabel->Automatic,PlotRange->All,FrameLabel->{"Subscript[Log, 10]v\[Chi]","Subscript[Log, 10]Subscript[b, E]"},PlotLabel->"Subscript[P, cap,Nuc] - Subscript[Log, 10]"<>ToString[StringForm["\[Kappa] = ``",Log10@\[Kappa]]]<>", Subscript[Log, 10]Subscript[m, Subscript[p, D]]"<>ToString[StringForm[" = `` [eV]",Log10@mpD]]<>ToString[StringForm[", v0 = `` [m/s]",v0]]],ContourPlot[Pfunc[v\[Chi],bE],{v\[Chi],Pfunc["Domain"][[1,1]],Pfunc["Domain"][[1,2]]},{bE,Pfunc["Domain"][[2,1]],Pfunc["Domain"][[2,2]]},ContourShading->False,ContourLabels->True]}]];
,{p,Length[PDict]}]

]*)


(* ::Subsubsection::Closed:: *)
(*Plot KE change*)


(* ::Input::Initialization:: *)
Clear[PlotKEchange]
PlotKEchange[KEchangeMassDicts_]:=Module[{AllbutMassFixed,\[Kappa],v0,fD,meDovermpD,plotname},
Do[
AllbutMassFixed=Table[KEchangeMassDicts[[i,j]],{i,Length[KEchangeMassDicts]}];
\[Kappa]=Union["\[Kappa]"/.AllbutMassFixed][[1]];

If[\[Kappa]<10^-7,
v0=Union["v0"/.AllbutMassFixed][[1]];
fD=Union["fD"/.AllbutMassFixed][[1]];
meDovermpD=Union["meD"/.AllbutMassFixed][[1]]/Union["mpD"/.AllbutMassFixed][[1]];

(*plotname=("Subscript[\[CapitalDelta]v, e]/Subscript[Overscript[v, _], e]Sqrt[(Subscript[\[Alpha], D]Subscript[f, D])/(\[Alpha]10^-2)] - Subscript[Log, 10]"<>ToString[StringForm["\[Kappa] = ``",Log10@\[Kappa]]]<>"\n Subscript[v, 0]"<>ToString[StringForm["= ``",N@(v0/1000)]]<>"[kms^-1], Subscript[f, D] = "<>ToString[Log10@fD]<>", Subscript[Log, 10]FractionBox[SubscriptBox[m,SubscriptBox[e,D]],Subscript[m, Subscript[p, D]]\]="<>ToString[Log10@meDovermpD]);*)
plotname=("\!\(\*FractionBox[SubscriptBox[\(\[CapitalDelta]v\), \(e\)], SubscriptBox[OverscriptBox[\(v\), \(_\)], \(e\)]]\)\!\(\*SqrtBox[FractionBox[\(\[Alpha]\\\ 0.05\), \(\*SubscriptBox[\(\[Alpha]\), \(D\)] \*SubscriptBox[\(f\), \(D\)]\)]]\) - \!\(\*SubscriptBox[\(Log\), \(10\)]\)"<>ToString[StringForm["\[Kappa] = ``",Log10@\[Kappa]]]<>", \!\(\*SubscriptBox[\(v\), \(0\)]\)"<>ToString[StringForm["= ``",N@(v0/1000)]]<>"[\!\(\*SuperscriptBox[\(kms\), \(-1\)]\)]"<>", \!\(\*SubscriptBox[\(Log\), \(10\)]\)\!\(\*FractionBox[SubscriptBox[\(m\), SubscriptBox[\(e\), \(D\)]], SubscriptBox[\(m\), SubscriptBox[\(p\), \(D\)]]]\)="<>ToString[Log10@meDovermpD]);

Print[ListLinePlot[{Log10@("mpD" ("c")^2/("JpereV")/.Constants`SIConstRepl),Log10@"\[CapitalDelta]eKE/vpbar"}/.AllbutMassFixed,AxesLabel->{"\!\(\*SubscriptBox[\(Log\), \(10\)]\)mpD [eV]"},PlotLabel->plotname]];

];
,{j,Dimensions[KEchangeMassDicts][[2]]}];
]


(* ::Input::Initialization:: *)
Clear[PlotKEchange\[Kappa]Contour]
PlotKEchange\[Kappa]Contour[KEchangeMassDicts_]:=Module[{AllbutMassand\[Kappa]Fixed,v0,meDovermpD,plotname},
AllbutMassand\[Kappa]Fixed=Gather[Flatten[KEchangeMassDicts],(#1["v0"]==#2["v0"]&&#1["meD"]/#1["mpD"]==#2["meD"]/#2["mpD"])&];
Do[
(*AllbutMassFixed=[[i]];*)

(*\[Kappa]=Union["\[Kappa]"/.AllbutMassFixed][[1]];*)

v0=Union["v0"/.AllbutMassand\[Kappa]Fixed[[j]]][[1]];
(*fD=Union["fD"/.AllbutMassFixed][[1]];*)
meDovermpD=Union[("meD")/("mpD")/.AllbutMassand\[Kappa]Fixed[[j]]][[1]];
(*Print[meDovermpD];
Print[v0];*)
(*plotname=("Subscript[\[CapitalDelta]v, e]/Subscript[Overscript[v, _], e]Sqrt[(Subscript[\[Alpha], D]Subscript[f, D])/(\[Alpha]10^-2)] - Subscript[Log, 10]"<>ToString[StringForm["\[Kappa] = ``",Log10@\[Kappa]]]<>"\n Subscript[v, 0]"<>ToString[StringForm["= ``",N@(v0/1000)]]<>"[kms^-1], Subscript[f, D] = "<>ToString[Log10@fD]<>", Subscript[Log, 10]FractionBox[SubscriptBox[m,SubscriptBox[e,D]],Subscript[m, Subscript[p, D]]\]="<>ToString[Log10@meDovermpD]);*)
plotname=("\!\(\*SubscriptBox[\(Log\), \(10\)]\)\!\(\*FractionBox[SubscriptBox[\(\[CapitalDelta]v\), \(e\)], SubscriptBox[OverscriptBox[\(v\), \(_\)], \(e\)]]\)\!\(\*SqrtBox[FractionBox[\(\[Alpha]\\\ 0.05\), \(\*SubscriptBox[\(\[Alpha]\), \(D\)] \*SubscriptBox[\(f\), \(D\)]\)]]\) - "<>" \!\(\*SubscriptBox[\(v\), \(0\)]\)"<>ToString[StringForm["= ``",N@(v0/1000)]]<>"[\!\(\*SuperscriptBox[\(kms\), \(-1\)]\)]"<>",  \!\(\*SubscriptBox[\(Log\), \(10\)]\)\!\(\*FractionBox[SubscriptBox[\(m\), SubscriptBox[\(e\), \(D\)]], SubscriptBox[\(m\), SubscriptBox[\(p\), \(D\)]]]\)="<>ToString[Log10@meDovermpD]);

(*Print[{Log10@("mpD"("c")^2/("JpereV")/.Constants`SIConstRepl),Log10@"\[Kappa]",Log10@"\[CapitalDelta]eKE/vpbar"}/.AllbutMassand\[Kappa]Fixed[[j]]];*)

(*Print[ListContourPlot[{Log10@("mpD"("c")^2/("JpereV")/.Constants`SIConstRepl),Log10@"\[Kappa]",Log10@"\[CapitalDelta]eKE/vpbar"}/.AllbutMassand\[Kappa]Fixed[[j]],ContourShading->False,ContourLabels->True,FrameLabel->{Style["Subscript[Log, 10]mpD [eV]",Black,14],Style["\[Kappa] [ ]",Black,14]},PlotLabel->Style[plotname,Black,12],FrameStyle->Directive[Black,Thickness[0.003]],Contours->Join[Table[i,{i,-2,10,1}],Table[i,{i,-5,-20,-5}]]]];*)
(*Print[ListContourPlot[{Log10@("mpD"("c")^2/("JpereV")/.Constants`SIConstRepl),Log10@"\[Kappa]",Log10@"\[CapitalDelta]eKE/vpbar"}/.AllbutMassand\[Kappa]Fixed[[j]],ContourShading->False,ContourLabels->True,FrameLabel->{Style["Subscript[Log, 10]mpD [eV]",Black,14],Style["\[Kappa] [ ]",Black,14]},PlotLabel->Style[plotname,Black,12],FrameStyle->Directive[Black,Thickness[0.003]],Contours->{0},ContourStyle->Red]];*)

Print[Show[{ListContourPlot[{Log10@("mpD" ("c")^2/("JpereV")/.Constants`SIConstRepl),Log10@"\[Kappa]",Log10@"\[CapitalDelta]eKE/vpbar"}/.AllbutMassand\[Kappa]Fixed[[j]],ContourShading->False,ContourLabels->True,FrameLabel->{Style["\!\(\*SubscriptBox[\(Log\), \(10\)]\)mpD [eV]",Black,14],Style["\[Kappa] [ ]",Black,14]},PlotLabel->Style[plotname,Black,12],FrameStyle->Directive[Black,Thickness[0.003]],Contours->Join[Table[i,{i,-10,10,1}],Table[i,{i,-10,-20,-5}]]],ListContourPlot[{Log10@("mpD" ("c")^2/("JpereV")/.Constants`SIConstRepl),Log10@"\[Kappa]",Log10@"\[CapitalDelta]eKE/vpbar"}/.AllbutMassand\[Kappa]Fixed[[j]],ContourShading->False,ContourLabels->True,FrameLabel->{Style["\!\(\*SubscriptBox[\(Log\), \(10\)]\)mpD [eV]",Black,14],Style["\[Kappa] [ ]",Black,14]},PlotLabel->Style[plotname,Black,12],FrameStyle->Directive[Black,Thickness[0.003]],Contours->{0},ContourStyle->Directive[Red,Thick]]}]];

Print[Show[{ListDensityPlot[{Log10@("mpD" ("c")^2/("JpereV")/.Constants`SIConstRepl),Log10@"\[Kappa]",Log10@"\[CapitalDelta]eKE/vpbar"}/.AllbutMassand\[Kappa]Fixed[[j]],FrameLabel->{"\!\(\*SubscriptBox[\(Log\), \(10\)]\)mpD [eV]","\[Kappa] [ ]"},PlotLabel->plotname,PlotLegends->Automatic],ListContourPlot[{Log10@("mpD" ("c")^2/("JpereV")/.Constants`SIConstRepl),Log10@"\[Kappa]",Log10@"\[CapitalDelta]eKE/vpbar"}/.AllbutMassand\[Kappa]Fixed[[j]],ContourShading->False,ContourLabels->True,Contours->Join[Table[i,{i,-2,10,1}],Table[i,{i,-5,-20,-5}]]],ListContourPlot[{Log10@("mpD" ("c")^2/("JpereV")/.Constants`SIConstRepl),Log10@"\[Kappa]",Log10@"\[CapitalDelta]eKE/vpbar"}/.AllbutMassand\[Kappa]Fixed[[j]],ContourShading->False,ContourLabels->False,Contours->{0},ContourStyle->Directive[Red,Thick]]}]];



plotname=("\!\(\*SubscriptBox[\(Log\), \(10\)]\)\!\(\*FractionBox[SubscriptBox[\(\[CapitalDelta]v\), \(p\)], SubscriptBox[OverscriptBox[\(v\), \(_\)], \(p\)]]\)\!\(\*SqrtBox[FractionBox[\(\[Alpha]\\\ 0.05\), \(\*SubscriptBox[\(\[Alpha]\), \(D\)] \*SubscriptBox[\(f\), \(D\)]\)]]\) - "<>" \!\(\*SubscriptBox[\(v\), \(0\)]\)"<>ToString[StringForm["= ``",N@(v0/1000)]]<>"[\!\(\*SuperscriptBox[\(kms\), \(-1\)]\)]"<>",  \!\(\*SubscriptBox[\(Log\), \(10\)]\)\!\(\*FractionBox[SubscriptBox[\(m\), SubscriptBox[\(e\), \(D\)]], SubscriptBox[\(m\), SubscriptBox[\(p\), \(D\)]]]\)="<>ToString[Log10@meDovermpD]);
Print[Show[{ListDensityPlot[{Log10@("mpD" ("c")^2/("JpereV")/.Constants`SIConstRepl),Log10@"\[Kappa]",Log10@"\[CapitalDelta]pKE/vpbar"}/.AllbutMassand\[Kappa]Fixed[[j]],FrameLabel->{"\!\(\*SubscriptBox[\(Log\), \(10\)]\)mpD [eV]","\[Kappa] [ ]"},PlotLabel->plotname,PlotLegends->Automatic],ListContourPlot[{Log10@("mpD" ("c")^2/("JpereV")/.Constants`SIConstRepl),Log10@"\[Kappa]",Log10@"\[CapitalDelta]pKE/vpbar"}/.AllbutMassand\[Kappa]Fixed[[j]],ContourShading->False,ContourLabels->True,Contours->Join[Table[i,{i,-2,10,1}],Table[i,{i,-5,-20,-5}]]],ListContourPlot[{Log10@("mpD" ("c")^2/("JpereV")/.Constants`SIConstRepl),Log10@"\[Kappa]",Log10@"\[CapitalDelta]pKE/vpbar"}/.AllbutMassand\[Kappa]Fixed[[j]],ContourShading->False,ContourLabels->False,Contours->{0},ContourStyle->Directive[Red,Thick]]}]]

,{j,Length@AllbutMassand\[Kappa]Fixed}];
(*,{j,1}];*)
]


(* ::Subsubsection::Closed:: *)
(*Plot Nc*)


(* ::Input::Initialization:: *)
Clear[PlotNc\[Kappa]Contour]
PlotNc\[Kappa]Contour[KEchangeMassDicts_]:=Module[{AllbutMassand\[Kappa]Fixed,v0,meDovermpD,plotname},
AllbutMassand\[Kappa]Fixed=Gather[Flatten[KEchangeMassDicts],(#1["v0"]==#2["v0"]&&#1["meD"]/#1["mpD"]==#2["meD"]/#2["mpD"])&];
Do[


v0=Union["v0"/.AllbutMassand\[Kappa]Fixed[[j]]][[1]];

meDovermpD=Union[("meD")/("mpD")/.AllbutMassand\[Kappa]Fixed[[j]]][[1]];

plotname=("\!\(\*SubscriptBox[\(Log\), \(\(10\)\(\\\ \)\)]\)\!\(\*SubscriptBox[\(N\), \(c\)]\)\!\(\*FractionBox[\(0.05\), SubscriptBox[\(f\), \(D\)]]\) [ ]- "<>" \!\(\*SubscriptBox[\(v\), \(0\)]\)"<>ToString[StringForm["= ``",N@(v0/1000)]]<>"[\!\(\*SuperscriptBox[\(kms\), \(-1\)]\)]"<>",  \!\(\*SubscriptBox[\(Log\), \(10\)]\)\!\(\*FractionBox[SubscriptBox[\(m\), SubscriptBox[\(e\), \(D\)]], SubscriptBox[\(m\), SubscriptBox[\(p\), \(D\)]]]\)="<>ToString[Log10@meDovermpD]);

Do[With[{k=i,l=j},AssociateTo[AllbutMassand\[Kappa]Fixed[[l,k]],"Nc"->(AllbutMassand\[Kappa]Fixed[[l,k]]["Nc"]/.AllbutMassand\[Kappa]Fixed[[l,k]])]],{i,Length[AllbutMassand\[Kappa]Fixed[[j]]]}];

Print[Show[{ListDensityPlot[{Log10@("mpD" ("c")^2/("JpereV")/.Constants`SIConstRepl),Log10@"\[Kappa]",Log10@"Nc"}/.AllbutMassand\[Kappa]Fixed[[j]],FrameLabel->{"\!\(\*SubscriptBox[\(Log\), \(10\)]\)mpD [eV]","\[Kappa] [ ]"},PlotLabel->plotname,PlotLegends->Automatic],ListContourPlot[{Log10@("mpD" ("c")^2/("JpereV")/.Constants`SIConstRepl),Log10@"\[Kappa]",Log10@"Nc"}/.AllbutMassand\[Kappa]Fixed[[j]],ContourShading->False,ContourLabels->True,Contours->Join[Table[i,{i,10,40,2}],Table[i,{i,15,-30,-5}]]]}]];

,{j,Length@AllbutMassand\[Kappa]Fixed}];
]


(* ::Subsubsection::Closed:: *)
(*Plot Capture Rate*)


(* ::Input::Initialization:: *)
Clear[PlotCaptureRate\[Kappa]Contour]
PlotCaptureRate\[Kappa]Contour[KEchangeMassDicts_]:=Module[{AllbutMassand\[Kappa]Fixed,v0,meDovermpD,plotname,VEarth,list},
AllbutMassand\[Kappa]Fixed=Gather[Flatten[KEchangeMassDicts],(#1["v0"]==#2["v0"]&&#1["meD"]/#1["mpD"]==#2["meD"]/#2["mpD"])&];
Do[
(*AllbutMassFixed=[[i]];*)

(*\[Kappa]=Union["\[Kappa]"/.AllbutMassFixed][[1]];*)

v0=Union["v0"/.AllbutMassand\[Kappa]Fixed[[j]]][[1]];
(*fD=Union["fD"/.AllbutMassFixed][[1]];*)
meDovermpD=Union["meD"/"mpD"/.AllbutMassand\[Kappa]Fixed[[j]]][[1]];

plotname=("\!\(\*SubscriptBox[\(Log\), \(10\)]\)\!\(\*SubscriptBox[\(\[CapitalGamma]\), \(cap\)]\)\!\(\*SubscriptBox[\(n\), \(aDM\)]\)\!\(\*FractionBox[\(0.05\), SubscriptBox[\(f\), \(D\)]]\) [\!\(\*SuperscriptBox[\(m\), \(-3\)]\)\!\(\*SuperscriptBox[\(s\), \(-1\)]\)] - "<>" \!\(\*SubscriptBox[\(v\), \(0\)]\)"<>ToString[StringForm["= ``",N@(v0/1000)]]<>"[\!\(\*SuperscriptBox[\(kms\), \(-1\)]\)]"<>",  \!\(\*SubscriptBox[\(Log\), \(10\)]\)\!\(\*FractionBox[SubscriptBox[\(m\), SubscriptBox[\(e\), \(D\)]], SubscriptBox[\(m\), SubscriptBox[\(p\), \(D\)]]]\)="<>ToString[Log10@meDovermpD]);

VEarth=((4\[Pi])/3)("rE")^3/.Constants`EarthRepl;
(*Print[Log10@VEarth];*)

Do[With[{k=i,l=j},AssociateTo[AllbutMassand\[Kappa]Fixed[[l,k]],"dNcMaxdt"->(AllbutMassand\[Kappa]Fixed[[l,k]]["dNcMaxdt"]/.AllbutMassand\[Kappa]Fixed[[l,k]])]],{i,Length[AllbutMassand\[Kappa]Fixed[[j]]]}];

(*Print[Show[{ListDensityPlot[{Log10@("mpD"(("c")^2/"JpereV")/.Constants`SIConstRepl),Log10@"\[Kappa]",Log10@"dNcMaxdt"}/.AllbutMassand\[Kappa]Fixed[[j]],FrameLabel->{"Subscript[Log, 10]mpD [eV]","\[Kappa] [ ]"},PlotLabel->plotname,PlotLegends->Automatic],ListContourPlot[{Log10@("mpD"(("c")^2/"JpereV")/.Constants`SIConstRepl),Log10@"\[Kappa]",Log10@"dNcMaxdt"}/.AllbutMassand\[Kappa]Fixed[[j]],ContourShading->False,ContourLabels->True,Contours->15]}]];*)



Print[Show[{ListDensityPlot[{Log10@("mpD"(("c")^2/"JpereV")/.Constants`SIConstRepl),Log10@"\[Kappa]",Log10@(("dNcMaxdt")/VEarth)}/.AllbutMassand\[Kappa]Fixed[[j]],FrameLabel->{"\!\(\*SubscriptBox[\(Log\), \(10\)]\)mpD [eV]","\[Kappa] [ ]"},PlotLabel->plotname,PlotLegends->Automatic],ListContourPlot[{Log10@("mpD"(("c")^2/"JpereV")/.Constants`SIConstRepl),Log10@"\[Kappa]",Log10@(("dNcMaxdt")/VEarth)}/.AllbutMassand\[Kappa]Fixed[[j]],ContourShading->False,ContourLabels->True,(*Contours->Join[Table[i,{i,-21,-1,2}],Table[i,{i,-25,-71,-5}]]*)Contours->Join[Table[i,{i,-15,5,2}],Table[i,{i,-15,-65,-5}]]]}]];

(*Print[ListDensityPlot[{Log10@("mpD"(("c")^2/"JpereV")/.Constants`SIConstRepl),Log10@"\[Kappa]",Log10["dNcMaxdt"/VEarth]}/.AllbutMassand\[Kappa]Fixed[[j]],FrameLabel->{"Subscript[Log, 10]mpD [eV]","\[Kappa] [ ]"},PlotLabel->plotname,PlotLegends->Automatic]];*)

(*list ={Log10@("mpD"(("c")^2/"JpereV")/.Constants`SIConstRepl),Log10@"\[Kappa]",Log10["dNcMaxdt"/VEarth]}/.AllbutMassand\[Kappa]Fixed[[j]]/.AllbutMassand\[Kappa]Fixed[[j]];
Print[Show[{ListDensityPlot[list,FrameLabel->{"Subscript[Log, 10]mpD [eV]","\[Kappa] [ ]"},PlotLabel->plotname,PlotLegends->Automatic],ListContourPlot[list,ContourLabels->All,ContourShading->None,Contours->5,ContourStyle->Black]}]];*)

(*Print[ListContourPlot[{Log10@("mpD"(("c")^2/"JpereV")/.Constants`SIConstRepl),Log10@"\[Kappa]",Log10["dNcMaxdt"/VEarth]}/.AllbutMassand\[Kappa]Fixed[[j]]/.AllbutMassand\[Kappa]Fixed[[j]],ContourShading->False,ContourLabels->True]];*)
(*plotname=("Subscript[Log, 10]Max(Subscript[P, cap,Nuc])[ ] - "<>" Subscript[v, 0]"<>ToString[StringForm["= ``",N@(v0/1000)]]<>"[kms^-1]");
Print[Show[{ListDensityPlot[Table[{Log10@("mpD"(("c")^2/"JpereV")/.Constants`SIConstRepl),Log10@"\[Kappa]",Max[("PfNuc"/.AllbutMassand\[Kappa]Fixed[[j,i]])["ValuesOnGrid"]]}/.AllbutMassand\[Kappa]Fixed[[j,i]],{i,Length[AllbutMassand\[Kappa]Fixed[[j]]]}],FrameLabel->{"Subscript[Log, 10]mpD [eV]","\[Kappa] [ ]"},PlotLabel->plotname,PlotLegends->Automatic],ListContourPlot[Table[{Log10@("mpD"(("c")^2/"JpereV")/.Constants`SIConstRepl),Log10@"\[Kappa]",Max[("PfNuc"/.AllbutMassand\[Kappa]Fixed[[j,i]])["ValuesOnGrid"]]}/.AllbutMassand\[Kappa]Fixed[[j,i]],{i,Length[AllbutMassand\[Kappa]Fixed[[j]]]}],ContourShading->False,ContourLabels->True,Contours->Join[Table[i,{i,-20,-0,2}],Table[i,{i,-25,-71,-5}]]]}]];*)

(*plotname=("Subscript[Log, 10]Max(Subscript[P, cap,e])[ ] - "<>" Subscript[v, 0]"<>ToString[StringForm["= ``",N@(v0/1000)]]<>"[kms^-1]");
Print[Show[{ListDensityPlot[Table[{Log10@("mpD"(("c")^2/"JpereV")/.Constants`SIConstRepl),Log10@"\[Kappa]",Max[("PfE"/.AllbutMassand\[Kappa]Fixed[[j,i]])["ValuesOnGrid"]]}/.AllbutMassand\[Kappa]Fixed[[j,i]],{i,Length[AllbutMassand\[Kappa]Fixed[[j]]]}],FrameLabel->{"Subscript[Log, 10]mpD [eV]","\[Kappa] [ ]"},PlotLabel->plotname,PlotLegends->Automatic],ListContourPlot[Table[{Log10@("mpD"(("c")^2/"JpereV")/.Constants`SIConstRepl),Log10@"\[Kappa]",Max[("PfE"/.AllbutMassand\[Kappa]Fixed[[j,i]])["ValuesOnGrid"]]}/.AllbutMassand\[Kappa]Fixed[[j,i]],{i,Length[AllbutMassand\[Kappa]Fixed[[j]]]}],ContourShading->False,ContourLabels->True,Contours->Join[Table[i,{i,-20,0,2}],Table[i,{i,-25,-71,-5}]]]}]];*)

(*plotname=("Subscript[Log, 10]Max(Subscript[P, cap])[ ] - "<>" Subscript[v, 0]"<>ToString[StringForm["= ``",N@(v0/1000)]]<>"[kms^-1]");
Print[Show[{ListDensityPlot[Table[{Log10@("mpD"(("c")^2/"JpereV")/.Constants`SIConstRepl),Log10@"\[Kappa]",Max[("Pf"/.AllbutMassand\[Kappa]Fixed[[j,i]])["ValuesOnGrid"]]}/.AllbutMassand\[Kappa]Fixed[[j,i]],{i,Length[AllbutMassand\[Kappa]Fixed[[j]]]}],FrameLabel->{"Subscript[Log, 10]mpD [eV]","\[Kappa] [ ]"},PlotLabel->plotname,PlotLegends->Automatic],ListContourPlot[Table[{Log10@("mpD"(("c")^2/"JpereV")/.Constants`SIConstRepl),Log10@"\[Kappa]",Max[("Pf"/.AllbutMassand\[Kappa]Fixed[[j,i]])["ValuesOnGrid"]]}/.AllbutMassand\[Kappa]Fixed[[j,i]],{i,Length[AllbutMassand\[Kappa]Fixed[[j]]]}],ContourShading->False,ContourLabels->True,Contours->Join[Table[i,{i,-20,0,2}],Table[i,{i,-25,-71,-5}]]]}]];*)

,{j,Length@AllbutMassand\[Kappa]Fixed}];
(*,{j,1}];*)
]


(* ::Subsubsection::Closed:: *)
(*Plot Evaporation Rate*)


(* ::Input::Initialization:: *)
Clear[PlotEvaporationRate\[Kappa]Contour]
PlotEvaporationRate\[Kappa]Contour[KEchangeMassDicts_]:=Module[{AllbutMassand\[Kappa]Fixed,v0,meDovermpD,plotname,VEarth},
AllbutMassand\[Kappa]Fixed=Gather[Flatten[KEchangeMassDicts],(#1["v0"]==#2["v0"]&&#1["meD"]/#1["mpD"]==#2["meD"]/#2["mpD"])&];
Do[
(*AllbutMassFixed=[[i]];*)

(*\[Kappa]=Union["\[Kappa]"/.AllbutMassFixed][[1]];*)

v0=Union["v0"/.AllbutMassand\[Kappa]Fixed[[j]]][[1]];
(*fD=Union["fD"/.AllbutMassFixed][[1]];*)
meDovermpD=Union["meD"/"mpD"/.AllbutMassand\[Kappa]Fixed[[j]]][[1]];

(*plotname=("Subscript[\[CapitalGamma], evap] [s^-1]- "<>" Subscript[v, 0]"<>ToString[StringForm["= ``",N@(v0/1000)]]<>"[kms^-1]"<>",  Subscript[Log, 10]Subscript[m, Subscript[e, D]]/Subscript[m, Subscript[p, D]]="<>ToString[Log10@meDovermpD]);*)

VEarth=((4\[Pi])/3)("rE")^3/.Constants`EarthRepl;

(*Print[ListDensityPlot[{Log10@("mpD"(("c")^2/"JpereV")/.Constants`SIConstRepl),Log10@"\[Kappa]",Log10["\[CapitalGamma]total"("\[Kappa]")^2]}/.AllbutMassand\[Kappa]Fixed[[j]]/.AllbutMassand\[Kappa]Fixed[[j]],FrameLabel->{"Subscript[Log, 10]mpD [eV]","\[Kappa] [ ]"},PlotLabel->plotname,PlotLegends->Automatic]];*)
Do[With[{k=i,l=j},AssociateTo[AllbutMassand\[Kappa]Fixed[[l,k]],"\[CapitalGamma]total"->(AllbutMassand\[Kappa]Fixed[[l,k]]["\[CapitalGamma]total"]/.AllbutMassand\[Kappa]Fixed[[l,k]])]],{i,Length[AllbutMassand\[Kappa]Fixed[[j]]]}];

plotname=("\!\(\*SubscriptBox[\(\[CapitalGamma]\), \(evap\)]\) [\!\(\*SuperscriptBox[\(s\), \(-1\)]\)]");

Print[Show[{ListDensityPlot[{Log10@("mpD"(("c")^2/"JpereV")/.Constants`SIConstRepl),Log10@"\[Kappa]",Log10["\[CapitalGamma]total" ("\[Kappa]")^2]}/.AllbutMassand\[Kappa]Fixed[[j]],FrameLabel->{"\!\(\*SubscriptBox[\(Log\), \(10\)]\)mpD [eV]","\[Kappa] [ ]"},PlotLabel->plotname,PlotLegends->Automatic],ListContourPlot[{Log10@("mpD"(("c")^2/"JpereV")/.Constants`SIConstRepl),Log10@"\[Kappa]",Log10["\[CapitalGamma]total" ("\[Kappa]")^2]}/.AllbutMassand\[Kappa]Fixed[[j]],ContourShading->False,ContourLabels->True,Contours->{4,2,0,-2,-4,-6,-8,-10,-12,-14,-20,-30,-40,-50}]}]];

plotname=("\!\(\*SubscriptBox[\(\[CapitalGamma]\), \(evap, e\)]\) [\!\(\*SuperscriptBox[\(s\), \(-1\)]\)]");

Print[Show[{ListDensityPlot[{Log10@("mpD"(("c")^2/"JpereV")/.Constants`SIConstRepl),Log10@"\[Kappa]",Log10["\[CapitalGamma]totale" ("\[Kappa]")^2]}/.AllbutMassand\[Kappa]Fixed[[j]],FrameLabel->{"\!\(\*SubscriptBox[\(Log\), \(10\)]\)mpD [eV]","\[Kappa] [ ]"},PlotLabel->plotname,PlotLegends->Automatic],ListContourPlot[{Log10@("mpD"(("c")^2/"JpereV")/.Constants`SIConstRepl),Log10@"\[Kappa]",Log10["\[CapitalGamma]totale" ("\[Kappa]")^2]}/.AllbutMassand\[Kappa]Fixed[[j]],ContourShading->False,ContourLabels->True,Contours->{4,2,0,-2,-4,-6,-8,-10,-12,-14,-20,-30,-40,-50}]}]];

plotname=("\!\(\*SubscriptBox[\(\[CapitalGamma]\), \(evap, Nuc\)]\) [\!\(\*SuperscriptBox[\(s\), \(-1\)]\)]");

Print[Show[{ListDensityPlot[{Log10@("mpD"(("c")^2/"JpereV")/.Constants`SIConstRepl),Log10@"\[Kappa]",Log10["\[CapitalGamma]totalNuc" ("\[Kappa]")^2]}/.AllbutMassand\[Kappa]Fixed[[j]],FrameLabel->{"\!\(\*SubscriptBox[\(Log\), \(10\)]\)mpD [eV]","\[Kappa] [ ]"},PlotLabel->plotname,PlotLegends->Automatic],ListContourPlot[{Log10@("mpD"(("c")^2/"JpereV")/.Constants`SIConstRepl),Log10@"\[Kappa]",Log10["\[CapitalGamma]totalNuc" ("\[Kappa]")^2]}/.AllbutMassand\[Kappa]Fixed[[j]],ContourShading->False,ContourLabels->True,Contours->{4,2,0,-2,-4,-6,-8,-10,-12,-14,-20,-30,-40,-50}]}]];

,{j,Length@AllbutMassand\[Kappa]Fixed}];
(*,{j,1}];*)
]


(* ::Input::Initialization:: *)
(*Clear[PlotNc\[Kappa]Contour]
PlotNc\[Kappa]Contour[KEchangeMassDicts_]:=Module[{AllbutMassand\[Kappa]Fixed,v0,meDovermpD,plotname,VEarth},
AllbutMassand\[Kappa]Fixed=Gather[Flatten[KEchangeMassDicts],(#1["v0"]==#2["v0"]&&#1["meD"]/#1["mpD"]==#2["meD"]/#2["mpD"])&];
Do[
(*AllbutMassFixed=[[i]];*)

(*\[Kappa]=Union["\[Kappa]"/.AllbutMassFixed][[1]];*)

v0=Union["v0"/.AllbutMassand\[Kappa]Fixed[[j]]][[1]];
(*fD=Union["fD"/.AllbutMassFixed][[1]];*)
meDovermpD=Union[("meD")/("mpD")/.AllbutMassand\[Kappa]Fixed[[j]]][[1]];

plotname=("Subscript[N, c]0.05/Subscript[f, D] []- "<>" Subscript[v, 0]"<>ToString[StringForm["= ``",N@(v0/1000)]]<>"[kms^-1]"<>",  Subscript[Log, 10]Subscript[m, Subscript[e, D]]/Subscript[m, Subscript[p, D]]="<>ToString[Log10@meDovermpD]);

VEarth=(4\[Pi])/3("rE")^3/.Constants`EarthRepl;
Print[VEarth];
Print[ListDensityPlot[{Log10@("mpD"("c")^2/("JpereV")/.Constants`SIConstRepl),Log10@"\[Kappa]",Log10["Nc"]}/.AllbutMassand\[Kappa]Fixed[[j]]/.AllbutMassand\[Kappa]Fixed[[j]],FrameLabel->{"Subscript[Log, 10]mpD [eV]","\[Kappa] [ ]"},PlotLabel->plotname,PlotLegends->Automatic]];

,{j,Length@AllbutMassand\[Kappa]Fixed}];
(*,{j,1}];*)
]*)


(* ::Subsubsection::Closed:: *)
(*Plot Energy Fraction of Energy Lost*)


(* ::Input::Initialization:: *)
Clear[PlotELfrac]
PlotELfrac[PDict_,Electronic\[Sigma]Dicsts_,Nuclear\[Sigma]Dicts_]:=Module[{Pfunc,Pfs,PDom,\[Kappa],mpD,v0},
mpD=Union["mpD"/.PDict][[1]] ("c")^2/("JpereV")/.Constants`SIConstRepl;
v0=Union["v0"/.PDict][[1]];
Print[mpD];
Print[v0];
Pfs=InterpolatePc[PDict];
(*Total*)
Do[
PDom=("PfNuc"/.Pfs[[p]])["Domain"];
\[Kappa]="\[Kappa]"/.Pfs[[p]];

(*Electronic*)
(*Pfunc=ELAveragede[mpD(("c")^2/("JpereV"))^-1/.Constants`SIConstRepl,#1,\[Kappa],#2]&;*)
Pfunc=\[Kappa]^2 ((GeteELFunction[Electronic\[Sigma]Dicsts]["dEdlofr"][#1,#2])/(EK[mpD,#2](("JpereV")/("c")^2/.Constants`SIConstRepl)))&;
Print[Show[{DensityPlot[Log10@Pfunc[bE,10^v\[Chi]],{v\[Chi],PDom[[1,1]],PDom[[1,2]]},{bE,10^PDom[[2,1]],10^PDom[[2,2]]},PlotLegends->Automatic,AxesLabel->Automatic,PlotRange->All,FrameLabel->{"\!\(\*SubscriptBox[\(Log\), \(10\)]\)v\[Chi] [\!\(\*SuperscriptBox[\(ms\), \(-1\)]\)]","\!\(\*SubscriptBox[\(b\), \(E\)]\) [m]"},PlotLabel->"\!\(\*SubscriptBox[\(Log\), \(10\)]\)\!\(\*FractionBox[SubscriptBox[\(\[CapitalDelta]E\), \(soft, e\)], SubscriptBox[\(E\), \(\[Chi]\)]]\) - \!\(\*SubscriptBox[\(Log\), \(10\)]\)"<>ToString[StringForm["\[Kappa] = ``",Log10@\[Kappa]]]<>", \!\(\*SubscriptBox[\(Log\), \(10\)]\)\!\(\*SubscriptBox[\(m\), SubscriptBox[\(p\), \(D\)]]\)"<>ToString[StringForm[" = `` [eV]",Log10@mpD]]<>ToString[StringForm[", v0 = `` [m/s]",v0]]],ContourPlot[Log10@Pfunc[bE,10^v\[Chi]],{v\[Chi],PDom[[1,1]],PDom[[1,2]]},{bE,10^PDom[[2,1]],10^PDom[[2,2]]},ContourShading->False,ContourLabels->True]}]];


(*Nuclear*)
(*Pfunc=ELAveragedNuc0[mpD(("c")^2/("JpereV"))^-1/.Constants`SIConstRepl,#1,\[Kappa],#2]&;*)
Pfunc=\[Kappa]^2 ((GetNucELFunction[Nuclear\[Sigma]Dicts]["dEdlofr"][#1,#2])/(EK[mpD,#2](("JpereV")/("c")^2/.Constants`SIConstRepl)))&;
Print[Show[{DensityPlot[Log10@Pfunc[bE,10^v\[Chi]],{v\[Chi],PDom[[1,1]],PDom[[1,2]]},{bE,10^PDom[[2,1]],10^PDom[[2,2]]},PlotLegends->Automatic,AxesLabel->Automatic,PlotRange->All,FrameLabel->{"\!\(\*SubscriptBox[\(Log\), \(10\)]\)v\[Chi] [\!\(\*SuperscriptBox[\(ms\), \(-1\)]\)]","\!\(\*SubscriptBox[\(b\), \(E\)]\) [m]"},PlotLabel->"\!\(\*SubscriptBox[\(Log\), \(10\)]\)\!\(\*FractionBox[SubscriptBox[\(\[CapitalDelta]E\), \(soft, Nuc\)], SubscriptBox[\(E\), \(\[Chi]\)]]\) - \!\(\*SubscriptBox[\(Log\), \(10\)]\)"<>ToString[StringForm["\[Kappa] = ``",Log10@\[Kappa]]]<>", \!\(\*SubscriptBox[\(Log\), \(10\)]\)\!\(\*SubscriptBox[\(m\), SubscriptBox[\(p\), \(D\)]]\)"<>ToString[StringForm[" = `` [eV]",Log10@mpD]]<>ToString[StringForm[", v0 = `` [m/s]",v0]]],ContourPlot[Log10@Pfunc[bE,10^v\[Chi]],{v\[Chi],PDom[[1,1]],PDom[[1,2]]},{bE,10^PDom[[2,1]],10^PDom[[2,2]]},ContourShading->False,ContourLabels->True]}]];
,{p,Length[Pfs]}]
]


(* ::Input::Initialization:: *)
Clear[PlotELNucoverE]
PlotELNucoverE[PDict_,Electronic\[Sigma]Dicsts_,Nuclear\[Sigma]Dicts_]:=Module[{PfuncE,PfuncNuc,Pfs,PDom,\[Kappa],mpD,v0},
mpD=Union["mpD"/.PDict][[1]] ("c")^2/("JpereV")/.Constants`SIConstRepl;
v0=Union["v0"/.PDict][[1]];
Print[mpD];
Print[v0];
Pfs=InterpolatePc[PDict];
(*Total*)
Do[
PDom=("PfNuc"/.Pfs[[p]])["Domain"];
\[Kappa]="\[Kappa]"/.Pfs[[p]];

(*PfuncE=ELAveragede[mpD(("c")^2/("JpereV"))^-1/.Constants`SIConstRepl,#1,\[Kappa],#2]&;
PfuncNuc=ELAveragedNuc0[mpD(("c")^2/("JpereV"))^-1/.Constants`SIConstRepl,#1,\[Kappa],#2]&;*)
PfuncE=\[Kappa]^2 (GeteELFunction[Electronic\[Sigma]Dicsts]["dEdlofr"][#1,#2])&;
PfuncNuc=\[Kappa]^2 (GetNucELFunction[Nuclear\[Sigma]Dicts]["dEdlofr"][#1,#2])&;

(*Print[PfuncE[10^6,10^4]];
Print[PfuncNuc[10^6,10^4]];*)
(*Print[Show[{DensityPlot[-Log10[PfuncNuc[10^v\[Chi],bE]/PfuncE[10^v\[Chi],bE]],{v\[Chi],PDom[[1,1]],PDom[[1,2]]},{bE,10^PDom[[2,1]],10^PDom[[2,2]]},PlotLegends->Automatic,AxesLabel->Automatic,PlotRange->All,FrameLabel->{"Subscript[Log, 10]v\[Chi] [ms^-1]","Subscript[b, E] [m]"},PlotLabel->"Subscript[Log, 10]Subscript[\[CapitalDelta]E, soft,Nuc]/Subscript[\[CapitalDelta]E, soft,e] - Subscript[Log, 10]"<>", Subscript[Log, 10]Subscript[m, Subscript[p, D]]"<>ToString[StringForm[" = `` [eV]",Log10@mpD]]<>ToString[StringForm[", v0 = `` [km/s]",v0/1000]]],ContourPlot[-Log10[PfuncNuc[10^v\[Chi],bE]/PfuncE[10^v\[Chi],bE]],{v\[Chi],PDom[[1,1]],PDom[[1,2]]},{bE,10^PDom[[2,1]],10^PDom[[2,2]]},ContourShading->False,ContourLabels->True]}]];*)
Print[Show[{DensityPlot[-Log10[PfuncNuc[bE,10^v\[Chi]]/PfuncE[bE,10^v\[Chi]]],{v\[Chi],PDom[[1,1]],PDom[[1,2]]},{bE,10^PDom[[2,1]],10^PDom[[2,2]]},PlotLegends->Automatic,AxesLabel->Automatic,PlotRange->All,FrameLabel->{"\!\(\*SubscriptBox[\(Log\), \(10\)]\)v\[Chi] [\!\(\*SuperscriptBox[\(ms\), \(-1\)]\)]","\!\(\*SubscriptBox[\(b\), \(E\)]\) [m]"},PlotLabel->"\!\(\*SubscriptBox[\(Log\), \(10\)]\)\!\(\*FractionBox[SubscriptBox[\(\[CapitalDelta]E\), \(soft, Nuc\)], SubscriptBox[\(\[CapitalDelta]E\), \(soft, e\)]]\) - \!\(\*SubscriptBox[\(Log\), \(10\)]\)"<>", \!\(\*SubscriptBox[\(Log\), \(10\)]\)\!\(\*SubscriptBox[\(m\), SubscriptBox[\(p\), \(D\)]]\)"<>ToString[StringForm[" = `` [eV]",Log10@mpD]]<>ToString[StringForm[", v0 = `` [km/s]",v0/1000]]],ContourPlot[-Log10[PfuncNuc[bE,10^v\[Chi]]/PfuncE[bE,10^v\[Chi]]],{v\[Chi],PDom[[1,1]],PDom[[1,2]]},{bE,10^PDom[[2,1]],10^PDom[[2,2]]},ContourShading->False,ContourLabels->True]}]];

,{p,1}](*all \[Kappa] give the same plot*)
]


(* ::Input:: *)
(**)


(* ::Subsubsection::Closed:: *)
(*Plot Capture Rate - compare Nuc E Hard Soft*)


(* ::Input::Initialization:: *)
Clear[PlotCaptureRateComparison]
PlotCaptureRateComparison[KEchangeMassDicts_]:=Module[{AllbutMassand\[Kappa]Fixed,v0,meDovermpD,plotname,VEarth,list,analysiskeys},
AllbutMassand\[Kappa]Fixed=Gather[Flatten[KEchangeMassDicts],(#1["v0"]==#2["v0"]&&#1["meD"]/#1["mpD"]==#2["meD"]/#2["mpD"])&];
Do[

v0=Union["v0"/.AllbutMassand\[Kappa]Fixed[[j]]][[1]];

meDovermpD=Union["meD"/"mpD"/.AllbutMassand\[Kappa]Fixed[[j]]][[1]];

VEarth=((4\[Pi])/3)("rE")^3/.Constants`EarthRepl;

plotname=("\!\(\*SubscriptBox[\(Log\), \(10\)]\)\!\(\*FractionBox[SubscriptBox[\(\[CapitalGamma]\), \(cap, E\)], SubscriptBox[\(\[CapitalGamma]\), \(cap, Nuc\)]]\) [ ] - "<>" \!\(\*SubscriptBox[\(v\), \(0\)]\)"<>ToString[StringForm["= ``",N@(v0/1000)]]<>"[\!\(\*SuperscriptBox[\(kms\), \(-1\)]\)]"<>",  \!\(\*SubscriptBox[\(Log\), \(10\)]\)\!\(\*FractionBox[SubscriptBox[\(m\), SubscriptBox[\(e\), \(D\)]], SubscriptBox[\(m\), SubscriptBox[\(p\), \(D\)]]]\)="<>ToString[Log10@meDovermpD]);

Print[Show[{ListDensityPlot[{Log10@("mpD"(("c")^2/"JpereV")/.Constants`SIConstRepl),Log10@"\[Kappa]",Round@Log10@(("dNcMaxdt_PfE")/("dNcMaxdt_PfNuc"))}/.AllbutMassand\[Kappa]Fixed[[j]],FrameLabel->{"\!\(\*SubscriptBox[\(Log\), \(10\)]\)mpD [eV]","\[Kappa] [ ]"},PlotLabel->plotname,PlotLegends->Automatic,PlotRange->All],ListContourPlot[{Log10@("mpD"(("c")^2/"JpereV")/.Constants`SIConstRepl),Log10@"\[Kappa]",Round@Log10@(("dNcMaxdt_PfE")/("dNcMaxdt_PfNuc"))}/.AllbutMassand\[Kappa]Fixed[[j]],ContourShading->False,ContourLabels->True,PlotRange->All,Contours->Join[Table[i,{i,-10,0,1}],Table[i,{i,-20,-70,-10}]]],ListContourPlot[{Log10@("mpD" ("c")^2/("JpereV")/.Constants`SIConstRepl),Log10@"\[Kappa]",Round@Log10@(("dNcMaxdt_PfE")/("dNcMaxdt_PfNuc"))}/.AllbutMassand\[Kappa]Fixed[[j]]/."nD"->1,ContourShading->False,ContourLabels->True,Contours->{0},(*ContourStyle->Directive[Red,Thick]*)ContourStyle->Directive[Red,AbsoluteThickness[5]]],
ListContourPlot[{Log10@("mpD" ("c")^2/("JpereV")/.Constants`SIConstRepl),Log10@"\[Kappa]",Round@Log10@(("dNcMaxdt_PfE")/("dNcMaxdt_PfNuc"))}/.AllbutMassand\[Kappa]Fixed[[j]]/."nD"->1,ContourShading->False,ContourLabels->True,Contours->{-1},(*ContourStyle->Directive[Red,Thick]*)ContourStyle->Directive[Black,Dashed,AbsoluteThickness[5]]]}]];


plotname=("\!\(\*SubscriptBox[\(Log\), \(10\)]\)\!\(\*FractionBox[SubscriptBox[\(\[CapitalGamma]\), \(cap, hard, E\)], SubscriptBox[\(\[CapitalGamma]\), \(cap, hard, Nuc\)]]\) [ ] - "<>" \!\(\*SubscriptBox[\(v\), \(0\)]\)"<>ToString[StringForm["= ``",N@(v0/1000)]]<>"[\!\(\*SuperscriptBox[\(kms\), \(-1\)]\)]"<>",  \!\(\*SubscriptBox[\(Log\), \(10\)]\)\!\(\*FractionBox[SubscriptBox[\(m\), SubscriptBox[\(e\), \(D\)]], SubscriptBox[\(m\), SubscriptBox[\(p\), \(D\)]]]\)="<>ToString[Log10@meDovermpD]);

analysiskeys={"PfE","PfNuc","int\[Lambda]invE","int\[Lambda]invNuc"};

Print[Show[{ListDensityPlot[{Log10@("mpD"(("c")^2/"JpereV")/.Constants`SIConstRepl),Log10@"\[Kappa]",Log10@(("dNcMaxdt_int\[Lambda]invE")/("dNcMaxdt_int\[Lambda]invNuc"))}/.AllbutMassand\[Kappa]Fixed[[j]],FrameLabel->{"\!\(\*SubscriptBox[\(Log\), \(10\)]\)mpD [eV]","\[Kappa] [ ]"},PlotLabel->plotname,PlotLegends->Automatic,PlotRange->All],ListContourPlot[{Log10@("mpD"(("c")^2/"JpereV")/.Constants`SIConstRepl),Log10@"\[Kappa]",Log10@(("dNcMaxdt_int\[Lambda]invE")/("dNcMaxdt_int\[Lambda]invNuc"))}/.AllbutMassand\[Kappa]Fixed[[j]],ContourShading->False,ContourLabels->True],ListContourPlot[{Log10@("mpD" ("c")^2/("JpereV")/.Constants`SIConstRepl),Log10@"\[Kappa]",Log10@(("dNcMaxdt_int\[Lambda]invE")/("dNcMaxdt_int\[Lambda]invNuc"))}/.AllbutMassand\[Kappa]Fixed[[j]]/."nD"->1,ContourShading->False,ContourLabels->False,Contours->{0},(*ContourStyle->Directive[Red,Thick]*)ContourStyle->Directive[Red,AbsoluteThickness[5]]],ListContourPlot[{Log10@("mpD" ("c")^2/("JpereV")/.Constants`SIConstRepl),Log10@"\[Kappa]",Log10@(("dNcMaxdt_int\[Lambda]invE")/("dNcMaxdt_int\[Lambda]invNuc"))}/.AllbutMassand\[Kappa]Fixed[[j]]/."nD"->1,ContourShading->False,ContourLabels->True,Contours->{-1},(*ContourStyle->Directive[Red,Thick]*)ContourStyle->Directive[Black,Dashed,AbsoluteThickness[5]]]}]];


plotname=("\!\(\*SubscriptBox[\(Log\), \(10\)]\)\!\(\*FractionBox[SubscriptBox[\(\[CapitalGamma]\), \(cap, soft\)], SubscriptBox[\(\[CapitalGamma]\), \(cap, hard\)]]\) [] - "<>" \!\(\*SubscriptBox[\(v\), \(0\)]\)"<>ToString[StringForm["= ``",N@(v0/1000)]]<>"[\!\(\*SuperscriptBox[\(kms\), \(-1\)]\)]"<>",  \!\(\*SubscriptBox[\(Log\), \(10\)]\)\!\(\*FractionBox[SubscriptBox[\(m\), SubscriptBox[\(e\), \(D\)]], SubscriptBox[\(m\), SubscriptBox[\(p\), \(D\)]]]\)="<>ToString[Log10@meDovermpD]);

Print[Show[{ListDensityPlot[{Log10@("mpD"(("c")^2/"JpereV")/.Constants`SIConstRepl),Log10@"\[Kappa]",Log10@Max[("dNcMaxdt_PfE"+"dNcMaxdt_PfNuc"-"dNcMaxdt_int\[Lambda]invE"-"dNcMaxdt_int\[Lambda]invNuc"),10^-100]-Log10@Max[("dNcMaxdt_int\[Lambda]invE"+"dNcMaxdt_int\[Lambda]invNuc"),10^-100]}/.AllbutMassand\[Kappa]Fixed[[j]]/."nD"->1,FrameLabel->{"\!\(\*SubscriptBox[\(Log\), \(10\)]\)mpD [eV]","\[Kappa] [ ]"},PlotLabel->plotname,PlotLegends->Automatic,PlotRange->All],ListContourPlot[{Log10@("mpD"(("c")^2/"JpereV")/.Constants`SIConstRepl),Log10@"\[Kappa]",Log10@Max[("dNcMaxdt_PfE"+"dNcMaxdt_PfNuc"-"dNcMaxdt_int\[Lambda]invE"-"dNcMaxdt_int\[Lambda]invNuc"),10^-100]-Log10@Max[("dNcMaxdt_int\[Lambda]invE"+"dNcMaxdt_int\[Lambda]invNuc"),10^-100]}/.AllbutMassand\[Kappa]Fixed[[j]]/."nD"->1,ContourShading->False,ContourLabels->True(*,Contours->Join[Table[i,{i,-21,-1,2}],Table[i,{i,-25,-71,-5}]]*),PlotRange->All],ListContourPlot[{Log10@("mpD" ("c")^2/("JpereV")/.Constants`SIConstRepl),Log10@"\[Kappa]",Log10@Max[("dNcMaxdt_PfE"+"dNcMaxdt_PfNuc"-"dNcMaxdt_int\[Lambda]invE"-"dNcMaxdt_int\[Lambda]invNuc"),10^-100]-Log10@Max[("dNcMaxdt_int\[Lambda]invE"+"dNcMaxdt_int\[Lambda]invNuc"),10^-100]}/.AllbutMassand\[Kappa]Fixed[[j]]/."nD"->1,ContourShading->False,ContourLabels->False,Contours->{0},(*ContourStyle->Directive[Red,Thick]*)ContourStyle->Directive[Red,AbsoluteThickness[5]]]}]];

plotname=("\!\(\*SubscriptBox[\(Log\), \(10\)]\)\!\(\*FractionBox[SubscriptBox[\(\[CapitalGamma]\), \(cap, soft, E\)], SubscriptBox[\(\[CapitalGamma]\), \(cap, hard, E\)]]\) [] - "<>" \!\(\*SubscriptBox[\(v\), \(0\)]\)"<>ToString[StringForm["= ``",N@(v0/1000)]]<>"[\!\(\*SuperscriptBox[\(kms\), \(-1\)]\)]"<>",  \!\(\*SubscriptBox[\(Log\), \(10\)]\)\!\(\*FractionBox[SubscriptBox[\(m\), SubscriptBox[\(e\), \(D\)]], SubscriptBox[\(m\), SubscriptBox[\(p\), \(D\)]]]\)="<>ToString[Log10@meDovermpD]);

Print[Show[{ListDensityPlot[{Log10@("mpD"(("c")^2/"JpereV")/.Constants`SIConstRepl),Log10@"\[Kappa]",Log10@Max[("dNcMaxdt_PfE"-"dNcMaxdt_int\[Lambda]invE"),10^-100]-Log10@Max[("dNcMaxdt_int\[Lambda]invE"),10^-100]}/.AllbutMassand\[Kappa]Fixed[[j]]/."nD"->1,FrameLabel->{"\!\(\*SubscriptBox[\(Log\), \(10\)]\)mpD [eV]","\[Kappa] [ ]"},PlotLabel->plotname,PlotLegends->Automatic,PlotRange->All],ListContourPlot[{Log10@("mpD"(("c")^2/"JpereV")/.Constants`SIConstRepl),Log10@"\[Kappa]",Log10@Max[("dNcMaxdt_PfE"-"dNcMaxdt_int\[Lambda]invE"),10^-100]-Log10@Max[("dNcMaxdt_int\[Lambda]invE"),10^-100]}/.AllbutMassand\[Kappa]Fixed[[j]]/."nD"->1,ContourShading->False,ContourLabels->True(*,Contours->Join[Table[i,{i,-21,-1,2}],Table[i,{i,-25,-71,-5}]]*),PlotRange->All],ListContourPlot[{Log10@("mpD" ("c")^2/("JpereV")/.Constants`SIConstRepl),Log10@"\[Kappa]",Log10@Max[("dNcMaxdt_PfE"-"dNcMaxdt_int\[Lambda]invE"),10^-100]-Log10@Max[("dNcMaxdt_int\[Lambda]invE"),10^-100]}/.AllbutMassand\[Kappa]Fixed[[j]]/."nD"->1,ContourShading->False,ContourLabels->False,Contours->{0},(*ContourStyle->Directive[Red,Thick]*)ContourStyle->Directive[Red,AbsoluteThickness[5]]]}]];

plotname=("\!\(\*SubscriptBox[\(Log\), \(10\)]\)\!\(\*FractionBox[SubscriptBox[\(\[CapitalGamma]\), \(cap, soft, Nuc\)], SubscriptBox[\(\[CapitalGamma]\), \(cap, hard, Nuc\)]]\) [] - "<>" \!\(\*SubscriptBox[\(v\), \(0\)]\)"<>ToString[StringForm["= ``",N@(v0/1000)]]<>"[\!\(\*SuperscriptBox[\(kms\), \(-1\)]\)]"<>",  \!\(\*SubscriptBox[\(Log\), \(10\)]\)\!\(\*FractionBox[SubscriptBox[\(m\), SubscriptBox[\(e\), \(D\)]], SubscriptBox[\(m\), SubscriptBox[\(p\), \(D\)]]]\)="<>ToString[Log10@meDovermpD]);

Print[Show[{ListDensityPlot[{Log10@("mpD"(("c")^2/"JpereV")/.Constants`SIConstRepl),Log10@"\[Kappa]",Log10@Max[("dNcMaxdt_PfNuc"-"dNcMaxdt_int\[Lambda]invNuc"),10^-100]-Log10@Max[("dNcMaxdt_int\[Lambda]invNuc"),10^-100]}/.AllbutMassand\[Kappa]Fixed[[j]]/."nD"->1,FrameLabel->{"\!\(\*SubscriptBox[\(Log\), \(10\)]\)mpD [eV]","\[Kappa] [ ]"},PlotLabel->plotname,PlotLegends->Automatic,PlotRange->All],ListContourPlot[{Log10@("mpD"(("c")^2/"JpereV")/.Constants`SIConstRepl),Log10@"\[Kappa]",Log10@Max[("dNcMaxdt_PfNuc"-"dNcMaxdt_int\[Lambda]invNuc"),10^-100]-Log10@Max[("dNcMaxdt_int\[Lambda]invNuc"),10^-100]}/.AllbutMassand\[Kappa]Fixed[[j]]/."nD"->1,ContourShading->False,ContourLabels->True(*,Contours->Join[Table[i,{i,-21,-1,2}],Table[i,{i,-25,-71,-5}]]*),PlotRange->All],ListContourPlot[{Log10@("mpD" ("c")^2/("JpereV")/.Constants`SIConstRepl),Log10@"\[Kappa]",Log10@Max[("dNcMaxdt_PfNuc"-"dNcMaxdt_int\[Lambda]invNuc"),10^-100]-Log10@Max[("dNcMaxdt_int\[Lambda]invNuc"),10^-100]}/.AllbutMassand\[Kappa]Fixed[[j]]/."nD"->1,ContourShading->False,ContourLabels->False,Contours->{0},(*ContourStyle->Directive[Red,Thick]*)ContourStyle->Directive[Red,AbsoluteThickness[5]]]}]];

plotname=("\!\(\*SubscriptBox[\(Log\), \(10\)]\)\!\(\*SubscriptBox[\(\[CapitalGamma]\), \(cap, E\)]\) [\!\(\*SuperscriptBox[\(m\), \(-3\)]\)\!\(\*SuperscriptBox[\(s\), \(-1\)]\)] - "<>" \!\(\*SubscriptBox[\(v\), \(0\)]\)"<>ToString[StringForm["= ``",N@(v0/1000)]]<>"[\!\(\*SuperscriptBox[\(kms\), \(-1\)]\)]"<>",  \!\(\*SubscriptBox[\(Log\), \(10\)]\)\!\(\*FractionBox[SubscriptBox[\(m\), SubscriptBox[\(e\), \(D\)]], SubscriptBox[\(m\), SubscriptBox[\(p\), \(D\)]]]\)="<>ToString[Log10@meDovermpD]);

Do[With[{k=i,l=j},AssociateTo[AllbutMassand\[Kappa]Fixed[[l,k]],"dNcMaxdt_PfE"->(AllbutMassand\[Kappa]Fixed[[l,k]]["dNcMaxdt_PfE"]/.AllbutMassand\[Kappa]Fixed[[l,k]])]],{i,Length[AllbutMassand\[Kappa]Fixed[[j]]]}];

Print[Show[{ListDensityPlot[{Log10@("mpD"(("c")^2/"JpereV")/.Constants`SIConstRepl),Log10@"\[Kappa]",Log10@(("dNcMaxdt_PfE")/VEarth)}/.AllbutMassand\[Kappa]Fixed[[j]],FrameLabel->{"\!\(\*SubscriptBox[\(Log\), \(10\)]\)mpD [eV]","\[Kappa] [ ]"},PlotLabel->plotname,PlotLegends->Automatic],ListContourPlot[{Log10@("mpD"(("c")^2/"JpereV")/.Constants`SIConstRepl),Log10@"\[Kappa]",Log10@(("dNcMaxdt_PfE")/VEarth)}/.AllbutMassand\[Kappa]Fixed[[j]],ContourShading->False,ContourLabels->True,Contours->Join[Table[i,{i,-18,2,2}],Table[i,{i,-18,-118,-5}]]]}]];

plotname=("\!\(\*SubscriptBox[\(Log\), \(10\)]\)\!\(\*SubscriptBox[\(\[CapitalGamma]\), \(cap, Nuc\)]\) [\!\(\*SuperscriptBox[\(m\), \(-3\)]\)\!\(\*SuperscriptBox[\(s\), \(-1\)]\)] - "<>" \!\(\*SubscriptBox[\(v\), \(0\)]\)"<>ToString[StringForm["= ``",N@(v0/1000)]]<>"[\!\(\*SuperscriptBox[\(kms\), \(-1\)]\)]"<>",  \!\(\*SubscriptBox[\(Log\), \(10\)]\)\!\(\*FractionBox[SubscriptBox[\(m\), SubscriptBox[\(e\), \(D\)]], SubscriptBox[\(m\), SubscriptBox[\(p\), \(D\)]]]\)="<>ToString[Log10@meDovermpD]);

Do[With[{k=i,l=j},AssociateTo[AllbutMassand\[Kappa]Fixed[[l,k]],"dNcMaxdt_PfNuc"->(AllbutMassand\[Kappa]Fixed[[l,k]]["dNcMaxdt_PfNuc"]/.AllbutMassand\[Kappa]Fixed[[l,k]])]],{i,Length[AllbutMassand\[Kappa]Fixed[[j]]]}];

Print[Show[{ListDensityPlot[{Log10@("mpD"(("c")^2/"JpereV")/.Constants`SIConstRepl),Log10@"\[Kappa]",Log10@(("dNcMaxdt_PfNuc")/VEarth)}/.AllbutMassand\[Kappa]Fixed[[j]],FrameLabel->{"\!\(\*SubscriptBox[\(Log\), \(10\)]\)mpD [eV]","\[Kappa] [ ]"},PlotLabel->plotname,PlotLegends->Automatic],ListContourPlot[{Log10@("mpD"(("c")^2/"JpereV")/.Constants`SIConstRepl),Log10@"\[Kappa]",Log10@(("dNcMaxdt_PfNuc")/VEarth)}/.AllbutMassand\[Kappa]Fixed[[j]],ContourShading->False,ContourLabels->True,Contours->Join[Table[i,{i,-18,2,2}],Table[i,{i,-18,-118,-5}]]]}]];


plotname=("\!\(\*SubscriptBox[\(Log\), \(10\)]\)\!\(\*SubscriptBox[\(\[CapitalGamma]\), \(cap, Hard, Nuc\)]\) [\!\(\*SuperscriptBox[\(m\), \(-3\)]\)\!\(\*SuperscriptBox[\(s\), \(-1\)]\)] - "<>" \!\(\*SubscriptBox[\(v\), \(0\)]\)"<>ToString[StringForm["= ``",N@(v0/1000)]]<>"[\!\(\*SuperscriptBox[\(kms\), \(-1\)]\)]"<>",  \!\(\*SubscriptBox[\(Log\), \(10\)]\)\!\(\*FractionBox[SubscriptBox[\(m\), SubscriptBox[\(e\), \(D\)]], SubscriptBox[\(m\), SubscriptBox[\(p\), \(D\)]]]\)="<>ToString[Log10@meDovermpD]);

Do[With[{k=i,l=j},AssociateTo[AllbutMassand\[Kappa]Fixed[[l,k]],"dNcMaxdt_int\[Lambda]invNuc"->(AllbutMassand\[Kappa]Fixed[[l,k]]["dNcMaxdt_int\[Lambda]invNuc"]/.AllbutMassand\[Kappa]Fixed[[l,k]])]],{i,Length[AllbutMassand\[Kappa]Fixed[[j]]]}];

Print[Show[{ListDensityPlot[{Log10@("mpD"(("c")^2/"JpereV")/.Constants`SIConstRepl),Log10@"\[Kappa]",Log10@(("dNcMaxdt_int\[Lambda]invNuc")/VEarth)}/.AllbutMassand\[Kappa]Fixed[[j]],FrameLabel->{"\!\(\*SubscriptBox[\(Log\), \(10\)]\)mpD [eV]","\[Kappa] [ ]"},PlotLabel->plotname,PlotLegends->Automatic,PlotRange->All],ListContourPlot[{Log10@("mpD"(("c")^2/"JpereV")/.Constants`SIConstRepl),Log10@"\[Kappa]",Log10@(("dNcMaxdt_int\[Lambda]invNuc")/VEarth)}/.AllbutMassand\[Kappa]Fixed[[j]],ContourShading->False,ContourLabels->True,Contours->Join[Table[i,{i,-18,2,2}],Table[i,{i,-18,-118,-5}]],PlotRange->All]}]];

plotname=("\!\(\*SubscriptBox[\(Log\), \(10\)]\)\!\(\*SubscriptBox[\(\[CapitalGamma]\), \(cap, Hard, e\)]\) [\!\(\*SuperscriptBox[\(m\), \(-3\)]\)\!\(\*SuperscriptBox[\(s\), \(-1\)]\)] - "<>" \!\(\*SubscriptBox[\(v\), \(0\)]\)"<>ToString[StringForm["= ``",N@(v0/1000)]]<>"[\!\(\*SuperscriptBox[\(kms\), \(-1\)]\)]"<>",  \!\(\*SubscriptBox[\(Log\), \(10\)]\)\!\(\*FractionBox[SubscriptBox[\(m\), SubscriptBox[\(e\), \(D\)]], SubscriptBox[\(m\), SubscriptBox[\(p\), \(D\)]]]\)="<>ToString[Log10@meDovermpD]);

Do[With[{k=i,l=j},AssociateTo[AllbutMassand\[Kappa]Fixed[[l,k]],"dNcMaxdt_int\[Lambda]invE"->(AllbutMassand\[Kappa]Fixed[[l,k]]["dNcMaxdt_int\[Lambda]invE"]/.AllbutMassand\[Kappa]Fixed[[l,k]])]],{i,Length[AllbutMassand\[Kappa]Fixed[[j]]]}];

Print[Show[{ListDensityPlot[{Log10@("mpD"(("c")^2/"JpereV")/.Constants`SIConstRepl),Log10@"\[Kappa]",Log10@(("dNcMaxdt_int\[Lambda]invE")/VEarth)}/.AllbutMassand\[Kappa]Fixed[[j]],FrameLabel->{"\!\(\*SubscriptBox[\(Log\), \(10\)]\)mpD [eV]","\[Kappa] [ ]"},PlotLabel->plotname,PlotLegends->Automatic,PlotRange->All],ListContourPlot[{Log10@("mpD"(("c")^2/"JpereV")/.Constants`SIConstRepl),Log10@"\[Kappa]",Log10@(("dNcMaxdt_int\[Lambda]invE")/VEarth)}/.AllbutMassand\[Kappa]Fixed[[j]],ContourShading->False,ContourLabels->True,Contours->Join[Table[i,{i,-18,2,2}],Table[i,{i,-18,-118,-5}]],PlotRange->All]}]];

plotname=("\!\(\*SubscriptBox[\(Log\), \(10\)]\)\!\(\*SubscriptBox[\(\[CapitalGamma]\), \(cap, soft, Nuc\)]\) [\!\(\*SuperscriptBox[\(m\), \(-3\)]\)\!\(\*SuperscriptBox[\(s\), \(-1\)]\)] - "<>" \!\(\*SubscriptBox[\(v\), \(0\)]\)"<>ToString[StringForm["= ``",N@(v0/1000)]]<>"[\!\(\*SuperscriptBox[\(kms\), \(-1\)]\)]"<>",  \!\(\*SubscriptBox[\(Log\), \(10\)]\)\!\(\*FractionBox[SubscriptBox[\(m\), SubscriptBox[\(e\), \(D\)]], SubscriptBox[\(m\), SubscriptBox[\(p\), \(D\)]]]\)="<>ToString[Log10@meDovermpD]);

Print[Show[{ListDensityPlot[{Log10@("mpD"(("c")^2/"JpereV")/.Constants`SIConstRepl),Log10@"\[Kappa]",Log10@Max[("dNcMaxdt_PfNuc"-"dNcMaxdt_int\[Lambda]invNuc")/VEarth,10^-100]}/.AllbutMassand\[Kappa]Fixed[[j]],FrameLabel->{"\!\(\*SubscriptBox[\(Log\), \(10\)]\)mpD [eV]","\[Kappa] [ ]"},PlotLabel->plotname,PlotLegends->Automatic,PlotRange->All],ListContourPlot[{Log10@("mpD"(("c")^2/"JpereV")/.Constants`SIConstRepl),Log10@"\[Kappa]",Log10@Max[("dNcMaxdt_PfNuc"-"dNcMaxdt_int\[Lambda]invNuc"),10^-100]}/.AllbutMassand\[Kappa]Fixed[[j]],ContourShading->False,ContourLabels->True(*,Contours->Join[Table[i,{i,-21,-1,2}],Table[i,{i,-25,-71,-5}]]*),PlotRange->All],ListContourPlot[{Log10@("mpD" ("c")^2/("JpereV")/.Constants`SIConstRepl),Log10@"\[Kappa]",Log10@Max[("dNcMaxdt_PfNuc"-"dNcMaxdt_int\[Lambda]invNuc")/VEarth,10^-100]}/.AllbutMassand\[Kappa]Fixed[[j]],ContourShading->False,ContourLabels->False,Contours->{0},(*ContourStyle->Directive[Red,Thick]*)ContourStyle->Directive[Red,AbsoluteThickness[5]]]}]];

plotname=("\!\(\*SubscriptBox[\(Log\), \(10\)]\)\!\(\*SubscriptBox[\(\[CapitalGamma]\), \(cap, soft, e\)]\) [\!\(\*SuperscriptBox[\(m\), \(-3\)]\)\!\(\*SuperscriptBox[\(s\), \(-1\)]\)] - "<>" \!\(\*SubscriptBox[\(v\), \(0\)]\)"<>ToString[StringForm["= ``",N@(v0/1000)]]<>"[\!\(\*SuperscriptBox[\(kms\), \(-1\)]\)]"<>",  \!\(\*SubscriptBox[\(Log\), \(10\)]\)\!\(\*FractionBox[SubscriptBox[\(m\), SubscriptBox[\(e\), \(D\)]], SubscriptBox[\(m\), SubscriptBox[\(p\), \(D\)]]]\)="<>ToString[Log10@meDovermpD]);

Print[Show[{ListDensityPlot[{Log10@("mpD"(("c")^2/"JpereV")/.Constants`SIConstRepl),Log10@"\[Kappa]",Log10@Max[("dNcMaxdt_PfE"-"dNcMaxdt_int\[Lambda]invE")/VEarth,10^-100]}/.AllbutMassand\[Kappa]Fixed[[j]],FrameLabel->{"\!\(\*SubscriptBox[\(Log\), \(10\)]\)mpD [eV]","\[Kappa] [ ]"},PlotLabel->plotname,PlotLegends->Automatic,PlotRange->All],ListContourPlot[{Log10@("mpD"(("c")^2/"JpereV")/.Constants`SIConstRepl),Log10@"\[Kappa]",Log10@Max[("dNcMaxdt_PfE"-"dNcMaxdt_int\[Lambda]invE")/VEarth,10^-100]}/.AllbutMassand\[Kappa]Fixed[[j]],ContourShading->False,ContourLabels->True(*,Contours->Join[Table[i,{i,-21,-1,2}],Table[i,{i,-25,-71,-5}]]*),PlotRange->All],ListContourPlot[{Log10@("mpD" ("c")^2/("JpereV")/.Constants`SIConstRepl),Log10@"\[Kappa]",Log10@Max[("dNcMaxdt_PfE"-"dNcMaxdt_int\[Lambda]invE")/VEarth,10^-100]}/.AllbutMassand\[Kappa]Fixed[[j]],ContourShading->False,ContourLabels->False,Contours->{0},(*ContourStyle->Directive[Red,Thick]*)ContourStyle->Directive[Red,AbsoluteThickness[5]]]}]];

,{j,Length@AllbutMassand\[Kappa]Fixed}];
(*,{j,1}];*)
]


(* ::Input:: *)
(**)


(* ::Subsection:: *)
(*MC Processes*)


(* ::Subsubsection::Closed:: *)
(*Gravitational Potential*)


Vgrav[r_]:=Piecewise[{{-(("G" "ME")/r),r>"rE"/. Constants`EarthRepl},{-(("ME" "G")/(2("rE")^3))( 3("rE")^2-r^2),r<"rE"/. Constants`EarthRepl}},-(("G" "ME")/("rE"/. Constants`EarthRepl))]/.Constants`SIConstRepl/.Constants`EarthRepl


(* ::Subsubsection::Closed:: *)
(*EL and \[CapitalDelta]r total*)


(* ::Input::Initialization:: *)
dc[bE_] := 2 Sqrt[("rcore")^2-bE^2]/.Constants`EarthRepl (*path length in the core for given impact parameter*)
dE[bE_] := 2 Sqrt[("rE")^2-bE^2]/.Constants`EarthRepl (*path length through the earth*)
dM[bE_]:= If[(("rcore")^2/.Constants`EarthRepl)>bE^2,dE[bE] - dc[bE],Re[dE[bE]]](*path length through the Mantle, iff condition is satisfied does it go through the core*)


(* ::Input::Initialization:: *)
Clear[GetTotalELFunction]
GetTotalELFunction[Electronic\[Sigma]Dicsts_,Nuclear\[Sigma]Dicts_]:=Module[{Region\[Sigma]Dicts,RegionELfuncs,RegionELDomains,CrustIndex,CoreIndex,rc,rE,dEdlofr,\[CapitalDelta]EthroughEarth},
Region\[Sigma]Dicts=Gather[Flatten[Join[Electronic\[Sigma]Dicsts,Nuclear\[Sigma]Dicts]],#1[["regionindex"]]==#2[["regionindex"]]&];

CrustIndex=FirstPosition["regionindex"/.Region\[Sigma]Dicts,1][[1]];
CoreIndex = FirstPosition["regionindex"/.Region\[Sigma]Dicts,2][[1]];

rc = "rcore"/.Constants`EarthRepl;
rE = "rE"/.Constants`EarthRepl;

(*Returns: dEdl[r,v\[Chi]]*)
dEdlofr=Quiet@Re@Total@Piecewise[{{10^"dEdlf"[ Log10[#2]]HeavisideTheta[Log10[#2]-"dEdlf"["Domain"][[1,1]]]/.Region\[Sigma]Dicts[[CrustIndex,;;]],rE> #1 >=  rc},{10^"dEdlf"[ Log10[#2]]HeavisideTheta[Log10[#2]-"dEdlf"["Domain"][[1,1]]]/.Region\[Sigma]Dicts[[CoreIndex,;;]], rc>#1}},0]&;

\[CapitalDelta]EthroughEarth=Quiet@Re@Piecewise[{{(dM[#1] dEdlofr[1.1rc,#2] + dc[#1] dEdlofr[0.9rc,#2]),rc^2>#1^2},{dM[#1] dEdlofr[1.1rc,#2],rE^2>#1^2>=rc^2}},0]&;

<|"dEdlofr"->dEdlofr,"\[CapitalDelta]EthroughEarth"->\[CapitalDelta]EthroughEarth|>
]


(* ::Input::Initialization:: *)
Clear[GeteELFunction]
GeteELFunction[Electronic\[Sigma]Dicsts_]:=Module[{Region\[Sigma]Dicts,RegionELfuncs,RegionELDomains,CrustIndex,CoreIndex,rc,rE,dEdlofr,\[CapitalDelta]EthroughEarth},
Region\[Sigma]Dicts=Gather[Flatten[Join[Electronic\[Sigma]Dicsts]],#1[["regionindex"]]==#2[["regionindex"]]&];

CrustIndex=FirstPosition["regionindex"/.Region\[Sigma]Dicts,1][[1]];
CoreIndex = FirstPosition["regionindex"/.Region\[Sigma]Dicts,2][[1]];

rc = "rcore"/.Constants`EarthRepl;
rE = "rE"/.Constants`EarthRepl;

(*Returns: dEdl[r,v\[Chi]]*)
dEdlofr=Quiet@Re@Total@Piecewise[{{10^"dEdlf"[ Log10[#2]]HeavisideTheta[Log10[#2]-"dEdlf"["Domain"][[1,1]]]/.Region\[Sigma]Dicts[[CrustIndex,;;]],rE> #1 >=  rc},{10^"dEdlf"[ Log10[#2]]HeavisideTheta[Log10[#2]-"dEdlf"["Domain"][[1,1]]]/.Region\[Sigma]Dicts[[CoreIndex,;;]], rc>#1}},0]&;


(*Returns: \[CapitalDelta]E[bE,v\[Chi]]*)
\[CapitalDelta]EthroughEarth=Quiet@Re@Piecewise[{{(dM[#1] dEdlofr[1.1rc,#2] + dc[#1] dEdlofr[0.9rc,#2]),rc^2>#1^2},{dM[#1] dEdlofr[1.1rc,#2],rE^2>#1^2>=rc^2}},0]&;

<|"dEdlofr"->dEdlofr,"\[CapitalDelta]EthroughEarth"->\[CapitalDelta]EthroughEarth|>
]


(* ::Input::Initialization:: *)
Clear[GetNucELFunction]
GetNucELFunction[Nuclear\[Sigma]Dicts_]:=Module[{Region\[Sigma]Dicts,RegionELfuncs,RegionELDomains,CrustIndex,CoreIndex,rc,rE,dEdlofr,\[CapitalDelta]EthroughEarth},
Region\[Sigma]Dicts=Gather[Flatten[Join[Nuclear\[Sigma]Dicts]],#1[["regionindex"]]==#2[["regionindex"]]&];

CrustIndex=FirstPosition["regionindex"/.Region\[Sigma]Dicts,1][[1]];
CoreIndex = FirstPosition["regionindex"/.Region\[Sigma]Dicts,2][[1]];

rc = "rcore"/.Constants`EarthRepl;
rE = "rE"/.Constants`EarthRepl;

(*Returns: dEdl[r,v\[Chi]]*)
dEdlofr=Quiet@Re@Total@Piecewise[{{10^"dEdlf"[ Log10[#2]]HeavisideTheta[Log10[#2]-"dEdlf"["Domain"][[1,1]]]/.Region\[Sigma]Dicts[[CrustIndex,;;]],rE> #1 >=  rc},{10^"dEdlf"[ Log10[#2]]HeavisideTheta[Log10[#2]-"dEdlf"["Domain"][[1,1]]]/.Region\[Sigma]Dicts[[CoreIndex,;;]], rc>#1}},0]&;


(*Returns: \[CapitalDelta]E[bE,v\[Chi]]*)
\[CapitalDelta]EthroughEarth=Quiet@Re@Piecewise[{{(dM[#1] dEdlofr[1.1rc,#2] + dc[#1] dEdlofr[0.9rc,#2]),rc^2>#1^2},{dM[#1] dEdlofr[1.1rc,#2],rE^2>#1^2>=rc^2}},0]&;

<|"dEdlofr"->dEdlofr,"\[CapitalDelta]EthroughEarth"->\[CapitalDelta]EthroughEarth|>
]


(* ::Input::Initialization:: *)
EK[m\[Chi]_,v\[Chi]_] :=1/2 m\[Chi] v\[Chi]^2
\[CapitalDelta]rTotby\[Kappa]squared[m\[Chi]_,v\[Chi]_,bE_,\[CapitalDelta]Efunc_,y_:10^-2]:=Re[dE[bE] y EK[m\[Chi],v\[Chi]]/\[CapitalDelta]Efunc[bE,v\[Chi]] ]


(* ::Subsubsection::Closed:: *)
(*Weak Capture Regime Trajectory - Linear*)


(* ::Input::Initialization:: *)
(*Clear[GetWeakCaptureTrajectory]
GetWeakCaptureTrajectory[m\[Chi]_,v\[Chi]E_,bE_,\[Kappa]_,ELfunc_]:=Module[{RE,lsol,ldom,voft},
(*this version crashes*)
(*use get GetWeakCaptureProbability* with this and GetCaptureRateOldTraj*)
RE = dE[bE];(*path length through the Earth*)
lsol=l/.NDSolve[{D[l[t],{t,2}]==-1/m\[Chi]\[Kappa]^2("JpereV"/.Constants`SIConstRepl)ELfunc[Sqrt[bE^2+l[t]^2],l'[t]],l'[0]==v\[Chi]E,l[0]==-RE/2,WhenEvent[l[t]==RE/2,"StopIntegration"]},l,{t,0,10^4}][[1]];(*Solution for trajectory of particle along a straight line through the earth parallel to it's velocity as a function of time.*)
ldom=lsol["Domain"][[1]];
(*Print["endpoint vs max"];
Print[lsol[ldom[[2]]]];
Print[RE/2];*)
voft=lsol';
<|"l(t)"->lsol,"v(t)"->voft,"domain"->ldom,"RE"->RE,"\[Kappa]"->\[Kappa],"bE"->bE,"m\[Chi]"->m\[Chi]|>
]*)


(* ::Input::Initialization:: *)
Clear[GetWeakCaptureTrajectoryinv]
GetWeakCaptureTrajectoryinv[m\[Chi]_,v\[Chi]E_,bE_,\[Kappa]_,ELfunc_]:=Module[{RE,vsol,ldom,voft},
(*use get GetWeakCaptureProbability*fromv with this and GetCaptureRate*)
(*first order solution for trajectory*)
RE = dE[bE];(*path length through the Earth*)
vsol=v/.NDSolve[{D[v[l],l]==-(1/(m\[Chi] v[l])) \[Kappa]^2 ("JpereV"/.Constants`SIConstRepl)ELfunc[Sqrt[bE^2+l^2],v[l]],v[0]==v\[Chi]E,WhenEvent[v[l]<="vesc"/.Constants`EarthRepl,"StopIntegration"]},v,{l,-(RE/2),RE/2}][[1]];
(*vsol=v/.NDSolve[{D[v[l],l]==-((\[Kappa]^2("JpereV"/.Constants`SIConstRepl))/m\[Chi])ELfunc[Sqrt[bE^2+l^2],v[l]],v[0]==v\[Chi]E},v,{l,-(RE/2),RE/2}][[1]];*)(*Solution for trajectory of particle along a straight line through the earth parallel to it's velocity as a function of time.*)
ldom=vsol["Domain"][[1]];

(*Print[vsol["Domain"]];*)

<|"v(l)"->vsol,"domain"->ldom,"RE"->RE,"\[Kappa]"->\[Kappa],"bE"->bE,"m\[Chi]"->m\[Chi]|>
]


(* ::Input::Initialization:: *)
Clear[GetWeakCaptureTrajectoryinvFE]
GetWeakCaptureTrajectoryinvFE[m\[Chi]_,v\[Chi]E_,bE_,\[Kappa]_,ELfunc_,N_:100]:=Module[{RE,h,vofls,vsol,ldom,voft},
(*use get GetWeakCaptureProbability*fromv with this and GetCaptureRate*)
(*Forward Euler solution for trajectory*)

RE = dE[bE];(*path length through the Earth*)
h=RE/N;
vofls=Table[{h s-RE/2,0},{s,N}];
vofls[[1,2]]=v\[Chi]E;

Do[
(*vofls[[s+1,2]]=vofls[[s,2]]-h(\[Kappa]^2("JpereV"/.Constants`SIConstRepl))/(m\[Chi] vofls[[s,2]])ELfunc[Sqrt[bE^2+vofls[[s,1]]^2],vofls[[s,2]]];*)(*I don't think this JpereV should be there.*)
vofls[[s+1,2]]=vofls[[s,2]]-h \[Kappa]^2/(m\[Chi] vofls[[s,2]]) ELfunc[Sqrt[bE^2+vofls[[s,1]]^2],vofls[[s,2]]];
If[vofls[[s+1,2]]<="vesc"/.Constants`EarthRepl,Break[]];
,{s,N-1}]; (*forward euler for dv/dl=-1/(Subscript[m, \[Chi]]v)dE/dl - (ie. from energy conservation / Raleigh-Lagrange Formalism)*)

vsol=Interpolation[vofls]; (*Solution for trajectory of particle along a straight line through the earth parallel to it's velocity as a function of length.*)
ldom=vsol["Domain"][[1]];

(*Print[vsol["Domain"]];*)

<|"v(l)"->vsol,"domain"->ldom,"RE"->RE,"\[Kappa]"->\[Kappa],"bE"->bE,"m\[Chi]"->m\[Chi]|>
]


(* ::Subsubsection::Closed:: *)
(*Weak Capture Regime Trajectory - w Potential - Forward Euler*)


(* ::Input::Initialization:: *)
Clear[GetWeakCaptureTrajectoryinvFEPotential]
GetWeakCaptureTrajectoryinvFEPotential[m\[Chi]_,v\[Chi]E_,bE_,\[Kappa]_,ELfunc_,N_:50,Vofr_:Vgrav]:=Module[{rE,\[CapitalDelta]l,E0,L0,steps,traj,ti,l,r,v,vr,\[CapitalDelta]E,L,Log\[CapitalDelta]L,vsol,ldom,voft},
(*Forward Euler solution for trajectory*)

rE="rE"/.Constants`EarthRepl;
\[CapitalDelta]l=rE/N; (*\[CapitalDelta]l*)

E0 = m\[Chi]/2 v\[Chi]E^2 +m\[Chi] Vofr[rE];
L0 = m\[Chi] v\[Chi]E bE;

traj={<|"l"->0,"r"->rE,"v"->v\[Chi]E,"vr"->-v\[Chi]E Sqrt[1 - (bE/rE)^2],"\[CapitalDelta]E"->0,"V"->m\[Chi] Vofr[rE],"L"->L0,"Log\[CapitalDelta]L"->0|>};
(*initialize list of points along trajectory at Earth's surface
	l - Path length
	r - Distance from Earth's centre
	v\[Chi] - aDM speed
	vr - radial aDM speed (initially pointing inward)
	\[CapitalDelta]E - Energy lost from dissipation
	L - Angular Momentum
	Log\[CapitalDelta]L - Log of time dependent part of L 
*)

steps = 0;
Do[
steps++;
ti = traj[[-1]];(*current trajectory point*)

l = ti["l"]+\[CapitalDelta]l;
r = ti["r"]+\[CapitalDelta]l ti["vr"]/ti["v"];

If[E0- m\[Chi] Vofr[r]-ti["\[CapitalDelta]E"]< 0,Break[]];
v = Sqrt[(2/m\[Chi]) (E0 - m\[Chi] Vofr[r] - ti["\[CapitalDelta]E"]) ];
\[CapitalDelta]E = ti["\[CapitalDelta]E"] + \[CapitalDelta]l \[Kappa]^2 ELfunc[ti["r"],ti["v"]];
vr = ti["vr"] + \[CapitalDelta]l ((ti["v"]^2-ti["vr"]^2)/(ti["r"]ti["v"]) - (m\[Chi] Vofr'[r])/(m\[Chi] ti["v"])-ti["vr"]/(m\[Chi] ti["v"]^2) \[Kappa]^2 ELfunc[ti["r"],ti["v"]]);(*assumes ELfunc > 0, so we include - for dissipation*)
L=m\[Chi] ti["r"]Sqrt[ti["v"]^2-ti["vr"]^2];
Log\[CapitalDelta]L = ti["Log\[CapitalDelta]L"]+\[CapitalDelta]l   (\[Kappa]^2 ELfunc[ti["r"],ti["v"]])/(m\[Chi] ti["v"]^2);


AppendTo[traj,<|"l"->l,"r"->r,"v"->v,"vr"->vr,"\[CapitalDelta]E"->\[CapitalDelta]E,"V"->m\[Chi] Vofr[r],"L"->L,"Log\[CapitalDelta]L"->Log\[CapitalDelta]L|>];

(*If[Abs@vr<="vesc"/.Constants`EarthRepl,Break[]];*)
If[r>rE,Break[]];
,{s,3N}];(*will break once reaches surface (<~ 2 N) this is just a buffer*)

vsol=Interpolation[{"l","v"}/.traj]; 
ldom=vsol["Domain"][[1]];

<|"v(l)"->vsol,"domain"->ldom,"traj"->traj,"\[Kappa]"->\[Kappa],"bE"->bE,"m\[Chi]"->m\[Chi],"v\[Chi]E"->v\[Chi]E,"E0"->E0,"L0"->L0,"steps"->steps |>
]


(* ::Subsubsection::Closed:: *)
(*Weak Capture Regime Trajectory - w Potential - Runga Kutta (ish)*)


(* ::Input::Initialization:: *)
Clear[GetWeakCaptureTrajectoryinvRK2wPotential]
GetWeakCaptureTrajectoryinvRK2wPotential[m\[Chi]_,v\[Chi]E_,bE_,\[Kappa]_,ELfunc_,N_:70,Vofr_:Vgrav]:=Module[{rE,\[CapitalDelta]l,E0,L0,steps,traj,ti,l,r,v,vr,\[CapitalDelta]E,L,Log\[CapitalDelta]L,Lcomped,vsol,ldom,lsol,rdom,voft,rfinal(*,vfinal,vesc*),lMC1,lMC2,lMC},
(*slightly better than FW for conserved quantities, not quite an actual RK2*)
(*
ELfunc - dE/dl
Vofr - (V(r))/m\[Chi]
*)

rE="rE"/.Constants`EarthRepl;
\[CapitalDelta]l=rE/N; (*\[CapitalDelta]l*)

E0 = m\[Chi]/2 v\[Chi]E^2 +m\[Chi] Vofr[rE];
L0 = m\[Chi] v\[Chi]E bE;

traj={<|"l"->0,"r"->rE,"v"->v\[Chi]E,"vr"->-v\[Chi]E Sqrt[1 - (bE/rE)^2],"\[CapitalDelta]E"->0,"V"->m\[Chi] Vofr[rE],"L"->L0,"Log\[CapitalDelta]L"->0(*,"Lcomped"->L0*)|>};
(*initialize list of points along trajectory at Earth's surface
	l - Path length
	r - Distance from Earth's centre
	v\[Chi] - aDM speed
	vr - radial aDM speed (initially pointing inward)
	\[CapitalDelta]E - Energy lost from dissipation
	L - Angular Momentum
	Log\[CapitalDelta]L - Log of time dependent part of L 
*)

steps = 0;
Do[
steps++;
ti = traj[[-1]];(*current trajectory point*)

l = ti["l"]+\[CapitalDelta]l;
r = ti["r"]+\[CapitalDelta]l ti["vr"]/ti["v"];

If[E0- m\[Chi] Vofr[r]-ti["\[CapitalDelta]E"]< 0,Break[]]; (*numerical, weaker than the below - will be satisfied when v<0, above will be when v<vesc*)
v = \[Sqrt](2/m\[Chi] (E0 - m\[Chi] Vofr[r] - ti["\[CapitalDelta]E"]) );(*this choice of incrementation, preserves energy but wildly violates L conservation / behaviour*)
(*v =ti["v"]+\[CapitalDelta]l (-(( m\[Chi] ti["vr"])/(m\[Chi] ti["v"]^2))Vofr'[r] - 1/(m\[Chi] ti["v"])\[Kappa]^2ELfunc[r,ti["v"]]);*)
\[CapitalDelta]E = ti["\[CapitalDelta]E"] + \[CapitalDelta]l \[Kappa]^2 ELfunc[r,v];


vr = ti["vr"] + \[CapitalDelta]l ((v^2-ti["vr"]^2)/(r v) - (m\[Chi] Vofr'[r])/(m\[Chi] v)-ti["vr"]/(m\[Chi] v^2) \[Kappa]^2 ELfunc[r,v]);(*assumes ELfunc > 0, so we include - for dissipation*)
(*vr = ti["vr"] + \[CapitalDelta]l (L^2/(m\[Chi]^2 r^3 v) - (m\[Chi] Vofr'[r])/(m\[Chi] v)-ti["vr"]/(m\[Chi] v^2)\[Kappa]^2ELfunc[r,v]);*)

L=m\[Chi] r Sqrt[v^2-vr^2];(*this choice of incrementation, preserves energy but wildly violates L conservation / behaviour*)
(*L=m\[Chi] ti["r"]Sqrt[ti["v"]^2-ti["vr"]^2];*)
(*Lcomped = ti["Lcomped"] (1-\[CapitalDelta]l (\[Kappa]^2ELfunc[r,v])/(m\[Chi] v^2));*)
(*L = ti["L"] (1-\[CapitalDelta]l (\[Kappa]^2ELfunc[r,v])/(m\[Chi] v^2));*)

(*Log\[CapitalDelta]L = ti["Log\[CapitalDelta]L"]+\[CapitalDelta]l   (\[Kappa]^2ELfunc[r,v])/(m\[Chi] v^2);*)
Log\[CapitalDelta]L = ti["Log\[CapitalDelta]L"]+ (\[CapitalDelta]l \[Kappa]^2 ELfunc[r,v])/(m\[Chi] v^2);


AppendTo[traj,<|"l"->l,"r"->r,"v"->v,"vr"->vr,"\[CapitalDelta]E"->\[CapitalDelta]E,"V"->m\[Chi] Vofr[r],"L"->L,"Log\[CapitalDelta]L"->Log\[CapitalDelta]L(*,"Lcomped"->Lcomped*)|>];

If[m\[Chi]/2 v^2+ m\[Chi] Vofr[r]< 0,Break[]];(*particle is captured, ie. bound to the earth, placed after the append statement so that at least one additional point is added. This is needed for the captured condition used further down the pipeline in the case where initial conditions already imply the particle is bound.*)

(*If[Abs@vr<="vesc"/.Constants`EarthRepl,Break[]];*)
If[r>rE,Break[]];
,{s,3N}];(*will break once reaches surface (<~ 2 N) this is just a buffer*)

If[Max[Abs[(("m\[Chi]")/2 ("v")^2+"V"+"\[CapitalDelta]E")/E0-1]/.traj]>0.20,Print["Energy conservation has been violated by more than 20% in Forward Euler: ",Max@Abs[(("m\[Chi]")/2 ("v")^2+"V"+"\[CapitalDelta]E")/E0-1]/.traj]];
(*If[Max[Abs[("L"Exp["Log\[CapitalDelta]L"])/L0-1]/.traj]>0.20,Print["Angular Momentum conservation has been violated by more than 20% in Forward Euler: ",Max[Abs[("L"Exp["Log\[CapitalDelta]L"])/L0-1]/.traj]]];*)

(*vfinal =  traj[[-1]]["v"];
vesc = Sqrt[(2 Abs@traj[[-1]]["V"])/m\[Chi]]; (*vesc(V(Subscript[r, final]))*)*)

vsol=Quiet@Interpolation[{"l","v"}/.traj]; 
ldom=vsol["Domain"][[1]];

lsol = Quiet@Interpolation[{"l","r"}/.traj];
rdom = lsol["Domain"][[1]];

rfinal = "r"/.traj[[-1]]; (*this is the capture condition - if we stop before we reach the surface, we are captured, otherwise we have left the Earth with E>0 so motion is unbounded. *)

(*values of l where particle crosses mantle-core boundary*)
If[rfinal>"rE"/.Constants`EarthRepl,
(*passes completely through the earth, so approximately symmetric about the l domain. *)
If[Min["r"/.traj]<"rcore"/.Constants`EarthRepl,
(*particle passes through the core*)
(*lMC1 = l/.Check[FindRoot[lsol[l]-"rcore"/.Constants`EarthRepl,{l,rdom[[2]]/4,rdom[[1]],rdom[[2]]/2}],Print[{m\[Chi],v\[Chi]E,bE,\[Kappa]}];(*{l->"broken"}*)FindRoot[lsol[l]-"rcore"/.Constants`EarthRepl,{l,rdom[[2]]/4,rdom[[1]],rdom[[2]]/2}]];
lMC2 = l/.Check[FindRoot[lsol[l]-"rcore"/.Constants`EarthRepl,{l,(3rdom[[2]])/4,rdom[[2]]/2,rdom[[2]]}],Print[{m\[Chi],v\[Chi]E,bE,\[Kappa]}];(*{l->"broken"}*)FindRoot[lsol[l]-"rcore"/.Constants`EarthRepl,{l,(3rdom[[2]])/4,rdom[[2]]/2,rdom[[2]]}]];
lMC=If[lMC1=="broken"||lMC2=="broken",{"Cap"},{lMC1,lMC2}];,*)
lMC1 = l/.Quiet@FindRoot[lsol[l]-"rcore"/.Constants`EarthRepl,{l,rdom[[2]]/4,rdom[[1]],rdom[[2]]/2}];
lMC2 = l/.Quiet@FindRoot[lsol[l]-"rcore"/.Constants`EarthRepl,{l,(3rdom[[2]])/4,rdom[[2]]/2,rdom[[2]]}];
lMC={lMC1,lMC2};,
(*doesn't pass through the core (No Core)*)
lMC={"NC"}
],
(*It's captured somewhere, so we won't need to distinguish between crust and core for hard capture.*)
lMC = {"Cap"}
];

(*,"vfinal"->vfinal,"vesc"->vesc*)
<|"v(l)"->vsol,"domain"->ldom,"r(l)"->lsol,"rdom"->rdom,"traj"->traj,"\[Kappa]"->\[Kappa],"bE"->bE,"m\[Chi]"->m\[Chi],"v\[Chi]E"->v\[Chi]E,"E0"->E0,"L0"->L0,"steps"->steps,"rfinal"->rfinal,"lMCs"->lMC|>
]


(* ::Subsubsection::Closed:: *)
(*unit test w FE*)


(* ::Input::Initialization:: *)
(*Clear[GetWeakCaptureTrajectoryinvRK2wPotential]
GetWeakCaptureTrajectoryinvRK2wPotential[m\[Chi]_,v\[Chi]E_,bE_,\[Kappa]_,ELfunc_,N_:30,Vofr_:Vgrav]:=Module[{rE,\[CapitalDelta]l,E0,L0,steps,traj,ti,l,r,v,vr,\[CapitalDelta]E,L,Log\[CapitalDelta]L,Lcomped,vsol,ldom,lsol,rdom,voft,rfinal(*,vfinal,vesc*),lMC1,lMC2,lMC},
(*slightly better than FW for conserved quantities, not quite an actual RK2*)
(*
ELfunc - dE/dl
Vofr - (V(r))/m\[Chi]
*)

rE="rE"/.Constants`EarthRepl;
\[CapitalDelta]l=rE/N; (*\[CapitalDelta]l*)

E0 = m\[Chi]/2v\[Chi]E^2 +m\[Chi] Vofr[rE];
L0 = m\[Chi] v\[Chi]E bE;

traj={<|"l"->0,"r"->rE,"v"->v\[Chi]E,"vr"->-v\[Chi]E Sqrt[1 - (bE/rE)^2],"\[CapitalDelta]E"->0,"V"->m\[Chi] Vofr[rE],"L"->L0,"Log\[CapitalDelta]L"->0,"Lcomped"->L0|>};
(*initialize list of points along trajectory at Earth's surface
	l - Path length
	r - Distance from Earth's centre
	v\[Chi] - aDM speed
	vr - radial aDM speed (initially pointing inward)
	\[CapitalDelta]E - Energy lost from dissipation
	L - Angular Momentum
	Log\[CapitalDelta]L - Log of time dependent part of L 
*)

steps = 0;
Do[
steps++;
ti = traj[[-1]];(*current trajectory point*)

l = ti["l"]+\[CapitalDelta]l;
r = ti["r"]+\[CapitalDelta]l ti["vr"]/ti["v"];

If[E0- m\[Chi] Vofr[ ti["r"]]-ti["\[CapitalDelta]E"]< 0,Break[]]; (*numerical, weaker than the below - will be satisfied when v<0, above will be when v<vesc*)
v = \[Sqrt](2/m\[Chi](E0 - m\[Chi] Vofr[ti["r"]] - ti["\[CapitalDelta]E"]) );
\[CapitalDelta]E = ti["\[CapitalDelta]E"] + \[CapitalDelta]l \[Kappa]^2ELfunc[ti["r"],ti["v"]];
vr = ti["vr"] + \[CapitalDelta]l ((ti["v"]^2-ti["vr"]^2)/(ti["r"]ti["v"]) - (m\[Chi] Vofr'[ti["r"]])/(m\[Chi] ti["v"])-ti["vr"]/(m\[Chi] ti["v"]^2)\[Kappa]^2ELfunc[ti["r"],ti["v"]]);(*assumes ELfunc > 0, so we include - for dissipation*)
L=m\[Chi] ti["r"]Sqrt[ti["v"]^2-ti["vr"]^2];
(*Log\[CapitalDelta]L = ti["Log\[CapitalDelta]L"]+\[CapitalDelta]l   (\[Kappa]^2ELfunc[r,v])/(m\[Chi] v^2);*)
Log\[CapitalDelta]L = ti["Log\[CapitalDelta]L"]+ (\[CapitalDelta]l \[Kappa]^2ELfunc[ti["r"],ti["v"]])/(m\[Chi] ti["v"]^2);
Lcomped = ti["Lcomped"] (1-\[CapitalDelta]l (\[Kappa]^2ELfunc[ti["r"],ti["v"]])/(m\[Chi] ti["v"]^2));

AppendTo[traj,<|"l"->l,"r"->r,"v"->v,"vr"->vr,"\[CapitalDelta]E"->\[CapitalDelta]E,"V"->m\[Chi] Vofr[r],"L"->L,"Log\[CapitalDelta]L"->Log\[CapitalDelta]L,"Lcomped"->Lcomped|>];

If[m\[Chi]/2v^2+ m\[Chi] Vofr[r]< 0,Break[]];(*particle is captured, ie. bound to the earth, placed after the append statement so that at least one additional point is added. This is needed for the captured condition used further down the pipeline in the case where initial conditions already imply the particle is bound.*)

(*If[Abs@vr<="vesc"/.Constants`EarthRepl,Break[]];*)
If[r>rE,Break[]];
,{s,3N}];(*will break once reaches surface (<~ 2 N) this is just a buffer*)

If[Max[Abs[(("m\[Chi]")/2("v")^2+"V"+"\[CapitalDelta]E")/E0-1]/.traj]>0.20,Print["Energy conservation has been violated by more than 20% in Forward Euler: ",Max@Abs[(("m\[Chi]")/2("v")^2+"V"+"\[CapitalDelta]E")/E0-1]/.traj]];
(*If[Max[Abs[("L"Exp["Log\[CapitalDelta]L"])/L0-1]/.traj]>0.20,Print["Angular Momentum conservation has been violated by more than 20% in Forward Euler: ",Max[Abs[("L"Exp["Log\[CapitalDelta]L"])/L0-1]/.traj]]];*)

(*vfinal =  traj[[-1]]["v"];
vesc = Sqrt[(2 Abs@traj[[-1]]["V"])/m\[Chi]]; (*vesc(V(Subscript[r, final]))*)*)

vsol=Quiet@Interpolation[{"l","v"}/.traj]; 
ldom=vsol["Domain"][[1]];

lsol = Quiet@Interpolation[{"l","r"}/.traj];
rdom = lsol["Domain"][[1]];

rfinal = "r"/.traj[[-1]]; (*this is the capture condition - if we stop before we reach the surface, we are captured, otherwise we have left the Earth with E>0 so motion is unbounded. *)

(*values of l where particle crosses mantle-core boundary*)
If[rfinal>"rE"/.Constants`EarthRepl,
(*passes completely through the earth, so approximately symmetric about the l domain. *)
If[Min["r"/.traj]<"rcore"/.Constants`EarthRepl,
(*particle passes through the core*)
lMC1 = l/.Check[FindRoot[lsol[l]-"rcore"/.Constants`EarthRepl,{l,rdom[[2]]/4,rdom[[1]],rdom[[2]]/2}],Print[{m\[Chi],v\[Chi]E,bE,\[Kappa]}];(*{l->"broken"}*)FindRoot[lsol[l]-"rcore"/.Constants`EarthRepl,{l,rdom[[2]]/4,rdom[[1]],rdom[[2]]/2}]];
lMC2 = l/.Check[FindRoot[lsol[l]-"rcore"/.Constants`EarthRepl,{l,(3rdom[[2]])/4,rdom[[2]]/2,rdom[[2]]}],Print[{m\[Chi],v\[Chi]E,bE,\[Kappa]}];(*{l->"broken"}*)FindRoot[lsol[l]-"rcore"/.Constants`EarthRepl,{l,(3rdom[[2]])/4,rdom[[2]]/2,rdom[[2]]}]];
lMC=If[lMC1=="broken"||lMC2=="broken",{"Cap"},{lMC1,lMC2}];,
(*doesn't pass through the core (No Core)*)
lMC={"NC"}
],
(*It's captured somewhere, so we won't need to distinguish between crust and core for hard capture.*)
lMC = {"Cap"}
];

(*,"vfinal"->vfinal,"vesc"->vesc*)
<|"v(l)"->vsol,"domain"->ldom,"r(l)"->lsol,"rdom"->rdom,"traj"->traj,"\[Kappa]"->\[Kappa],"bE"->bE,"m\[Chi]"->m\[Chi],"v\[Chi]E"->v\[Chi]E,"E0"->E0,"L0"->L0,"steps"->steps,"rfinal"->rfinal,"lMCs"->lMC|>
]*)


(* ::Subsection:: *)
(*Weak Capture Probabilities*)


(* ::Subsubsection::Closed:: *)
(*Old Weak capture regime probability*)


(* ::Input::Initialization:: *)
(*Clear[GetWeakCaptureProbability]
GetWeakCaptureProbability[lDict_,\[Sigma]Dict_]:= Module[{v\[Chi]min,v\[Chi]raw,v\[Chi],\[Chi]traj,bE,m\[Chi],\[Omega]min,params,tdom,\[Kappa],regionindex,\[Omega]capture,\[Delta]\[Xi],\[Xi]of\[Omega]andt,\[Omega]maxoft,d\[Sigma]dERoftand\[Omega],\[Sigma]oftand\[Omega]min,\[Lambda]invoft,integralof\[Lambda]invwrtr,SurvivalProb,CaptureProb},
\[Omega]min="\[Omega]min"/.\[Sigma]Dict;
params=\[Sigma]Dict["params"];
\[Delta]\[Xi]=(10^"\[Xi]"/.\[Sigma]Dict)[[1]];
m\[Chi]=("m\[Chi]"/.lDict);
tdom = "domain"/.lDict;
\[Kappa]=("\[Kappa]"/.lDict);
regionindex = "regionindex"/.\[Sigma]Dict;

v\[Chi]min=3Sqrt[("\[HBar]" \[Omega]min)/(2 m\[Chi])]/.params; (*minimum allowable velocity given the lower bound on energy transfer, times an O(1) fudge factor for numerical stability (3) *)
v\[Chi]raw=("v(t)"/.lDict);
v\[Chi]=If[v\[Chi]raw[#]>v\[Chi]min,("v(t)"/.lDict)[#],0]&;
\[Chi]traj=("l(t)"/.lDict);
bE=("bE"/.lDict);(*impact parameter*)

\[Xi]of\[Omega]andt=Abs[EnergyLoss`\[Xi]of\[Omega]andv\[Chi][#1,\[Omega]min,("v(t)"/.lDict)[#2],m\[Chi],params]-\[Delta]\[Xi]]&;
\[Omega]maxoft=EnergyLoss`\[Omega]maxofv\[Chi][("v(t)"/.lDict)[#1],m\[Chi],params]&;

(*d\[Sigma]dERoftand\[Omega]=10^Re@("d\[Sigma]f"/.\[Sigma]Dict)[Log10@v\[Chi][#1],Log10@\[Xi]of\[Omega]andt[#2,#1]]&; (*not currently used*)*)

\[Sigma]oftand\[Omega]min=10^(Re[("\[Sigma]of\[Xi]f"/.\[Sigma]Dict)[Log10@v\[Chi][#1],Log10@\[Xi]of\[Omega]andt[#2,#1]]])&;

(*Post \[Omega]capture*)
\[Omega]capture = 1/(2 ("\[HBar]"/.params))m\[Chi] (v\[Chi][#]^2 - ("vesc"/.Constants`EarthRepl)^2)&;


\[Lambda]invoft=If[v\[Chi]raw[#]>v\[Chi]min && \[Omega]maxoft[#]>\[Omega]capture[#],\[Kappa]^2("ne"/.params)\[Sigma]oftand\[Omega]min[#,\[Omega]capture[#]],0]&;
(*\[Lambda]invoft=If[v\[Chi]raw[#]>v\[Chi]min ,\[Kappa]^2("ne"/.params)\[Sigma]oftand\[Omega]min[#,\[Omega]capture[#]],0]&;*)

(*Print[tdom];*)

(*integralof\[Lambda]invwrtr=Total@Table[\[Lambda]invoft[t]v\[Chi][t],{t,Round[tdom[[1]]],Round[tdom[[2]]],1}];*)
integralof\[Lambda]invwrtr=Total@Table[\[Lambda]invoft[t]v\[Chi][t]Region\[CapitalTheta][ltor[\[Chi]traj[t],bE],regionindex],{t,Round[tdom[[1]]],Round[tdom[[2]]],1}];

integralof\[Lambda]invwrtr
]*)


(* ::Input:: *)
(*Image[CompressedData["*)
(*1:eJztXW9QFEfap973w75fqLsv+8niw1vUW3xI+cEvlnVeVYq6KqmUV1fkUvFM*)
(*DsXoK3duRcMletQe5R05IcFsOJcQscgiiQJxQwicMWgWgWWRvyFAMAiCBGRP*)
(*RERENLq5qTv3Zqane3pmev7s7iyzS/pXG7Ps9j7T/XRPP7/u5+ln/nfva7/O*)
(*+6+UlJQ3/of959evFGS+/vorzhd+yv7x4qE3Dvzu0P7/f+7Q4f2/2//6lr3/*)
(*zX74f+x/t9ny3PswBQUFBQUFBQUFBQUFBQUFBQUFBQUFBQUFBQUFBQUFBQUF*)
(*BQUFBQUFBQUFBQUFBQUFBQUFBQUFBQUFBQUFBQUFBQUFBQUFBQUFBQUFBQUF*)
(*BQUFBQUFBQUFBQUFBQUFBQUFBQUFBQUFBQUFRbKCebA0PTTY5usDr+7xVatr*)
(*REERNeaHhJE8FrS6KusIVKuJglsjPikC38zefcjELjg00+/z9U+ZIYo1Kg/v*)
(*TszcMkMSh9Wpbtjar4N4/aAyRm5J/yZhJJbqhO7PTgWT2TCa3ACporun1ko1*)
(*zN2uk1W5vzi+E3sdqZ1RlAvdm/5u9oEBgYYbEro3OzOj3kxsiK6tRkxFsmsj*)
(*NP9NwFfvLvh0TL2Ids3NwJN7NydmjV9iwC2M5AZ//OoU1pwcSXMjbw9InSfr*)
(*XNPmWB0kolaNq5QziECfrL2+H5LJiVildwIVxZ5mX/+EQpS2ZCvmpY78FAJS*)
(*n9njGVyMnmwwA0VpgqQDvtgaNXb2+U12Gyspv8PQlR8G2VkGoH+CzL2mPVth*)
(*Q237LmLVg8qAlyIrR1omQtwJ/HnbhlT291s901EJsBpxaYBU0WulGmakqlIg*)
(*KttPvnnss/ff5V51nYtYmcWA+3TedrbMB95JAyINNGTR76391V/Lf370o7Pq*)
(*zcSG6FpqxGQksTaY4IUD6Tate91QzWPDfGdz4QtlKixaDWvFWDQmR5m+QjOX*)
(*nFtS1TpP1rkmzLE6SFytGlApszhYd/BZu03ypc3+q/LhZVFOxCoVf5CaVTmq*)
(*YbATYF6CSrJn5hVwyBW1kZbfsawvgIjVi/uglGzvHVMqqD98JTcGROoWZ2tQ*)
(*Rlskat9Yfk1xLXipsU8LBLywUaqmAq11lxbQtZPUCMWnAUDxG18Aqq0IxDZk*)
(*jGKwIoufhbLqLi2plZnx5oGZKhLGotmQm2fLWTOnY+nYRQ8QkZdpT97Bkrza*)
(*uOXNtvGm4NncAnc9cWFqqOaxYbL2A/V9PzWsLWMBnQumR9jR2NzIGljPy+m4*)
(*fVV0HupcXAh4Hx/GkrhaNaDS6cYd/F5A6oZtuezHjl9C5doya5BiI1bp6tRl*)
(*T3FBLrCeaSVDqjVMgHkJGmnx2sxia36a0ppHiOXxJndBgbtpPFrSI6+gzvBd*)
(*7jjC950t/ZdOTzO/xdLscYL+tGV6JnHSIiWK9sJe+KWcsShrEXMXUcZCAlBv*)
(*fCYoVUz6joBZqLB7RbVQVIxFsyGRWTqg7iQdLMmrDVDz3zRpzF6UsYidS+7o*)
(*O76DyJ7a7al604a+QHOQuFo1pAEmeOVLzP/BjLrgarpoICqBmOjhtzMM3l/W*)
(*zUtEU4ysUmbNrPCR4GIDjivoegl8E5S6XQgONML6xIgs5NI7vl1gFtuPq0ud*)
(*rslkbw3bpqIumSsL0q+NrlHxC9nWlugZio6xGGkOUsxZRwaQlOE4q+kP1PZU*)
(*QoTuzyInGMm/aKzb0BWRsP6JBek149UAAMMTFPNgaaynD0XJYq+hCTHOhFld*)
(*mOoDn/dMLTxQ8W4ixlKquNPDqxNdQGxb2W/BTHWyzNunH52r3pCH418N+K6w*)
(*L3/RO8DSeYq+BJ9wr66ph0SBEcwMTx8/fERsKvNohfQFauPgt4otpuAg+Gpg*)
(*ZB5+9OT2cCerz+9WZaKYf3zF66Sjf+57mRQrtREbTK95JFgZHwLDrPHYSWBb*)
(*Xz3WJg71rmkFwWYH/E2ud3x9gaGOd9RsK/NoYWwECOkbW5L2IxgM3E30JHgt*)
(*wBYY54cMc3eUFdt5LfhEUUsD1vBa+UZuDfkyF2Ggv+aLM2NJAq1Gp4HJU5uB*)
(*apXujEgFGr+/jJbkQ8HMCmflQR5Iyk+FT2xHOnAHL4u0HV5xA4PgQCOoCBGD*)
(*puHyLMyNY0s/cPGWrIwKJFKBD8qW7WV/vDxcd3BbRnp6+jZna+u5goKCM4GW*)
(*fDv7bX6HopZb9+3bKGlmdIzFSHM0PItKscxiV4nMvcXe9nXXpFZ/eVi230qS*)
(*Z6zbwtMtCt+ozDsahwYoVKjN/2+Pev9UkSONkhVfWY0BvtiTqd7K/WXSb8sO*)
(*l/XNosbOjQlzlNf7KijwhwsY8wFnAdDWCvmlujpTbwhakpNfr7TcJAo0OjM8*)
(*/eHe9V6/v/fGimxuePr49kjAHxiZ//7f8t9crRbCeI5+Mi/5ghmuzOZbmt3U*)
(*w8+w4421jiyh7TmH/XNi0WV/aTn3+YsNfqVnzTJtxAzTax4J0CYA+ZXnw3f7*)
(*lgZ9pbvKSCVx2xq6cemzw9ulBbZXVbb+A44WMOA/OF3bkCcUcL3T8q3X4RI6*)
(*fV/LiGxgGbGGq6NX4HaA5YwlCbQanQbQ1LzrgnwdZTVjATsJKdA4mwLtPRZR*)
(*B7BcWlqa4N51bBN8R7YdjZDbITcXCvzQYiw2m00hC+53oBASJAi59ORRJKsX*)
(*dqUAD9ZyB/RnIWS4RsOBIzbwRta8rZ7euuf4d4JnKDbGotUcmf8vBY+Ikbn4*)
(*USN4UVhkEe6phI52Lkj6eUdBgeP5ZwSG8JOf7cHkGes2UCx1QwbvG8WDmeLX*)
(*AIIKNe6mpQH3i4CZnMg/7H3/3XoQO8e+HK83cBGzNYOs0WWu+ZxZkKUc/Ywt*)
(*VpAtTA55pQOCRW1v0Jq4hBlpMVADAnHPOYHt/sUJ55tCaK4iOtdQQxb9588f*)
(*/5h9NeYVA+tWmfcR+IR7VfeTJRqeQ5hHQZaZ+P2BwVlxnQfoip/EZDjMtR8l*)
(*TdfMlaa9/OeHqoWP+8tP8CqtzefUe6KqVyg4XlvNccjtpxuukRZR1mkjVphe*)
(*80gw33kRDDOX4wQYlnsd59DYA0MdYKWzEZrCskPcffFZ6cHyXLlthd3E29NS*)
(*VsJRz17hNil3t4M1CbCtrpys4zl7akv2uOBCoOLPzmp+MLgqOtVVZMAaWs5Y*)
(*kkCrUWkAxYwSYjgsZiz4Mte0viQMpFBvoWDUnqu7JS+HRYUIXi9ybbTCT3BZ*)
(*p8DCW12WfhwLT0g4J95QCWeWM13cxgCz2LSL60hHKyMoTvy9yFimYXcDyxwb*)
(*YzHUHL0wEKa30C7RchjrENFTOVuTCUap6O26JbAv0ZMXNtxt031+6cZd3BtA*)
(*UqFqD8808KuSHMeFgduolnB1/1o7XPLPN70GZozK04NwO4e5AVc0Jyqv8L8d*)
(*7RTmqDfBpHF85+4z4sT1bueo9NJJFMfCrEz386Rl5Pbjp2jbhaUr1+/98JT4*)
(*i9XWv7gUrWP8peDDmiaoWabzsneUUylPXY6/WjXOvl9q51eOWdVeIl0JJ38c*)
(*i4k1jwp6EReTtbtlJjJMiLhAvHT3519DN4TI7YVtNDjOs71tS+GVljrsuoJA*)
(*Z/2c5OLJxlgQElerUWiAmTy1VWGsoxZoMmMRjwzHEhIrAxxIMCqhvgRGH0vj*)
(*VdGZIkcr6kW07iacAjDEWLDINmQI5c44fcbCl2B1B/x5HEXhMVBkF8z3qCuD*)
(*c4wE5Nfi9A0NLN/hMTIWI83RM/g8/WKx2dmIxaXAYB70IygG2wgEW00ywRF1*)
(*Gx7IgsJVJATIvAaQVKjWw73NHLVQnugRJo3qBiBz+rITzCHSSFqm1QsmCmBn*)
(*RWjFsSAkE2MJc6TlBk9SAiPB2+Bd1+jCYzJd4csrtlPYacbNz7o5f+khBCSD*)
(*vtjdMiJMzvisroDV2ogWQjQjti1LgPWMBfSFvKfktnW6vhoIeasFdxowbcdA*)
(*sZMfDofROBfuEWEfsu7SA1GgvA7rlbFYqNWINSBuZ8sOl0Qp8I43mzMWf+zS*)
(*Scti+E5kglc8xQXuprEYUqXIQI4X4aKlhpdJ5YzOFoYYC/adqiE0zlj4/28+*)
(*BWZeQO/4NT3YfMCoA85YxL/shb2tsTEWI83RMfhGw0UQe7Vt+j13NKrZ8/tN*)
(*NtgMbHAY7DYm2OpUBrLEswEkFar0sDDDEE70gLu+vPor/i/o7nnjzA3p71WY*)
(*yXpkLOHw03/eB1TFL+62aAGuKHe3XOX/Xrn0cQ6+JSWDwGfKcjn3vctZe0Nr*)
(*JrJeG5GjpziDTzek4cMEsJyxzDeeIX0rt63+UvBn1cfSVS4S7m4Po3EuiBJu*)
(*JSDhx8VYrNRqZBpY7i0S5nzVPCSRU6ALvGeCixB46dysRkHrzwqhqAR3PfEo*)
(*SQIzFugVAp6Sre+NLMz2ecCp9c3vfj3RtC9NFvkjZSxh0TO0b5fapdacsUii*)
(*dkjhIosNzysJQeqWkt7IiSaKiUmx2TdxUTF4uErcGqCspkoPC5ODklcwfWU4*)
(*l4CMRT4R/bgYC3ee6rt+QFiuXL/3Tx3CEg7PfVIDXGm13I4C9BNBAqOUDleR*)
(*LmeV8oyRFImgjUgB75ifPPchadUqwnLGovKtmm2VD2DKWBJOqxFogJn0ZEK6*)
(*0qq6gRG5SpdbHSoTvxTWMxa9a1vGWKCXQV3XKPIWD7y1Ze7ItiNDnlWDzz4y*)
(*xiKGXthUL2UaY4EBKCqKBNtyKbqJ9wQunJZfH6h3I6Y5r9jMM9JtYpVOEU59*)
(*xakBpGqq9PD4mSruTj/WJz9WCzwagss4HB6+IJz9kZZk2htyovcKzTU4TGYs*)
(*wbqKONrof38/DwJwu7o0Qm4leNDzThbUz4Put3idyE8PQYie+j0XVSgNBsu0*)
(*EcMjQrjjmPVHuMWrNB+2AoZrHn1lkOuBaFuRu1PyrUDjRds6UiUc5i1rldwW*)
(*MFpJ4r9IHsayHrVqVAPMZI1wKFUny2vEjIUPpUhJ21HerKPaiE43kyxT1Eh4*)
(*xqJxdgsCnm6uuxlaHAOJ6zi/WWjmEvu+2POl3IkmZyxC4gARcWQsYuAQ2UmO*)
(*arLx7WFZL4dm0KOVIC9wtOrcsEa6DdVVlCZuOMapAcRqqtxNAuXIOtM0jbV2*)
(*6dvTebIADHgmVxLxsnyp0IXtIWAwxFjQYgpMQXowMi201Qp5PCq/1RcYEWOR*)
(*HA2CfwT6p3VIC9PzN/70xG8v9ILgwKyPW4kPUUIntrjXmb+TSQ0Ga7Qhxvup*)
(*ePf1ARZBWOgbCYZqHltl4LahnGwDoAGMgs+Zu20gHB2zrShUKaewW7wtlrrf*)
(*ysK9gcnFWNapVo1pQC94BUN0kbcwskK/pOHTzSlpRwImPYLIVMaCJZBTJn4T*)
(*k4xFZOKZjnyh0albnPWaGeRYbujpk7K50P2JMUVyMyVjwRwjeL3EJ1PB5ohp*)
(*00aUx6gMMRb0jS395fJmZQI2caPIZn/2IEjfW+928A9XQhdAJ9pSt+Ry2yuC*)
(*GGWqNiPdJj5SgePWvmbPQTyXSnwaIIP23cSMfwhiLbZXlXra23xd3vfrQRqE*)
(*HIdvHLtnJ+vF84YnarvafO3VB4WpRjzdLJY2xFjQNm/Oztraz6PPIIdJbPkI*)
(*5O7IPP6Jp8O0nGlPHy+Mdkm3VdCJ55EgObccxGjLIZ7UHeL1vPdvw6QEg8Kp*)
(*q7zS7sZCrlhRo95BXmu0gYdS4WmYIoFpNY+tMmiIZlWUVHcpcp3N//0wYOOu*)
(*/L/6vvC1SzMRoXO4wlE79pW7/xMvO4YbzhcKzFN2DtdkxoInFNWaQQ0LXOda*)
(*NaAB0RuUsumwV545NAaVYno1Lx8LHiUrP8IRLUxlLJphl4psdAZNPHqUAg6F*)
(*4pdR9jabPT0DAETQKTdnCIxFPDOESyeHJctrEGFzVgdKNsljXCWlGFmyNwGs*)
(*/a8Sn9cRbNplVxThi6W/jOWGM9Rt4iajKOSAq2BzPBsghd7dxMx2l+10YfMG*)
(*N584/tR+VZ6Kc7lH8SxmruSx7lnlrGaMsYQffFO9R3bpaDLIYVi9+l5Zhck5*)
(*0/61OssHr/ReW5TE2qLDQ1cXf9D4OToYTtqM4iCcJQcUUYhO1Hq6AQ+LtIF7*)
(*h6P0tJtX89gqs9p/8pQ8ayKePEey6yUM9cZjMtvKFhuulqdV5Ox1WYss15nJ*)
(*jEU7DJ+cx8KYeV2fWjWgAS2TFJtKOZjNWJhRt2gqTIpKgpna9J4+Z6icmEBO*)
(*M+gSZYfD8sBpnZQWnMtuTBjxOYScGyh32wZkfPmUaMWey/JNFngt6aWmW4rl*)
(*0sUnISoh1iDy5jAPp/qbPcWiMGUpvMXuevlzwJlJL2Bxtk0vQRliEjns8LvR*)
(*7uVqBC5X7GnunwlxvKQmjg2QwcjdxDya6Gk/y2dN8dR2CfmuSXhyb6rvc5+H*)
(*L3m2YXA4qHLh+cE6kIOl+brWdYVLd3nfjymDnFzi+FD32U/NzJnGrM4OjZKO*)
(*BrGkZejGfZ0Q3JXey0LrPANz8i+X/S4u23DOzkYhse28v0h0HjHjoyonhizT*)
(*hrj+iPJBrGbWPMbKMKvjI1/UnifmOuPw5Pbwl/xo97SDm2Ku7QtSZqHQvTEk*)
(*p6VJniJeSJkoDGwhZxGQcP1z4pg3YA017YFiDo/MvK5HrRrQgJZJilWl8cjS*)
(*D1KLcNDOFUCxroHCW2WPakNPmDBrB26tEM+jAWsK0xti/ZMQYW5PSR5+sCfj*)
(*cr7nq/5DRU5Wg59IWSzTBjykSUyrZQRm1jzmyiQmTD/aE+HZ3nWoVYtVGg/G*)
(*glKRSrNuUPy4gLYGpZG3aNmhEzKYeKCMRQ1WMxa1xLYwbQuIGvqwOYrEgFFA*)
(*Rxt83ES9O1eIw4o68BadmdjXFP0z3MyrTGLCGvO6rrVqLWMJzfe9vdXoppXe*)
(*vMSFgNa7HTAdrc6ZJop1jjuNO1A0ci6Iu232OOHokD6oOikA7qbUDSD+SDt3*)
(*UULDrIbMnnsJiEjn8/pZxVgWhr0nuL1r7xXlfBO60c5tm589T3qwL8BaawOL*)
(*m0h9Zo/Gkzd1gQU4Rms3zKtMYgLvXOCJhx2dUdwTgRzUubgQ8J6g93WtVWtU*)
(*GpYGHGme6jE8L4nxNjb7syVd5qW8pUhKMIutR4jPbZY8bzl5II0ns9IHEiPM*)
(*aogsZjFJNbLW2uDjJnSDpgwiNHOphD/gFiVjMbUyCQiNGNCI9KURnkuQs661*)
(*ao1K0Q9SN2w7WDeoSS4Mz0tcvE2xp9nUTCwUSY7Q/Yl+H4pHTurhEbq/gOFu*)
(*1DvxlsOshjAP764DjVBtrGNIO1eCiMiErHOjlrMOkPgqpXciBQUFBQUFBQUF*)
(*BQUFBQUFBQUFBQUFBQUFRRT4D6SaJeA=*)
(*"], "Byte", ColorSpace -> "RGB", ImageResolution -> {191.9986, 191.9986}, Interleaving -> True]*)


(* ::Text:: *)
(*Image[CompressedData["*)
(*1:eJzt3W+MVted4PnSrlbPvill86JeWaXRqjTiBUIrpJGFlhcW88LIshXBZJYo*)
(*lhCe9iYzRnKW7mkGIY9aTAvW8iDtQNIz7MzToqUEJkzjge5lYrvIhkw5bgxb*)
(*bQzTYBoaQ7pC2RRO4RTupOIH49rnuf/POb/z77n3qSqK70fuVqh77znnnnvu*)
(*uff8nnvP/Z9/5//4R9/+74aGhv75/9j9f//on+za8Hu/9092f/1/6v7jf/vO*)
(*P3/pn37nW//7M9/5/W/902/93rrf+e+7f7zb/b8f/Q9DQ73/vQAAAAAAAAAA*)
(*AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA*)
(*AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA*)
(*AAAAAAAAAAAAAAAAAAAAAAA8Nm6/P5545/rcUhcFDvlxGn//dv3E5q6/0+RR*)
(*L8pWeARbU1EnpSaqGgAAAABWjs79u1fPpQOmc1fv3u/US23+3q3rU+bgMR+d*)
(*iUOybOEgBp1ndgwl1rdvBG7R3YGLE3l13JtvvESWPKVKWyLmSPrc1Y8HXBP5*)
(*cRracaZ+Yjfa66OP+s3uKXDuutj6i7IVPOm6Eiv1zrubDUYo1CiOfjIVdVIy*)
(*qtqIA01cvOXrEML2NVDTVVISokODjxL1cyJ1Zi5PdOt9upnTbXl1LH1Impf3*)
(*2lD22d0WG9tTFa1eyiNf+JP/dsfcwNuWykurVK7YNqmt/wjGaAEAAJa3+Ztv*)
(*7n5qpKUOmlojT+0+08+I4c7EHzz9xLBl8Hjr8Ab76HdiZ1KG7acbGWMpYsIj*)
(*UnV0a+Nfnx/gbaiz0paIMJLO2sWbNwcVJFnS8Ejn/J7RdIPhl8bNg335+K7U*)
(*tzeM+NP1JJam+P3Na5OG1sTe5tQojl7EOxPfy/bi62tsVW3GgZIDP/bcgbdn*)
(*xFMzZF8DDaRKSnKbbqzFRWU6vO7l/3xFPpFuH9uU9j+tf/DvrtTKezl2LHF6*)
(*3fG6YUtzznRmJtsvrB7WK3j1C4cvBTfHO8c2ZQ19xxm9lV/avypdtvVUmZ58*)
(*mihtqVuuIy+b15KvHbwwWyYe2ya19R/ZAwsAALAszZ7ZOdbK79vWbt7eGzlt*)
(*e3qs/yGKe1Ca3VOKKbuW1RQcHhGqY/vmJ3sDjIHehvbxoMPAFWVa8/VkPL19*)
(*czEAaW1oX2s+hrWwxOGRuTdezAcym47dsa8XlG5IYo3urZpodsx2fW/CtiP2*)
(*zPMlIxu+vSvtDp7ID/zwxsPCgQ+tuODSD6gXqBy6fN9Kxy8PIkMl0+JEWpuP*)
(*mFtjO8/MmltcOZjHrkb2nG8m72XUsQTrRT2eH6sGF6S96Fw7vLFon088vW1X*)
(*2WdH7XURHzHqvKhEpXVXT7Q01pifdEVbuvH6ltFKsbY/V1xbNhwuClZELKuJ*)
(*pP9bPAmK9dMw7aN4YAEAAJar2VNbs58pN+yfVH4a7ty//uaZfkYMj3J4pKiO*)
(*odEdp6vV0a2NyUv1Bn5Oy3EUU5SpPCDzVw5tyNrLi28M4mGaJX65ZvaDEwd2*)
(*7Tpw4gNhyBqdrj+xAYZH/Cn6wyOV3asc+E3HhGfKwiouuPSDD48s5nkmnEid*)
(*mdM7sgduRvecNwNOnamftffu2tv+2VTNKORy7FjC3Bl/uQxUjwzb9qJzaX8W*)
(*Shp98YT6VNv89Lunz0d02kVdrdp/SSlJHjdRY3/VE8160nUP5FuVC2tZWiHu*)
(*FZSgUd5H7sACAAAsX+/tS+/RW1tPBY1rHC9RF29Ef3979iDyqu3fN9+Qzl6g*)
(*2XD4lpl6ekv45KFroXkqK1WmCjFf7BdGfJ2pv8wL9+6Hv07+VPxoO/LKWf+w*)
(*xJOnOs3K/L10J9R5TAIrzTLLqDCPaUimFeXL+sJ0M8KorrvTZ3Zkw+SdE33l*)
(*6Z4fID9O6QBlfjpZdeLilDynxfy9j/PpcqTpMYodSB6J79y/bilY7HSlziFn*)
(*QGLFcXvt2SwW8Oxr5urFWn9pDJFnP5hQG65ZhQ2HR3rlOfJM+uc1B/P3PaLn*)
(*ebWey2FVUl9ctCCo69FX1b1z/aJ0Ii10zr6SvqNViTTGzv5rLWBwxxKyr2li*)
(*vW3SEzKbYiaZHkWeH8VdcdlepoXo3J+ynORJb9wae749OdOxhrfvvL4l7Y/W*)
(*7L8UFElydXlF/69kM3dqq5h5bDQjde3Qk1KspZ8ECY8AAAA06/yekeBoQGfm*)
(*7X3r1Le7u7euR4qX5x1vUFdvLbPV8n8md8rZzXp2G1qNYHjyTFP8kfF2t/l+*)
(*t3F7PXum+PX2xVPZ+LN4w1z79dAUkmdeIesP/ezUS9UnxIc3FqsFVpplWCcM*)
(*cPP1Nhy+dOHgxkrVtdbuOatEwOavHFGfW+9NhrCvOrmEGB4R/hqap/mofDr6*)
(*EY/TocnTOyvrtsZeOlWNEsyd/6PNxlwDQ8Orv1Ndqyzqm1fa1YJVjoDlKDjH*)
(*Jc4xdkBijkkLKqsXb6wYZ2ceqJCDmoMKj5Q/ohexzYiKk2aHGF79wus3Yqqk*)
(*vuDwSFDXU+yZ3qyrW+2ckE+ksvstI43hs/966jO4YwnZ12wY/mo77zKHWi+e*)
(*utDekD/b8Q8OViITARWX7+WOE1p/oZ/kl36WP3hhC48UbfKZI/4Imr/LK64A*)
(*ldhFERDW8+4vPKLGbGslSHgEAACgUfk8qcKdmq6IJrRGnuq9Rb0tjw4UL1Fr*)
(*b0QPaa/353MgqOGR9HY0mwtPi5wE5JlIbiSHn1iVvnRerqW8/6HeXs+e3bM2*)
(*T6qcRKO41fbf5obkmd8Ht1q9ZcOrN2/ftS0fN+SrBVZafHhkZHS0lWX6zbX5*)
(*+0L73stX61zLBzfJKpVpRUZ3Tug7ID89UoaQwvIs41GtseeULKvHIN+lpNKS*)
(*416+r18+tZDl2RoZezKdLKcyLUplrWIHktTkI6AcBcd0pVXOMXZAYsUcr8Ua*)
(*xaQF1YkLyscL1Lki89+4LUHNgYVHyjlJ8vcCgiuuPPjFNAzp/Bv56oFVUlvo*)
(*e1FhXY+6ZjJTUTlPy9//2u9m57AtPGJOaRE6+6+3PoM7lpB9TcvZO4laa79Z*)
(*nNq98/ibaeplfCeo4kJP8gpbk8weRvQHtMO6vCJiVbkk5n8ysugrPFIEPYUd*)
(*JTwCAACwpIrfyoq7MO3h7vyZ9nygVh3Izp99ZVQdLGXcI5AsJJMuywuQJpFt*)
(*mBUmOM8b7/5U/Z5o58Kr+n5Vb6/L++Sh0R1npIcI/PfaIXlWfsAtJjLpXDuU*)
(*/VH7tdM3bIsNj1R3bvbEN7Rt88cPKvvf+fmR7GMZxVhNmjKhWKsyNg/Jsxzq*)
(*bz1R/FpbDvKK2ih/Pm9tOJT94jx7ervxkNOdSz/THuo39zLuCNhqNPRgmPyJ*)
(*+dbIAyHKRC/5SWMZSA4wPGIb5vv2pBgSlsc0MT/97gW9Ahdr7hHlXRPt/Z3w*)
(*7u7G4aQrUSduypt10cCs9WafEcn9/lZEffraasi+FvHP7adntVdNsj3InicK*)
(*rLjQk9xfU+Z7L7av3oZ1eZUnevKQZH66mReFPsIjIb0P4REAAIAlYg5FtIe7*)
(*87/nv9E9ufv1yq1nPk2Adn/muSFPs0julZM1V61ald17quGRqDzViUCK9+3L*)
(*GU6K2+vyCXHj6ytFwfXZT2zceZYjsVcvlPkUr56rny9uPDyy/lC5c8Vv0tmP*)
(*8PlTQ1/95h9VqrbYgzy1Iq1sGoiT7ZfzZy+UuFJInpZfeYvayH+sLXbpGyfK*)
(*DDqnt1sOTGXykaJplD9lRx0BW42GHgxT/fBI+RpNER8p4kzWtwmWXXikeCor*)
(*ZCC3WOERQyW/4K4nHUub7zhlhy3vChoOj0TVp6+thuxrlkZ27maFTs/Y7B9Z*)
(*4oEVF3mSO2rK3DvLV28Du7zK03F6aKj+TKqVR2vkD38RHgEAAFhS5m17NqjN*)
(*n8jO/h71IntYeKSXdHKbv+nIsd79aG9dJTwSnmdn6vRucyIQbS3lee6er3zr*)
(*TX3ehjJL/6+AAXnK9VCMAtTaaTw84tqFwIkeLIdgeN3u01NiaMeeZ7GKHpLQ*)
(*hz3y6LiIaJR/nb1gTL2glz/uCDiyt+7NwMMjxYMCeXxE/7c10UGERxxvMrj2*)
(*JPwFCE9CDbB/2Ld8fyd2QiDh1cR0L/I3y2wniOPlRlcTi6pPT1sN2lf1rUfl*)
(*+2JKeCS04oJPcr1C7e97FX/P3yrKXtIy+hVZJUv17cE8FKW931ZJMzSaUb7Q*)
(*qT6z2HeChEcAAAAaVjyarP9gl9/n6uERZT4A8UX2Be8NeXZzv+NMcufZvQdN*)
(*S7Hp2J30jjDbKjTP28c2FV9/XJtORCG9uG8+PWL+gleMm+Xv6pTC8pTrwfLT*)
(*75KER4whojpQNA7B3vZJ6Xs0MeER3+/n8sipmBEgf/Hq/B5t6oVd0vQXUUfA*)
(*kb11bwYfHimfFkmG0Hl0xDWX8uCnZrXGA+R882VBD2Utgw/7Bnd3SqCgKutJ*)
(*tOmorVOzCnEOVzmj6jM0POLa1+jwiKfiwk5ycZ/1vSjeD7XEXfV+xdPl9VRf*)
(*2MkvCVJ0JCaaUb7QqX0wvt8EE4RHAAAAmlV8q9PyAEh2d1aMi8xvEQqK30Sd*)
(*N+Q7TvXuQpNHltPkt546VQ2PBOZZZlZ5sUMYElTnHinG1vrPeOXcrM7vIATm*)
(*KQ9MzK9VhFSamFgZzokMjxQfLDIfGBfz9IxUA9azBeKMzziIIyfjkyn57+fJ*)
(*ZAg5c9OoI2DP3r7DdcMj+X44s8snIOm2ybwinPkOKjxSdBbSGena12IIG/J1*)
(*kaAqSczfPFd8YzZYyKEL7u6y5x2Ml7SykFYRw5JPkCLwJZ2GrnJG1aenYwna*)
(*19DwSGjFhZ3k4ib6XpSVuOPMnLBBvn5gl5coOqtNx360pxKblMvkj2Z0rh3O*)
(*vs9jfW4kLsEc4REAAICGVSb5qz5LoYVHihHa0JpXL2gPD8zfvKnfo+e3svLD*)
(*31namzZt0j5Z88wzz1RvaMPyFN7aKJ9ilsMj3fvVS/uztLVHSIrqMJ8t6dyf*)
(*npmLylMa4cye2prVt/5yhLvSzPhCZ+Z08fGK2PBIedQr86Rmuzlzc8r55RpB*)
(*yHrF2HrN/vIjoMVxKGtDGDmVR6uYjNT8/bxz7diWUX3TuCNgyd61w3XDI0Z0*)
(*SJTX3ZodO56xllzPtuHwSPnztz4O9SVWGcIaw8POzIz1DTf317TKGKdlHgeL*)
(*oEMX2t1l72G0Nh35eaUEsxdeXa+0VvEEKaclFkviKmdUffo6lpB9DQ2PhFZc*)
(*2Eku7oVZHUWeamVo4ZHALi9VRGqSebGsbTEsmuGbcCQ6wRLhEQAAgMZVvrc6*)
(*8tS2A0d7M9Yd3aPOPaKt9XL7ZLLWgfRjksYtXHFr3xp7/uBJ/RMCxc+ZldvO*)
(*4vfQ6DzLgMbolm5eldlD1Ztp7fa6HO1pN62VL2a2xp7bneTay3T1cOXjCGF5*)
(*VqY5+P1eOpX1jM+D+iqtsnx44/95cvzo7mqWseGRyu73vnK572g68ere5KOk*)
(*wodfGgiPLMxN7Cy/69vL8WT7nwkfVy5GTmu+nrTFo/uKT35Wvr9ZDIlaa3ce*)
(*VVZSihF4BCofu8hnkcxno60egvKrTsWMjuXnT4rvngQmlh2HfJ6D3oQuR8f1*)
(*tIoCFg0u2RnHizVlFVqORVDxihF1unvdais+zzq2szIGDd7X8m207o6mnczJ*)
(*9u7eSWUUM7BKqtNciG8+WIRFtkK7u2JQn+1V7xxKW5jSrehzHB/dV3xUdnjj*)
(*4Ur3E9bEourT27EE7GtweCSw4sJOcuUjNHn7EqpD7czSLvtk+3fUuUcCuzy9*)
(*zlKWh2ECohmVXNf+/rFxjd6qCY8AAAAsA71pRteJk1wO7/6LylqnXqqOQYeK*)
(*u+B/f1lPce78vrX6qmakonrbWbx3rs6SF5Bn+eRysXjspf27nrRkKt0uqz87*)
(*WqqjNfLCn92JytMyUWF3lHLkijF/h6/SevGF3WpdjG45dmi7XmehIY2F2QsH*)
(*NwoHfXj1747nu9loeCR5wON582gObzx4wfzR11hr3b63Z5Qo1k4tqW46r35D*)
(*K0bgEXDOKCm1W4EvTz2xIu/Xi4dezLTKqiueF7B/z1erQsuxCCqevKfDq1/Q*)
(*Gm7Evs5fOSIc/NbYv37PLGJQlVTimDEDxMAHf4K7O2k94wyX66mbltqoQ5tY*)
(*T0R9+joW/76Gh0fCKi7wJHfP9FqtjtkLbaE2etlu/MGtylreLk/M2vYOU0A0*)
(*wz0jrDwRE+ERAACApde5f/3cyfbebKa6A0fHz129a77VPz99ceLogco65jyd*)
(*YnrKhIbFJ1//zU+K287OpR9Ks+SF5dnLK12jN3fozfneiOCwJdPqNLI3fpSX*)
(*718eff8zI8nxLNduohMXp7VMA/IsPyu7dc/erPgTF6cccyU4Ki1dPnP5rWTx*)
(*3vZbl3sDiXyvyjrLP9xgVKMvv2Qv1cKFphWRZ+f+VHE0exUnzBwxf6/3seT8*)
(*iCeHXJxfotsu0sroVWtSG/nxrBSjc//u1XPjRZ2Kh7Isv0Rot5Iiz8DE1B0p*)
(*Gri1Gi+8+vfS8d6mY775JpwDq6DiaXvarTX5GMTua29P82PROxeMIxFVJWXM*)
(*KGhKJK3I4qFwFMPR3ZV9gaVN6/XUTUve+bAmppQvqD59HYtnX7PyZxtlhUyL*)
(*k/1DS89TcUXQ4NkdZbWZletsXmZ1KI0mOdNvmcfL0+UVymtD9TKlCohmOA+o*)
(*sQuERwAAALBihf5QDfjkbxP5XqzpCZ175FFXzCITMkMplpEBf51oscRGMxpP*)
(*kPAIAAAAHhmER9CI4j2wgEdHFlZ8eCSZkeLogWyOj8iZWbEcEB5pJkHCIwAA*)
(*AHhkEB5BLelEncUsvKGBgHRgNfzEqsQ3f3hr0OVcXJVpIczZUPAoWFHhkfRE*)
(*e2K4etKt2vsX/s0Lt374zXSraiLp/xYrqFh/LJkwmcsLAAAAHgGER1BLdV5H*)
(*bQbb0M1WYONLZqRwT3yE5W1FhUdEUTvmmIJWTEdbf8Wd4QAAAFiJomaBBHTZ*)
(*vI5722/ZJ90UzN/7uEKYYBlYUsWEpSHzRy9f6ommiIrcde7fjUpHW58zHAAA*)
(*AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA*)
(*AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA*)
(*AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA*)
(*AAAAAACAFWF++uLExMXp+aUuR7jO/evnxs9dv99Z6oIAAAAAAICV4M7rW1pD*)
(*Xa0X35hb6rIEem/faK/EQ2sOXlnqogAAAAAAAF3n/t2r58Z7zl29a3m4Yf7e*)
(*rYsf3F6U4oRkdevwhiTWMLTp2J1FKZRk7vo7Sa395VTH/PP7+h5M7EwCOkMj*)
(*e87XyLN3rC5ff1RCQj2PXokBAAAAAI+b2QvtF1YPD1UNr37h9Rvaamd29Jbs*)
(*OLMIJQrMqjP1s/beve2fTS3hqyo32uuHhEdY0j+bezD7wYkDuw6c+GC2bo7r*)
(*2/rxWb4evRIDAAAAAB4znWvtDckTDcNPPL1tV8/2zWtHWsLQfvmFR5aDdOjf*)
(*arXU4b8tPNJYjo9SsOHRKzEAAAAA4PHSObOjFxxZ/90ryvym89PvXsgGs7ff*)
(*H0+99mwvZvHsa+MV7+gvTLjf0UnSyt44mb+XrDhxcapYLTSrYj1LIdR3W3rz*)
(*t2oZKbt671ayXGW8FmOTDv23v/LKyNDIK2c76p+L8Ej+Co4j+WSn0l1J5m9N*)
(*qvBeeVSKJL6/fVU36VXbv+9MMKve7m7fuidPXVseKrEyPcdKSkc75I2XGAAA*)
(*AACAwcjm79h+2vp+Svokh0X1iYAbP3r5qZGW+o7OxoMXZvW0tp6anTr10li5*)
(*Zmts55nZmKy09YTHErLoxGlLRpnO1KnvaG8V5YKf+8jjIFcOrqm+YKOFR/JX*)
(*cBzJJzu1vn3hwsGNlSINbzx8rSMmYS9vZ+r07nXV3WqNvXRKfQFp9sLBr2nH*)
(*Sq9Mz7FKy+Q85I2WGAAAAACAAeoN67tGX7QNRy8f35X6erLimq/vqvjeRDkp*)
(*ajKeHn7iyc3b01d00sBD9ZGKdJWR0dHeiHp4dW/FbemoOJ2qNDSrYr10RWt4*)
(*ZHh4OM+oKE45J2rn0v4km9bIU9sqBU6zPX45sP6KOEhSkc8cua39OXVn4ntZ*)
(*ib+9YWTIHh5pDQ+38hJtyyIP2cSzRRJpCiMbvl2tnbK8s2d2jJZ7VezW6M6J*)
(*4hmb2VNb0/epkkNQZJQkWVSz51gFHPIGSwwAAAAAwGDNXTqYTj4yNLzu5fZb*)
(*F6dtLzb4JgS58e5PryvvXqSj8OpHWrLHPobX7Xt7pnij5sgzRogjeO6R7JEL*)
(*S3ikO+Reu+ds9ixDFg0pVz6/pztob209UZRkYX7yD1fFfie4jIMku5t/sdc6*)
(*94h1Qf5IzOiWY9nzInkYQ1vZPZNHEqUZ3XG63Kt8epni+z53jm3qRYFevVAc*)
(*6c7Mn35jSPvccMixCjjkjZQYAAAAAIDB68xMHinfkhhevXnfmzfNIElQzEKZ*)
(*ySOdc6KyRZLEqlcvKM+pJM+CVB9DaSw8MrL9tPFqT5lq8pFd5eGWhUv7Vw1F*)
(*jsor4Y7O2d4EJDvOzC30Hx5Zf+haR19Z2z9nsCHZA22WD/1ApC9UFQ+69Myd*)
(*2jqkhTVCj5XnkDdSYgAAAAAAFs389MUT+57LJpoYffENbe5MX8zCNp+FHh4J*)
(*+IJJQ+ERbYGW6twbL/ZK2xp7bt/R7oD8ZPvl9OURJW7gpYQ7yhds+g2PqAvi*)
(*wyOu6VuKtLMXqnpPC53s7vnR7KCroaKQYxVwyJspMQAAAAAAi2z+ypEtyVwQ*)
(*W15XnhJwxyxuH9ukfh9YmmdjWYVHFjpv/4uvaOPx1tjzx65FTQmqhjuSt096*)
(*76gsaXhEm+XDnO3j53/yD/VIxPC6fWdnjZScxyrokDdUYgAAAAAAFt21Q0+a*)
(*Y3VnzCJ9YaO19VR1iG2EApZVeCT51symP3n/4lvtvb2x+IGj566Ln/510r/g*)
(*23skZeSVs9eWKDySvDG0av8lV5GTMo6+8tOr544e6O34XnHGGe+xCjvkjZQY*)
(*AAAAAICBmp0WpmLtXHi1N/ODNsNmMpWp7dWTdARcHeV2Zt7es7bVX3jEmVVV*)
(*jfBIssqTh675svDQ4wGd83tGh0ZeOXJoYOGRdGZV9U0YfeH2N2fsgZ4kJ8en*)
(*nCtruY5V2CFvpMQAAAAAAAxUdxCcfFD1wNFsTsyjB7anM1G0NhxWh8blXB27*)
(*kykret65nn3kJZ3bc6i19p/1lh09sK2Y6bWf8Igzq4Xb7xczeL727FB1Xs/3*)
(*83hKQHgke/hh7Lntvecn8lzOXb0b9wSJGe7ovWDTWrOmOrXo3PV31DlHn33N*)
(*2Kvg8Eg5dci24qiVez43sXN0SK87NafkeY3e5skjM9niiYu37imBMu+xCjvk*)
(*jZQYAAAAAICBeu/frB7Wp6GQZqLo6Vw7vFFbuRw+Z99iLbXGnj/y3W/1Fx5x*)
(*Z2WfzLPIKuTlmtmzr6wRUxneeDh8/hHhaZAsuDOkPaki0vcqJDyyMHtm55g2*)
(*IWplu861Y8ncMdacOteObBqRCtMa23mmOOr+YxV0yBspMQAAAAAAg9a5P3Vx*)
(*4mQ6/0bvOYqJi1P2Byjmb57LV92lf+O1M3NZm8Yj+Q5sZX5N6buwVtasklRE*)
(*RVZ3Jr5n5pRul6+TjtizZx/G01k4duVPQGw4fCuojHlO6hyinUs/VMuTriTR*)
(*90pNSNyNNIuZyxNZmaVJTDv3ryuVV00lCz2NPp19sCdZa/vmNEzW2jlRLY/v*)
(*WPkPeRMlBgAAAAAAA5FOb/IP/+Tn2t/TWWlX8FSh6VwfX/kXb2shsM7p7b34*)
(*yNZTvNECAAAAAMBjIg2P9KbDKCYeOdnemz09Mrrn/IqdJzQNjySTrlQnnEmf*)
(*HmltOuadDhcAAAAAAKwUt994SZ8OI59A41j41COPoLlL/06f1yUlTjgDAAAA*)
(*AABWtvnpctaVZN6Vc1fvmd85XoF6M32MlzOBHDjqnHAGAAAAAAAAAAAAAAAA*)
(*AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA*)
(*AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA*)
(*AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA*)
(*AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAWMbm5+fn5uZ+*)
(*+ctffhyju353q+62S118AAAAAACA/n3xxRf37t2bzt2JUWw1OzvbTWepdwXQ*)
(*3Ppx+7XE8UnfqpPHu6u1f3yrsqF/I1/WRXpAPWr7TP7ZfOui0caK6GGQSets*)
(*IBUW0G9rJ9IiWZpcl0fuALAU6t/ICxajO03y4KZiKf32t7/96KOP0vjG/Pz8*)
(*w4cPozbvrt/dqrttN4VuOt3UHCt/snPow68Offh0+zNh4ZlfdBftPJP/88bH*)
(*T9vWfKwlp4ztxNTOJ0Y6cTVAeATLGeGR5Wdxq6uIxDzq0RjCI4vtsQqPpHdJ*)
(*tkbgXtpgvsrpKld+dj9nrKBuqpALLm1QXdO33EzL3lik5dUdsSTt7L285Quo*)
(*TF/5bZUdsn2N+uunvXGf3xjCI+jX3bt3p6en/+7v/q5mOt0Uuul0U3Osk4VH*)
(*vjr08z++YSysHx7ppSClvOh+suPDr67/+G8GknbRw0unZv1uczBDrj40UpBa*)
(*143Bh0eWTWU/8h7LmiQ8svSE2hpUW9QiCo/VcDfGQO6G4624LumR2SHjaS69*)
(*5Pk9VMN7o+WrdQ3CGFldw9/xOtfwNfuI08I5NJeXd1OvpG0Nn6gRA2UFd/n8*)
(*lekrf0Rl27fvr/76aW/c5zcm5Mgtk4uGivDI0vrVr341PT3d/f+Lk1oSHln/*)
(*86e7/3/HJ/pCwiNB0utKW/7tjW7TTKTv/oXwyKPjsaxJwiNLT+hhFu2wcPck*)
(*WiZ3uiuuS3pUdsg8/mLJG98dPd9bk5O3zBXKLM0z2tP1us933/4E769zYB6w*)
(*PKE/HWbumr437vL5KtNXvojKtuxf3fqLbG/c5zeG8Aj68lGiwQTTV2xsS5Pw*)
(*yI5P/qb9cyUSkiI8EiQ7ZeSHk+k2zUQIjzwGHsuaJDyy9BYrPCIelWV5T7fk*)
(*lkmtrLgu6VHZoeUSHjEpnYW0uqvvrXszEbi/6WqTzuCDa3l1vcrumf2kvj+x*)
(*x8NybyeXL7yyrftXt/76CY9wn98EwiOI9+DBg+np6U8//bTBNLupddPsdDri*)
(*0iw8srDw2R+v//CrQ7/4SXWhOzyS/LO3NFkt+698BKV4bSf/rxKdSKMx5t+r*)
(*KX+1fOunKKSSb55IVuZeAKSSXVnsavH0iVbSvTb/nrl2+I9Wr/pXf/D/umu4*)
(*OGWkjtPdbRZdQPUVTO0KVqVfXoS/Wxar5cjzzdfKFlhfVq1REK0qhPRdr8iq*)
(*faTUYypbW6L65WLHZcu1jyGv8Rq0V2Mdbx1bnhkVG4W8B9KdjVza/OJnPK7q*)
(*bjHiKpYdb7AmxbqQm4jWlp2V4HltWj5yUs1pLTJfJfLkdFags9FKxRXWdDRF*)
(*7YBkG0u3SOrfLO3If3iN/CyTXtjvRKz7K/Qw/fVc9jZlrmOW2vewfVxLFevM*)
(*lpX6N0uGlrZfblu3j7HteMBZrP5BPjp9VJelGVjHaLacxTNcOgghNRPQeVjP*)
(*emcfK2Qh32CIqbkusMF9r31H7SVf6vCINWYjp+EbKgV12979VVqF5VrhWq5m*)
(*p79u43p6JHqIL1aIrXyhlW3fv9r11194hPt8a29jVIitU/FcELQuuLLI1Tn3*)
(*1Z3Wql4sql//+teNzDpS9emnn3700Uezs7Pi0krkIQ1KVF+xCQiPPL3+5+Um*)
(*yfpKkEF6eiSJYxR/NMIyWiJFqEQNj3TzVTLq/tEIv1QfhpGeHklfLMr/aO5+*)
(*dHhE6jgDus3e83rFGukTfPpFS7wBK/5o5KolUpzj2hi7S+tejN5Q35XIgpiM*)
(*/sWdqzs8oq2t152eWN4XR/7gk25mr16JcEjUHt63z/ZGId0VKBXjrpVkD6uJ*)
(*SyuZLSb6UDdTk0Vd6GUzq0uIMdkqwdjb4452Xxw5qeakS3Lb+JMZH7LXY3yj*)
(*7a7hqlJHUxTOmONR4RG9HYV0IWZ+UltyxzFd++u4xdb/Zj0McptS2X/qdI4S*)
(*41qq8IeIYyRlaG/7lW1r9jG2Hfedxbb0qido39UlVZh8QXE14bzZO1uOcvrG*)
(*hUccSZv847uAGwx3hQZcnj3nibmjtpIvfnhEaefy2raOyD849+XuD974Q5a+*)
(*5cqKYiBQvSg6rm1eUgms5QusbNf+1a6/fsMj3Od7OibvfV7ABUE8vM5Lf3x3*)
(*Wrd6saju378/PT09Pz+v/PUf/+Psv77+ODc3103zk0+UcX9BeTBDjyoEhEfU*)
(*mIMY61DDI0LApFeGPFn1QZFKMbTwiO9NGT0dMzwiBEz6fBVI6RnF0Y272xQ6*)
(*feOGzNetV/8kXQm1DkvI17Nf/RQkJFHPCq7wiNB/Vv9kvcWPC4+IBe7jXsS5*)
(*pdCE7I3CfW321Ip56bKUV1sv+lA3U5PiAVMrQFrFXQmOw+NqnlLNCZdkuSar*)
(*J6ejHvtqtCrPTa9tRVVYeCTo7lQ9UNaIh3BrHrbL+qpB4ZGQw+DL3zYccIdH*)
(*olpqzWMkZOhp+2oF9NvH6ELPYqWPEhOrUV3+8EhIwEs4wNXzO6LdhnQenrBT*)
(*YHjEdYMReyiFy7N/f3urBYxMA9qQxB1/chROCgyFhkcCDnQe2M6JQ07HciF4*)
(*I7Re/5VDTtxYw3IH4CifkJIeLLeWL6SyPfvXTP2Fj3W5z3f/yZZeWRItEuG8*)
(*IIQcHSHRiO60fvViUQ0iPNJNrZtmN2UxRy2MoMY3Qp4eUV9ISQIsrvCI9BRH*)
(*kmlaBjlAIbxc45sCJUmzkpGRbzUmo+ygPgGLn9Qz2u70xKiy0Eu4eiupZ9IG*)
(*f2bHJdzW+O5qpHFHTEFE3vCIq8tT/iXmpP1aEPxbkHWXLJcBZzq+IamwyBh8*)
(*OBuFloFn5OLt9QNaTPyhbqQmLUv7uu/RmoZYcmfrtP4AZ71JMv/qqcf+Gq2j*)
(*nH1GR0LDI74bhWoeMadM1I2r42yw7Y2vOQfWuHkzapu5z5ZqSEvt/xjZMpQL*)
(*6Dy+UX2MLvQsto091E36rC5veCSkb7beNCvda+C5GtB5xIagxSwc15KahzLs*)
(*PBEysV9ynBmLrP2EP86j5BYTHgnq/Tz5eZa77xz9y4XdEMIL7rhyaPmlhe7y*)
(*+Ss7bv/i628h5PwxcuA+Xy69e209l4ALQtA9gBmxCu9OG6heLKolD4+o75gE*)
(*zj1SkYRHKvENPdyhzPVhzkyiR1fEQjqCGNo0I47wiDK9iW1mkkDGKVPtOAO6*)
(*TWnE5eit7PcJyVqWXkXoNuWzXIvJ91+QwMpy5+oIj+g/Hmj3S3L3FhsesVZU*)
(*HyEQKX2pXAGNwh668NSKd4wg72Efh7qJmgwvm3nz7a0E5Q+eMop7JGzha3S+*)
(*euyv0RbrGDvaT1O0LjSHz86bZb2J+Ed63oFqhby/RkpyWb3NOTw8U02qt7F3*)
(*hBfVUmseI+vgV2j7Rno1+piAHddzEH/6M5KsU13eri+kb7ZWuhbaqbQlu4DO*)
(*wz9kDQiPOK4lgYcy7PJsK6KwC7aSW1aP50woW+gZr1b+GndYnGXy3RvZAn2+*)
(*uJT/IqGuIe1CQHjREikSTyxX+TyV3c/+xdSfVvigQ8l9vlgQTdB9XkjhrQnZ*)
(*Lv1x3WkT1YtFJYdH6okMj1RfsRlQeMT8grCyeT/hEW2aVv/TI30+KCISTpmy*)
(*4xxQt2k/Q/vuNvP+whp6jiyITK8sT66e8EjUvbC0T76tlml4pJp98KjVtkJ/*)
(*IQifpQ6P+H8VMgcyixMecdRjP43WuD3V7wEWLzziPJk9B6VSVP9hsO2v8E8p*)
(*a29z7qtr824Y31IHER4plwj1WE2v/z5GzC2+hzELWae6Fik8UvmTei6Ylkl4*)
(*xJVC1OXZVkTLsNRa182FR4Qs5NBIkbc8gBYCif2U0bddsdwetUprx7dcTt53*)
(*k+E9nkb5bZXpL5+7svvbv/D6U3e43/AI9/nh5WokPOK+9PcRHqlZvVhU6dSs*)
(*n30W+/yCSze1bprdlMWl0lwfxZylAwqP2KcN0TcXCynka+5FUHgk/kERkXjK*)
(*5JHsAXWb9lsay2Jvt2nuRc2CyLxjmajwSMylPKjMRp9quY9zpeNdJhQ6aKQg*)
(*XrN6W1ZW9vbfljGCu8X0caibqEnbUuWvQW3ZU1Drb076igHhETPj8q8hzS+u*)
(*0QrL9HsAy5bOShL2VUtLWMN9MvsOitKiYxpEn+ER76iwj7sg93bxLbXmMfLu*)
(*RqXtC+n13cfo+juLhULWqS5veCSkb46LYZSjJ9EyCY9EVWhseMRSTEcc2JJc*)
(*GSuLGi0Lx9taZdbKaig6Eh4esRXOF14OaAzOAx8VHnFXpq98oZVt295fPv/2*)
(*cd089/lBQu7zAi4IZtl9l/4+wiM1qxeLqtPpLNWHfVXp52N2/OLpBsIjyjrm*)
(*N2UU1vlMnOERcytj+tbkcznV51KkDxn3yXLKpB3n8ePtBrpNs5vw3HRJ1wFn*)
(*t2luZVz9YgsiErpvV66O8Ih6b+/LqfK3uPGpeHQ9faujEw34iSooPJKt9WNj*)
(*XU+tBAyqyr9pZYo61M3UpHgLpu6i9QY29IpeLZfr8hcaHrH+wqT/w1MWNdHg*)
(*IIdahuimaF2m1WlImE09ekE3+OkvozGnll7nlvBIVM/VZ3jEM1zto6XWO0b+*)
(*3dDbvhQ9iO9jxJLGn8VCIWtUl3Vk4LgamVvFxTACYr91wyPu4+u9ljgrNPLy*)
(*bGXumf1aYEutj/CInm8fUdFlEx1pIDyiriFlFnNbE/uLibF+WGWHls0oX8D2*)
(*kb089/lhAu7zAi4IZkl9l/7I7rSB6sVi+/jjj2dmZr788stGUuum002tm6Zt*)
(*BUt4JIsnfFj36REh3JHGJbRHSsp/pvlqX8/xhEf0WEeyU2p4RHguxfwCjvp1*)
(*4IU+PuwrLFDPpz66TVu/rl9r1JGLOWhzd5t6b5FvU6sgAinK7sjVFR6RRqOT*)
(*x7W19cp2jDQt+2jmEnB/pG/TTVeNNDjqNSw8Uu6O+IOGtVbEIZy/xUQf6mZq*)
(*0txJufqMNFyVMHncUfv2IxccHtGjCEo5PPUY3WiFUaZ0LMOaYlkz0jjI2DHx*)
(*7sHRhTjyU7ZwtYiw/RX6tpieK+zGWfkgdJ9jVd/pWusYCRk62r71iEb3MeKO*)
(*R5zFrhO0/+qyDQ1cJ6d4j+646a6cXAveJlE7PBI4JHZfS1wVGnd5jiiFWPLe*)
(*HyMG3H4BFaxRe5fIgL5aW2pTEC8EzuVS4hHXW61rFQakemTMHJ66yhdQmb7y*)
(*B1S2Y/v69ae1N18UgPv80Js/731ewAVBqG/PpT+6O61dvVhs6Xd4f/WrXzWS*)
(*WjedbmrdNG0rWMMjRZChVniknBXEfJbDOhtquUm2VcjUrHlIZCiNvZhv8ZSz*)
(*wlayq271ofleT83wSHFG1eo2Kx2BOYwr2O5r82VCx2C5US+2aaQgYo0IIWBL*)
(*ru7wiL6j5nFQCtddFnhDaeyJupNB/aNWMvPu0HnsAn7XEa4olrz1K5W0jafF*)
(*CLXgO9RN1GReF/pxFFYJrwRtN4xtLUcuMDySt2FHLXnqMbLRKrnldeVobvbT*)
(*T/h5r/J3df/lduTtQhz5LQj3VX3sr31Mo+doPwyh4RF7tVrXj2qpYkHDj5EY*)
(*HrHmZG3f0X2MbcdDz2LPCdpEdaV/tXa14ta2SlLDI7a87dUSlrQ9CUdWsTcY*)
(*cquo7I7r8hy6o7aSx464vYS7B5G1MqTS2A+JER5Rc3G2M/dhXpCOmme5dhIF*)
(*dDzSkNiWQGBl+srvq2zH9vXrT2tvfYdHuM+3V4l8bAIuCGoaSri2/Juyt311*)
(*p7WqF0vg7t27jslUw6UTvXZTa6RUS8gRw0EgzmvEWRYtJvD+GyvDsmhzAJoU*)
(*2IsPODwCKBpvb8sAV1CsbJ9//vmdO3emp6c/+eST3/zmNw8ePIjavLt+d6vu*)
(*tt0Uuul0UxtQORdLk9OoPq78T00AVcujxXCD+zhZiferwOOO8AiWoRV4uVke*)
(*d23AID18+DCdUjX1UYxiq24K3XSWeldqs3zOBhHoNRFnmbQYbnAfIyvwdhVA*)
(*/pC66+ROn3Jv+PQPyBePp4G0t6W2TO7agMH7/PPPP/vss3v37t2N0V2/u9Uj*)
(*+tCIMk3rQj4PCY+OxNBnTXK8Pw4sLOcWQ3jksUFwBFipjNkZohYPLF88nlZC*)
(*s1i+d20ABkGbtbWhb+8+XvRJuR7tqwAGb9m2GMIjj4GQD9YAAIDUsr1rAwAA*)
(*AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA*)
(*AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA*)
(*AAAAAAAAAAAAWNE++c2XfzH9xX/66wd//N863f++f6Vz6sMHH/zyi87DpS4Z*)
(*AAAAAADAgN3/vBcYOXSx83/95efmf8f+uvPBLCESAAAAAACwYt381cM/+5sH*)
(*YmCk+O+7Fz7/6dQXv/zNl0tdWAAAAAAAgIb97dzDw38lPzRi/nfqRmf+ARES*)
(*YOW49eP2a6+9dnxyqcsBoIbJ45zH8KPDB5ZCeuYN/uRL8mk0j+Ta0v7xrQaT*)
(*XBzNV4XfY3shXmFXlvkHX/6XD+XnRv78xoPffrHw7kdfaH+fvPPFUpcaDeud*)
(*zu42XfTrFZYtrIklnUZG6GiVLIQE6i1Xd0Dq533p67ugr+Tc3lV/0jIpGVv9*)
(*+bYPSd/RAsRFEe2hXuXICanHL2D7vhtPE5Wb12LOcZdhOw4hjTNfT0y+1tmj*)
(*NnvfPth0U3Fu5aoirZpd6VirwHcInMsjqsDfmw5I+F1ZMx1+xAFztlvpnA4s*)
(*rKtDrl2+wHPWKH9DfULoOS/Wn9Zejc3p8Onw7cvp8IOqILSDqOzS4kQYCI8U*)
(*CI8sohUWHvn/Pn5oe1Dk3Y96YZC7v/5S+/vhv+pM3WcekpUiMJ4d1snYE1Mv*)
(*DeaFQu19zb643nItv/SqZt7MONIX/jh5vExB6xbM7WM7ab2GfPXn296/XCqh*)
(*53gG7k+TlZPdkAh3y47tPQfXWz5xhyIqN+zgWSs7tHz5zZq5TL1cmxfvkOV1*)
(*b42y3XPFHdUqqq46ebyavXD+KktsN7L29LW/CFf5sCpYrF8HZUF3ZQ11+L42*)
(*XfecDiisq0OuXb6IDtdRflt63vRj+iQh/+7mel3IA1A6fDp8Ovx+OvzYO7IB*)
(*RkfMoi5FTEDXRBtqAOGRgVmqql2cVvWbB19+94LrVZrxW71P2Aiv2Hz4YBGK*)
(*h8EqQt/ttj/k52uV7sQsd6e+cEmRSu3lk5Nq0bUNfNt7egLhwqdvH3lWa5v7*)
(*6s+3ffRyb+MI358mK6cslZGefXvPwQ0on1SK8Mr1HzxnZQeUr/ylLUlAjOtF*)
(*3x3re1DzMpj8SvpjS0JmCULatznisVWBN31jBXHAU+eEWQwh7baZDj+sTdc9*)
(*p12FdR+QuuWL6HAt5ZfWCu4TYvqkoPz9v+rR4dPh0+EHd/jRd2SDHDMTHnEg*)
(*PDIwKzs88v6M/uJM4H//4dLnn/6WGUgecWVDCzibQ8Ij1sSsv1PlKfquZnWX*)
(*y+V13S8IV1NramJe2j5HndUhl2rXDsY/OmLwNY7g/WmwctLVJuW7x5hQUbVM*)
(*IeXz75BjecDBc1V2SPnKf0prS7VTzce3vInLYJGHeKjEJhZ3yjmrwJu+f5ug*)
(*8EhwbzogAaO8Zjr8qDbtLbLnnI5vHbXLF97hWsvv3NaTfkSfFJZ/YDycDj94*)
(*hxzL6fCreazMDj/2jsxaqEYQHnEgPDIwKzs88uYtOTzy5zcevDfzxS/uP+z+*)
(*92/flyMkfzvH+zUrhv9sjngwUEjMe4kc9NMjvlL6to+/Gde2iag/oezhtxie*)
(*JWHLA9YO35/GKiffJPp3G/fBDSifsSSqcqMOnuVuObh8QlWICehjMffxrXu1*)
(*V7JIflVU0xL3xnntdZyQljtfd/reHxMjqiBw1fLnT/Pn8aJsxc/MYpJKEt0E*)
(*ggtZt8P3t+m657SnsJ6bpLrlCz1ng8of3ycEn/OB+cfdyNLhG0vo8OOs+A6/*)
(*nzuyKi0MY7kMlPlkCZg1YKSsHWP3BUTZ2n28LQHOavHV8J9QKnOHzR1yX9Py*)
(*UuRrVc9pS6rerm9JLsT2WhB30d4OnO3HUl1R5VnQD2d5iXrN+LN2Yuj5l+vJ*)
(*7ca9R65W5WzKaoqhw7D/fN36wZr/OtWbeOQX960zk1z+JeGRFSPwbO53dCBv*)
(*rF6PfJfnusudmfu2l85yvYcKuBsPqz/L3aCn/vx/D12us98tB+1PQ5VTdruh*)
(*9WMtgy8yZq7j2xvX8riDF/hjYkz55BvLMgHfcusqcuaWIYB+uiirRQ7I3Ach*)
(*9MdEbT3P4DniR4uwsUV3rco6yTbabUXvoXElQCsM+bRN3LdAMUV0nlMBbbru*)
(*Oe0ubEiHXKN8gedsYPnj+4TAcz4w/+bi4XT4Acvp8Fd+hx93jK2FSrsvY6yr*)
(*XReU64A1ZWFc7L6AaHmZFxhh57TwSLsyJa5RdLEG1VzE3XVe0/L9Mg63JxF3*)
(*eHOxL8TOWhB3UW4H3vYjV1dUefQkjRkX5R8M1LpQrrrp23C29N3HQ25VzqZs*)
(*HJ7jgWOfP70mf7OmiH6Yn60pV/iE79esGP7bJz1w5zjfzMTk5OUIhSX+18Ry*)
(*fU3x8ilun5ySx4+3tWur4+KqrxFef6EXe9vFuOnoiJx7RHtooHLMgY15ffIU*)
(*xn5w/eVzr+tZHnXwrEPBwPLZ7tYdd8O+5Qta3TlOLsdNo+sa6tjt6lqVo+w/*)
(*uaX68qSvtiP3Lxq++62IwahYaqFw0i2Alom1VfRRROc5FdCm657T7sJ6OuS6*)
(*5Qs5Z4PL30+fEHDOe/OvNNmY7p4O37cvnuV0+Cu/w487xpaNQuJgYV26JTzi*)
(*apjCQXLHEYTwiFQp2vHzBu6qfwq4pkn7JRdV7ajDL8eDvxB7akHcRSnF0Di0*)
(*r7qc5XG26cDwSHWN7LzSw2iedudsVe6m3NfdWOL/sXzSt/vf3Oe9qUUc8ZNr*)
(*93h6ZMWIbULpyWo5afoKj0hdqyuIGLvcV3jn9sIprZ+2Zi/umv7QUX/W3jbw*)
(*Yuw7lPG9RcAWzvZQu3K0AvjuWM3CBDWekPL1U7m175aj6m8gd8tCFvoFLkmg*)
(*tyT/H56Bo36R06/j3eXipIOVwjnv2eXfF+zpyz+72BqZ2CGopYu+IJs9jnTD*)
(*Xh1gxhy0mkXUzqno0UGtc9rMzdch1y2ff/vw8vfZJ3jO+bg+Megu2VM6cxU6*)
(*fDr8x7TDbyI8Yhu2K6uFddSWIaPjAiJmLl9U5OLaWr5rICulrw1kfe3S1xcI*)
(*20SGRwZ+IfbUgryLtnCfp/2EVFfAUXFG/r3hEV/0OO60MYrracpmhCbU27fl*)
(*h0N+8EGnu/S3Xyw4Zmf95DdMzbpi9HFDbz9rLO3fdTGR8nfcjEYvV8sW1PkY*)
(*v+UZO6BnUP29ofdn58XGVn/ikfDWn2f7iOV9b+LqRRdqVE7w7Z+1MEGNI6h8*)
(*fVVuxMFzZBFaf5brkauh9THSVk8I44c28Rzzy07O17IbUfdtjeduXsrdnn4f*)
(*p5LjqhtxjulV57qlqxbSUjcDC48sqOdUXJs2to88p23hEXuHXLN8vu1jyl+j*)
(*T7Ce8/30iWHjCnuR5bXo8OnwH8cOv58OxjvOk1Jxnzb2pHwXkGrVqJxNwF16*)
(*zw7KzaZoOkHXNEcbibmW1tq43oXYUwvhkauQ9hOw72Hlqf5B29YbHjFX8IdH*)
(*bMfDEnRzNuVyhbg7ng8/lacWSSceufGpdV7W//TXDx4SHVk5+rhbtm8TGh4s*)
(*13NcbBpZXslP7FG829tvgzz3h7G3VjG/SoT/chW+vMY2UUkHV46973N2d5XC*)
(*hDWOkMX9VW74wQvJw1f8oHC9ul7fjVveqfi7ZAvPnaG1pkIHgpX07QOG2Cpw*)
(*lsxYSft5c1mHR4xzKnYw2fc5LaTrabP1yufbPqb8TfQJ1Y367xPD4yN0+Gai*)
(*vkz1pOjw4z1CHX7dDsaV26KFR+IOXEPhkchTMCg8En8trbVxA+ER1wm/FOGR*)
(*wN5e74MHER5xHw9LeMTflB1BHocfXhVen7nxae/FmfFbD7r//eADYfrWS3eZ*)
(*eGQlGXB4xN7JOKIPlRXqLs//YTszvNuLV0N3v+LtdWKiG576824furzORjFp*)
(*91U5ysaePq6yfUDjCCtfv5UbevCCMvGVXqoc252M6xbMc82x3pQf76eFiXxX*)
(*vbp3y5X0LVsMMDziG9QE3ZVF/gLc14qWbaLadECezmNmu6A4OuS65Yvd3lb+*)
(*RvqEgMUBbX5JwyN0+HT4Po9Uh99AB2hbXVltkOGRuMPWUHgk+oLqDY/0cy2t*)
(*tXG9C7G3tYaGR0LaT2h4JDRkkA6jXEejZnjEH+0VwyPBTTkyoPu3cw8P/5Ue*)
(*ALn7696jIXOff/nnN4TgyY9uPvj8C54dWUn6uFu2byItMU8T5S/WTRydecxy*)
(*Xxfg2z70AqSm6DwLY6Ij3vrzbh+2vM5WMWn3VTmFgP7cuMNwH9yg8vVfuWEH*)
(*LywXe/kqKftvfdW/+Jb7c6lEH9Pof90b5hrjqaDrvZq+tdY8Ayo5l4BjaK6i*)
(*VpvvrkzMPvuBZEDhEd9tap3oSHR4xNsh1y1f5PbRd8nR5es3Vhe+hpIVHT4d*)
(*viuXx77Dj+9gjPTlcqtrhYdH/IVR8oscJ/YbHlHKoA2uxRw817Sg8bb/Wlpr*)
(*45oXYk8tBIdHQtpPQHjEWx5HnpZgRo3wiO94WNcIbcrRdz7v3dFfsfnBB53/*)
(*OvWF+NzIn/xV5xf3mZR1hbH9Plc9zaqLnVdEuQGqbdhYR2/izf7bf5Hxpad3*)
(*IsYvAdWvRYk9gL/+/GMEe/15tw9YbmG523Ttj954Gqgcx254t/ccXF/55Fxj*)
(*lgccPEfuQeVzFsM4mbV1PMuzifiUlc08Jidvqf8jRvVzcd4j7G4illtf932p*)
(*fo+g/zuoCqql81yB1RaR3U9F3JVZbhtqhEdiO3x3m657TnsKa5TXd0WJL1/4*)
(*OWstf40+IeacF3JSm6unvYolo8Onw1dWpsOXftUO7SDkVYxkhSMRct5aDnDc*)
(*BST5S/CoPSA8Yq11rWZjrmnSfvVxLa21ce0LsbsWgsMjAe0nJDziLo96JupN*)
(*WI5V1AiP+I6HfLa7mrJyIpuRrZDbpTNT1g/4Vv/7t+9//jefEhtZeULull9T*)
(*OBqV9UJRTUTo8tU8bH1BP8uN0gesJV2Tio7PWK7lYLmX99Sf7wrrqT/v9tGB*)
(*U/tmnv1xNp7+Kkdb3Xl4he1dB9dXPmstxCz3HTxHMkHlq64rpl9tvNbWJy+v*)
(*LvKUv19qFpbb/bASWO6W/Zur6wgX9cAChJ1nyi4dn/T/6mMb5FRKFHyGN9Ph*)
(*O9p03XPaU9jqEvsRqVu+wHPWWv4afULMOS/lr7fXmG6fDt9WCzHL6fCdVkSH*)
(*H95BFEl6TgZpN0Jqv0wlXbuPC4inpcWHR8xSlWtZ99hzTZPH+9HX0lobN3Ih*)
(*ttdCRHjEyFuszaC+31YerYnI17hK46kbHvEcD3UFvdakpqztmB40DbwyXvnl*)
(*F8f+2voZ3+5/b916MP0ZsREAAAAAwID0+aMe0KzfftF7jOTQRf2dmh980Ln8*)
(*CXOxAgAAAAAGyfeOG7CYHn658NHfPfzwVw8vf/LF1P2Hn/6WWVgBAAAAAANH*)
(*dAQAAAAAADxOlGlaF6xzAgMAAAAAAKxY+ny4zDkCAAAAAAAAAAAAAAAAAAAA*)
(*AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA*)
(*AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA*)
(*AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAYPHdaK8fyu04*)
(*Yyw+s2Oowlwh3Xx9+4Z9O3OhkHU1bfXvQra2TAEAAAAAAGIlcYY8+KD8Q9IL*)
(*eVRiEsn669vtHUakQlnxjLnczC1L60bxD1vsw54pAAAAAABAvF7kohoP6cUe*)
(*fPGRSjAljU8Y4Q8jutFdQ368RIu1+MMj9kwBAAAAAAD6oYdH5EBGzvJ0iT88*)
(*Ikdd0jdk0r8nr+Ioj5J4Yh+ERwAAAABgCcxPX5wYHx9//3ZfW89df2dc8M71*)
(*uX6ymr939Vxv84mLU/c7tlVuJYmMj5+7eteyUsA6nft306y6ed26Ny+vMpWm*)
(*4ihNQDJFmT7O1pQqx7nj3jq+/b60PKPXdrFb565qBbbko+QWlVV9dVqMGmVI*)
(*5wsxAhnlNCJiQMKMVAjvzdifBhFmF9HmHpGeZyE8AgAAAACLqXP/+pv7Nq8e*)
(*to/TAhgzTYpjzYCsOteOPT/WqiQxvG7f2zPKiHf2QltdJVnp7GzkOp2p07vX*)
(*DSurtMaeP3atktf8lSMvrFZWGd548IKSUUgy6Xozk201MbVyjBIb++StY3WK*)
(*UU21tmfP7lOK3C3wkSvznnyU3IKzqqmRFqPMzNp2zz6iTBBSECMVysys8kQh*)
(*ladHLHEZbS1vpgAAAACAAbjxo5ezcXJrbGy0xtA2GeGNbPj2LtX3Ju5EZTV7*)
(*ZkeyaPiJp7ft2rU9HxeP7jlfDHdvH9vUKtbY1Vtp7UjvD60tr9+JWKdzfk8v*)
(*p9bI2s3bk3W2Pf1EktnovveyVT75Ly8MV1bZ9lSSxlBr66kyahGQTLLatfaG*)
(*dAQ//MSqrFDVysmXt0ae2lYt74bD2vsbzjq+fHyX4He/9ve7SY3sOW/klZS5*)
(*2K1Nx7JnM+5MfE9K54X/9Svd1TYduxOeVT1NtRgtVc/cI+JrL/5IhXXukUpe*)
(*6us1ep7GEsIjAAAAALBY7hzblDw70J6c6Th+3vaTf3SPzKpz9pWR3tB2x5ki*)
(*/jB/5VAyln/mSDZ2v3V4Q/ef3zihPAdy4dVV1fRC1pnY2U121asXqoPo2RPf*)
(*GFL2YvbsiTdvlm+edC7tXzNUhgiCk8liKKNbevsu1s3rW1q9/T5dLs6CGCOv*)
(*nC3+5K9jSRoqKupv4faRZ5LAS7t4vqUzczoJMVTzMqT74FzFyKqmhlqMypiK*)
(*JGgNX6RCnklEmr9Vztz6dArhEQAAAABYFHOXfpYP2QccHvFn1R0Dd8e1rRff*)
(*UKbkuP+Tl79aCUmkLyKoY99skKzNBOFeJymCNtpPIweuvbh26El1jZBk5t54*)
(*sbtf6w8Z79sUZTu9vbvBk4euqX9NSlx5FKOf8EiWSFnAtG7WHLyirHX93z/l*)
(*fuwj3Sdn5MPISl/ese6/ZUkzLabC/viGZxVnpMJ+XMzZWMX1LMUiPAIAAAAA*)
(*S2HQ4RFvVpf2957v2H66fIJiZrK9ZTSd3aG1cyL943v7krdZ1u7J55eYPbtn*)
(*rTZIDlgnfWBjaHRLPulG59qxJC8tciAVccPhW/kfApJJgx9bT82VU8VqE4im*)
(*T7sYtTd3aquSWR/hESPgk6ZZjcTM33xz59ps6o7KjimyyMeOM9JEu7asVJ9N*)
(*/uHTW8wZWXpzruzfuPk/Tnl2pE6LqcwR4pj6Q1xFmGKl+mVe7U/2gptp2/7u*)
(*zhQAAAAAMHD1wyOrtn9/PPtQzMfOz7fIWSV/zQeC5ah9eGxspPowQHc8nc7j*)
(*0Rp7/mA7nRa1tXaPMo1pwDqdqVMvpgPp4XW72wfT2T1HxRF8sUn6co3yDIU/*)
(*meSJk7/3zJZsjo9cdZLX5BUdZU6TXsozf6q+pBNXxwtiVEOZBbQMJmSTe1ge*)
(*H/E/VRMSQOm9w7NOq99ebGTjfm2yW0mtFgMAAAAAQLD64RHV8OoXKp9DCcgq*)
(*H+z2Ru1pmKE38cSFWfOrqZ2Zd/7gfyky+urXDl818wlZZ/769zd/pVhnbOfE*)
(*Pccudq4d6pVDmPTTnUz5MEA+VWwxG2r+MEv25s/Q8LqX2yfHx8dPtncXX2vR*)
(*wiPBdSx/bbYIj8zffDP74M7wut2npxwzz8yd2eF9dMT1YdsKNUISHBtZqN9i*)
(*AAAAAAAIUys8onztZNvTY9mTEuXnUAKySv66atPzxag9mxY1HewWzwJ0Zt5O*)
(*P0s7vHrz5jzQMPbSqalK1CJknfkrR5IhdWvkqc15MML8bm+e4LXDG4fVGU1D*)
(*k0mnJ9n6p9cr79N0Zt7crj4Sk779U5U90VE+AxFXx2lUQ5+XIwscbCmDCdns*)
(*HmkphadHxMlKQrKS5RGSmYjYyELNFgMAAAAAQLBa4RHD/JXvJs86iNNZyFkl*)
(*r5gMqaP2nnSGiWztuYmdyVdgXjyRfVKmeKeifKwjZJ0bh3tv37Q27M8yKl40*)
(*EaIN81falthISDLJzurzrkrxiPmb544eSGIfB45OXJyeT6YkcX0l11XHlqhG*)
(*Os2JGkxIqiyZk8R83iIk8hEQQFF0Zk59a6S1/g8nQ2MjC7VaDAAAAAAAEZoN*)
(*jzhfuJCzSj7kOvSVb/zHm8rrIlcOrikfBUjX0T+Pkn2sNpuLM2Sd83t6D29o*)
(*X2LJ3nFRow3ZNCbDGw8Ls5KEJJN+7aY3NatZB6v2XzLSLDI+tbU11NpxxvEl*)
(*XWsd26MaaYmf/MN3lU8Mpx/XMUMxDT86kkjeqflXB1/U5yFx6r/FAAAAAAAQ*)
(*wxse6dy/fm58/NxNz4SgqezJBnnYbMkqjWCMbD9dPlbQmTmxtVWJQKTjdX0W*)
(*03TDPMGQddISrNl/qTpAz+IalWjD7JmdvbdQrDO2BiWTRCS0R09mT/derrEG*)
(*Pzozk0lQxv1IhrWO0/iAuHFeuOpnhucvvLpGCChlQRN35MORlaScb6Q3p214*)
(*hKT/FgMAAAAAgNft98dzrz3bG4E++1r2z3euq4Pizvk92RdTR/e950pn/OiB*)
(*7dkkHK0Nh29Iq9iyyt6KyacoPdl+OZ1VonwnJi9Fa+SpdBbTbm7b0pXKiU4D*)
(*1rl9bFP6kZPVm/cdHU9nQ30umY+jHOrffiP9Js3Y77R/NK55/3ZwMnkIobtX*)
(*2w4cHS93qxoxmbv+Tprwyfbe7ZvXJtOK6C/zBNRxmlYa1dDCQ4X0faCh1thz*)
(*SZGP7ksLbL5UlEU+1NhPVFYafS5WX4SkkRYDAAAAAIBX+VkVnf7Gxp3Xt+QT*)
(*PWQvqLjTST+IEptVNgVqhT6j6uyFg9oaaW77Kl/tDVinOzR/aaylr9Iae74c*)
(*rdsLXD7HEJBMslfHtoy6Cqx/laY18pRWe2F13OOPamTPxFQTMqekDYl8BARQ*)
(*ykzFuVilr/0WGmoxAAAAAAB4XD6+y+J7E/rUDbMfnDiwa9eBN6W3a6rpHDg6*)
(*fq76mZbYrOanL77V3pumJCXUXePe1XMn01V27W2fFFcKWKdzf+riRDYbajYd*)
(*aliBd+06fjk4mXKvspX2tt8yVim+StOrvKv3xPeXAuq46+rJf9ld/h/Oup/n*)
(*6L0oZS9N1yfv/N+9hX/mevklLKssw5nJn1o+CjQz+a78ceIGWwwAAAAAAAAA*)
(*AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA*)
(*AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA*)
(*AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA*)
(*AACA5e//Bw+HPD0=*)
(*"], "Byte", ColorSpace -> "RGB", ImageResolution -> {191.9986, 191.9986}, Interleaving -> True]*)


(* ::Text:: *)
(*We run into this issue when we try to use NIntegrate for the survival probability integral*)


(* ::Input::Initialization:: *)
(*Clear[GetWeakCaptureProbabilityNuc]
GetWeakCaptureProbabilityNuc[lDict_,\[Sigma]Dict_]:= Module[{v\[Chi]min,v\[Chi]raw,v\[Chi],\[Chi]traj,bE,m\[Chi],\[Omega]min,params,NucleusParams,tdom,\[Kappa],regionindex,mN,\[Omega]capture,\[Delta]\[Xi],\[Xi]of\[Omega]andt,\[Omega]maxoft,d\[Sigma]dERoftand\[Omega],\[Sigma]oftand\[Omega]min,\[Lambda]invoft,integralof\[Lambda]invwrtr},
\[Omega]min="\[Omega]min"/.\[Sigma]Dict;
params=\[Sigma]Dict["params"];
NucleusParams=\[Sigma]Dict["NucleusParams"];
\[Delta]\[Xi]=(10^"\[Xi]"/.\[Sigma]Dict)[[1]];
m\[Chi]=("m\[Chi]"/.lDict);
tdom = "domain"/.lDict;
\[Kappa]=("\[Kappa]"/.lDict);
mN = "mN"/.\[Sigma]Dict;
regionindex = "regionindex"/.\[Sigma]Dict;

(*Print[mN];*)

v\[Chi]min=3Sqrt[("\[HBar]" \[Omega]min)/(2 m\[Chi])]/.params; (*minimum allowable velocity given the lower bound on energy transfer, times an O(1) fudge factor for numerical stability (3) *)
v\[Chi]raw=("v(t)"/.lDict);
v\[Chi]=If[v\[Chi]raw[#]>v\[Chi]min,("v(t)"/.lDict)[#],0]&;

\[Chi]traj=("l(t)"/.lDict);
bE=("bE"/.lDict);(*impact parameter*)

\[Xi]of\[Omega]andt=Abs[EnergyLoss`\[Xi]of\[Omega]Nuc[#1,\[Omega]min,("v(t)"/.lDict)[#2],m\[Chi],mN,params]-\[Delta]\[Xi]]&;
\[Omega]maxoft=EnergyLoss`\[Omega]maxNuc[("v(t)"/.lDict)[#1],m\[Chi],mN,params]&;

(*d\[Sigma]dERoftand\[Omega]=10^Re@("d\[Sigma]f"/.\[Sigma]Dict)[Log10@v\[Chi][#1],Log10@\[Xi]of\[Omega]andt[#2,#1]]&; (*not currently used*)*)

\[Sigma]oftand\[Omega]min=10^(Re[("\[Sigma]of\[Xi]f"/.\[Sigma]Dict)[Log10@v\[Chi][#1],Log10@\[Xi]of\[Omega]andt[#2,#1]]])&;


(*Post \[Omega]capture*)
\[Omega]capture = 1/(2 ("\[HBar]"/.params))m\[Chi] (v\[Chi][#]^2 - ("vesc"/.Constants`EarthRepl)^2)&;


\[Lambda]invoft=If[v\[Chi]raw[#]>v\[Chi]min && \[Omega]maxoft[#]>\[Omega]capture[#],\[Kappa]^2("nI"/.NucleusParams)\[Sigma]oftand\[Omega]min[#,\[Omega]capture[#]],0]&;
(*\[Lambda]invoft=If[v\[Chi]raw[#]>v\[Chi]min ,\[Kappa]^2("nI"/.NucleusParams)\[Sigma]oftand\[Omega]min[#,\[Omega]capture[#]],0]&;*)



integralof\[Lambda]invwrtr=Total@Table[\[Lambda]invoft[t]v\[Chi][t]Region\[CapitalTheta][ltor[\[Chi]traj[t],bE],regionindex],{t,Round[tdom[[1]]],Round[tdom[[2]]],1}];

integralof\[Lambda]invwrtr
]*)


(* ::Subsubsection::Closed:: *)
(*Weak Capture Regime Capture Probability - e*)


(* ::Input::Initialization:: *)
Clear[GetWeakCaptureProbabilityfromv]
GetWeakCaptureProbabilityfromv[vDict_,\[Sigma]Dict_]:= Module[{v\[Chi]min,v\[Chi]raw,v\[Chi],\[Chi]traj,bE,m\[Chi],\[Omega]min,params,ldom,\[Kappa],regionindex,\[Sigma]dom,\[Omega]capture,\[Delta]\[Xi],\[Xi]of\[Omega]andl,\[Omega]maxofl,d\[Sigma]dERoftand\[Omega],\[Sigma]oftand\[Omega]min,rCore,rE,lcore,lcrust,\[Lambda]invoft,\[Lambda]integrand,integralof\[Lambda]invwrtr},
\[Omega]min="\[Omega]min"/.\[Sigma]Dict;
params=\[Sigma]Dict["params"];
\[Delta]\[Xi]=(10^"\[Xi]"/.\[Sigma]Dict)[[1]];
m\[Chi]=("m\[Chi]"/.vDict);
ldom = "domain"/.vDict;
\[Kappa]=("\[Kappa]"/.vDict);
regionindex = "regionindex"/.\[Sigma]Dict; (*1 in crust, 2 in core*)
\[Sigma]dom=("\[Sigma]of\[Xi]f"/.\[Sigma]Dict)["Domain"];

v\[Chi]min=Max[3Sqrt[("\[HBar]" \[Omega]min)/(2 m\[Chi])]/.params,"vesc"/.Constants`EarthRepl];(*minimum allowable velocity given the lower bound on energy transfer, times an O(1) fudge factor for numerical stability (3) *)
v\[Chi]raw=("v(l)"/.vDict);
v\[Chi]=If[v\[Chi]raw[#]>v\[Chi]min,("v(l)"/.vDict)[#],0]&;
(*\[Chi]traj=("l(t)"/.vDict);*)
bE=("bE"/.vDict);(*impact parameter*)

(*\[Xi]of\[Omega]andl=Abs[EnergyLoss`\[Xi]of\[Omega]andv\[Chi][#1,\[Omega]min,("v(l)"/.vDict)[#2],m\[Chi],params]-\[Delta]\[Xi]]&;*)
(*\[Xi]of\[Omega]andl=Abs[EnergyLoss`\[Xi]of\[Omega]andv\[Chi][#1,\[Omega]min,("v(l)"/.vDict)[#2],m\[Chi],params]]&;*)
\[Xi]of\[Omega]andl=EnergyLoss`\[Xi]of\[Omega]andv\[Chi][#1,\[Omega]min,("v(l)"/.vDict)[#2],m\[Chi],params]&;
\[Omega]maxofl=EnergyLoss`\[Omega]maxofv\[Chi][("v(l)"/.vDict)[#1],m\[Chi],params]&;

(*d\[Sigma]dERoftand\[Omega]=10^Re@("d\[Sigma]f"/.\[Sigma]Dict)[Log10@v\[Chi][#1],Log10@\[Xi]of\[Omega]andt[#2,#1]]&; (*not currently used*)*)

\[Sigma]oftand\[Omega]min=10^(Re[("\[Sigma]of\[Xi]f"/.\[Sigma]Dict)[Log10@v\[Chi][#1],Log10@\[Xi]of\[Omega]andl[#2,#1]]])&;

(*Post \[Omega]capture*)
\[Omega]capture = 1/(2 ("\[HBar]"/.params)) m\[Chi] (v\[Chi][#]^2 - ("vesc"/.Constants`EarthRepl)^2)&;


\[Lambda]invoft=If[v\[Chi]raw[#]>v\[Chi]min && \[Omega]maxofl[#]>\[Omega]capture[#]&&(\[Xi]of\[Omega]andl[\[Omega]capture[#],#])>10^\[Sigma]dom[[2,1]],\[Kappa]^2 ("ne"/.params)\[Sigma]oftand\[Omega]min[#,\[Omega]capture[#]],0]&;

rCore="rcore"/.Constants`EarthRepl;
rE = "rE"/.Constants`EarthRepl;

\[Lambda]integrand[l_?NumericQ]:=\[Lambda]invoft[l];

If[regionindex==1&&bE<rCore,
(*if material is in the crust, and trajectory passes through the core, divide the integration region by hand*)
lcore = Sqrt[rCore^2-bE^2];
lcrust = Sqrt[rE^2-bE^2]-lcore/2;
If[ldom[[1]]+lcrust +lcore>ldom[[2]],Print["Integration domain is fucky!"]];

integralof\[Lambda]invwrtr=NIntegrate[\[Lambda]integrand[l],{l,ldom[[1]],ldom[[1]]+lcrust},PrecisionGoal->3,MaxRecursion->8,AccuracyGoal->3] + NIntegrate[\[Lambda]integrand[l],{l,ldom[[1]]+lcrust +lcore,ldom[[2]]},PrecisionGoal->3,MaxRecursion->8,AccuracyGoal->3];

Return[integralof\[Lambda]invwrtr,Module];
];

If[regionindex==2,
(*if the material is in the core, restrict the integration region to only in the core*)

If[bE>rCore,
(*the material is in the core but the particle doesn't pass through the core. So return 0, this material doesn't contribute to the capture rate*)
Return[0,Module];
];

lcore = Sqrt[rCore^2-bE^2];
lcrust = Sqrt[rE^2-bE^2]-lcore/2;
If[ldom[[1]]+lcrust +lcore>ldom[[2]],Print["Integration domain is fucky!"]];

integralof\[Lambda]invwrtr=NIntegrate[\[Lambda]integrand[l],{l,ldom[[1]]+lcrust,ldom[[1]]+lcrust+lcore},PrecisionGoal->3,MaxRecursion->8,AccuracyGoal->3];

Return[integralof\[Lambda]invwrtr,Module];
];

(*If both of these are false, then the material is in the crust, and the trajectory does not pass through the core. So integration region is continuous*)


integralof\[Lambda]invwrtr=NIntegrate[\[Lambda]integrand[l],{l,ldom[[1]],ldom[[2]]},PrecisionGoal->3,MaxRecursion->8,AccuracyGoal->3];

integralof\[Lambda]invwrtr
]


(* ::Subsubsection::Closed:: *)
(*Weak Capture Regime Capture Probability - Nuc*)


(* ::Input::Initialization:: *)
Clear[GetWeakCaptureProbabilityNucfromv]
GetWeakCaptureProbabilityNucfromv[vDict_,\[Sigma]Dict_]:= Module[{v\[Chi]min,v\[Chi]raw,v\[Chi],\[Chi]traj,bE,m\[Chi],\[Omega]min,params,NucleusParams,ldom,\[Kappa],regionindex,\[Sigma]dom,mN,\[Omega]capture,\[Delta]\[Xi],\[Xi]of\[Omega]andl,\[Omega]maxofl,d\[Sigma]dERoftand\[Omega],\[Sigma]oftand\[Omega]min,\[Lambda]invoft,\[Lambda]integrand,rCore,rE,lcore,lcrust,integralof\[Lambda]invwrtr},
\[Omega]min="\[Omega]min"/.\[Sigma]Dict;
params=\[Sigma]Dict["params"];
NucleusParams=\[Sigma]Dict["NucleusParams"];
\[Delta]\[Xi]=(10^"\[Xi]"/.\[Sigma]Dict)[[1]];
m\[Chi]=("m\[Chi]"/.vDict);
ldom = "domain"/.vDict;
\[Kappa]=("\[Kappa]"/.vDict);
mN = "mN"/.\[Sigma]Dict;
regionindex = "regionindex"/.\[Sigma]Dict;
\[Sigma]dom=("\[Sigma]of\[Xi]f"/.\[Sigma]Dict)["Domain"];

(*Print[mN];*)

v\[Chi]min=Max[3Sqrt[("\[HBar]" \[Omega]min)/(2 m\[Chi])]/.params,"vesc"/.Constants`EarthRepl];(*minimum allowable velocity given the lower bound on energy transfer, times an O(1) fudge factor for numerical stability (3) *)
v\[Chi]raw=("v(l)"/.vDict);
v\[Chi]=If[v\[Chi]raw[#]>v\[Chi]min,("v(l)"/.vDict)[#],0]&;

(*\[Chi]traj=("l(t)"/.vDict);*)
bE=("bE"/.vDict);(*impact parameter*)

(*\[Xi]of\[Omega]andl=Abs[EnergyLoss`\[Xi]of\[Omega]Nuc[#1,\[Omega]min,("v(l)"/.vDict)[#2],m\[Chi],mN,params]-\[Delta]\[Xi]]&;*)
(*\[Xi]of\[Omega]andl=Abs[EnergyLoss`\[Xi]of\[Omega]Nuc[#1,\[Omega]min,("v(l)"/.vDict)[#2],m\[Chi],mN,params]]&;*)
\[Xi]of\[Omega]andl=EnergyLoss`\[Xi]of\[Omega]Nuc[#1,\[Omega]min,("v(l)"/.vDict)[#2],m\[Chi],mN,params]&;
\[Omega]maxofl=EnergyLoss`\[Omega]maxNuc[("v(l)"/.vDict)[#1],m\[Chi],mN,params]&;

(*d\[Sigma]dERoftand\[Omega]=10^Re@("d\[Sigma]f"/.\[Sigma]Dict)[Log10@v\[Chi][#1],Log10@\[Xi]of\[Omega]andt[#2,#1]]&; (*not currently used*)*)

\[Sigma]oftand\[Omega]min=10^(Re[("\[Sigma]of\[Xi]f"/.\[Sigma]Dict)[Log10@v\[Chi][#1],Log10@\[Xi]of\[Omega]andl[#2,#1]]])&;


(*Post \[Omega]capture*)
\[Omega]capture = 1/(2 ("\[HBar]"/.params)) m\[Chi] (v\[Chi][#]^2 - ("vesc"/.Constants`EarthRepl)^2)&;


(*\[Lambda]invoft=If[v\[Chi]raw[#]>v\[Chi]min && \[Omega]maxofl[#]>\[Omega]capture[#]&&Log10@(\[Xi]of\[Omega]andl[\[Omega]capture[#],#])>\[Sigma]dom[[2,1]],\[Kappa]^2("nI"/.NucleusParams)\[Sigma]oftand\[Omega]min[#,\[Omega]capture[#]],0]&;*)
\[Lambda]invoft=If[v\[Chi]raw[#]>v\[Chi]min && \[Omega]maxofl[#]>\[Omega]capture[#]&&(\[Xi]of\[Omega]andl[\[Omega]capture[#],#])>10^\[Sigma]dom[[2,1]],\[Kappa]^2 ("nI"/.NucleusParams)\[Sigma]oftand\[Omega]min[#,\[Omega]capture[#]],0]&;

rCore="rcore"/.Constants`EarthRepl;
rE = "rE"/.Constants`EarthRepl;

\[Lambda]integrand[l_?NumericQ]:=\[Lambda]invoft[l];

If[regionindex==1&&bE<rCore,
(*if material is in the crust, and trajectory passes through the core, divide the integration region by hand*)
lcore = Sqrt[rCore^2-bE^2];
lcrust = Sqrt[rE^2-bE^2]-lcore/2;
If[ldom[[1]]+lcrust +lcore>ldom[[2]],Print["Integration domain is fucky!"]];

integralof\[Lambda]invwrtr=NIntegrate[\[Lambda]integrand[l],{l,ldom[[1]],ldom[[1]]+lcrust},PrecisionGoal->3,MaxRecursion->8,AccuracyGoal->3] + NIntegrate[\[Lambda]integrand[l],{l,ldom[[1]]+lcrust +lcore,ldom[[2]]},PrecisionGoal->3,MaxRecursion->8,AccuracyGoal->3];

Return[integralof\[Lambda]invwrtr,Module];
];

If[regionindex==2,
(*if the material is in the core, restrict the integration region to only in the core*)
If[bE>rCore,
(*the material is in the core but the particle doesn't pass through the core. So return 0, this material doesn't contribute to the capture rate*)
Return[0,Module];
];

lcore = Sqrt[rCore^2-bE^2];
lcrust = Sqrt[rE^2-bE^2]-lcore/2;
If[ldom[[1]]+lcrust +lcore>ldom[[2]],Print["Integration domain is fucky!"]];

integralof\[Lambda]invwrtr=NIntegrate[\[Lambda]integrand[l],{l,ldom[[1]]+lcrust,ldom[[1]]+lcrust+lcore},PrecisionGoal->3,MaxRecursion->8,AccuracyGoal->3];

Return[integralof\[Lambda]invwrtr,Module];
];

(*If both of these are false, then the material is in the crust, and the trajectory does not pass through the core. So integration region is continuous*)


(*THESE NEED TO BE REWRITTEN. THEY ASSUME A LINEAR TRAJECTORY. Maybe include the l values where you pass through the mantle / core boundary*)

integralof\[Lambda]invwrtr=NIntegrate[\[Lambda]integrand[l],{l,ldom[[1]],ldom[[2]]},PrecisionGoal->3,MaxRecursion->8,AccuracyGoal->3];

integralof\[Lambda]invwrtr
]


(* ::Subsubsection::Closed:: *)
(*Weak Capture Regime Capture Interaction Length - e*)


(* ::Input::Initialization:: *)
Clear[GetWeakCaptureInteractionLength]
GetWeakCaptureInteractionLength[\[Sigma]Dict_]:= Module[{v\[Chi]min,m\[Chi],\[Omega]min,params,regionindex,\[Sigma]dom,\[Omega]capture,\[Delta]\[Xi],\[Xi]of\[Omega]andv\[Chi]t,\[Omega]maxofv\[Chi]t,d\[Sigma]dERoftand\[Omega],\[Sigma]ofv\[Chi]and\[Omega]min,\[Lambda]invof\[Chi]},
\[Omega]min="\[Omega]min"/.\[Sigma]Dict;
params=\[Sigma]Dict["params"];
\[Delta]\[Xi]=(10^"\[Xi]"/.\[Sigma]Dict)[[1]];
m\[Chi]=("m\[Chi]"/.\[Sigma]Dict);
regionindex = "regionindex"/.\[Sigma]Dict;
\[Sigma]dom=("\[Sigma]of\[Xi]f"/.\[Sigma]Dict)["Domain"];

v\[Chi]min=Max[3Sqrt[("\[HBar]" \[Omega]min)/(2 m\[Chi])]/.params,"vesc"/.Constants`EarthRepl]; (*restrict minimum velocity according to minimum velocity of interpolation (see energy loss package), and the escape velocity*)

\[Xi]of\[Omega]andv\[Chi]t=EnergyLoss`\[Xi]of\[Omega]andv\[Chi][#1,\[Omega]min,#2,m\[Chi],params]&;
\[Omega]maxofv\[Chi]t=EnergyLoss`\[Omega]maxofv\[Chi][#,m\[Chi],params]&;

\[Sigma]ofv\[Chi]and\[Omega]min=10^(Re[("\[Sigma]of\[Xi]f"/.\[Sigma]Dict)[Log10@#1,Log10@(\[Xi]of\[Omega]andv\[Chi]t[#2,#1])]])&;

\[Omega]capture = 1/(2 ("\[HBar]"/.params)) m\[Chi] (#^2 - ("vesc"/.Constants`EarthRepl)^2)&; (*energy transfer needed to capture the particle*)

(*Print[("\[Sigma]of\[Xi]f"/.\[Sigma]Dict)["Domain"]];*)
(*Print[("\[Sigma]of\[Xi]f"/.\[Sigma]Dict)["Domain"]];
Print@Plot[Log10@(\[Xi]of\[Omega]andv\[Chi]t[\[Omega]capture[#],#])&[v\[Chi]],{v\[Chi],"vesc"/.EarthRepl,5+"vesc"/.EarthRepl}];
Print[Plot[Log10@(\[Xi]of\[Omega]andv\[Chi]t[\[Omega]capture[#],#])-\[Sigma]dom[[2,1]]&[v\[Chi]],{v\[Chi],"vesc"/.EarthRepl,5+"vesc"/.EarthRepl}]];*)

(*With[{v\[Chi]="vesc"/.EarthRepl},
Print["v\[Chi]min"];
Print[v\[Chi]min];
Print["\[Omega]capture and \[Omega]maxofv\[Chi]"];
Print[\[Omega]capture[v\[Chi]]];
Print[\[Omega]maxofv\[Chi]t[v\[Chi]]];
Print["\[Lambda]^-1"];
Print[("ne"/.params)\[Sigma]ofv\[Chi]and\[Omega]min[v\[Chi],\[Omega]capture[v\[Chi]]]];
Print[Log10@\[Xi]of\[Omega]andv\[Chi]t[\[Omega]capture[v\[Chi]],v\[Chi]]];
Print[Log10@\[Xi]of\[Omega]andv\[Chi]t[\[Omega]maxofv\[Chi]t[v\[Chi]],v\[Chi]]];];*)

(*\[Lambda]invof\[Chi]=If[#>v\[Chi]min && \[Omega]maxofv\[Chi]t[#]>\[Omega]capture[#],("ne"/.params)\[Sigma]ofv\[Chi]and\[Omega]min[#,\[Omega]capture[#]],0]&;*)
(*\[Lambda]invof\[Chi]=If[#>v\[Chi]min && \[Omega]maxofv\[Chi]t[#]>\[Omega]capture[#]&&Log10@(\[Xi]of\[Omega]andv\[Chi]t[\[Omega]capture[#],#])>\[Sigma]dom[[2,1]],("ne"/.params)\[Sigma]ofv\[Chi]and\[Omega]min[#,\[Omega]capture[#]],0]&; (*must have v > minimum velocity allowed, must have maximum possible energy transfer greater than energy needed to capture, and must have that energy transfer be in the domain of the interpolation function*)*)
\[Lambda]invof\[Chi]=If[#>v\[Chi]min && \[Omega]maxofv\[Chi]t[#]>\[Omega]capture[#]&&(\[Xi]of\[Omega]andv\[Chi]t[\[Omega]capture[#],#])>10^\[Sigma]dom[[2,1]],("ne"/.params)\[Sigma]ofv\[Chi]and\[Omega]min[#,\[Omega]capture[#]],0]&; (*must have v > minimum velocity allowed, must have maximum possible energy transfer greater than energy needed to capture, and must have that energy transfer be in the domain of the interpolation function*)


(*Print@Plot[\[Lambda]invof\[Chi][v\[Chi]],{v\[Chi],"vesc"/.EarthRepl,5+"vesc"/.EarthRepl}];*)

\[Lambda]invof\[Chi]
]


(* ::Input::Initialization:: *)
Clear[GetTotaleWeakCaptureInteractionLength]
GetTotaleWeakCaptureInteractionLength[\[Sigma]Dicts_]:=Module[{},
Total@Table[GetWeakCaptureInteractionLength[\[Sigma]Dicts[[i]]][#],{i,Length[\[Sigma]Dicts]}]&
]


(* ::Input:: *)
(*?\[Xi]of\[Omega]andv\[Chi]*)


(* ::Input:: *)
(*?\[Omega]maxofv\[Chi]*)


(* ::Subsubsection::Closed:: *)
(*Weak Capture Regime Capture Interaction Length - Nuc*)


(* ::Input::Initialization:: *)
Clear[GetWeakCaptureInteractionLengthNuc]
GetWeakCaptureInteractionLengthNuc[\[Sigma]Dict_]:= Module[{v\[Chi]min,m\[Chi],mN,\[Omega]min,params,NucleusParams,regionindex,\[Sigma]dom,\[Omega]capture,\[Delta]\[Xi],\[Xi]of\[Omega]andv\[Chi]t,\[Omega]maxofv\[Chi]t,d\[Sigma]dERoftand\[Omega],\[Sigma]ofv\[Chi]and\[Omega]min,\[Lambda]invof\[Chi]},
\[Omega]min="\[Omega]min"/.\[Sigma]Dict;
params=\[Sigma]Dict["params"];
NucleusParams=\[Sigma]Dict["NucleusParams"];
\[Delta]\[Xi]=(10^"\[Xi]"/.\[Sigma]Dict)[[1]];
m\[Chi]=("m\[Chi]"/.\[Sigma]Dict);
mN = "mN"/.\[Sigma]Dict;
regionindex = "regionindex"/.\[Sigma]Dict;
\[Sigma]dom=("\[Sigma]of\[Xi]f"/.\[Sigma]Dict)["Domain"];

(*v\[Chi]min=3Sqrt[("\[HBar]" \[Omega]min)/(2 m\[Chi])]/.params; *)
v\[Chi]min=Max[3Sqrt[("\[HBar]" \[Omega]min)/(2 m\[Chi])]/.params,"vesc"/.Constants`EarthRepl];(*minimum allowable velocity given the lower bound on energy transfer, times an O(1) fudge factor for numerical stability (3) *)
(*v\[Chi]=If[#>v\[Chi]min,#,0]&;*)
(*\[Chi]traj=("l(t)"/.vDict);*)

(*\[Xi]of\[Omega]andv\[Chi]t=Abs[EnergyLoss`\[Xi]of\[Omega]Nuc[#1,\[Omega]min,#2,m\[Chi],mN,params]-\[Delta]\[Xi]]&;*)
(*\[Xi]of\[Omega]andv\[Chi]t=Abs[EnergyLoss`\[Xi]of\[Omega]Nuc[#1,\[Omega]min,#2,m\[Chi],mN,params]]&;*)
\[Xi]of\[Omega]andv\[Chi]t=EnergyLoss`\[Xi]of\[Omega]Nuc[#1,\[Omega]min,#2,m\[Chi],mN,params]&;
\[Omega]maxofv\[Chi]t=EnergyLoss`\[Omega]maxNuc[#,m\[Chi],mN,params]&;

(*d\[Sigma]dERoftand\[Omega]=10^Re@("d\[Sigma]f"/.\[Sigma]Dict)[Log10@v\[Chi][#1],Log10@\[Xi]of\[Omega]andt[#2,#1]]&; (*not currently used*)*)

(*\[Sigma]ofv\[Chi]and\[Omega]min=10^(Re[("\[Sigma]of\[Xi]f"/.\[Sigma]Dict)[Log10@#1,Log10@(Min[\[Xi]of\[Omega]andv\[Chi]t[#2,#1],1-\[Delta]\[Xi]])]])&;*)
\[Sigma]ofv\[Chi]and\[Omega]min=10^(Re[("\[Sigma]of\[Xi]f"/.\[Sigma]Dict)[Log10@#1,Log10@(\[Xi]of\[Omega]andv\[Chi]t[#2,#1])]])&;

(*Post \[Omega]capture*)
\[Omega]capture = 1/(2 ("\[HBar]"/.params)) m\[Chi] (#^2 - ("vesc"/.Constants`EarthRepl)^2)&; (*energy transfer needed to capture the particle*)


(*Print[("\[Sigma]of\[Xi]f"/.\[Sigma]Dict)["Domain"]];
Print@Plot[Log10@(\[Xi]of\[Omega]andv\[Chi]t[\[Omega]capture[#],#])&[v\[Chi]],{v\[Chi],"vesc"/.EarthRepl,5+"vesc"/.EarthRepl}];
Print[Plot[Log10@(\[Xi]of\[Omega]andv\[Chi]t[\[Omega]capture[#],#])-\[Sigma]dom[[2,1]]&[v\[Chi]],{v\[Chi],"vesc"/.EarthRepl,5+"vesc"/.EarthRepl}]];

With[{v\[Chi]=1.0001 "vesc"/.EarthRepl},
Print["v\[Chi]min"];
Print[v\[Chi]-v\[Chi]min];
Print["\[Omega]maxofv\[Chi]-\[Omega]capture"];
Print[\[Omega]maxofv\[Chi]t[v\[Chi]]-\[Omega]capture[v\[Chi]]];
Print["\[Lambda]^-1"];
Print[("nI"/.NucleusParams)\[Sigma]ofv\[Chi]and\[Omega]min[v\[Chi],\[Omega]capture[v\[Chi]]]];
Print["\[Xi]of\[Omega]andv\[Chi]s"];
Print[Log10@\[Xi]of\[Omega]andv\[Chi]t[\[Omega]capture[v\[Chi]],v\[Chi]]];
Print[Log10@\[Xi]of\[Omega]andv\[Chi]t[\[Omega]maxofv\[Chi]t[v\[Chi]],v\[Chi]]]];*)

\[Lambda]invof\[Chi]=If[#>v\[Chi]min && \[Omega]maxofv\[Chi]t[#]>\[Omega]capture[#]&&(\[Xi]of\[Omega]andv\[Chi]t[\[Omega]capture[#],#])>10^\[Sigma]dom[[2,1]],("nI"/.NucleusParams)\[Sigma]ofv\[Chi]and\[Omega]min[#,\[Omega]capture[#]],0]&;
(*\[Lambda]invof\[Chi]=If[#>v\[Chi]min ,("ne"/.params)\[Sigma]ofv\[Chi]and\[Omega]min[#,\[Omega]capture[#]],0]&; *)(*if v is greater than the lowest allowable velocity (set by \[Omega]min - the lowest allowable energy transfer), and the maximum possible energy transfer (set by kinetic energy) is greater than energy needed to capture (always true),
get the INVERSE interaction length / cross-section to transfer at least \[Omega]capture to the system and have energy less than the escape velocity*)

\[Lambda]invof\[Chi]
]


(* ::Input::Initialization:: *)
Clear[GetTotalNucWeakCaptureInteractionLength]
GetTotalNucWeakCaptureInteractionLength[\[Sigma]Dicts_]:=Module[{},
Total@Table[GetWeakCaptureInteractionLengthNuc[\[Sigma]Dicts[[i]]][#],{i,Length[\[Sigma]Dicts]}]&
]


(* ::Subsubsection::Closed:: *)
(*Weak Capture Regime Capture Probability - w Potential - e*)


(* ::Input::Initialization:: *)
Clear[GetWeakCaptureProbabilitywPot]
GetWeakCaptureProbabilitywPot[vDict_,\[Sigma]Dict_]:= Module[{v\[Chi]min,v\[Chi]raw,v\[Chi],\[Chi]traj,bE,m\[Chi],\[Omega]min,params,ldom,\[Kappa],regionindex,\[Sigma]dom,\[Omega]capture,\[Delta]\[Xi],\[Xi]of\[Omega]andl,\[Omega]maxofl,d\[Sigma]dERoftand\[Omega],\[Sigma]oftand\[Omega]min,rCore,rE,lcore,lcrust,\[Lambda]invoft,\[Lambda]integrand,integralof\[Lambda]invwrtr},
\[Omega]min="\[Omega]min"/.\[Sigma]Dict;
params=\[Sigma]Dict["params"];
\[Delta]\[Xi]=(10^"\[Xi]"/.\[Sigma]Dict)[[1]];
m\[Chi]=("m\[Chi]"/.vDict);
ldom = "domain"/.vDict;
\[Kappa]=("\[Kappa]"/.vDict);
regionindex = "regionindex"/.\[Sigma]Dict; (*1 in crust, 2 in core*)
\[Sigma]dom=("\[Sigma]of\[Xi]f"/.\[Sigma]Dict)["Domain"];

v\[Chi]min=Max[3Sqrt[("\[HBar]" \[Omega]min)/(2 m\[Chi])]/.params,"vesc"/.Constants`EarthRepl];(*minimum allowable velocity given the lower bound on energy transfer, times an O(1) fudge factor for numerical stability (3) *)
v\[Chi]raw=("v(l)"/.vDict);
v\[Chi]=If[v\[Chi]raw[#]>v\[Chi]min,("v(l)"/.vDict)[#],0]&;
(*\[Chi]traj=("l(t)"/.vDict);*)
bE=("bE"/.vDict);(*impact parameter*)

(*\[Xi]of\[Omega]andl=Abs[EnergyLoss`\[Xi]of\[Omega]andv\[Chi][#1,\[Omega]min,("v(l)"/.vDict)[#2],m\[Chi],params]-\[Delta]\[Xi]]&;*)
(*\[Xi]of\[Omega]andl=Abs[EnergyLoss`\[Xi]of\[Omega]andv\[Chi][#1,\[Omega]min,("v(l)"/.vDict)[#2],m\[Chi],params]]&;*)
\[Xi]of\[Omega]andl=EnergyLoss`\[Xi]of\[Omega]andv\[Chi][#1,\[Omega]min,("v(l)"/.vDict)[#2],m\[Chi],params]&;
\[Omega]maxofl=EnergyLoss`\[Omega]maxofv\[Chi][("v(l)"/.vDict)[#1],m\[Chi],params]&;

(*d\[Sigma]dERoftand\[Omega]=10^Re@("d\[Sigma]f"/.\[Sigma]Dict)[Log10@v\[Chi][#1],Log10@\[Xi]of\[Omega]andt[#2,#1]]&; (*not currently used*)*)

\[Sigma]oftand\[Omega]min=10^(Re[("\[Sigma]of\[Xi]f"/.\[Sigma]Dict)[Log10@v\[Chi][#1],Log10@\[Xi]of\[Omega]andl[#2,#1]]])&;

(*Post \[Omega]capture*)
\[Omega]capture = 1/(2 ("\[HBar]"/.params)) m\[Chi] (v\[Chi][#]^2 - ("vesc"/.Constants`EarthRepl)^2)&;


\[Lambda]invoft=If[v\[Chi]raw[#]>v\[Chi]min && \[Omega]maxofl[#]>\[Omega]capture[#]&&(\[Xi]of\[Omega]andl[\[Omega]capture[#],#])>10^\[Sigma]dom[[2,1]],\[Kappa]^2 ("ne"/.params)\[Sigma]oftand\[Omega]min[#,\[Omega]capture[#]],0]&;

rCore="rcore"/.Constants`EarthRepl;
rE = "rE"/.Constants`EarthRepl;

\[Lambda]integrand[l_?NumericQ]:=\[Lambda]invoft[l];


If[vDict["lMCs"]=={ "NC"},
(*Doesn't pass through the core, just integrate*)
If[regionindex==2,
(*material is in the core, return 0*)
Return[0,Module],
(*material is in the mantle, just integrate*)
integralof\[Lambda]invwrtr=NIntegrate[\[Lambda]integrand[l],{l,ldom[[1]],ldom[[2]]},PrecisionGoal->3,MaxRecursion->8,AccuracyGoal->3];
],
(*Passes through the core*)
If[regionindex==2,
(*material is in the core, integrate between the two Mantle-Core boundary points*)
integralof\[Lambda]invwrtr=NIntegrate[\[Lambda]integrand[l],{l,vDict["lMCs"][[1]],vDict["lMCs"][[2]]},PrecisionGoal->3,MaxRecursion->8,AccuracyGoal->3];,
(*material is in the mantle, integrate over the regions before and after the core, then add*)
integralof\[Lambda]invwrtr=NIntegrate[\[Lambda]integrand[l],{l,ldom[[1]],vDict["lMCs"][[1]]},PrecisionGoal->3,MaxRecursion->8,AccuracyGoal->3]+NIntegrate[\[Lambda]integrand[l],{l,vDict["lMCs"][[2]],ldom[[2]]},PrecisionGoal->3,MaxRecursion->8,AccuracyGoal->3];
]
];

integralof\[Lambda]invwrtr
]


(* ::Subsubsection::Closed:: *)
(*unit test*)


(* ::Input:: *)
(*temptrajdict = GetWeakCaptureTrajectoryinvRK2wPotential[10^9 ("JpereV")/("c")^2/.SIConstRepl,5 10^4,1 10^6,10^-10,"dEdlofr"/.GetTotalELFunction[Electronic\[Sigma]Dicts[[4]],Nuclear\[Sigma]Dicts[[4]]]];*)


(* ::Input:: *)
(*GetWeakCaptureProbabilitywPot[temptrajdict,Electronic\[Sigma]Dicts[[4,1]]]*)


(* ::Input:: *)
(*"materialname"/.Electronic\[Sigma]Dicts[[4,1]]*)


(* ::Input:: *)
(*Gather[temptrajdict["traj"],Sign@#1["vr"]==Sign@#2["vr"]&];*)
(*Table[Interpolation[{"l","r"}/.%[[i]]],{i,2}]*)
(**)


(* ::Input:: *)
(*InverseFunction[#^2&]*)


(* ::Input:: *)
(*temptrajdict["lMCs"]*)


(* ::Subsubsection::Closed:: *)
(*Weak Capture Regime Capture Probability - w Potential - Nuc*)


(* ::Input::Initialization:: *)
Clear[GetWeakCaptureProbabilityNucwPot]
GetWeakCaptureProbabilityNucwPot[vDict_,\[Sigma]Dict_]:= Module[{v\[Chi]min,v\[Chi]raw,v\[Chi],\[Chi]traj,bE,m\[Chi],\[Omega]min,params,NucleusParams,ldom,\[Kappa],regionindex,\[Sigma]dom,mN,\[Omega]capture,\[Delta]\[Xi],\[Xi]of\[Omega]andl,\[Omega]maxofl,d\[Sigma]dERoftand\[Omega],\[Sigma]oftand\[Omega]min,\[Lambda]invoft,\[Lambda]integrand,rCore,rE,lcore,lcrust,integralof\[Lambda]invwrtr},
\[Omega]min="\[Omega]min"/.\[Sigma]Dict;
params=\[Sigma]Dict["params"];
NucleusParams=\[Sigma]Dict["NucleusParams"];
\[Delta]\[Xi]=(10^"\[Xi]"/.\[Sigma]Dict)[[1]];
m\[Chi]=("m\[Chi]"/.vDict);
ldom = "domain"/.vDict;
\[Kappa]=("\[Kappa]"/.vDict);
mN = "mN"/.\[Sigma]Dict;
regionindex = "regionindex"/.\[Sigma]Dict;
\[Sigma]dom=("\[Sigma]of\[Xi]f"/.\[Sigma]Dict)["Domain"];

(*Print[mN];*)

v\[Chi]min=Max[3Sqrt[("\[HBar]" \[Omega]min)/(2 m\[Chi])]/.params,"vesc"/.Constants`EarthRepl];(*minimum allowable velocity given the lower bound on energy transfer, times an O(1) fudge factor for numerical stability (3) *)
v\[Chi]raw=("v(l)"/.vDict);
v\[Chi]=If[v\[Chi]raw[#]>v\[Chi]min,("v(l)"/.vDict)[#],0]&;

(*\[Chi]traj=("l(t)"/.vDict);*)
bE=("bE"/.vDict);(*impact parameter*)

(*\[Xi]of\[Omega]andl=Abs[EnergyLoss`\[Xi]of\[Omega]Nuc[#1,\[Omega]min,("v(l)"/.vDict)[#2],m\[Chi],mN,params]-\[Delta]\[Xi]]&;*)
(*\[Xi]of\[Omega]andl=Abs[EnergyLoss`\[Xi]of\[Omega]Nuc[#1,\[Omega]min,("v(l)"/.vDict)[#2],m\[Chi],mN,params]]&;*)
\[Xi]of\[Omega]andl=EnergyLoss`\[Xi]of\[Omega]Nuc[#1,\[Omega]min,("v(l)"/.vDict)[#2],m\[Chi],mN,params]&;
\[Omega]maxofl=EnergyLoss`\[Omega]maxNuc[("v(l)"/.vDict)[#1],m\[Chi],mN,params]&;

(*d\[Sigma]dERoftand\[Omega]=10^Re@("d\[Sigma]f"/.\[Sigma]Dict)[Log10@v\[Chi][#1],Log10@\[Xi]of\[Omega]andt[#2,#1]]&; (*not currently used*)*)

\[Sigma]oftand\[Omega]min=10^(Re[("\[Sigma]of\[Xi]f"/.\[Sigma]Dict)[Log10@v\[Chi][#1],Log10@\[Xi]of\[Omega]andl[#2,#1]]])&;


(*Post \[Omega]capture*)
\[Omega]capture = 1/(2 ("\[HBar]"/.params)) m\[Chi] (v\[Chi][#]^2 - ("vesc"/.Constants`EarthRepl)^2)&;


(*\[Lambda]invoft=If[v\[Chi]raw[#]>v\[Chi]min && \[Omega]maxofl[#]>\[Omega]capture[#]&&Log10@(\[Xi]of\[Omega]andl[\[Omega]capture[#],#])>\[Sigma]dom[[2,1]],\[Kappa]^2("nI"/.NucleusParams)\[Sigma]oftand\[Omega]min[#,\[Omega]capture[#]],0]&;*)
\[Lambda]invoft=If[v\[Chi]raw[#]>v\[Chi]min && \[Omega]maxofl[#]>\[Omega]capture[#]&&(\[Xi]of\[Omega]andl[\[Omega]capture[#],#])>10^\[Sigma]dom[[2,1]],\[Kappa]^2 ("nI"/.NucleusParams)\[Sigma]oftand\[Omega]min[#,\[Omega]capture[#]],0]&;

rCore="rcore"/.Constants`EarthRepl;
rE = "rE"/.Constants`EarthRepl;

\[Lambda]integrand[l_?NumericQ]:=\[Lambda]invoft[l];

If[vDict["lMCs"]=={ "NC"},
(*Doesn't pass through the core, just integrate*)
If[regionindex==2,
(*material is in the core, return 0*)
Return[0,Module],
(*material is in the mantle, just integrate*)
integralof\[Lambda]invwrtr=NIntegrate[\[Lambda]integrand[l],{l,ldom[[1]],ldom[[2]]},PrecisionGoal->3,MaxRecursion->8,AccuracyGoal->3];
],
(*Passes through the core*)
If[regionindex==2,
(*material is in the core, integrate between the two Mantle-Core boundary points*)
integralof\[Lambda]invwrtr=NIntegrate[\[Lambda]integrand[l],{l,vDict["lMCs"][[1]],vDict["lMCs"][[2]]},PrecisionGoal->3,MaxRecursion->8,AccuracyGoal->3];,
(*material is in the mantle, integrate over the regions before and after the core, then add*)
integralof\[Lambda]invwrtr=NIntegrate[\[Lambda]integrand[l],{l,ldom[[1]],vDict["lMCs"][[1]]},PrecisionGoal->3,MaxRecursion->8,AccuracyGoal->3]+NIntegrate[\[Lambda]integrand[l],{l,vDict["lMCs"][[2]],ldom[[2]]},PrecisionGoal->3,MaxRecursion->8,AccuracyGoal->3];
]
];

integralof\[Lambda]invwrtr
]


(* ::Subsection::Closed:: *)
(*Initial Conditions*)


(* ::Subsubsection::Closed:: *)
(*vmax*)


(* ::Input::Initialization:: *)
(*Clear[Getv\[Chi]Max]
Getv\[Chi]Max[Pcut_,\[Beta]D_,m\[Chi]_]:=Module[{oneminusCDFvmax,vesc},
vesc=("vesc"/.Constants`EarthRepl);
oneminusCDFvmax = \[ExponentialE]^(1/2 m\[Chi] vesc^2 \[Beta]D) (1+2 \[ExponentialE]^(-(1/2) m\[Chi] vmax^2 \[Beta]D) Sqrt[(m\[Chi] \[Beta]D)/(2 \[Pi])]vmax- Erf[Sqrt[(m\[Chi] \[Beta]D)/(2 \[Pi])]vmax] );(*Integral of MB dist, from Subscript[v, max] to \[Infinity] - so probability of v >= Subscript[v, max]*)
vmax/.NSolve[{Pcut==oneminusCDFvmax&&0<vmax<10^7},vmax][[1]](*make sure it's non-relativistic*)
]*)


(* ::Input::Initialization:: *)
Clear[Getv\[Chi]Max]
Getv\[Chi]Max[Pcut_,\[Beta]D_,m\[Chi]_]:=Module[{oneminusCDFvmax,vesc,vmaxsol,vmaxbound=0.2 "c"/.SIConstRepl},
vesc=("vesc"/.Constants`EarthRepl);
oneminusCDFvmax = E^(1/2 m\[Chi] vesc^2 \[Beta]D) (1+2 E^(-(1/2) m\[Chi] vmax^2 \[Beta]D) Sqrt[(m\[Chi] \[Beta]D)/(2 \[Pi])]vmax- Erf[Sqrt[(m\[Chi] \[Beta]D)/(2 \[Pi])]vmax] );

(*Integral of MB dist, from Subscript[v, max] to \[Infinity] - so probability of v >= Subscript[v, max]*)
vmaxsol=NSolve[{Pcut==oneminusCDFvmax&&0<vmax<vmaxbound},vmax];

If[(vmaxsol)!={},vmax/.vmaxsol[[1]],vmaxbound](*make sure it's non-relativistic and solve for velocity corresponding to the Pcut*)
]


(* ::Input:: *)
(*Integrate[ E^(-(1/2) m (v^2-vesc^2) \[Beta]) Sqrt[2/\[Pi]] v^2 (m \[Beta])^(3/2) ,{v,vmax,\[Infinity]}]*)
(*%//Normal//PowerExpand//Simplify*)


(* ::Subsubsection::Closed:: *)
(*Get v\[Chi] grid at earth*)


(* ::Input::Initialization:: *)
(*Clear[Getv\[Chi]atEarth]
Getv\[Chi]atEarth[m\[Chi]_,\[Beta]D_,v\[Chi]Max_,V_:Vgrav,nbs_:10,nv\[Chi]s_:40,lowerbrat_:10^-6,lowervrat_:10^-6]:=Module[{rE,vmin,v\[Chi]s,bs,bEarth,vrEarth,v\[Chi],\[Chi]speedEarth},

    rE = "rE"/.Constants`EarthRepl;
    vmin = "vesc"/.Constants`EarthRepl;
   
v\[Chi]s = vmin+ 10^Subdivide[ Log10[lowervrat(v\[Chi]Max-vmin)],Log10@(v\[Chi]Max-vmin),nv\[Chi]s];

bs = Subdivide[lowerbrat rE,rE,nbs];

Table[<|"\[Chi]speedEarth"->v\[Chi]s[[i]],"bEarth"->bs[[j]]|>,{i,nv\[Chi]s},{j,nbs}]
]*)


(* ::Input::Initialization:: *)
(*Clear[Getv\[Chi]atEarth]
Getv\[Chi]atEarth[m\[Chi]_,\[Beta]D_,v\[Chi]Max_,V_:Vgrav,nbs_:10,nv\[Chi]s_:40,lowerbrat_:10^-6,lowervrat_:10^-8]:=Module[{rE,vmin,v\[Chi]s,bs,bEarth,vrEarth,v\[Chi],\[Chi]speedEarth},

    rE = "rE"/.Constants`EarthRepl;
    vmin = "vesc"/.Constants`EarthRepl;
   
v\[Chi]s = vmin+ 10^Subdivide[ Log10[lowervrat(v\[Chi]Max-vmin)],Log10@(v\[Chi]Max-vmin),nv\[Chi]s];

bs = Subdivide[lowerbrat rE,rE,nbs];

Table[<|"\[Chi]speedEarth"->v\[Chi]s[[i]],"bEarth"->bs[[j]]|>,{i,nv\[Chi]s},{j,nbs}]
]*)


(* ::Input::Initialization:: *)
Clear[Getv\[Chi]atEarth]
Getv\[Chi]atEarth[v\[Chi]Max_,nbs_:10,nv\[Chi]s_:40,lowerbrat_:10^-6,lowervrat_:10^-8]:=Module[{rE,vmin,v\[Chi]s,bs,bEarth,vrEarth,v\[Chi],\[Chi]speedEarth},

    rE = "rE"/.Constants`EarthRepl;
    vmin = "vesc"/.Constants`EarthRepl;
   
v\[Chi]s = vmin+ 10^Subdivide[ Log10[lowervrat(v\[Chi]Max-vmin)],Log10@(v\[Chi]Max-vmin),nv\[Chi]s];

bs = Subdivide[lowerbrat rE,rE,nbs];

Table[<|"\[Chi]speedEarth"->v\[Chi]s[[i]],"bEarth"->bs[[j]]|>,{i,nv\[Chi]s},{j,nbs}]
]


(* ::Subsection::Closed:: *)
(*Get Capture Probability*)


(* ::Subsubsection::Closed:: *)
(*Get Capture Rate old*)


(* ::Input::Initialization:: *)
(*Clear[GetCaptureRateOldTraj]
GetCaptureRateOldTraj[meDovermpD_,\[Kappa]_,v0_,\[Sigma]DictsElectronic_,\[Sigma]DictsNuclear_]:= Module[{mpD,meD,\[Beta]D,v0pD,v\[Chi]Max,ICdict,ELfuncTotalby\[Kappa]squared,ELthroughEarthby\[Kappa]squared,\[CapitalDelta]rby\[Kappa]squared,RE,TrajDict,vfinal,captured,integralof\[Lambda]invwrtrE ,integralof\[Lambda]invwrtrNuc,SurvivalProb,CaptureProb,PDictList={},TrajDictinv},

(*logname = NotebookDirectory[]<>"capturelog.txt";
logstrm = OpenAppend[logname];
Print[logname];*)

(*Write[logstrm,"init"];*)

mpD=Join["m\[Chi]"/.\[Sigma]DictsElectronic,"m\[Chi]"/.\[Sigma]DictsNuclear];
mpD=If[Length@Union[mpD]==1,mpD[[1]],Print["m\[Chi] must be the same in all processes:",mpD];Return[<||>,Module]];
meD = meDovermpD mpD;

(*Print@v0to\[Beta]D[v0,meD,mpD];
*)
\[Beta]D ="\[Beta]D"/.v0to\[Beta]D[ v0,meD,mpD];

v\[Chi]Max=Getv\[Chi]Max[10^-6,\[Beta]D,mpD];
(*Print[v\[Chi]Max];*)

(*Write[logstrm,"post init"];*)

(*Print[v0];
Print[meDovermpD];*)

ICdict=Flatten@Getv\[Chi]atEarth[mpD,\[Beta]D,v\[Chi]Max];(*dictionary of initial conditions at earth*)

(*Print[ICdict];*)

(*Write[logstrm,"post ICdict"];*)

(*Print["\[Chi]speedEarth"/.ICdict];*)

{ELfuncTotalby\[Kappa]squared,ELthroughEarthby\[Kappa]squared} = {"dEdlofr","\[CapitalDelta]EthroughEarth"}/.GetTotalELFunction[\[Sigma]DictsElectronic,\[Sigma]DictsNuclear];

(*Write[logstrm,"post EL func"];*)

(*Print[ICdict];*)
(*Print[ELfuncTotalby\[Kappa]squared[10^6,10^6]];
Print[ELthroughEarthby\[Kappa]squared[10^6,10^6]];*)

(*Print[ICdict];*)
(*Print@fTot0ELCore[Log10[mpD],Log10[ "\[Chi]speedEarth"]/.ICdict];*)
(*Print@fTot0ELCore[Log10[mpD],1+Log10[ "\[Chi]speedEarth"]/.ICdict];*)
(*Print@fTot0ELCore[Log10[mpD],2+Log10[ "\[Chi]speedEarth"]/.ICdict];
(*Print@dM["bEarth"/.ICdict];*)
Print@dc["bEarth"/.ICdict];
Print[("rcore"/.Constants`EarthRepl)>("bEarth"/.ICdict)^2];*)
(*"JpereV"\[CapitalDelta]ETot0[mpD,"\[Chi]speedEarth"/.ICdict,\[Kappa],"bEarth"/.ICdict,Constants`SIConstRepl]/.Constants`SIConstRepl*)
Monitor[
Do[
(*\[CapitalDelta]rby\[Kappa]squared=\[CapitalDelta]rAveragedTot0[mpD,"\[Chi]speedEarth"/.ICdict[[n]],1,"bEarth"/.ICdict[[n]]];*)
\[CapitalDelta]rby\[Kappa]squared=\[CapitalDelta]rTotby\[Kappa]squared[mpD,"\[Chi]speedEarth"/.ICdict[[n]],"bEarth"/.ICdict[[n]],ELthroughEarthby\[Kappa]squared];
RE=dE["bEarth"/.ICdict[[n]]];

(*Write[logstrm,"post delta r"];*)

(*Print["\[CapitalDelta]r: ",\[CapitalDelta]rby\[Kappa]squared];*)
(*Print[\[CapitalDelta]rby\[Kappa]squared];*)

(*Print@\[CapitalDelta]rby\[Kappa]squared;*)
(*Print@(\[CapitalDelta]rby\[Kappa]squared/RE);*)
(*\[Kappa]StrongWeakBoundary=\[Kappa]t/.Solve[\[CapitalDelta]rby\[Kappa]squared /\[Kappa]t^2==RE/100,\[Kappa]t][[2]];*)
(*Print[ELEarthTot[r,mpD,"\[Chi]speedEarth"/.ICdict]/.r->1/2"rE"/.Constants`EarthRepl];*)
(*Plot[(ELEarthTot[#,mpD,"\[Chi]speedEarth"/.ICdict])&[r],{r,1,"rE"/.Constants`EarthRepl}];*)

(*If[\[CapitalDelta]rby\[Kappa]squared /\[Kappa]^2>RE/100,*)
(*Print["LHS and RHS of \[CapitalDelta]r inequality"];
Print[\[CapitalDelta]rby\[Kappa]squared /\[Kappa]^2];
Print[RE/100000];*)

If[\[CapitalDelta]rby\[Kappa]squared /\[Kappa]^2>RE/100,
(*on average will loose 100% of kinetic energy as it travels through the earth, if speed is fixed*)

TrajDict=GetWeakCaptureTrajectory[mpD,"\[Chi]speedEarth"/.ICdict[[n]],"bEarth"/.ICdict[[n]],\[Kappa],ELfuncTotalby\[Kappa]squared ];
vfinal= ("v(t)"/.TrajDict)[TrajDict["domain"][[2]]];
captured=vfinal<"vesc"/.Constants`EarthRepl;

(*Print[TrajDict];
Print@Plot[("v(t)"/.TrajDict)[t],{t,TrajDict["domain"][[1]],TrajDict["domain"][[2]]},PlotLabel->"v(t)"];
*)




(*TrajDict=GetWeakCaptureTrajectoryinvFE[mpD,"\[Chi]speedEarth"/.ICdict[[n]],"bEarth"/.ICdict[[n]],\[Kappa],ELfuncTotalby\[Kappa]squared ];*)
(*TrajDict=GetWeakCaptureTrajectoryinv[mpD,"\[Chi]speedEarth"/.ICdict[[n]],"bEarth"/.ICdict[[n]],\[Kappa],ELfuncTotalby\[Kappa]squared ];*)
(*vfinal= ("v(l)"/.TrajDict)[TrajDict["domain"][[2]]];
captured=vfinal<"vesc"/.Constants`EarthRepl;*)
(*Write[logstrm,"post cap traj"];*)

(*Print[TrajDictinv];
Print@Plot[("v(l)"/.TrajDictinv)[l],{l,TrajDictinv["domain"][[1]],TrajDictinv["domain"][[2]]},PlotLabel->"v(l)"];
*)


If[!captured,
(*allow continuous slowing down, escapes with v>vesc, check if can be captured through hard scatter*)

(*Write[logstrm,"not continuous captured"];*)

(*integralof\[Lambda]invwrtrE = EchoTiming[Total@Table[GetWeakCaptureProbability[TrajDict,\[Sigma]DictsElectronic[[i]]],{i,Length[\[Sigma]DictsElectronic]}],"E oldalg: "];

integralof\[Lambda]invwrtrNuc = EchoTiming[Total@Table[GetWeakCaptureProbabilityNuc[TrajDict,\[Sigma]DictsNuclear[[i]]],{i,Length[\[Sigma]DictsNuclear]}],"Nuc oldalg: "];*)(*uses l(t) solution instead of v(l)*)
integralof\[Lambda]invwrtrE = Total@Table[GetWeakCaptureProbability[TrajDict,\[Sigma]DictsElectronic[[i]]],{i,Length[\[Sigma]DictsElectronic]}];

integralof\[Lambda]invwrtrNuc = Total@Table[GetWeakCaptureProbabilityNuc[TrajDict,\[Sigma]DictsNuclear[[i]]],{i,Length[\[Sigma]DictsNuclear]}];
(*Print["E capture integrals of \[Lambda]"];
Print[integralof\[Lambda]invwrtrE];

Print@Total@Table[GetWeakCaptureProbabilityfromv[TrajDictinv,\[Sigma]DictsElectronic[[i]]],{i,Length[\[Sigma]DictsElectronic]}];

Print["Nuc capture integrals of \[Lambda]"];
Print[integralof\[Lambda]invwrtrNuc];
Print@Total@Table[GetWeakCaptureProbabilityNucfromv[TrajDictinv,\[Sigma]DictsNuclear[[i]]],{i,Length[\[Sigma]DictsNuclear]}];*)

(*integralof\[Lambda]invwrtrE = Total@Table[GetWeakCaptureProbabilityfromv[TrajDict,\[Sigma]DictsElectronic[[i]]],{i,Length[\[Sigma]DictsElectronic]}];

integralof\[Lambda]invwrtrNuc = Total@Table[GetWeakCaptureProbabilityNucfromv[TrajDict,\[Sigma]DictsNuclear[[i]]],{i,Length[\[Sigma]DictsNuclear]}];*)

(*Write[logstrm,"post integrals of \[Lambda]"];

Quiet@Clear[v\[Chi]min,v\[Chi]raw,v\[Chi],\[Chi]traj,bE,m\[Chi],\[Omega]min,params,tdom,regionindex,\[Omega]capture,\[Delta]\[Xi],\[Xi]of\[Omega]andt,\[Omega]maxoft,d\[Sigma]dERoftand\[Omega],\[Sigma]oftand\[Omega]min,\[Lambda]invoft,integralof\[Lambda]invwrtr];
Quiet@Clear[v\[Chi]min,v\[Chi]raw,v\[Chi],\[Chi]traj,bE,m\[Chi],\[Omega]min,params,NucleusParams,tdom,regionindex,mN,\[Omega]capture,\[Delta]\[Xi],\[Xi]of\[Omega]andt,\[Omega]maxoft,d\[Sigma]dERoftand\[Omega],\[Sigma]oftand\[Omega]min,\[Lambda]invoft,integralof\[Lambda]invwrtr];*)

(*put get NcMax here*)

SurvivalProb =E^-(integralof\[Lambda]invwrtrE + integralof\[Lambda]invwrtrNuc);
CaptureProb=N[1-SurvivalProb,5];
AppendTo[PDictList,<|"Psurv"->N[SurvivalProb,5],"Pcap"->Max[CaptureProb,10^-50],"CSDcap"->captured,"v\[Chi]E"->"\[Chi]speedEarth"/.ICdict[[n]],"bEarth"->("bEarth"/.ICdict[[n]]),"vfinal"->vfinal,"\[Kappa]"->\[Kappa],"mpD"->mpD,"meD"->meD,"v\[Chi]Max"->v\[Chi]Max,"\[Beta]D"->\[Beta]D,"v0"->v0,"y"->\[CapitalDelta]rby\[Kappa]squared /(\[Kappa]^2 RE),"RE"->RE,"PcapE"->Max[N[1-E^-integralof\[Lambda]invwrtrE,5],10^-100],"PcapNuc"->Max[N[1-E^-integralof\[Lambda]invwrtrNuc,5],10^-100],"int\[Lambda]invE"->Log10@integralof\[Lambda]invwrtrE,"int\[Lambda]invNuc"->Log10@integralof\[Lambda]invwrtrNuc|>];,


(*allowed continuous slowing down, escapes with v<vesc, it's captured*)
(*Write[logstrm,"soft captured"];*)
AppendTo[PDictList,<|"Psurv"->0,"Pcap"->1,"CSDcap"->True,"v\[Chi]E"->"\[Chi]speedEarth"/.ICdict[[n]],"vfinal"->vfinal,"bEarth"->("bEarth"/.ICdict[[n]]),"\[Kappa]"->\[Kappa],"mpD"->mpD,"meD"->meD,"v\[Chi]Max"->v\[Chi]Max,"\[Beta]D"->\[Beta]D,"v0"->v0,"y"->\[CapitalDelta]rby\[Kappa]squared /(\[Kappa]^2 RE),"RE"->RE,"PcapE"->1,"PcapNuc"->1,"int\[Lambda]invE"->100,"int\[Lambda]invNuc"->100|>]
];
,
(*allowed call it captured since on average KE is lost as it traverses the earth*)
(*Write[logstrm,"soft captured (delta r only)"];*)
AppendTo[PDictList,<|"Psurv"->0,"Pcap"->1,"CSDcap"->True,"v\[Chi]E"->"\[Chi]speedEarth"/.ICdict[[n]],"vfinal"->0,"bEarth"->("bEarth"/.ICdict[[n]]),"\[Kappa]"->\[Kappa],"mpD"->mpD,"meD"->meD,"v\[Chi]Max"->v\[Chi]Max,"\[Beta]D"->\[Beta]D,"v0"->v0,"y"->\[CapitalDelta]rby\[Kappa]squared /(\[Kappa]^2 RE),"RE"->RE,"PcapE"->1,"PcapNuc"->1,"int\[Lambda]invE"->100,"int\[Lambda]invNuc"->100|>]
];
(*Print[PDictList[[-1]]];*)
,{n,Length[ICdict]}
];
,n];

(*Close[logstrm];*)

PDictList
]*)


(* ::Subsubsection::Closed:: *)
(*Get Probability*)


(* ::Input::Initialization:: *)
Clear[GetCaptureProbability]
GetCaptureProbability[meDovermpD_,\[Kappa]_,v0_,\[Sigma]DictsElectronic_,\[Sigma]DictsNuclear_,mpDcap_:True]:= Module[{m\[Chi],mpD,meD,\[Beta]D,v0pD,v\[Chi]Max,ICdict,ELfuncTotalby\[Kappa]squared,ELthroughEarthby\[Kappa]squared,\[CapitalDelta]rby\[Kappa]squared,RE,TrajDict,vfinal,captured,integralof\[Lambda]invwrtrE ,integralof\[Lambda]invwrtrNuc,SurvivalProb,CaptureProb,PDictList={},TrajDictinv,MaxBoltzSuppress=10^-12},

(*Get the capture probability for given microphysics parameters

	mpDcap - [boole] True=> compute capture of dark protons, False => compute capture of dark electrons*)

m\[Chi]=Join["m\[Chi]"/.\[Sigma]DictsElectronic,"m\[Chi]"/.\[Sigma]DictsNuclear];
m\[Chi]=If[Length@Union[m\[Chi]]==1,m\[Chi][[1]],Print["m\[Chi] must be the same in all processes:",m\[Chi]];Return[<||>,Module]];

If[mpDcap,
mpD = m\[Chi];meD = meDovermpD m\[Chi];,
meD = m\[Chi];mpD =  m\[Chi]/meDovermpD;
];

(*what we can do to change from p->e is introduce a flag for mpD or meD. Then here we use the \[Sigma]Dict for e capture. Then all we need to do is relabel the variables so mpD -> meD, (ie. meD = m\[Chi], mpD = m\[Chi]/meDovermpD, with m\[Chi] the mass used in the \[Sigma]Dicts). We also need to pass this to v\[Chi] at Earth and v\[Chi]Max (can just pass meD instead of mpD). This will also be passed to Get weak capture trajectory*)

\[Beta]D ="\[Beta]D"/.v0to\[Beta]D[ v0,meD,mpD];(*for fixed m\[Chi], this is the only thing that will change for different values of the mpDcap flag*)

v\[Chi]Max=Getv\[Chi]Max[MaxBoltzSuppress,\[Beta]D,m\[Chi]]; (*from the above comment, the only difference between pD and eD capture will be velocity range we look at*)

 ICdict=Flatten@Getv\[Chi]atEarth[v\[Chi]Max];(*dictionary of initial conditions at earth*) 

{ELfuncTotalby\[Kappa]squared,ELthroughEarthby\[Kappa]squared} = {"dEdlofr","\[CapitalDelta]EthroughEarth"}/.GetTotalELFunction[\[Sigma]DictsElectronic,\[Sigma]DictsNuclear];

(*Print@ELthroughEarthby\[Kappa]squared;*)

(*Monitor[*)

Do[
\[CapitalDelta]rby\[Kappa]squared=\[CapitalDelta]rTotby\[Kappa]squared[m\[Chi],"\[Chi]speedEarth"/.ICdict[[n]],"bEarth"/.ICdict[[n]],ELthroughEarthby\[Kappa]squared]; 


RE=dE["bEarth"/.ICdict[[n]]];

If[\[CapitalDelta]rby\[Kappa]squared /\[Kappa]^2>RE/100,
(*on average will loose 100% of kinetic energy as it travels through the earth, if speed is fixed*)

TrajDict=GetWeakCaptureTrajectoryinvRK2wPotential[m\[Chi],"\[Chi]speedEarth"/.ICdict[[n]],"bEarth"/.ICdict[[n]],\[Kappa],ELfuncTotalby\[Kappa]squared ]; 
(*TrajDict=GetWeakCaptureTrajectoryinvFE[m\[Chi],"\[Chi]speedEarth"/.ICdict[[n]],"bEarth"/.ICdict[[n]],\[Kappa],ELfuncTotalby\[Kappa]squared ];*) 
(*TrajDict=GetWeakCaptureTrajectoryinv[m\[Chi],"\[Chi]speedEarth"/.ICdict[[n]],"bEarth"/.ICdict[[n]],\[Kappa],ELfuncTotalby\[Kappa]squared ]; *)

(*vfinal= ("v(l)"/.TrajDict)[TrajDict["domain"][[2]]];*)
(*captured=vfinal<"vesc"/.Constants`EarthRepl;*)
(*captured=vfinal<"vesc"/.TrajDict;*)
vfinal= "vfinal"/.TrajDict;
captured = ("rfinal"/.TrajDict)<("rE"/.Constants`EarthRepl);(*stop condition is bounded motion (E<0) so if we exit the Earth the motion is unbounded, otherwise its bounded.*)

If[!captured,
(*allow continuous slowing down, escapes with v>vesc, check if can be captured through hard scatter*)

(*integralof\[Lambda]invwrtrE = Total@Table[GetWeakCaptureProbabilityfromv[TrajDict,\[Sigma]DictsElectronic[[i]]],{i,Length[\[Sigma]DictsElectronic]}];
integralof\[Lambda]invwrtrNuc = Total@Table[GetWeakCaptureProbabilityNucfromv[TrajDict,\[Sigma]DictsNuclear[[i]]],{i,Length[\[Sigma]DictsNuclear]}];*)
integralof\[Lambda]invwrtrE = Total@Table[GetWeakCaptureProbabilitywPot[TrajDict,\[Sigma]DictsElectronic[[i]]],{i,Length[\[Sigma]DictsElectronic]}];
integralof\[Lambda]invwrtrNuc = Total@Table[GetWeakCaptureProbabilityNucwPot[TrajDict,\[Sigma]DictsNuclear[[i]]],{i,Length[\[Sigma]DictsNuclear]}];

SurvivalProb =If[integralof\[Lambda]invwrtrE + integralof\[Lambda]invwrtrNuc>250,0,N[E^-(integralof\[Lambda]invwrtrE + integralof\[Lambda]invwrtrNuc),5]];
CaptureProb=If[integralof\[Lambda]invwrtrE + integralof\[Lambda]invwrtrNuc>250,1,If[N[integralof\[Lambda]invwrtrE + integralof\[Lambda]invwrtrNuc,5]<0.001,N[integralof\[Lambda]invwrtrE + integralof\[Lambda]invwrtrNuc,5],N[1-E^-(integralof\[Lambda]invwrtrE + integralof\[Lambda]invwrtrNuc),5]]];

AppendTo[PDictList,<|"Psurv"->N[SurvivalProb,5],"Pcap"->Max[CaptureProb,10^-100],"CSDcap"->Boole@captured,"v\[Chi]E"->"\[Chi]speedEarth"/.ICdict[[n]],"bEarth"->("bEarth"/.ICdict[[n]]),"vfinal"->vfinal,"\[Kappa]"->\[Kappa],"mpD"->mpD,"meD"->meD,"m\[Chi]"->m\[Chi],"mpDcap"->mpDcap,"v\[Chi]Max"->v\[Chi]Max,"\[Beta]D"->\[Beta]D,"v0"->v0,"y"->\[CapitalDelta]rby\[Kappa]squared /(\[Kappa]^2 RE),"RE"->RE,"PcapE"->Max[N[integralof\[Lambda]invwrtrE,5],10^-100],"PcapNuc"->Max[N[integralof\[Lambda]invwrtrNuc,5],10^-100],"int\[Lambda]invE"->integralof\[Lambda]invwrtrE,"int\[Lambda]invNuc"->integralof\[Lambda]invwrtrNuc|>];,


(*allowed continuous slowing down, v<vesc, it's captured*)
AppendTo[PDictList,<|"Psurv"->0,"Pcap"->1,"CSDcap"->Boole@True,"v\[Chi]E"->"\[Chi]speedEarth"/.ICdict[[n]],"vfinal"->vfinal,"bEarth"->("bEarth"/.ICdict[[n]]),"\[Kappa]"->\[Kappa],"mpD"->mpD,"meD"->meD,"m\[Chi]"->m\[Chi],"mpDcap"->mpDcap,"v\[Chi]Max"->v\[Chi]Max,"\[Beta]D"->\[Beta]D,"v0"->v0,"y"->\[CapitalDelta]rby\[Kappa]squared /(\[Kappa]^2 RE),"RE"->RE,"PcapE"->1,"PcapNuc"->1,"int\[Lambda]invE"->10^-100,"int\[Lambda]invNuc"->10^-100|>](*captured by CSD so take inverse capture length to Subscript[r, E] so that it's integral over the trajectory is ~1, If we are CSD captured, we want to set the inverse hard capture interaction length to be small so that we can use it to distinguish between hard and soft capture*)
];
,
(*call it captured since on average more than 100 percent of KE is lost as it traverses the earth*)
AppendTo[PDictList,<|"Psurv"->0,"Pcap"->1,"CSDcap"->Boole@True,"v\[Chi]E"->"\[Chi]speedEarth"/.ICdict[[n]],"vfinal"->0,"bEarth"->("bEarth"/.ICdict[[n]]),"\[Kappa]"->\[Kappa],"mpD"->mpD,"meD"->meD,"m\[Chi]"->m\[Chi],"mpDcap"->mpDcap,"v\[Chi]Max"->v\[Chi]Max,"\[Beta]D"->\[Beta]D,"v0"->v0,"y"->\[CapitalDelta]rby\[Kappa]squared /(\[Kappa]^2 RE),"RE"->RE,"PcapE"->1,"PcapNuc"->1,"int\[Lambda]invE"->10^-100,"int\[Lambda]invNuc"->10^-100|>]
];
,{n,Length[ICdict]}
];
(*,n];*)

PDictList
]


(* ::Input::Initialization:: *)
Clear[IterateGetCaptureProbability]
IterateGetCaptureProbability[meDovermpD_,\[Kappa]s_,v0_,\[Sigma]DictsElectronic_,\[Sigma]DictsNuclear_,mpDcap_:True]:= Module[{PDictList={}},

(*Monitor[*)
Do[(*Print["\[Kappa] is: ",N@\[Kappa]s[[\[Kappa]]]];*)
PDictList=Join[PDictList,GetCaptureProbability[meDovermpD,\[Kappa]s[[\[Kappa]]],v0,\[Sigma]DictsElectronic,\[Sigma]DictsNuclear,mpDcap]]; (*flag for p->e*)
,{\[Kappa],Length[\[Kappa]s]}];
(*,\[Kappa]];*)

PDictList
]


(* ::Subsubsection::Closed:: *)
(*Get Capture Probability - mpD only*)


(* ::Input::Initialization:: *)
Clear[GetCaptureRate]
GetCaptureRate[meDovermpD_,\[Kappa]_,v0_,\[Sigma]DictsElectronic_,\[Sigma]DictsNuclear_,mpDcap_:True]:= Module[{mpD,meD,\[Beta]D,v0pD,v\[Chi]Max,ICdict,ELfuncTotalby\[Kappa]squared,ELthroughEarthby\[Kappa]squared,\[CapitalDelta]rby\[Kappa]squared,RE,TrajDict,vfinal,captured,integralof\[Lambda]invwrtrE ,integralof\[Lambda]invwrtrNuc,SurvivalProb,CaptureProb,PDictList={},TrajDictinv,MaxBoltzSuppress=10^-12},

mpD=Join["m\[Chi]"/.\[Sigma]DictsElectronic,"m\[Chi]"/.\[Sigma]DictsNuclear];
mpD=If[Length@Union[mpD]==1,mpD[[1]],Print["m\[Chi] must be the same in all processes:",mpD];Return[<||>,Module]];
meD = meDovermpD mpD;(*flag for p->e*)

(*what we can do to change from p->e is introduce a flag for mpD or meD. Then here we use the \[Sigma]Dict for e capture. Then all we need to do is relabel the variables so mpD -> meD, (ie. meD = m\[Chi], mpD = m\[Chi]/meDovermpD, with m\[Chi] the mass used in the \[Sigma]Dicts). We also need to pass this to v\[Chi] at Earth and v\[Chi]Max (can just pass meD instead of mpD). This will also be passed to Get weak capture trajectory*)

\[Beta]D ="\[Beta]D"/.v0to\[Beta]D[ v0,meD,mpD];

v\[Chi]Max=Getv\[Chi]Max[MaxBoltzSuppress,\[Beta]D,mpD]; (*flag for p->e*)

(*ICdict=Flatten@Getv\[Chi]atEarth[mpD,\[Beta]D,v\[Chi]Max];*)
ICdict=Flatten@Getv\[Chi]atEarth[v\[Chi]Max];(*dictionary of initial conditions at earth*) (*flag for p->e*)

{ELfuncTotalby\[Kappa]squared,ELthroughEarthby\[Kappa]squared} = {"dEdlofr","\[CapitalDelta]EthroughEarth"}/.GetTotalELFunction[\[Sigma]DictsElectronic,\[Sigma]DictsNuclear];

Monitor[
Do[
\[CapitalDelta]rby\[Kappa]squared=\[CapitalDelta]rTotby\[Kappa]squared[mpD,"\[Chi]speedEarth"/.ICdict[[n]],"bEarth"/.ICdict[[n]],ELthroughEarthby\[Kappa]squared]; (*flag for p->e*)
RE=dE["bEarth"/.ICdict[[n]]];

If[\[CapitalDelta]rby\[Kappa]squared /\[Kappa]^2>RE/100,
(*on average will loose 100% of kinetic energy as it travels through the earth, if speed is fixed*)

TrajDict=GetWeakCaptureTrajectoryinvFE[mpD,"\[Chi]speedEarth"/.ICdict[[n]],"bEarth"/.ICdict[[n]],\[Kappa],ELfuncTotalby\[Kappa]squared ]; (*flag for p->e*)

vfinal= ("v(l)"/.TrajDict)[TrajDict["domain"][[2]]];
captured=vfinal<"vesc"/.Constants`EarthRepl;

If[!captured,
(*allow continuous slowing down, escapes with v>vesc, check if can be captured through hard scatter*)

integralof\[Lambda]invwrtrE = Total@Table[GetWeakCaptureProbabilityfromv[TrajDict,\[Sigma]DictsElectronic[[i]]],{i,Length[\[Sigma]DictsElectronic]}];

integralof\[Lambda]invwrtrNuc = Total@Table[GetWeakCaptureProbabilityNucfromv[TrajDict,\[Sigma]DictsNuclear[[i]]],{i,Length[\[Sigma]DictsNuclear]}];

SurvivalProb =E^-(integralof\[Lambda]invwrtrE + integralof\[Lambda]invwrtrNuc);
CaptureProb=N[integralof\[Lambda]invwrtrE + integralof\[Lambda]invwrtrNuc,5];
AppendTo[PDictList,<|"Psurv"->N[SurvivalProb,5],"Pcap"->Max[CaptureProb,10^-100],"CSDcap"->Boole@captured,"v\[Chi]E"->"\[Chi]speedEarth"/.ICdict[[n]],"bEarth"->("bEarth"/.ICdict[[n]]),"vfinal"->vfinal,"\[Kappa]"->\[Kappa],"mpD"->mpD,"meD"->meD,"v\[Chi]Max"->v\[Chi]Max,"\[Beta]D"->\[Beta]D,"v0"->v0,"y"->\[CapitalDelta]rby\[Kappa]squared /(\[Kappa]^2 RE),"RE"->RE,"PcapE"->Max[N[integralof\[Lambda]invwrtrE,5],10^-100],"PcapNuc"->Max[N[integralof\[Lambda]invwrtrNuc,5],10^-100],"int\[Lambda]invE"->integralof\[Lambda]invwrtrE,"int\[Lambda]invNuc"->integralof\[Lambda]invwrtrNuc|>];,


(*allowed continuous slowing down, escapes with v<vesc, it's captured*)
AppendTo[PDictList,<|"Psurv"->0,"Pcap"->1,"CSDcap"->Boole@True,"v\[Chi]E"->"\[Chi]speedEarth"/.ICdict[[n]],"vfinal"->vfinal,"bEarth"->("bEarth"/.ICdict[[n]]),"\[Kappa]"->\[Kappa],"mpD"->mpD,"meD"->meD,"v\[Chi]Max"->v\[Chi]Max,"\[Beta]D"->\[Beta]D,"v0"->v0,"y"->\[CapitalDelta]rby\[Kappa]squared /(\[Kappa]^2 RE),"RE"->RE,"PcapE"->1,"PcapNuc"->1,"int\[Lambda]invE"->10^-100,"int\[Lambda]invNuc"->10^-100|>](*captured by CSD so take inverse capture length to Subscript[r, E] so that it's integral over the trajectory is ~1, If we are CSD captured, we want to set the inverse hard capture interaction length to be small so that we can use it to distinguish between hard and soft capture*)
];
,
(*call it captured since on average more than 100 percent of KE is lost as it traverses the earth*)
AppendTo[PDictList,<|"Psurv"->0,"Pcap"->1,"CSDcap"->Boole@True,"v\[Chi]E"->"\[Chi]speedEarth"/.ICdict[[n]],"vfinal"->0,"bEarth"->("bEarth"/.ICdict[[n]]),"\[Kappa]"->\[Kappa],"mpD"->mpD,"meD"->meD,"v\[Chi]Max"->v\[Chi]Max,"\[Beta]D"->\[Beta]D,"v0"->v0,"y"->\[CapitalDelta]rby\[Kappa]squared /(\[Kappa]^2 RE),"RE"->RE,"PcapE"->1,"PcapNuc"->1,"int\[Lambda]invE"->10^-100,"int\[Lambda]invNuc"->10^-100|>]
];
,{n,Length[ICdict]}
];
,n];

PDictList
]


(* ::Input::Initialization:: *)
Clear[IterateGetCaptureRate]
IterateGetCaptureRate[meDovermpD_,\[Kappa]s_,v0_,\[Sigma]DictsElectronic_,\[Sigma]DictsNuclear_]:= Module[{PDictList={}},

Monitor[
Do[(*Print["\[Kappa] is: ",N@\[Kappa]s[[\[Kappa]]]];*)
PDictList=Join[PDictList,GetCaptureRate[meDovermpD,\[Kappa]s[[\[Kappa]]],v0,\[Sigma]DictsElectronic,\[Sigma]DictsNuclear]]; (*flag for p->e*)
,{\[Kappa],Length[\[Kappa]s]}];
,\[Kappa]];

PDictList
]


(* ::Subsubsection::Closed:: *)
(*Interpolate Capture Probability distributions*)


(* ::Input::Initialization:: *)
Clear[InterpolatePc]
InterpolatePc[PDictList_]:=Module[{Intertable,Interf,\[Kappa]GatheredDicts,mpD,meD,mpDcap,m\[Chi],v\[Chi]Max,\[Beta]D,v0},

{mpD,meD,m\[Chi],mpDcap,v\[Chi]Max,\[Beta]D,v0}=Union[{"mpD","meD","m\[Chi]","mpDcap","v\[Chi]Max","\[Beta]D","v0"}/.PDictList][[1]];
(*Print[Union[{"mpD","meD","v\[Chi]Max","\[Beta]D","v0"}/.PDictList][[1]]];*)
\[Kappa]GatheredDicts=Gather[PDictList,(#1["\[Kappa]"]==#2["\[Kappa]"])&];
Interf=Table[<|"Pf"->Interpolation[Log10@{{"v\[Chi]E","bEarth"},"Pcap"}/.\[Kappa]GatheredDicts[[i]],InterpolationOrder->1],"\[Kappa]"->(Union["\[Kappa]"/.\[Kappa]GatheredDicts[[i]]][[1]]),"mpD"->mpD,"meD"->meD,"m\[Chi]"->m\[Chi],"mpDcap"->mpDcap,"v\[Chi]Max"->v\[Chi]Max,"\[Beta]D"->\[Beta]D,"v0"->v0,"PfE"->Interpolation[Log10@{{"v\[Chi]E","bEarth"},"PcapE"}/.\[Kappa]GatheredDicts[[i]],InterpolationOrder->1],"PfNuc"->Interpolation[Log10@{{"v\[Chi]E","bEarth"},"PcapNuc"}/.\[Kappa]GatheredDicts[[i]],InterpolationOrder->1],"y"->Interpolation[Log10@{{"v\[Chi]E","bEarth"},"y"}/.\[Kappa]GatheredDicts[[i]],InterpolationOrder->1],"CSDcap"->Interpolation[{Log10@{"v\[Chi]E","bEarth"},1-2"CSDcap"}/.\[Kappa]GatheredDicts[[i]],InterpolationOrder->1],"int\[Lambda]invE"->Interpolation[Log10@{{"v\[Chi]E","bEarth"},Max["int\[Lambda]invE",10^-100]}/.\[Kappa]GatheredDicts[[i]],InterpolationOrder->1],"int\[Lambda]invNuc"->Interpolation[Log10@{{"v\[Chi]E","bEarth"},Max["int\[Lambda]invNuc",10^-100]}/.\[Kappa]GatheredDicts[[i]],InterpolationOrder->1]|>,{i,Length[\[Kappa]GatheredDicts]}];
Interf
]


(* ::Input::Initialization:: *)
(*Clear[InterpolatePc]
InterpolatePc[PDictList_]:=Module[{Intertable,Interf,\[Kappa]GatheredDicts,mpD,meD,m\[Chi],v\[Chi]Max,\[Beta]D,v0},

{mpD,meD,m\[Chi],v\[Chi]Max,\[Beta]D,v0}=Union[{"mpD","meD","m\[Chi]","v\[Chi]Max","\[Beta]D","v0"}/.PDictList][[1]];
(*Print[Union[{"mpD","meD","v\[Chi]Max","\[Beta]D","v0"}/.PDictList][[1]]];*)
\[Kappa]GatheredDicts=Gather[PDictList,(#1["\[Kappa]"]==#2["\[Kappa]"])&];
Interf=Table[<|"Pf"->Interpolation[Log10@{{"v\[Chi]E","bEarth"},"Pcap"}/.\[Kappa]GatheredDicts[[i]],InterpolationOrder->1],"\[Kappa]"->(Union["\[Kappa]"/.\[Kappa]GatheredDicts[[i]]][[1]]),"mpD"->mpD,"meD"->meD,"m\[Chi]"->m\[Chi],"v\[Chi]Max"->v\[Chi]Max,"\[Beta]D"->\[Beta]D,"v0"->v0,"PfE"->Interpolation[Log10@{{"v\[Chi]E","bEarth"},"PcapE"}/.\[Kappa]GatheredDicts[[i]],InterpolationOrder->1],"PfNuc"->Interpolation[Log10@{{"v\[Chi]E","bEarth"},"PcapNuc"}/.\[Kappa]GatheredDicts[[i]],InterpolationOrder->1],"y"->Interpolation[Log10@{{"v\[Chi]E","bEarth"},"y"}/.\[Kappa]GatheredDicts[[i]],InterpolationOrder->1],"CSDcap"->Interpolation[{Log10@{"v\[Chi]E","bEarth"},1-2"CSDcap"}/.\[Kappa]GatheredDicts[[i]],InterpolationOrder->1],"int\[Lambda]invE"->Interpolation[Log10@{{"v\[Chi]E","bEarth"},Max["int\[Lambda]invE",10^-100]}/.\[Kappa]GatheredDicts[[i]],InterpolationOrder->1],"int\[Lambda]invNuc"->Interpolation[Log10@{{"v\[Chi]E","bEarth"},Max["int\[Lambda]invNuc",10^-100]}/.\[Kappa]GatheredDicts[[i]],InterpolationOrder->1]|>,{i,Length[\[Kappa]GatheredDicts]}];
Interf
]*)


(* ::Subsection:: *)
(*Get Total Capture Rate*)


(* ::Subsubsection::Closed:: *)
(*aDM number density*)


(* ::Input::Initialization:: *)
nD0\[CapitalLambda]CDM[fD_,mDMtot_] := Module[{\[CapitalLambda]CDMrepl,h,\[Rho]critbaumanncm,\[Rho]critbaumannm, n0,mSMtot,\[Rho]DM0},
(*mDMtot - [kg] sum of masses of DM species *)

\[CapitalLambda]CDMrepl= {"\[CapitalOmega]k" ->0.00,"\[CapitalOmega]r"->9.4 10^-5, "\[CapitalOmega]m"->0.32,"\[CapitalOmega]\[CapitalLambda]"->0.68,"\[CapitalOmega]b"->0.05,"\[CapitalOmega]c"->0.27}; (*[] \[CapitalLambda]CDM density parameters*)
h= 0.67; (*baumann 1.2.64*)

\[Rho]critbaumanncm = 1.1 10^-5 h^2; (*[protons cm^-3]*) 
\[Rho]critbaumannm = \[Rho]critbaumanncm  (10^6); (*[protons m^-3] ((10^2cm)/m)^3*) 

n0 = "\[CapitalOmega]b" \[Rho]critbaumannm /.\[CapitalLambda]CDMrepl;  (*[m^-3] proton / electron number density today*)
mSMtot= "m" + 2000 "m"/.Constants`SIConstRepl;

(*Print[n0];*)

\[Rho]DM0= "\[CapitalOmega]c"/"\[CapitalOmega]b" n0 mSMtot/.\[CapitalLambda]CDMrepl;
\[Rho]DM0/mDMtot fD
]


(* ::Input::Initialization:: *)
naDM[fD_,mDMtot_] := Module[{\[Rho]DM},
(*fD - [ ] ionized mass fraction of local DM density in the galaxy
mDMtot - [kg] sum of masses of DM species *)

(*Print[n0];*)

\[Rho]DM = 0.3 GeV cm^-3/.cm->10^-2/.GeV->10^9 ("JpereV")/("c")^2/.Constants`SIConstRepl;  (*[kg m^-3] Local DM number density today*)

\[Rho]DM/mDMtot fD
]


(* ::Subsubsection::Closed:: *)
(*Flux and NcMax defn*)


(* ::Input::Initialization:: *)
Clear[GetNcMax]
GetNcMax[InterPDicts_,keys_:{"Pf"},mpDcap_:True]:=Module[{tempInterPDicts,Pfs,rE,vesc,tEarth,AEarth,m\[Chi],\[Beta]D,Pfdom,fintegrand,intregboundaries,dNcMaxdt},
(*Get maximum possible accumulation over the lifetime of Earth, with \[Kappa]=1

InterPDicts - Probability Associations (output of IterateGetCaptureRate)
keys - a list of Association Keys for Probability Interpolations (function of v\[Chi] and bE). If more than one key is given, they will be summed when computing the capture rate. Default value corresponds to the Total capture probability.*)

tempInterPDicts = InterPDicts;

If[AnyTrue[Table[MemberQ[Flatten@Keys@InterPDicts,keys[[k]]],{k,Length[keys]}],!#&],Print["One of these keys is not a valid association for the interpolation dictionary"];Return[0,Module]];
Pfs=keys/.InterPDicts;

rE="rE"/.Constants`EarthRepl;
vesc = "vesc"/.Constants`EarthRepl;
tEarth = 4.543 10^9 (356  24 3600);(*[s] age of the Earth*)
AEarth = 4 \[Pi] rE^2; (*[(m^2)] Surface Area of the Earth*)

Do[
(*m\[Chi]="mpD"/.InterPDicts[[i]]; (*flag for p->e*)*)
m\[Chi]="m\[Chi]"/.InterPDicts[[i]]; (*flag for p->e*)
\[Beta]D="\[Beta]D"/.InterPDicts[[i]];
Pfdom=Union@Table[Pfs[[i,k]]["Domain"],{k,Length[Pfs[[i]]]}];

If[(Dimensions@Pfdom)[[1]]!=1,Print["The interpolation functions corresponding to these keys have inequivalent domains, so will not be integrated over."];Return[0,Module]];

Pfdom = Pfdom[[1]];
(*Print[Pfdom];
Print[Dimensions@Pfs];*)
(*fintegrand[bE_?NumericQ,v\[Chi]_?NumericQ]:=bE^2\[ExponentialE]^(-(1/2) m\[Chi] (v\[Chi]^2-vesc^2) \[Beta]D)  v\[Chi]^3/Sqrt[1 + bE^2/rE^2]10^(Pfs[[i]][Log10@v\[Chi],Log10@bE]);*)
fintegrand[bE_?NumericQ,v\[Chi]_?NumericQ]:=bE^2 E^(-(1/2) m\[Chi] (v\[Chi]^2-vesc^2) \[Beta]D)  v\[Chi]^3/Sqrt[1 + bE^2/rE^2] Sum[10^(Pfs[[i,k]][Log10@v\[Chi],Log10@bE]),{k,Length[Pfs[[i]]]}];(*[(s^-1)] goes as <Subscript[v, \[Chi]]> Subscript[n, aDM] Subscript[A, Earth] so this rate over Subscript[V, Earth] is <Subscript[v, \[Chi]]> Subscript[n, aDM] Subscript[R, Earth]^-1 which is what was used in the appendix*)


intregboundaries=If[mpDcap,Table[10^c,{c,-1,4}],Table[10^c,{c,-3,4}]]; (*THIS MAY GO OUT OF RANGE, CHANGE ACCORDING TO INTERPOLATION PRESCRIPTION*)
(*intregboundaries=If[mpDcap,Table[10^c,{c,1,4}],Table[10^c,{c,-3,4}]];*)
(*dNcMaxdt=AEarth(3 "nD" Sqrt[2/\[Pi]]  (m\[Chi] \[Beta]D)^(3/2))/rE^3;Quiet@(NIntegrate[fintegrand[bE,v\[Chi]],{v\[Chi],10^Pfdom[[1,1]],10^Pfdom[[1,1]]+intregboundaries[[1]]},
{bE,10^Pfdom[[2,1]],10^Pfdom[[2,2]]},PrecisionGoal->3,MaxRecursion->8,AccuracyGoal->3]+Sum[With[{ct=c},NIntegrate[fintegrand[bE,v\[Chi]],{v\[Chi],10^Pfdom[[1,1]]+intregboundaries[[ct-1]],10^Pfdom[[1,1]]+intregboundaries[[ct]]},
{bE,10^Pfdom[[2,1]],10^Pfdom[[2,2]]},PrecisionGoal->3,MaxRecursion->8,AccuracyGoal->3]],{c,2,Length[intregboundaries]}]+NIntegrate[fintegrand[bE,v\[Chi]],{v\[Chi],10^Pfdom[[1,1]]+intregboundaries[[-1]],10^Pfdom[[1,2]]},
{bE,10^Pfdom[[2,1]],10^Pfdom[[2,2]]},PrecisionGoal->3,MaxRecursion->8,AccuracyGoal->3]);   (*we break up the integration region because Nuclear capture probability is EXTREMELY sharply peaked near v\[Chi]=vesc. Mathematica will ignore this contribution for most of the parameter space if we don't divide up the integration region like this. *)*)
dNcMaxdt=AEarth (3 "nD" Sqrt[2/\[Pi]]  (m\[Chi] \[Beta]D)^(3/2))/rE^3 Quiet@(NIntegrate[fintegrand[bE,v\[Chi]],{v\[Chi],10^Pfdom[[1,1]],10^Pfdom[[1,1]]+intregboundaries[[1]]},
{bE,10^Pfdom[[2,1]],10^Pfdom[[2,2]]},PrecisionGoal->3,MaxRecursion->8,AccuracyGoal->3]+Sum[NIntegrate[fintegrand[bE,v\[Chi]],{v\[Chi],10^Pfdom[[1,1]]+intregboundaries[[c-1]],10^Pfdom[[1,1]]+intregboundaries[[c]]},
{bE,10^Pfdom[[2,1]],10^Pfdom[[2,2]]},PrecisionGoal->3,MaxRecursion->8,AccuracyGoal->3],{c,2,Length[intregboundaries]}]+NIntegrate[fintegrand[bE,v\[Chi]],{v\[Chi],10^Pfdom[[1,1]]+intregboundaries[[-1]],10^Pfdom[[1,2]]},
{bE,10^Pfdom[[2,1]],10^Pfdom[[2,2]]},PrecisionGoal->3,MaxRecursion->8,AccuracyGoal->3]);   (*we break up the integration region because Nuclear capture probability is EXTREMELY sharply peaked near v\[Chi]=vesc. Mathematica will ignore this contribution for most of the parameter space if we don't divide up the integration region like this. *)

If[keys=={"Pf"},
AppendTo[tempInterPDicts[[i]],<|"dNcMaxdt"->dNcMaxdt|>];
AppendTo[tempInterPDicts[[i]],<|"NcMax"->dNcMaxdt tEarth|>];
,
AppendTo[tempInterPDicts[[i]],<|"dNcMaxdt"<>StringJoin[Table["_"<>keys[[k]],{k,Length[keys]}]]->dNcMaxdt|>];
AppendTo[tempInterPDicts[[i]],<|"NcMax"<>StringJoin[Table["_"<>keys[[k]],{k,Length[keys]}]]->dNcMaxdt tEarth|>];
];

,{i,Length[Pfs]}];

tempInterPDicts
]


(* ::Input:: *)
(**)


(* ::Subsubsection::Closed:: *)
(*Evaporation Rate*)


(* ::Input::Initialization:: *)
(*Clear[GetEvaporationRate]
GetEvaporationRate[\[Sigma]Dict_,mpDcap_:True]:= Module[{v\[Chi]min,v\[Chi]max,m\[Chi],\[Omega]min,params,regionindex,NucleusParams,mN,nT,\[Beta]E,\[Omega]capture,\[Delta]\[Xi],\[Xi]min,\[Xi]of\[Omega]andv\[Chi],\[Omega]maxofv\[Chi],d\[Sigma]dERofv\[Chi]and\[Omega],\[Sigma]ofv\[Chi]and\[Omega]min,\[CapitalOmega]evapperparticle,\[CapitalGamma]evappercapturedparticle,v\[Chi]st\[Omega]capis\[Omega]max,RegionVolume,\[CapitalOmega]evapperparticlewvT,EarthVolume,EvapRescale,fintegrand,intregboundaries,Nv\[Chi]=100},
\[Omega]min="\[Omega]min"/.\[Sigma]Dict;
params=\[Sigma]Dict["params"];
NucleusParams=\[Sigma]Dict["NucleusParams"];
\[Delta]\[Xi]=(10^"\[Xi]"/.\[Sigma]Dict)[[1]];
m\[Chi]=("m\[Chi]"/.\[Sigma]Dict);
mN = "mN"/.\[Sigma]Dict;
regionindex="regionindex"/.\[Sigma]Dict;
\[Beta]E = {"\[Beta]crust","\[Beta]core"}[[regionindex]]/.Constants`EarthRepl;(*Earth Temeperature*)
nT = "nI"/.NucleusParams;

\[Xi]of\[Omega]andv\[Chi]=EnergyLoss`\[Xi]of\[Omega]Nuc[#1,\[Omega]min,#2,m\[Chi],mN,params]&;

\[Sigma]ofv\[Chi]and\[Omega]min=10^(Re[("\[Sigma]of\[Xi]f"/.\[Sigma]Dict)[Log10@#1,Log10@\[Xi]of\[Omega]andv\[Chi][#2,#1]]])&;

v\[Chi]min = 10^(("\[Sigma]of\[Xi]f"/.\[Sigma]Dict)["Domain"][[1,1]]);
\[Xi]min = 10^(("\[Sigma]of\[Xi]f"/.\[Sigma]Dict)["Domain"][[2,1]]);

(*Post \[Omega]capture*)
\[Omega]maxofv\[Chi]=EnergyLoss`\[Omega]maxNuc[#1,m\[Chi],mN,params]&; (*maximum energy transferable in a collision*)

\[Omega]capture = 1/(2 ("\[HBar]"/.params))m\[Chi] (("vesc"/.Constants`EarthRepl)^2-#^2 )&;(*energy transfered to aDM such that it has v > Subscript[v, esc]*)
(*We are going outside the range of the interpolating function at times. *)

v\[Chi]st\[Omega]capis\[Omega]max=v\[Chi]/.NSolve[\[Omega]maxofv\[Chi][v\[Chi]]==\[Omega]capture[v\[Chi]],v\[Chi]][[2]];(*find the lowest velocity such that aDM can be evaporated through a single scatter*)

EvapRescale=1/(-If[x^2<500,\[ExponentialE]^(-(1/2)x^2) Sqrt[2/\[Pi]] x,0]+Erf[x/Sqrt[2]])/.x->Sqrt[m\[Chi]]"vesc" Sqrt[\[Beta]E]/.Constants`EarthRepl; (*account for the fact that the MB distribution is truncated, renormalize*)

(*fintegrand[vT_?NumericQ,v\[Chi]_?NumericQ]:=vT^2 E^(-\[Beta]E 1/2mN vT^2)((vT+v\[Chi])^3-((vT-v\[Chi])^2)^(3/2))/(3 vT v\[Chi])v\[Chi]^2E^(-\[Beta]E 1/2m\[Chi] v\[Chi]^2) If[\[Omega]maxofv\[Chi][v\[Chi]]>\[Omega]capture[v\[Chi]],1,0]\[Sigma]ofv\[Chi]and\[Omega]min[v\[Chi],\[Omega]capture[v\[Chi]]];*)
fintegrand[vT_?NumericQ,v\[Chi]_?NumericQ]:=vT^2 E^(-\[Beta]E 1/2mN vT^2)((vT+v\[Chi])^3-((vT-v\[Chi])^2)^(3/2))/(3 vT v\[Chi])v\[Chi]^2E^(-\[Beta]E 1/2m\[Chi] v\[Chi]^2) If[\[Omega]maxofv\[Chi][v\[Chi]]>\[Omega]capture[v\[Chi]]&&(\[Xi]of\[Omega]andv\[Chi][\[Omega]capture[v\[Chi]],v\[Chi]])>\[Xi]min,1,0]\[Sigma]ofv\[Chi]and\[Omega]min[v\[Chi],\[Omega]capture[v\[Chi]]];


(*Including the measures, This scales as \[Sigma] v\[Chi]^3 vT^3 Max (v\[Chi],vT)^ *)

(*\[CapitalOmega]evapperparticlewvT=EvapRescale  nT  ((Sqrt[m\[Chi] mN ]\[Beta]E)/(2 \[Pi]))^3 4 \[Pi] 2 \[Pi] NIntegrate[fintegrand[vT,v\[Chi]],{v\[Chi],v\[Chi]st\[Omega]capis\[Omega]max,"vesc"/.Constants`EarthRepl},{vT,0,\[Infinity]},PrecisionGoal->3,MaxRecursion->5,AccuracyGoal->3];*)

(*[(s^-1)] ignoring measures, This scales as nT v\[Chi]^-3 vT^-3 so total scaling is nT \[Sigma] <Max (v\[Chi],vT)>, so for light aDM compared to nuclei(v\[Chi]>>vT) this should scale identically to the appendix*)


(*first factor in table comes from integral of |vT - v\[Chi]| over cos\[Theta] *)


intregboundaries=If[mpDcap,Table[10^c,{c,-1,4}],Table[10^c,{c,-3,4}]];
(*intregboundaries=If[mpDcap,Table[(("vesc"/.Constants`EarthRepl)-v\[Chi]st\[Omega]capis\[Omega]max)10^c,{c,-5,0}],Table[(("vesc"/.Constants`EarthRepl)-v\[Chi]st\[Omega]capis\[Omega]max)10^c,{c,-7,0}]];*)

\[CapitalOmega]evapperparticlewvT=EvapRescale  nT  ((Sqrt[m\[Chi] mN ]\[Beta]E)/(2 \[Pi]))^3 4 \[Pi] 2 \[Pi]Quiet@(NIntegrate[fintegrand[vT,v\[Chi]],{v\[Chi],v\[Chi]st\[Omega]capis\[Omega]max,v\[Chi]st\[Omega]capis\[Omega]max + intregboundaries[[1]]},{vT,0,\[Infinity]},PrecisionGoal->3,MaxRecursion->5,AccuracyGoal->3]+Sum[NIntegrate[fintegrand[vT,v\[Chi]],{v\[Chi],v\[Chi]st\[Omega]capis\[Omega]max+ intregboundaries[[c-1]],v\[Chi]st\[Omega]capis\[Omega]max + intregboundaries[[c]]},{vT,0,\[Infinity]},PrecisionGoal->3,MaxRecursion->5,AccuracyGoal->3],{c,2,Length[intregboundaries]}]+NIntegrate[fintegrand[vT,v\[Chi]],{v\[Chi],v\[Chi]st\[Omega]capis\[Omega]max+ intregboundaries[[-1]],"vesc"/.Constants`EarthRepl},{vT,0,\[Infinity]},PrecisionGoal->3,MaxRecursion->5,AccuracyGoal->3]); (*[(s^-1)] ignoring measures, This scales as nT v\[Chi]^-3 vT^-3 so total scaling is nT \[Sigma] <Max (v\[Chi],vT)>, so for light aDM compared to nuclei(v\[Chi]>>vT) this should scale identically to the appendix*)(*first factor in table comes from integral of |vT - v\[Chi]| over cos\[Theta] *)

RegionVolume =  4\[Pi] Integrate[r^2Region\[CapitalTheta][r,regionindex],{r,0,"rE"/.Constants`EarthRepl}]; (*since this is done independent of \[Kappa], the right place to integrate over r with MFP is further down pipeline*)

EarthVolume =  (4\[Pi])/3 ("rE"/.Constants`EarthRepl)^3;
\[CapitalGamma]evappercapturedparticle=\[CapitalOmega]evapperparticlewvT RegionVolume/EarthVolume; (*\[CapitalGamma]\[Kappa]^2 corresponds to Revp in appendix*)(*Assumes homogeneously distributed aDM in the Earth*)

\[CapitalGamma]evappercapturedparticle
]*)


(* ::Input::Initialization:: *)
Clear[GetEvaporationRate]
GetEvaporationRate[\[Sigma]Dict_,mpDcap_:True]:= Module[{v\[Chi]min,v\[Chi]max,m\[Chi],\[Omega]min,params,regionindex,NucleusParams,mN,nT,\[Beta]E,\[Omega]capture,\[Delta]\[Xi],\[Xi]min,\[Xi]of\[Omega]andv\[Chi],\[Omega]maxofv\[Chi],d\[Sigma]dERofv\[Chi]and\[Omega],\[Sigma]ofv\[Chi]and\[Omega]min,\[CapitalOmega]evapperparticle,\[CapitalGamma]evappercapturedparticle,v\[Chi]st\[Omega]capis\[Omega]max,RegionVolume,\[CapitalOmega]evapperparticlewvT,EarthVolume,EvapRescale,fintegrand,intregboundaries,Nv\[Chi]=100},
\[Omega]min="\[Omega]min"/.\[Sigma]Dict;
params=\[Sigma]Dict["params"];
NucleusParams=\[Sigma]Dict["NucleusParams"];
\[Delta]\[Xi]=(10^"\[Xi]"/.\[Sigma]Dict)[[1]];
m\[Chi]=("m\[Chi]"/.\[Sigma]Dict);
mN = "mN"/.\[Sigma]Dict;
regionindex="regionindex"/.\[Sigma]Dict;
\[Beta]E = {"\[Beta]crust","\[Beta]core"}[[regionindex]]/.Constants`EarthRepl;(*Earth Temeperature*)
nT = "nI"/.NucleusParams;

\[Xi]of\[Omega]andv\[Chi]=EnergyLoss`\[Xi]of\[Omega]Nuc[#1,\[Omega]min,#2,m\[Chi],mN,params]&;

\[Sigma]ofv\[Chi]and\[Omega]min=10^(Re[("\[Sigma]of\[Xi]f"/.\[Sigma]Dict)[Log10@#1,Log10@\[Xi]of\[Omega]andv\[Chi][#2,#1]]])&;

v\[Chi]min = 10^(("\[Sigma]of\[Xi]f"/.\[Sigma]Dict)["Domain"][[1,1]]);
\[Xi]min = 10^(("\[Sigma]of\[Xi]f"/.\[Sigma]Dict)["Domain"][[2,1]]);

(*Post \[Omega]capture*)
\[Omega]maxofv\[Chi]=EnergyLoss`\[Omega]maxNuc[#1,m\[Chi],mN,params]&; (*maximum energy transferable in a collision*)

\[Omega]capture = 1/(2 ("\[HBar]"/.params)) m\[Chi] (("vesc"/.Constants`EarthRepl)^2-#^2 )&;(*energy transfered to aDM such that it has v > Subscript[v, esc]*)
(*We are going outside the range of the interpolating function at times. *)

v\[Chi]st\[Omega]capis\[Omega]max=v\[Chi]/.NSolve[\[Omega]maxofv\[Chi][v\[Chi]]==\[Omega]capture[v\[Chi]],v\[Chi]][[2]];(*find the lowest velocity such that aDM can be evaporated through a single scatter*)

EvapRescale=1/(-If[x^2<500,E^(-(1/2) x^2) Sqrt[2/\[Pi]] x,0]+Erf[x/Sqrt[2]])/.x->Sqrt[m\[Chi]]"vesc" Sqrt[\[Beta]E]/.Constants`EarthRepl; (*account for the fact that the MB distribution is truncated, renormalize*)

(*fintegrand[vT_?NumericQ,v\[Chi]_?NumericQ]:=vT^2 E^(-\[Beta]E 1/2mN vT^2)((vT+v\[Chi])^3-((vT-v\[Chi])^2)^(3/2))/(3 vT v\[Chi])v\[Chi]^2E^(-\[Beta]E 1/2m\[Chi] v\[Chi]^2) If[\[Omega]maxofv\[Chi][v\[Chi]]>\[Omega]capture[v\[Chi]],1,0]\[Sigma]ofv\[Chi]and\[Omega]min[v\[Chi],\[Omega]capture[v\[Chi]]];*)
fintegrand[vT_?NumericQ,v\[Chi]_?NumericQ]:=vT^2 E^(-\[Beta]E 1/2 mN vT^2) ((vT+v\[Chi])^3-((vT-v\[Chi])^2)^(3/2))/(3 vT v\[Chi]) v\[Chi]^2 E^(-\[Beta]E 1/2 m\[Chi] v\[Chi]^2) If[\[Omega]maxofv\[Chi][v\[Chi]]>\[Omega]capture[v\[Chi]]&&(\[Xi]of\[Omega]andv\[Chi][\[Omega]capture[v\[Chi]],v\[Chi]])>\[Xi]min,1,0]\[Sigma]ofv\[Chi]and\[Omega]min[v\[Chi],\[Omega]capture[v\[Chi]]];


(*Including the measures, This scales as \[Sigma] v\[Chi]^3 vT^3 Max (v\[Chi],vT)^ *)

(*\[CapitalOmega]evapperparticlewvT=EvapRescale  nT  ((Sqrt[m\[Chi] mN ]\[Beta]E)/(2 \[Pi]))^3 4 \[Pi] 2 \[Pi] NIntegrate[fintegrand[vT,v\[Chi]],{v\[Chi],v\[Chi]st\[Omega]capis\[Omega]max,"vesc"/.Constants`EarthRepl},{vT,0,\[Infinity]},PrecisionGoal->3,MaxRecursion->5,AccuracyGoal->3];*)

(*[(s^-1)] ignoring measures, This scales as nT v\[Chi]^-3 vT^-3 so total scaling is nT \[Sigma] <Max (v\[Chi],vT)>, so for light aDM compared to nuclei(v\[Chi]>>vT) this should scale identically to the appendix*)


(*first factor in table comes from integral of |vT - v\[Chi]| over cos\[Theta] *)


intregboundaries=If[mpDcap,Table[10^c,{c,-1,4}],Table[10^c,{c,-3,4}]];
(*intregboundaries=If[mpDcap,Table[(("vesc"/.Constants`EarthRepl)-v\[Chi]st\[Omega]capis\[Omega]max)10^c,{c,-5,0}],Table[(("vesc"/.Constants`EarthRepl)-v\[Chi]st\[Omega]capis\[Omega]max)10^c,{c,-7,0}]];*)

\[CapitalOmega]evapperparticlewvT=EvapRescale  nT  ((Sqrt[m\[Chi] mN ]\[Beta]E)/(2 \[Pi]))^3 4 \[Pi] 2 \[Pi] Quiet@(NIntegrate[fintegrand[vT,v\[Chi]],{v\[Chi],v\[Chi]st\[Omega]capis\[Omega]max,v\[Chi]st\[Omega]capis\[Omega]max + intregboundaries[[1]]},{vT,0,\[Infinity]},PrecisionGoal->3,MaxRecursion->5,AccuracyGoal->3]+Sum[NIntegrate[fintegrand[vT,v\[Chi]],{v\[Chi],v\[Chi]st\[Omega]capis\[Omega]max+ intregboundaries[[c-1]],v\[Chi]st\[Omega]capis\[Omega]max + intregboundaries[[c]]},{vT,0,\[Infinity]},PrecisionGoal->3,MaxRecursion->5,AccuracyGoal->3],{c,2,Length[intregboundaries]}]+NIntegrate[fintegrand[vT,v\[Chi]],{v\[Chi],v\[Chi]st\[Omega]capis\[Omega]max+ intregboundaries[[-1]],"vesc"/.Constants`EarthRepl},{vT,0,\[Infinity]},PrecisionGoal->3,MaxRecursion->5,AccuracyGoal->3]); (*[(s^-1)] ignoring measures, This scales as nT v\[Chi]^-3 vT^-3 so total scaling is nT \[Sigma] <Max (v\[Chi],vT)>, so for light aDM compared to nuclei(v\[Chi]>>vT) this should scale identically to the appendix*)(*first factor in table comes from integral of |vT - v\[Chi]| over cos\[Theta] *)

RegionVolume =  4\[Pi] Integrate[r^2 Region\[CapitalTheta][r,regionindex],{r,0,"rE"/.Constants`EarthRepl}]; (*since this is done independent of \[Kappa], the right place to integrate over r with MFP is further down pipeline*)

EarthVolume =  (4\[Pi])/3 ("rE"/.Constants`EarthRepl)^3;
\[CapitalGamma]evappercapturedparticle=\[CapitalOmega]evapperparticlewvT RegionVolume/EarthVolume; (*\[CapitalGamma]\[Kappa]^2 corresponds to Revp in appendix*)(*Assumes homogeneously distributed aDM in the Earth*)

<|"\[CapitalGamma]"->\[CapitalGamma]evappercapturedparticle,"regionindex"->regionindex|>
]


(* ::Input::Initialization:: *)
Clear[IterateGetEvaporationRate]
IterateGetEvaporationRate[Nuclear\[Sigma]Dicts_]:=Module[{\[CapitalGamma]list={},\[CapitalGamma]Total},
Monitor[
Do[
AppendTo[\[CapitalGamma]list,GetEvaporationRate[Nuclear\[Sigma]Dicts[[n]]]],
{n,Length[Nuclear\[Sigma]Dicts]}
];
,n];
\[CapitalGamma]Total=Total[\[CapitalGamma]list[[;;,1]]];
<|"\[CapitalGamma]Total"->\[CapitalGamma]Total,"\[CapitalGamma]s"->\[CapitalGamma]list|>
]


(* ::Subsubsection::Closed:: *)
(*Evaporation Rate Electronic*)


(* ::Input::Initialization:: *)
Clear[GetEvaporationRateend]
GetEvaporationRateend[\[Sigma]Dict_,mpDcap_:True]:= Module[{v\[Chi]min,v\[Chi]max,m\[Chi],\[Omega]min,params,regionindex,mT,nT,vF,vesc,vTm,vTp,\[Beta]E,\[Omega]capture,\[Delta]\[Xi],\[Xi]min,\[Xi]of\[Omega]andv\[Chi]t,\[Omega]maxofv\[Chi]t,d\[Sigma]dERofv\[Chi]and\[Omega],\[Sigma]ofv\[Chi]and\[Omega]min,\[CapitalOmega]evapperparticle,\[CapitalGamma]evappercapturedparticle,v\[Chi]st\[Omega]capis\[Omega]max,RegionVolume,\[CapitalOmega]evapperparticlewvT,EarthVolume,EvapRescale,fintegrand,intregboundaries,PLOT=False},
(*\[Omega]min="\[Omega]min"/.\[Sigma]Dict;*)
params=\[Sigma]Dict["params"];
\[Delta]\[Xi]=(10^"\[Xi]"/.\[Sigma]Dict)[[1]];
m\[Chi]=("m\[Chi]"/.\[Sigma]Dict);
mT = "m"/.\[Sigma]Dict;
regionindex="regionindex"/.\[Sigma]Dict;
\[Beta]E = {"\[Beta]crust","\[Beta]core"}[[regionindex]]/.Constants`EarthRepl;(*Earth Temeperature*)
nT = "ne"/.params;
vF="vF"/.params;
vesc = ("vesc"/.Constants`EarthRepl);
(*\[Xi]of\[Omega]andv\[Chi]t=EnergyLoss`\[Xi]of\[Omega]andv\[Chi][#1,\[Omega]min,#2,m\[Chi],params]&;*)

(*v\[Chi]min = 10^(("\[Sigma]of\[Xi]f"/.\[Sigma]Dict)["Domain"][[1,1]]);*)
v\[Chi]min = Max[10^(("\[Sigma]of\[Xi]f"/.\[Sigma]Dict)["Domain"][[1,1]]),Re[Sqrt[vesc^2 - 8/(m\[Chi] \[Beta]E)]]];
\[Xi]min = 10^(("\[Sigma]of\[Xi]f"/.\[Sigma]Dict)["Domain"][[2,1]]);


\[Omega]min="\[Omega]min"/.\[Sigma]Dict;
\[Xi]of\[Omega]andv\[Chi]t=EnergyLoss`\[Xi]of\[Omega]andv\[Chi]nd[#1,\[Omega]min,#2,m\[Chi],params]&;

\[Sigma]ofv\[Chi]and\[Omega]min=10^(("\[Sigma]of\[Xi]f"/.\[Sigma]Dict)[Log10@#1,Log10@\[Xi]of\[Omega]andv\[Chi]t[#2,#1]])&;


(*Post \[Omega]capture*)
(*\[Omega]maxofv\[Chi]t=EnergyLoss`\[Omega]maxofv\[Chi][#1,m\[Chi],params]&; (*maximum energy transferable in a collision*)
*)
\[Omega]maxofv\[Chi]t=EnergyLoss`\[Omega]maxofv\[Chi]nd[#1,m\[Chi],params]&; (*maximum energy transferable in a collision*)

\[Omega]capture = 1/(2 ("\[HBar]"/.params)) m\[Chi] (vesc^2-#^2 )&;(*energy transfered to aDM such that it has v > Subscript[v, esc]*)
(*We are going outside the range of the interpolating function at times. *)

(*v\[Chi]st\[Omega]capis\[Omega]max=v\[Chi]/.NSolve[\[Omega]maxofv\[Chi]t[v\[Chi]]==\[Omega]capture[v\[Chi]],v\[Chi]][[2]];*)(*find the lowest velocity such that aDM can be evaporated through a single scatter*)


(*v\[Chi]st\[Omega]capis\[Omega]max=NSolve[\[Omega]maxofv\[Chi]t[v\[Chi]]==\[Omega]capture[v\[Chi]],v\[Chi]][[2]];
If[v\[Chi]st\[Omega]capis\[Omega]max==v\[Chi],v\[Chi]st\[Omega]capis\[Omega]max=v\[Chi]min,v\[Chi]st\[Omega]capis\[Omega]max=v\[Chi]/.v\[Chi]st\[Omega]capis\[Omega]max];*)
v\[Chi]st\[Omega]capis\[Omega]max=v\[Chi]min;

(*v\[Chi]st\[Omega]capis\[Omega]max=\[Sqrt](("vesc"/.Constants`EarthRepl)^2- (2 "\[HBar]")/m\[Chi] \[Omega]maxofv\[Chi]t[#]/.params)&;*)
(*v\[Chi]st\[Omega]capis\[Omega]max=v\[Chi]/.NSolve[v\[Chi]==\[Sqrt](("vesc"/.Constants`EarthRepl)^2- (2 ("\[HBar]"/.params))/m\[Chi] \[Omega]maxofv\[Chi]t[#])&[v\[Chi]],v\[Chi]];*)
(*Print[Show[{LogPlot[\[Omega]maxofv\[Chi]t[v\[Chi]],{v\[Chi],0,"vesc"/.Constants`EarthRepl},PlotStyle->Black,PlotRange->All],LogPlot[\[Omega]capture[v\[Chi]],{v\[Chi],0,"vesc"/.Constants`EarthRepl}]}]];*)
(*Print[v\[Chi]st\[Omega]capis\[Omega]max];
Print[\[Omega]maxofv\[Chi]t[10^4]];*)
(*Print[v\[Chi]st\[Omega]capis\[Omega]max];*)

EvapRescale=1/(-If[x^2<500,E^(-(1/2) x^2) Sqrt[2/\[Pi]] x,0]+Erf[x/Sqrt[2]])/.x->Sqrt[m\[Chi]]"vesc" Sqrt[\[Beta]E]/.Constants`EarthRepl; (*account for the fact that the aDM MB distribution is truncated, normalize to 1*)

(*Print[LogLogPlot[v\[Chi]^2E^(-\[Beta]E 1/2m\[Chi] v\[Chi]^2) ,{v\[Chi],0,"vesc"/.EarthRepl}]];*)

fintegrand[vT_?NumericQ,v\[Chi]_?NumericQ]:=vT^2 (fnd[vT/vF]/.params)(((vT+v\[Chi])^3-((vT-v\[Chi])^2)^(3/2))/(3 vT v\[Chi]))v\[Chi]^2 E^(-\[Beta]E 1/2 m\[Chi] v\[Chi]^2) (*If[\[Omega]maxofv\[Chi]t[v\[Chi]]>\[Omega]capture[v\[Chi]]&&(\[Xi]of\[Omega]andv\[Chi]t[\[Omega]capture[v\[Chi]],v\[Chi]])>\[Xi]min&&v\[Chi]>v\[Chi]min,1,0]*)
If[\[Omega]maxofv\[Chi]t[v\[Chi]]>\[Omega]capture[v\[Chi]]&&(\[Xi]of\[Omega]andv\[Chi]t[\[Omega]capture[v\[Chi]],v\[Chi]])>\[Xi]min,1,0]\[Sigma]ofv\[Chi]and\[Omega]min[v\[Chi],\[Omega]capture[v\[Chi]]];

(*intregboundaries=If[mpDcap,Table[10^c,{c,-1,4}],Table[10^c,{c,-3,4}]];*)
intregboundaries=If[mpDcap,Table[(("vesc"/.Constants`EarthRepl)-v\[Chi]st\[Omega]capis\[Omega]max)10^c,{c,-5,-1}],Table[(("vesc"/.Constants`EarthRepl)-v\[Chi]st\[Omega]capis\[Omega]max)10^c,{c,-7,-1}]];

{vTm,vTp}={vF Dielectrics`\[Zeta]m,vF Dielectrics`\[Zeta]p}/.params;



If[PLOT,
Print[DensityPlot[EvapRescale fintegrand[vT,v\[Chi]],{v\[Chi],v\[Chi]st\[Omega]capis\[Omega]max,"vesc"/.Constants`EarthRepl},{vT,vTm,vTp},PlotLegends->Automatic,PlotRange->All](*,PlotLabel->"fintegrand"*)];

Print[DensityPlot[Log10@\[Sigma]ofv\[Chi]and\[Omega]min[v\[Chi],EnergyLoss`\[Omega]of\[Xi]andv\[Chi]nd[\[Xi],\[Omega]min,v\[Chi],m\[Chi],params]],{v\[Chi],v\[Chi]st\[Omega]capis\[Omega]max,"vesc"/.Constants`EarthRepl},{\[Xi],10^-8,1-10^-8},PlotLegends->Automatic,PlotRange->All]];

Print[v\[Chi]st\[Omega]capis\[Omega]max];
Print[10^(("\[Sigma]of\[Xi]f"/.\[Sigma]Dict)["Domain"][[1,1]])];
Print[Plot[Log10[\[Sigma]ofv\[Chi]and\[Omega]min[#,\[Omega]capture[#]](*If[\[Omega]maxofv\[Chi]t[#]>\[Omega]capture[#]&&(\[Xi]of\[Omega]andv\[Chi]t[\[Omega]capture[#],#])>\[Xi]min,1,0]*)]&[10^v\[Chi]],{v\[Chi],Log10@v\[Chi]st\[Omega]capis\[Omega]max,Log10@"vesc"/.Constants`EarthRepl},PlotLegends->Automatic,PlotRange->{-30,0},PlotLabel->"\[Sigma](\!\(\*SubscriptBox[\(\[Omega]\), \(capture\)]\))"]];

Print[Plot[\[Xi]of\[Omega]andv\[Chi]t[\[Omega]capture[v\[Chi]],v\[Chi]],{v\[Chi],v\[Chi]st\[Omega]capis\[Omega]max,"vesc"/.Constants`EarthRepl},PlotLegends->Automatic,PlotRange->All,PlotLabel->"\[Xi](\!\(\*SubscriptBox[\(\[Omega]\), \(capture\)]\)(\!\(\*SubscriptBox[\(v\), \(\[Chi]\)]\)))"]];
(*Print[\[Xi]min];*)
Print[Plot[\[Omega]maxofv\[Chi]t[v\[Chi]],{v\[Chi],v\[Chi]st\[Omega]capis\[Omega]max,"vesc"/.Constants`EarthRepl},PlotLegends->Automatic,PlotRange->All,PlotLabel->"\!\(\*SubscriptBox[\(\[Omega]\), \(max\)]\)(\!\(\*SubscriptBox[\(v\), \(\[Chi]\)]\))"]];
Print[\[Omega]min];
Print[Plot[\[Omega]capture[v\[Chi]],{v\[Chi],v\[Chi]st\[Omega]capis\[Omega]max,"vesc"/.Constants`EarthRepl},PlotLegends->Automatic,PlotRange->All,PlotLabel->"\!\(\*SubscriptBox[\(\[Omega]\), \(capture\)]\)(\!\(\*SubscriptBox[\(v\), \(\[Chi]\)]\))"]];

Print["This is actually okay, \[Omega]min > \[Omega]capture, this is not a precision issue with \[Xi]."];
Print["The cross-section decreases with v\[Chi]"];
Print["Why are there bumps at high masses?"];
Print[{{v\[Chi]st\[Omega]capis\[Omega]max,v\[Chi]st\[Omega]capis\[Omega]max + intregboundaries[[1]]},Table[With[{ct=c},{v\[Chi]st\[Omega]capis\[Omega]max+ intregboundaries[[ct-1]],v\[Chi]st\[Omega]capis\[Omega]max + intregboundaries[[ct]]}],{c,2,Length[intregboundaries]}],{v\[Chi]st\[Omega]capis\[Omega]max+ intregboundaries[[-1]],"vesc"/.Constants`EarthRepl}}]
];

\[CapitalOmega]evapperparticlewvT=EvapRescale  nT (("D"/.params)/Sqrt[2 \[Pi]]) (Sqrt[m\[Chi] \[Beta]E]/vF)^3 (*Quiet@*)(NIntegrate[fintegrand[vT,v\[Chi]],{v\[Chi],v\[Chi]st\[Omega]capis\[Omega]max,v\[Chi]st\[Omega]capis\[Omega]max + intregboundaries[[1]]},{vT,vTm,vTp},PrecisionGoal->3,MaxRecursion->5,AccuracyGoal->3]+Sum[With[{ct=c},NIntegrate[fintegrand[vT,v\[Chi]],{v\[Chi],v\[Chi]st\[Omega]capis\[Omega]max+ intregboundaries[[ct-1]],v\[Chi]st\[Omega]capis\[Omega]max + intregboundaries[[ct]]},{vT,vTm,vTp},PrecisionGoal->3,MaxRecursion->5,AccuracyGoal->3]],{c,2,Length[intregboundaries]}]+NIntegrate[fintegrand[vT,v\[Chi]],{v\[Chi],v\[Chi]st\[Omega]capis\[Omega]max+ intregboundaries[[-1]],"vesc"/.Constants`EarthRepl},{vT,vTm,vTp},PrecisionGoal->3,MaxRecursion->5,AccuracyGoal->3]); (*[(s^-1)] ignoring measures, This scales as nT v\[Chi]^-3 vT^-3 so total scaling is nT \[Sigma] <Max (v\[Chi],vT)>, so for light aDM compared to nuclei(v\[Chi]>>vT) this should scale identically to the appendix*)(*first factor in table comes from integral of |vT - v\[Chi]| over cos\[Theta] *)

RegionVolume =  4\[Pi] Integrate[r^2 Region\[CapitalTheta][r,regionindex],{r,0,"rE"/.Constants`EarthRepl}];
EarthVolume =  (4\[Pi])/3 ("rE"/.Constants`EarthRepl)^3;
\[CapitalGamma]evappercapturedparticle=\[CapitalOmega]evapperparticlewvT RegionVolume/EarthVolume; (*\[CapitalGamma]\[Kappa]^2 corresponds to Revp in appendix*)(*Assumes homogeneously distributed aDM in the Earth*)

(*\[CapitalGamma]evappercapturedparticle*)
<|"\[CapitalGamma]"->\[CapitalGamma]evappercapturedparticle,"regionindex"->regionindex|>
]


(* ::Input::Initialization:: *)
Clear[IterateGetEvaporationRateend]
IterateGetEvaporationRateend[Electronic\[Sigma]Dicts_]:=Module[{\[CapitalGamma]list={},\[CapitalGamma]Total},
(*Monitor[*)
Do[
If[Electronic\[Sigma]Dicts[[n]]!=<|"NA"->Null|>,
AppendTo[\[CapitalGamma]list,GetEvaporationRateend[Electronic\[Sigma]Dicts[[n]]]]],
{n,Length[Electronic\[Sigma]Dicts]}
];
(*,n];*)
(*\[CapitalGamma]Total=Total[\[CapitalGamma]list];*)
\[CapitalGamma]Total=Total[\[CapitalGamma]list[[;;,1]]];
<|"\[CapitalGamma]Total"->\[CapitalGamma]Total,"\[CapitalGamma]s"->\[CapitalGamma]list|>
]


(* ::Subsubsection:: *)
(*Get Penetration Depth*)


(* ::Input::Initialization:: *)
Clear[GetPenetrationDepthfnofv]
GetPenetrationDepthfnofv[v\[Chi]_,Electronic\[Sigma]Dicts_,Nuclear\[Sigma]Dicts_,vesc_:"vesc"/.Constants`EarthRepl]:=Module[{m\[Chi],\[Beta]C,\[Beta]M,\[Lambda]inve,\[Lambda]invNuc,eregiontable,eregiondict,vdom,nT,Nucregiontable,Nucregiondict},

(* \[Sigma]Dicts should be a single list format (ie. for 1 mass)
	\[Kappa] = 1 *)

m\[Chi]=Union@Join["m\[Chi]"/.Electronic\[Sigma]Dicts,"m\[Chi]"/.Nuclear\[Sigma]Dicts];
If[Length[m\[Chi]]>1,Print["Masses in all \[Sigma]Dicts must be equal."];Return[<|"N/A"->Null|>,Module](*,m\[Chi]=m\[Chi][[1]]*)];

{\[Beta]C,\[Beta]M}={"\[Beta]core","\[Beta]crust"}/.Constants`EarthRepl;

(*vbar = <|2->Log10@Max[Sqrt[2/(m\[Chi] \[Beta]C)],vesc],1->Log10@Max[Sqrt[2/(m\[Chi] \[Beta]M)],vesc]|>;*)
(*vbar = <|2->Log10@Max[Sqrt[2/(m\[Chi] \[Beta]C)],0],1->Log10@Max[Sqrt[2/(m\[Chi] \[Beta]M)],0]|>;*)

eregiontable = Gather[Electronic\[Sigma]Dicts,#1["regionindex"]==#2["regionindex"]&];
eregiondict=Association@Table[Union["regionindex"/.eregiontable[[i]]][[1]]->eregiontable[[i]],{i,Length[eregiontable]}];

\[Lambda]inve =Association@Table[r->Sum[vdom=("\[Sigma]f"/.eregiondict[#2][[o]])["Domain"][[1]];
nT=("ne"/.eregiondict[#2][[o]]["params"]);If[#1>vdom[[1]]&&#1<vdom[[2]],nT 10^("\[Sigma]f"/.eregiondict[#2][[o]])[#1],0],{o,Length[eregiondict[#2]]}]&[Log10@v\[Chi],r],{r,2}];(*\[Lambda]inve[v\[Chi]][r] - r is the regionindex*)

Nucregiontable = Gather[Nuclear\[Sigma]Dicts,#1["regionindex"]==#2["regionindex"]&];
Nucregiondict=Association@Table[Union["regionindex"/.Nucregiontable[[i]]][[1]]->Nucregiontable[[i]],{i,Length[eregiontable]}];

\[Lambda]invNuc=Association@Table[r->Sum[vdom=("\[Sigma]f"/.Nucregiondict[#2][[o]])["Domain"][[1]];
nT=("nI"/.Nucregiondict[#2][[o]]["NucleusParams"]);If[#1>vdom[[1]]&&#1<vdom[[2]],nT 10^("\[Sigma]f"/.Nucregiondict[#2][[o]])[#1],0],{o,Length[Nucregiondict[#2]]}]&[Log10@v\[Chi],r],{r,2}];(*\[Lambda]invNuc[v\[Chi]][r] - r is the regionindex*)

(*Print[\[Lambda]invNuc];
Print[\[Lambda]inve];*)

<|"m\[Chi]"->m\[Chi][[1]],"\[Lambda]C"->(\[Lambda]inve[2]+\[Lambda]invNuc[2])^-1,"\[Lambda]M"->(\[Lambda]inve[1]+\[Lambda]invNuc[1])^-1,"v\[Chi]"->v\[Chi],"\[Lambda]inve"->\[Lambda]inve,"\[Lambda]invNuc"->\[Lambda]invNuc|>
]


(* ::Input::Initialization:: *)
Clear[GetPenetrationDepth]
GetPenetrationDepth[Electronic\[Sigma]Dicts_,Nuclear\[Sigma]Dicts_,vesc_:"vesc"/.Constants`EarthRepl]:=Module[{m\[Chi],\[Beta]C,\[Beta]M,vbar,\[Lambda]inve,\[Lambda]invNuc,eregiontable,eregiondict,vdom,nT,Nucregiontable,Nucregiondict},

(* *\[Sigma]Dicts should be a single list format (ie. for 1 mass) *)

m\[Chi]=Union@Join["m\[Chi]"/.Electronic\[Sigma]Dicts,"m\[Chi]"/.Nuclear\[Sigma]Dicts];
If[Length[m\[Chi]]>1,Print["Masses in all \[Sigma]Dicts must be equal."];Return[<|"N/A"->Null|>,Module](*,m\[Chi]=m\[Chi][[1]]*)];

{\[Beta]C,\[Beta]M}={"\[Beta]core","\[Beta]crust"}/.Constants`EarthRepl;

vbar = <|2->Log10@Max[Sqrt[2/(m\[Chi] \[Beta]C)],vesc],1->Log10@Max[Sqrt[2/(m\[Chi] \[Beta]M)],vesc]|>;
(*vbar = <|2->Log10@Max[Sqrt[2/(m\[Chi] \[Beta]C)],0],1->Log10@Max[Sqrt[2/(m\[Chi] \[Beta]M)],0]|>;*)

eregiontable = Gather[Electronic\[Sigma]Dicts,#1["regionindex"]==#2["regionindex"]&];
eregiondict=Association@Table[Union["regionindex"/.eregiontable[[i]]][[1]]->eregiontable[[i]],{i,Length[eregiontable]}];

\[Lambda]inve =Association@Table[r->Sum[vdom=("\[Sigma]f"/.eregiondict[#2][[o]])["Domain"][[1]];
nT=("ne"/.eregiondict[#2][[o]]["params"]);If[#1>vdom[[1]]&&#1<vdom[[2]],nT 10^("\[Sigma]f"/.eregiondict[#2][[o]])[#1],0],{o,Length[eregiondict[#2]]}]&[vbar[r],r],{r,2}];(*\[Lambda]inve[v\[Chi]][r] - r is the regionindex*)

Nucregiontable = Gather[Nuclear\[Sigma]Dicts,#1["regionindex"]==#2["regionindex"]&];
Nucregiondict=Association@Table[Union["regionindex"/.Nucregiontable[[i]]][[1]]->Nucregiontable[[i]],{i,Length[eregiontable]}];

\[Lambda]invNuc=Association@Table[r->Sum[vdom=("\[Sigma]f"/.Nucregiondict[#2][[o]])["Domain"][[1]];
nT=("nI"/.Nucregiondict[#2][[o]]["NucleusParams"]);If[#1>vdom[[1]]&&#1<vdom[[2]],nT 10^("\[Sigma]f"/.Nucregiondict[#2][[o]])[#1],0],{o,Length[Nucregiondict[#2]]}]&[vbar[r],r],{r,2}];(*\[Lambda]invNuc[v\[Chi]][r] - r is the regionindex*)

(*Print[\[Lambda]invNuc];
Print[\[Lambda]inve];*)

<|"m\[Chi]"->m\[Chi][[1]],"\[Lambda]C"->(\[Lambda]inve[2]+\[Lambda]invNuc[2])^-1,"\[Lambda]M"->(\[Lambda]inve[1]+\[Lambda]invNuc[1])^-1,"vbar"->vbar,"\[Lambda]inve"->\[Lambda]inve,"\[Lambda]invNuc"->\[Lambda]invNuc|>
]


(* ::Input::Initialization:: *)
Clear[EvaporationwMFP]
EvaporationwMFP[\[CapitalGamma]total_,\[Lambda]MFPs_,\[Kappa]_,regionindex_]:=Module[{rC,rE,\[Lambda]eff,EvaporationVolume,RegionVolume},

rC = "rcore"/.Constants`EarthRepl;
rE="rE"/.Constants`EarthRepl;

(*\[Lambda]eff = If[\[Lambda]MFPs["\[Lambda]M"]<rC,\[Lambda]MFPs["\[Lambda]M"],rC + \[Lambda]MFPs["\[Lambda]C"]];*)
(*\[Lambda]eff = If[\[Lambda]MFPs["\[Lambda]M"]/\[Kappa]^2<rC,\[Lambda]MFPs["\[Lambda]M"]/\[Kappa]^2,rC + \[Lambda]MFPs["\[Lambda]C"]/\[Kappa]^2];*)
\[Lambda]eff = If[\[Lambda]MFPs["\[Lambda]M"]/\[Kappa]^2<rE-rC,\[Lambda]MFPs["\[Lambda]M"]/\[Kappa]^2,rE-rC + \[Lambda]MFPs["\[Lambda]C"]/\[Kappa]^2];

(*EvaporationVolume =  4\[Pi] Integrate[r^2Region\[CapitalTheta][r,regionindex],{r,Max[rE-\[Lambda]eff/\[Kappa]^2,0],rE}];*)
EvaporationVolume =  4\[Pi] Integrate[r^2 Region\[CapitalTheta][r,regionindex],{r,Max[rE-\[Lambda]eff,0],rE}]; 

RegionVolume =  4\[Pi] Integrate[r^2 Region\[CapitalTheta][r,regionindex],{r,0,"rE"/.Constants`EarthRepl}]; 

(*\[CapitalGamma]total EvaporationVolume/RegionVolume*)
EvaporationVolume

]


(* ::Input::Initialization:: *)
Clear[EvaporationVolume]
EvaporationVolume[\[Lambda]MFPs_,\[Kappa]_,regionindex_:1]:=Module[{rC,rE,\[Lambda]eff,EvapVolume},

rC = "rcore"/.Constants`EarthRepl;
rE="rE"/.Constants`EarthRepl;

\[Lambda]eff = If[\[Lambda]MFPs["\[Lambda]M"]/\[Kappa]^2<rE-rC,\[Lambda]MFPs["\[Lambda]M"]/\[Kappa]^2,rE-rC + \[Lambda]MFPs["\[Lambda]C"]/\[Kappa]^2];

(*if \[Lambda]eff < rE-rC and regionindex = 2, 0 
if \[Lambda]eff < rE-rC and region index = 1, compute as below
if \[Lambda]eff > rE-rC and region index = 1, Mantle volume
if \[Lambda]eff > rE-rC and region index = 2, rE -> rC below,
if \[Lambda]eff > rE, region volume. *)

(*EvapVolume =  (4 \[Pi])/3rE^3 Piecewise[{{1,rE<\[Lambda]eff},{1-(1-\[Lambda]eff/rE)^3,rE>= \[Lambda]eff}}];*)
EvapVolume =  (4 \[Pi])/3 Piecewise[{{rE^3-rC^3,rE<\[Lambda]eff&&regionindex==1},{rE^3-(rE-\[Lambda]eff)^3,rE-rC>\[Lambda]eff&&regionindex==1},{rE^3-rC^3,rE-rC<\[Lambda]eff&&regionindex==1},
{rC^3,rE<\[Lambda]eff&&regionindex==2},{rE^3-(rE-\[Lambda]eff)^3- (rE^3-rC^3),rE-rC<\[Lambda]eff&&regionindex==2},{0,\[Lambda]eff < rE-rC&&regionindex==2}}];

<|"Vevap"->EvapVolume,"\[Lambda]eff"->\[Lambda]eff|>
]


(* ::Subsubsection::Closed:: *)
(*Rough Evaporation Rate*)


(* ::Input::Initialization:: *)
Clear[GetinvMFT]
GetinvMFT[\[Sigma]Dict_]:=Module[{params,NucleusParams,mN,m\[Chi],regionindex,\[Beta]E,nT,vbar,vesc,\[Sigma]f,v\[Chi]min,fintegrand},
(*Collision time from pg 55 of the appendix *)

NucleusParams=\[Sigma]Dict["NucleusParams"];
m\[Chi]=("m\[Chi]"/.\[Sigma]Dict);
regionindex="regionindex"/.\[Sigma]Dict;
\[Beta]E = {"\[Beta]crust","\[Beta]core"}[[regionindex]]/.Constants`EarthRepl;(*Earth Temeperature*)
nT = "nI"/.NucleusParams;
m\[Chi]=("m\[Chi]"/.\[Sigma]Dict);
vbar=  Sqrt[8/(\[Pi] m\[Chi] \[Beta]E)];
vesc ="vesc"/.Constants`EarthRepl;

\[Sigma]f="\[Sigma]f"/.\[Sigma]Dict;

(*vbar nT 10^\[Sigma]f[Log10@#]&*)
(*nT 10^\[Sigma]f[Log10@#]&*)

v\[Chi]min=\[Sigma]f["Domain"][[1,1]];

(*fintegrand[v\[Chi]_?NumberQ]:=4 \[Pi] ((m\[Chi] \[Beta]E)/(2 \[Pi]))^(3/2) v\[Chi]^2 E^(-\[Beta]E m\[Chi]/2 v\[Chi]^2)10^\[Sigma]f[Log10@v\[Chi]];

nT vbar NIntegrate[fintegrand[v\[Chi]],{v\[Chi],v\[Chi]min,"vesc"/.Constants`EarthRepl}]*)

(*Print[vbar];*)

(*fintegrand[v\[Chi]_?NumberQ]:=4 \[Pi] ((m\[Chi] \[Beta]E)/(2 \[Pi]))^(3/2) v\[Chi]^3 E^(-\[Beta]E m\[Chi]/2 v\[Chi]^2)10^\[Sigma]f[Log10@v\[Chi]];*)(*for MFT*)
(*fintegrand[v\[Chi]_?NumberQ]:=4 \[Pi] ((m\[Chi] \[Beta]E)/(2 \[Pi]))^(3/2) v\[Chi]^2 E^(-\[Beta]E m\[Chi]/2 v\[Chi]^2)10^\[Sigma]f[Log10@v\[Chi]];*)(*for MFP*)
fintegrand[v\[Chi]_?NumberQ]:=4 \[Pi] ((m\[Chi] \[Beta]E)/(2 \[Pi]))^(3/2) v\[Chi]^2 E^(-\[Beta]E m\[Chi]/2 v\[Chi]^2) (v\[Chi]^2/(v\[Chi]^2+vesc^2))^2 10^\[Sigma]f[Log10@v\[Chi]];(*for MFP with FF as in (A.9) of the appendix*)

nT  NIntegrate[fintegrand[v\[Chi]],{v\[Chi],v\[Chi]min,vesc}](*Thermally averaged interaction rate*)

]


(* ::Input::Initialization:: *)
GetTotalinvMFT[\[Sigma]Dicts_]:=Module[{},
Sum[GetinvMFT[\[Sigma]Dicts[[i]]],{i,Length[\[Sigma]Dicts]}]
]


(* ::Input::Initialization:: *)
GetPejection[\[Sigma]Dict_]:=Module[{m\[Chi],regionindex,\[Beta]E,vesc},
(*Eqn A.57 of the appendix*)
m\[Chi]=("m\[Chi]"/.\[Sigma]Dict);
regionindex="regionindex"/.\[Sigma]Dict;
\[Beta]E = {"\[Beta]crust","\[Beta]core"}[[regionindex]]/.Constants`EarthRepl;(*Earth Temeperature*)
vesc=("vesc"/.Constants`EarthRepl);

(1-Sqrt[2/\[Pi]]  (-(m\[Chi] \[Beta]E)^(1/2) E^(-(1/2) m\[Chi] vesc^2 \[Beta]E) vesc+Sqrt[\[Pi]/2] Erf[(Sqrt[m\[Chi]] vesc Sqrt[\[Beta]E])/Sqrt[2]]))
]


(* ::Input::Initialization:: *)
Clear[GetAppendixEvaporationRate]
GetAppendixEvaporationRate[\[Sigma]Dict_]:=Module[{regionindex,RegionVolume,EarthVolume},
(*Subscript[R, evp] from Eqn A.56 of the appendix*)
(*Print[\[Sigma]Dict];*)

regionindex="regionindex"/.\[Sigma]Dict;

RegionVolume =  4\[Pi] Integrate[r^2 Region\[CapitalTheta][r,regionindex],{r,0,"rE"/.Constants`EarthRepl}];
EarthVolume =  (4\[Pi])/3 ("rE"/.Constants`EarthRepl)^3;

(*Print[GetinvMFT[\[Sigma]Dict][10^4]];
*)

(*Print[GetPejection[\[Sigma]Dict]];
*)
(*GetinvMFT[\[Sigma]Dict][#]GetPejection[\[Sigma]Dict]RegionVolume/EarthVolume& (*function of v\[Chi]*)*)
GetinvMFT[\[Sigma]Dict]GetPejection[\[Sigma]Dict] RegionVolume/EarthVolume (*function of v\[Chi]*)
]


(* ::Input::Initialization:: *)
GetTotalAppendixEvaporationRate[\[Sigma]Dicts_]:=Module[{},
(*Sum[GetAppendixEvaporationRate[\[Sigma]Dicts[[i]]][#],{i,Length[\[Sigma]Dicts]}]&*)
Sum[GetAppendixEvaporationRate[\[Sigma]Dicts[[i]]],{i,Length[\[Sigma]Dicts]}]
]


(* ::Subsubsection::Closed:: *)
(*Get Subscript[N, c] in the presence of evaporation*)


(* ::Input::Initialization:: *)
(*Clear[GetNc]
GetNc[NcMaxDictList_,Nuc\[Sigma]Dict_,Electronic\[Sigma]Dicts_,keys_:{"Pf"}]:=Module[{indexforcapturethreshold,\[Kappa],m\[Chi],eind,Nucind,\[Lambda]MFPs,NcMax,dNcMaxdt,\[CapitalGamma]total,namesepchar,keysuffix,tE,timeenoughforequil,Nc,tempNcDict={},NcDicts={}},
\[CapitalGamma]total="\[CapitalGamma]Total"/.IterateGetEvaporationRate[Nuc\[Sigma]Dict];(*evaporation rate per captured particle assuming homogeneous distribution in the Earth - Conservative: it is dominated by the core, the lower g there will reduce the number density*)

keysuffix =If[keys=={"Pf"},"",StringJoin[Table["_"<>keys[[k]],{k,Length[keys]}]]];

Do[

NcMax="NcMax"<>keysuffix/.NcMaxDictList[[i]];
dNcMaxdt="dNcMaxdt"<>keysuffix/.NcMaxDictList[[i]];

tE=NcMax/dNcMaxdt;(*[s] age of the Earth*)

\[Kappa]="\[Kappa]"/.NcMaxDictList[[i]];


(*m\[Chi]=NcMaxDictList[[i]];*)
(*{eind, Nucind}= {FirstPosition["m\[Chi]"/.Electronic\[Sigma]Dicts,m\[Chi]][[1]],FirstPosition["m\[Chi]"/.Nuc\[Sigma]Dict,m\[Chi]][[1]]},*)
(*\[Lambda]MFPs = GetPenetrationDepth[Electronic\[Sigma]Dicts[[eind]],Nuc\[Sigma]Dict[[Nucind]]];*)
\[Lambda]MFPs = GetPenetrationDepth[Electronic\[Sigma]Dicts,Nuc\[Sigma]Dict]; 
(*need to pass region index, also need to be able to test this. So do a test call of this with a PDictList*)

timeenoughforequil=1/tE<(\[Kappa]^2 \[CapitalGamma]total);
Nc=If[timeenoughforequil,dNcMaxdt/(\[Kappa]^2 \[CapitalGamma]total),NcMax ];

tempNcDict = Append[NcMaxDictList[[i]] ,{"timeenoughforequil"<>keysuffix->timeenoughforequil,"Nc"<>keysuffix->Nc,"\[CapitalGamma]total"->\[CapitalGamma]total}];
AppendTo[NcDicts,tempNcDict];

,{i,Length[NcMaxDictList]}];
NcDicts
]*)


(* ::Subsubsection::Closed:: *)
(*Get Subscript[N, c] in the presence of evaporation w Penetration Depth*)


(* ::Input::Initialization:: *)
Clear[GetNc]
GetNc[NcMaxDictList_,Nuc\[Sigma]Dict_,Electronic\[Sigma]Dicts_,ElectronicEvapDict_,keys_:{"Pf"}]:=Module[{\[CapitalGamma]DictNuc,\[CapitalGamma]Dicte,\[CapitalGamma]total,\[CapitalGamma]totalNuc,\[CapitalGamma]totale,indexforcapturethreshold,\[Kappa],UseMFP=True,m\[Chi],eind,Nucind,\[Lambda]MFPs,\[CapitalGamma]totalMFP,regionindex,NcMax,dNcMaxdt,namesepchar,keysuffix,tE,timeenoughforequil,Nc,timeenoughforequilMFP, NcMFP,tempNcDict={},NcDicts={}},

\[CapitalGamma]DictNuc=IterateGetEvaporationRate[Nuc\[Sigma]Dict];
(*\[CapitalGamma]Dicte=IterateGetEvaporationRateend[Electronic\[Sigma]Dicts]; *)
\[CapitalGamma]Dicte=IterateGetEvaporationRateend[ElectronicEvapDict]; 

\[CapitalGamma]totalNuc="\[CapitalGamma]Total"/.\[CapitalGamma]DictNuc;
\[CapitalGamma]totale="\[CapitalGamma]Total"/.\[CapitalGamma]Dicte;(*evaporation rate per captured particle assuming homogeneous distribution in the Earth - Conservative: it is dominated by the core, the lower g there will reduce the number density*)

keysuffix =If[keys=={"Pf"},"",StringJoin[Table["_"<>keys[[k]],{k,Length[keys]}]]];

Do[

NcMax="NcMax"<>keysuffix/.NcMaxDictList[[i]];
dNcMaxdt="dNcMaxdt"<>keysuffix/.NcMaxDictList[[i]];

tE=NcMax/dNcMaxdt;(*[s] age of the Earth*)

\[Kappa]="\[Kappa]"/.NcMaxDictList[[i]];

If[UseMFP,
regionindex = "regionindex"/.NcMaxDictList[[i]];
\[Lambda]MFPs = GetPenetrationDepth[Electronic\[Sigma]Dicts,Nuc\[Sigma]Dict]; 
\[CapitalGamma]totalNuc =Total@Table[EvaporationwMFP[("\[CapitalGamma]s"/.\[CapitalGamma]DictNuc)[[i]]["\[CapitalGamma]"],\[Lambda]MFPs,\[Kappa],("\[CapitalGamma]s"/.\[CapitalGamma]DictNuc)[[i]]["regionindex"]],{i,Length[("\[CapitalGamma]s"/.\[CapitalGamma]DictNuc)]}];
\[CapitalGamma]totale=Total@Table[EvaporationwMFP[("\[CapitalGamma]s"/.\[CapitalGamma]Dicte)[[i]]["\[CapitalGamma]"],\[Lambda]MFPs,\[Kappa],("\[CapitalGamma]s"/.\[CapitalGamma]Dicte)[[i]]["regionindex"]],{i,Length[("\[CapitalGamma]s"/.\[CapitalGamma]Dicte)]}]
];

\[CapitalGamma]total = \[CapitalGamma]totalNuc+\[CapitalGamma]totale;

(*Print["\[Kappa] :",N@Log10@\[Kappa]];
Print[\[Kappa]^2\[CapitalGamma]total];*)
timeenoughforequil=1/tE<(\[Kappa]^2 \[CapitalGamma]total);
Nc=If[timeenoughforequil,dNcMaxdt/(\[Kappa]^2 \[CapitalGamma]total),NcMax ];

(*timeenoughforequilMFP=1/tE<(\[Kappa]^2 \[CapitalGamma]totalMFP);
NcMFP=If[timeenoughforequilMFP,dNcMaxdt/(\[Kappa]^2 \[CapitalGamma]totalMFP),NcMax ];*)

tempNcDict = Append[NcMaxDictList[[i]] ,{"timeenoughforequil"<>keysuffix->timeenoughforequil,"Nc"<>keysuffix->Nc,"\[CapitalGamma]total"->\[CapitalGamma]total,"\[CapitalGamma]totalNuc"->\[CapitalGamma]totalNuc,"\[CapitalGamma]totale"->\[CapitalGamma]totale}];
AppendTo[NcDicts,tempNcDict];

,{i,Length[NcMaxDictList]}];

NcDicts
]


(* ::Subsubsection::Closed:: *)
(*Scan parameter space - old*)


(* ::Input::Initialization:: *)
Clear[ScanGetCaptureRateandGetNc]
ScanGetCaptureRateandGetNc[meDratios_,\[Kappa]s_,v0s_,Electronic\[Sigma]Dicsts_,Nuclear\[Sigma]Dicts_]:=Module[{PDictList,InterPDicts,PDictLists={},ScanDictList={},meDratio,v0,mpD,\[Beta]D,NcMaxDict,NcDict,MemoryUsed={}},
mpD= "m\[Chi]"/.Electronic\[Sigma]Dicsts[[1]];

Print["mpD: ",mpD ( ("JpereV")/("c")^2)^-1/.Constants`SIConstRepl];
Do[
meDratio=("meD/mpD"/.meDratios[[m]]);
v0 = v0s[[n]];
\[Beta]D ="\[Beta]D"/.v0to\[Beta]D[v0,meDratio mpD,mpD];
Print["meD/mpD: ",meDratio];
Print["v0: ",v0];
Print["\[Beta]D: ",\[Beta]D];

PDictList=IterateGetCaptureRate[meDratio,\[Kappa]s,v0,Electronic\[Sigma]Dicsts,Nuclear\[Sigma]Dicts];(*uses forward Euler*)

PlotPDist[PDictList];

InterPDicts=InterpolatePc[PDictList];

NcMaxDict=GetNcMax[InterPDicts];

NcDict=GetNc[NcMaxDict,Nuclear\[Sigma]Dicts];

AppendTo[ScanDictList, NcDict];

,{m,Length[meDratios]},{n,Length[v0s]}];

ScanDictList
]


(* ::Subsubsection::Closed:: *)
(*Scan parameter space*)


(* ::Input::Initialization:: *)
Clear[ScanGetCapture]
ScanGetCapture[meDratios_,\[Kappa]s_,v0s_,Electronic\[Sigma]Dicts_,Nuclear\[Sigma]Dicts_,mpDcap_:True,ElectronicEvapDict_:<|"NA"->Null|>]:=Module[{PDictList,InterPDicts,PDictLists={},ScanDictList={},meDratio,v0,mpD,meD,m\[Chi],\[Beta]D,NcMaxDict,NcDict},
m\[Chi]= "m\[Chi]"/.Electronic\[Sigma]Dicts[[1]];

(*Print["m\[Chi]: ",m\[Chi]( ("JpereV")/("c")^2)^-1/.Constants`SIConstRepl];*)
Do[
meDratio=("meD/mpD"/.meDratios[[me]]);
If[mpDcap,
mpD = m\[Chi];meD = meDratio m\[Chi];,
meD = m\[Chi];mpD =  m\[Chi]/meDratio;
];
(*meDratio=("meD/mpD"/.meDratios[[me]]);*)
v0 = v0s[[vs]];
\[Beta]D ="\[Beta]D"/.v0to\[Beta]D[v0,meD,mpD];

(*If[mpDcap,
Print["meD: ",meD( ("JpereV")/("c")^2)^-1/.Constants`SIConstRepl];,
Print["mpD: ",mpD( ("JpereV")/("c")^2)^-1/.Constants`SIConstRepl];
];*)
(*Print["v0: ",v0];
Print["\[Beta]D: ",\[Beta]D];*)

PDictList=IterateGetCaptureProbability[meDratio,\[Kappa]s,v0,Electronic\[Sigma]Dicts,Nuclear\[Sigma]Dicts,mpDcap];(*uses Forward Euler*)

PlotPDist[PDictList];

InterPDicts=InterpolatePc[PDictList];

NcMaxDict=GetNcMax[InterPDicts];

NcDict=GetNc[NcMaxDict,Nuclear\[Sigma]Dicts,Electronic\[Sigma]Dicts,ElectronicEvapDict];(*need to pass both nuclear and electronic \[Sigma]dicts. They will be used both for penetration depth, and for nuclear and electronic evaporation.*)

AppendTo[ScanDictList, NcDict];

,{me,Length[meDratios]},{vs,Length[v0s]}];

(*ScanDictList=IterateKEchange[Flatten@ScanDictList,{0.05},"\[Alpha]"/.Constants`SIConstRepl];*)

ScanDictList
]


(* ::Subsection::Closed:: *)
(*Get Significant Capture Regimes*)


(* ::Subsubsection::Closed:: *)
(*KE Change *)


(* ::Input::Initialization:: *)
VCoulombE[N_,\[Alpha]D_]:=\[Alpha]D/("\[Alpha]") ("e")^2/(4 \[Pi] "\[Epsilon]0") N/("rE")/.Constants`EarthRepl/.Constants`SIConstRepl


(* ::Input::Initialization:: *)
Clear[KEchange]
KEchange[ScanDict_,fD_,\[Alpha]D_]:=Module[{mpD,meD,nD0,Nc,Ntot,pDvescchange,eDvescchange,pDKEfracchange,eDKEfracchange},
{mpD,meD}={"mpD","meD"}/.ScanDict;
(*Print[mpD," ",meD];*)
(*nD0=nD0\[CapitalLambda]CDM[fD,mpD +meD];*)
nD0=naDM[fD,mpD +meD];

Nc="Nc"/.ScanDict;
Ntot = Nc/."nD"->nD0;

pDvescchange=Sqrt[(2 VCoulombE[Ntot,\[Alpha]D])/mpD]; 
eDvescchange=Sqrt[(2 VCoulombE[Ntot,\[Alpha]D])/meD];
pDKEfracchange=pDvescchange/Sqrt[("vesc")^2+3/("\[Beta]D" mpD)]/.Constants`EarthRepl/.ScanDict;
eDKEfracchange=eDvescchange/Sqrt[("vesc")^2+3/("\[Beta]D" meD)]/.Constants`EarthRepl/.ScanDict;
(*Print[pDvescchange/Sqrt[("vesc")^2+3/("\[Beta]D" mpD)]/.Constants`EarthRepl/.ScanDict];*)
Append[ScanDict,{"pDvescchange"->pDvescchange,"eDvescchange"->eDvescchange,"\[CapitalDelta]pKE/vpbar"->pDKEfracchange,"\[CapitalDelta]eKE/vpbar"->eDKEfracchange,"Ntot"->Ntot,"nD"->nD0,"fD"->fD,"\[Alpha]D"->\[Alpha]D}]
]


(* ::Input::Initialization:: *)
Clear[IterateKEchange]
IterateKEchange[ScanDicts_,fDs_,\[Alpha]D_]:=Module[{KEchangeDicts={}},
Do[
AppendTo[KEchangeDicts,KEchange[ScanDicts[[i]],fDs[[j]],\[Alpha]D]];
,{i,Length[ScanDicts]},{j,Length[fDs]}];
KEchangeDicts
]


(* ::Chapter:: *)
(*End*)


End[];


EndPackage[];
