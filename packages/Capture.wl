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


(* ::Section:: *)
(*Definitions*)


(* ::Subsection:: *)
(*Utilities*)


(* ::Subsubsection::Closed:: *)
(*Profiler*)


(* ::Input:: *)
(*Clear[Profiler]*)
(*Profiler[fun_,inputs_List]:=TableForm[#&/@Join@@@({AbsoluteTiming[fun@@#],#}&/@inputs)]*)


(* ::Subsubsection:: *)
(*Convert from v0 to \[Beta]D*)


(* ::Input::Initialization:: *)
Clear[v0to\[Beta]D]
v0to\[Beta]D[v0_,meD_,mpD_]:=Module[{\[Beta]D,vpDmean,veDmean},
(*
v0 is the velocity corresponding to the temperature of the aDM particles. By equipartition: 3 Subscript[T, D] = 1/2(Subscript[m, Subscript[e, D]]Subscript[v^2, Subscript[e, D]]+Subscript[m, Subscript[p, D]]Subscript[v^2, Subscript[p, D]]) so we define 3 Subscript[T, D] = Overscript[m, _] Subscript[v^2, 0] 
*)
\[Beta]D = 6/((meD + mpD)v0^2);
vpDmean = Sqrt[3/(mpD \[Beta]D)];s
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


(* ::Subsubsection::Closed:: *)
(*Gather PDicts*)


(* ::Input::Initialization:: *)
sortbymeDmpDv0[pdictscans_]:=Table[Gather[Flatten@pdictscans[[i]],({"meD"/"mpD","v0"}/.#1)==({"meD"/"mpD","v0"}/.#2)&],{i,Length[pdictscans]}](*required formatting for rerunning analysis on existing PDicts*)


(* ::Subsection:: *)
(*Plotting*)


(* ::Subsubsection::Closed:: *)
(*Plot Subscript[P, cap]*)


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

ExportToDir[plot,ToString@StringForm["Pcap_total_Log10mpD_``_Log10\[Kappa]_``.png",Round@Log10@mpD,#]&@Round@Log10@\[Kappa],ToString@StringForm["Pcap_total_Log10mpD_``",Round@Log10@mpD]];

(*Electronic*)
Pfunc="PfE"/.Pfs[[p]];
plot=Show[{DensityPlot[Pfunc[v\[Chi],Log10@bE],{v\[Chi],Pfunc["Domain"][[1,1]],Pfunc["Domain"][[1,2]]},{bE,10^Pfunc["Domain"][[2,1]],10^Pfunc["Domain"][[2,2]]},PlotLegends->Automatic,AxesLabel->Automatic,PlotRange->All,FrameLabel->{"\!\(\*SubscriptBox[\(Log\), \(10\)]\)v\[Chi] [\!\(\*SuperscriptBox[\(ms\), \(-1\)]\)]","\!\(\*SubscriptBox[\(b\), \(E\)]\) [m]"},PlotLabel->"\!\(\*SubscriptBox[\(Log\), \(10\)]\)\!\(\*SubscriptBox[\(P\), \(cap, e\)]\) - \!\(\*SubscriptBox[\(Log\), \(10\)]\)"<>ToString[StringForm["\[Kappa] = ``",Log10@\[Kappa]]]<>", \!\(\*SubscriptBox[\(Log\), \(10\)]\)\!\(\*SubscriptBox[\(m\), SubscriptBox[\(p\), \(D\)]]\)"<>ToString[StringForm[" = `` [eV]",Log10@mpD]]<>ToString[StringForm[", v0 = `` [km/s]",v0/1000]]],ContourPlot[Pfunc[v\[Chi],Log10@bE],{v\[Chi],Pfunc["Domain"][[1,1]],Pfunc["Domain"][[1,2]]},{bE,10^Pfunc["Domain"][[2,1]],10^Pfunc["Domain"][[2,2]]},ContourShading->False,ContourLabels->True]}];

ExportToDir[plot,StringForm["Pcap_electronic_Log10mpD_``_Log10\[Kappa]_``.png",Round[Log10@mpD],#]&@Round@Log10@\[Kappa],ToString@StringForm["Pcap_electronic_Log10mpD_``",Round[Log10@mpD]]];

(*Nuclear*)
Pfunc="PfNuc"/.Pfs[[p]];
plot=Show[{DensityPlot[Pfunc[v\[Chi],Log10@bE],{v\[Chi],Pfunc["Domain"][[1,1]],Pfunc["Domain"][[1,2]]},{bE,10^Pfunc["Domain"][[2,1]],10^Pfunc["Domain"][[2,2]]},PlotLegends->Automatic,AxesLabel->Automatic,PlotRange->All,FrameLabel->{"\!\(\*SubscriptBox[\(Log\), \(10\)]\)v\[Chi] [\!\(\*SuperscriptBox[\(ms\), \(-1\)]\)]","\!\(\*SubscriptBox[\(b\), \(E\)]\) [m]"},PlotLabel->"\!\(\*SubscriptBox[\(Log\), \(10\)]\)\!\(\*SubscriptBox[\(P\), \(cap, Nuc\)]\) - \!\(\*SubscriptBox[\(Log\), \(10\)]\)"<>ToString[StringForm["\[Kappa] = ``",Log10@\[Kappa]]]<>", \!\(\*SubscriptBox[\(Log\), \(10\)]\)\!\(\*SubscriptBox[\(m\), SubscriptBox[\(p\), \(D\)]]\)"<>ToString[StringForm[" = `` [eV]",Log10@mpD]]<>ToString[StringForm[", v0 = `` [km/s]",v0/1000]]],ContourPlot[Pfunc[v\[Chi],Log10@bE],{v\[Chi],Pfunc["Domain"][[1,1]],Pfunc["Domain"][[1,2]]},{bE,10^Pfunc["Domain"][[2,1]],10^Pfunc["Domain"][[2,2]]},ContourShading->False,ContourLabels->True]}];

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

ExportToDir[plot,ToString@StringForm["Pcap_total_Log10mpD_``_Log10\[Kappa]_``.png",Round@Log10@mpD,#]&@Round@Log10@\[Kappa],ToString@StringForm["Pcap_total_Log10mpD_``",Round@Log10@mpD]];

(*Electronic*)
Pfunc="PfE"/.Pfs[[p]];
plot=Show[{DensityPlot[Pfunc[v\[Chi],Log10@bE],{v\[Chi],Pfunc["Domain"][[1,1]],Pfunc["Domain"][[1,2]]},{bE,10^Pfunc["Domain"][[2,1]],10^Pfunc["Domain"][[2,2]]},PlotLegends->Automatic,AxesLabel->Automatic,PlotRange->All,FrameLabel->{"\!\(\*SubscriptBox[\(Log\), \(10\)]\)v\[Chi] [\!\(\*SuperscriptBox[\(ms\), \(-1\)]\)]","\!\(\*SubscriptBox[\(b\), \(E\)]\) [m]"},PlotLabel->"\!\(\*SubscriptBox[\(Log\), \(10\)]\)\!\(\*SubscriptBox[\(P\), \(cap, e\)]\) - \!\(\*SubscriptBox[\(Log\), \(10\)]\)"<>ToString[StringForm["\[Kappa] = ``",Log10@\[Kappa]]]<>", \!\(\*SubscriptBox[\(Log\), \(10\)]\)\!\(\*SubscriptBox[\(m\), SubscriptBox[\(p\), \(D\)]]\)"<>ToString[StringForm[" = `` [eV]",Log10@mpD]]<>ToString[StringForm[", v0 = `` [km/s]",v0/1000]]],ContourPlot[Pfunc[v\[Chi],Log10@bE],{v\[Chi],Pfunc["Domain"][[1,1]],Pfunc["Domain"][[1,2]]},{bE,10^Pfunc["Domain"][[2,1]],10^Pfunc["Domain"][[2,2]]},ContourShading->False,ContourLabels->True]}];

ExportToDir[plot,StringForm["Pcap_electronic_Log10mpD_``_Log10\[Kappa]_``.png",Round[Log10@mpD],#]&@Round@Log10@\[Kappa],ToString@StringForm["Pcap_electronic_Log10mpD_``",Round[Log10@mpD]]];

(*Nuclear*)
Pfunc="PfNuc"/.Pfs[[p]];
plot=Show[{DensityPlot[Pfunc[v\[Chi],Log10@bE],{v\[Chi],Pfunc["Domain"][[1,1]],Pfunc["Domain"][[1,2]]},{bE,10^Pfunc["Domain"][[2,1]],10^Pfunc["Domain"][[2,2]]},PlotLegends->Automatic,AxesLabel->Automatic,PlotRange->All,FrameLabel->{"\!\(\*SubscriptBox[\(Log\), \(10\)]\)v\[Chi] [\!\(\*SuperscriptBox[\(ms\), \(-1\)]\)]","\!\(\*SubscriptBox[\(b\), \(E\)]\) [m]"},PlotLabel->"\!\(\*SubscriptBox[\(Log\), \(10\)]\)\!\(\*SubscriptBox[\(P\), \(cap, Nuc\)]\) - \!\(\*SubscriptBox[\(Log\), \(10\)]\)"<>ToString[StringForm["\[Kappa] = ``",Log10@\[Kappa]]]<>", \!\(\*SubscriptBox[\(Log\), \(10\)]\)\!\(\*SubscriptBox[\(m\), SubscriptBox[\(p\), \(D\)]]\)"<>ToString[StringForm[" = `` [eV]",Log10@mpD]]<>ToString[StringForm[", v0 = `` [km/s]",v0/1000]]],ContourPlot[Pfunc[v\[Chi],Log10@bE],{v\[Chi],Pfunc["Domain"][[1,1]],Pfunc["Domain"][[1,2]]},{bE,10^Pfunc["Domain"][[2,1]],10^Pfunc["Domain"][[2,2]]},ContourShading->False,ContourLabels->True]}];

ExportToDir[plot,ToString@StringForm["Pcap_nuclear_Log10mpD_``_Log10\[Kappa]_``.png",Round@Log10@mpD,#]&@Round@Log10@\[Kappa],ToString@StringForm["Pcap_nuclear_Log10mpD_``",Round@Log10@mpD]];


ClearAll[Pfunc,plot];
,{p,Length[Pfs]}];
Clear[Pfs];
]


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

v0=Union["v0"/.AllbutMassand\[Kappa]Fixed[[j]]][[1]];

meDovermpD=Union[("meD")/("mpD")/.AllbutMassand\[Kappa]Fixed[[j]]][[1]];

plotname=("\!\(\*SubscriptBox[\(Log\), \(10\)]\)\!\(\*FractionBox[SubscriptBox[\(\[CapitalDelta]v\), \(e\)], SubscriptBox[OverscriptBox[\(v\), \(_\)], \(e\)]]\)\!\(\*SqrtBox[FractionBox[\(\[Alpha]\\\ 0.05\), \(\*SubscriptBox[\(\[Alpha]\), \(D\)] \*SubscriptBox[\(f\), \(D\)]\)]]\) - "<>" \!\(\*SubscriptBox[\(v\), \(0\)]\)"<>ToString[StringForm["= ``",N@(v0/1000)]]<>"[\!\(\*SuperscriptBox[\(kms\), \(-1\)]\)]"<>",  \!\(\*SubscriptBox[\(Log\), \(10\)]\)\!\(\*FractionBox[SubscriptBox[\(m\), SubscriptBox[\(e\), \(D\)]], SubscriptBox[\(m\), SubscriptBox[\(p\), \(D\)]]]\)="<>ToString[Log10@meDovermpD]);

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
v0=Union["v0"/.AllbutMassand\[Kappa]Fixed[[j]]][[1]];

meDovermpD=Union["meD"/"mpD"/.AllbutMassand\[Kappa]Fixed[[j]]][[1]];

plotname=("\!\(\*SubscriptBox[\(Log\), \(10\)]\)\!\(\*SubscriptBox[\(\[CapitalGamma]\), \(cap\)]\)\!\(\*SubscriptBox[\(n\), \(aDM\)]\)\!\(\*FractionBox[\(0.05\), SubscriptBox[\(f\), \(D\)]]\) [\!\(\*SuperscriptBox[\(m\), \(-3\)]\)\!\(\*SuperscriptBox[\(s\), \(-1\)]\)] - "<>" \!\(\*SubscriptBox[\(v\), \(0\)]\)"<>ToString[StringForm["= ``",N@(v0/1000)]]<>"[\!\(\*SuperscriptBox[\(kms\), \(-1\)]\)]"<>",  \!\(\*SubscriptBox[\(Log\), \(10\)]\)\!\(\*FractionBox[SubscriptBox[\(m\), SubscriptBox[\(e\), \(D\)]], SubscriptBox[\(m\), SubscriptBox[\(p\), \(D\)]]]\)="<>ToString[Log10@meDovermpD]);

VEarth=((4\[Pi])/3)("rE")^3/.Constants`EarthRepl;

Do[With[{k=i,l=j},AssociateTo[AllbutMassand\[Kappa]Fixed[[l,k]],"dNcMaxdt"->(AllbutMassand\[Kappa]Fixed[[l,k]]["dNcMaxdt"]/.AllbutMassand\[Kappa]Fixed[[l,k]])]],{i,Length[AllbutMassand\[Kappa]Fixed[[j]]]}];

Print[Show[{ListDensityPlot[{Log10@("mpD"(("c")^2/"JpereV")/.Constants`SIConstRepl),Log10@"\[Kappa]",Log10@(("dNcMaxdt")/VEarth)}/.AllbutMassand\[Kappa]Fixed[[j]],FrameLabel->{"\!\(\*SubscriptBox[\(Log\), \(10\)]\)mpD [eV]","\[Kappa] [ ]"},PlotLabel->plotname,PlotLegends->Automatic],ListContourPlot[{Log10@("mpD"(("c")^2/"JpereV")/.Constants`SIConstRepl),Log10@"\[Kappa]",Log10@(("dNcMaxdt")/VEarth)}/.AllbutMassand\[Kappa]Fixed[[j]],ContourShading->False,ContourLabels->True,Contours->Join[Table[i,{i,-15,5,2}],Table[i,{i,-15,-65,-5}]]]}]];


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
v0=Union["v0"/.AllbutMassand\[Kappa]Fixed[[j]]][[1]];

meDovermpD=Union["meD"/"mpD"/.AllbutMassand\[Kappa]Fixed[[j]]][[1]];

VEarth=((4\[Pi])/3)("rE")^3/.Constants`EarthRepl;

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
Pfunc=\[Kappa]^2 ((GeteELFunction[Electronic\[Sigma]Dicsts]["dEdlofr"][#1,#2])/(EK[mpD,#2](("JpereV")/("c")^2/.Constants`SIConstRepl)))&;
Print[Show[{DensityPlot[Log10@Pfunc[bE,10^v\[Chi]],{v\[Chi],PDom[[1,1]],PDom[[1,2]]},{bE,10^PDom[[2,1]],10^PDom[[2,2]]},PlotLegends->Automatic,AxesLabel->Automatic,PlotRange->All,FrameLabel->{"\!\(\*SubscriptBox[\(Log\), \(10\)]\)v\[Chi] [\!\(\*SuperscriptBox[\(ms\), \(-1\)]\)]","\!\(\*SubscriptBox[\(b\), \(E\)]\) [m]"},PlotLabel->"\!\(\*SubscriptBox[\(Log\), \(10\)]\)\!\(\*FractionBox[SubscriptBox[\(\[CapitalDelta]E\), \(soft, e\)], SubscriptBox[\(E\), \(\[Chi]\)]]\) - \!\(\*SubscriptBox[\(Log\), \(10\)]\)"<>ToString[StringForm["\[Kappa] = ``",Log10@\[Kappa]]]<>", \!\(\*SubscriptBox[\(Log\), \(10\)]\)\!\(\*SubscriptBox[\(m\), SubscriptBox[\(p\), \(D\)]]\)"<>ToString[StringForm[" = `` [eV]",Log10@mpD]]<>ToString[StringForm[", v0 = `` [m/s]",v0]]],ContourPlot[Log10@Pfunc[bE,10^v\[Chi]],{v\[Chi],PDom[[1,1]],PDom[[1,2]]},{bE,10^PDom[[2,1]],10^PDom[[2,2]]},ContourShading->False,ContourLabels->True]}]];


(*Nuclear*)
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

PfuncE=\[Kappa]^2 (GeteELFunction[Electronic\[Sigma]Dicsts]["dEdlofr"][#1,#2])&;
PfuncNuc=\[Kappa]^2 (GetNucELFunction[Nuclear\[Sigma]Dicts]["dEdlofr"][#1,#2])&;

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

Print[Show[{ListDensityPlot[{Log10@("mpD"(("c")^2/"JpereV")/.Constants`SIConstRepl),Log10@"\[Kappa]",Log10@Max[(1/VEarth)("dNcMaxdt_PfNuc"-"dNcMaxdt_int\[Lambda]invNuc"),10^-100]}/.AllbutMassand\[Kappa]Fixed[[j]],FrameLabel->{"\!\(\*SubscriptBox[\(Log\), \(10\)]\)mpD [eV]","\[Kappa] [ ]"},PlotLabel->plotname,PlotLegends->Automatic,PlotRange->All],ListContourPlot[{Log10@("mpD"(("c")^2/"JpereV")/.Constants`SIConstRepl),Log10@"\[Kappa]",Log10@Max[(1/VEarth)("dNcMaxdt_PfNuc"-"dNcMaxdt_int\[Lambda]invNuc"),10^-100]}/.AllbutMassand\[Kappa]Fixed[[j]],ContourShading->False,ContourLabels->True(*,Contours->Join[Table[i,{i,-21,-1,2}],Table[i,{i,-25,-71,-5}]]*),PlotRange->All],ListContourPlot[{Log10@("mpD" ("c")^2/("JpereV")/.Constants`SIConstRepl),Log10@"\[Kappa]",Log10@Max[(1/VEarth)("dNcMaxdt_PfNuc"-"dNcMaxdt_int\[Lambda]invNuc"),10^-100]}/.AllbutMassand\[Kappa]Fixed[[j]],ContourShading->False,ContourLabels->False,Contours->{0},(*ContourStyle->Directive[Red,Thick]*)ContourStyle->Directive[Red,AbsoluteThickness[5]]]}]];

plotname=("\!\(\*SubscriptBox[\(Log\), \(10\)]\)\!\(\*SubscriptBox[\(\[CapitalGamma]\), \(cap, soft, e\)]\) [\!\(\*SuperscriptBox[\(m\), \(-3\)]\)\!\(\*SuperscriptBox[\(s\), \(-1\)]\)] - "<>" \!\(\*SubscriptBox[\(v\), \(0\)]\)"<>ToString[StringForm["= ``",N@(v0/1000)]]<>"[\!\(\*SuperscriptBox[\(kms\), \(-1\)]\)]"<>",  \!\(\*SubscriptBox[\(Log\), \(10\)]\)\!\(\*FractionBox[SubscriptBox[\(m\), SubscriptBox[\(e\), \(D\)]], SubscriptBox[\(m\), SubscriptBox[\(p\), \(D\)]]]\)="<>ToString[Log10@meDovermpD]);

Print[Show[{ListDensityPlot[{Log10@("mpD"(("c")^2/"JpereV")/.Constants`SIConstRepl),Log10@"\[Kappa]",Log10@Max[("dNcMaxdt_PfE"-"dNcMaxdt_int\[Lambda]invE")/VEarth,10^-100]}/.AllbutMassand\[Kappa]Fixed[[j]],FrameLabel->{"\!\(\*SubscriptBox[\(Log\), \(10\)]\)mpD [eV]","\[Kappa] [ ]"},PlotLabel->plotname,PlotLegends->Automatic,PlotRange->All],ListContourPlot[{Log10@("mpD"(("c")^2/"JpereV")/.Constants`SIConstRepl),Log10@"\[Kappa]",Log10@Max[("dNcMaxdt_PfE"-"dNcMaxdt_int\[Lambda]invE")/VEarth,10^-100]}/.AllbutMassand\[Kappa]Fixed[[j]],ContourShading->False,ContourLabels->True(*,Contours->Join[Table[i,{i,-21,-1,2}],Table[i,{i,-25,-71,-5}]]*),PlotRange->All],ListContourPlot[{Log10@("mpD" ("c")^2/("JpereV")/.Constants`SIConstRepl),Log10@"\[Kappa]",Log10@Max[("dNcMaxdt_PfE"-"dNcMaxdt_int\[Lambda]invE")/VEarth,10^-100]}/.AllbutMassand\[Kappa]Fixed[[j]],ContourShading->False,ContourLabels->False,Contours->{0},(*ContourStyle->Directive[Red,Thick]*)ContourStyle->Directive[Red,AbsoluteThickness[5]]]}]];

,{j,Length@AllbutMassand\[Kappa]Fixed}];
(*,{j,1}];*)
]


(* ::Input:: *)
(**)


(* ::Subsubsection::Closed:: *)
(*Plot trajectory*)


(* ::Input::Initialization:: *)
Clear[Plottraj]
Plottraj[ICtable_(*,Electronic\[Sigma]Dicts_,Nuclear\[Sigma]Dicts_*),N_:140]:=Module[{temptrajdict},

If[AssociationQ@ICtable,
temptrajdict=ICtable,
temptrajdict = (Print[#[[1]]];#[[2]])&[AbsoluteTiming@GetWeakCaptureTrajectoryinvRK2wPotential[Sequence@@ICtable(*,"dEdlofr"/.GetTotalELFunction[Electronic\[Sigma]Dicts,Nuclear\[Sigma]Dicts]*),N]];
];

Print["lMCs"/.temptrajdict];
Print["steps"/.temptrajdict];
Print[ListPlot[{"l","r"}/.temptrajdict["traj"],PlotLabel->"r(l)"]];
Print[ListPlot[{"l","v"}/.temptrajdict["traj"],PlotLabel->"v(l)"]];
Print[ListPlot[{"r","v"}/.temptrajdict["traj"],PlotLabel->"v(r)"]];
Print[ListPlot[{"r","vr"}/.temptrajdict["traj"],PlotLabel->"\!\(\*SubscriptBox[\(v\), \(r\)]\)(r)"]];
Print[ListPlot[{"l","vr"}/.temptrajdict["traj"],PlotLabel->"\!\(\*SubscriptBox[\(v\), \(r\)]\)(l)"]];
Print[ListPlot[{"l",Abs@"vr"}/.temptrajdict["traj"],PlotLabel->"|\!\(\*SubscriptBox[\(v\), \(r\)]\)(l)|",PlotRange->All]];
Print[ListPlot[{"l",("L")^2/(("m\[Chi]")^2 ("r")^2) 1/(("v")^2-("vr")^2)}/.temptrajdict["traj"]/.temptrajdict,PlotLabel->"\!\(\*SubscriptBox[\(v\), \(\[Perpendicular]\)]\)(L) / \!\(\*SubscriptBox[\(v\), \(\[Perpendicular]\)]\)(vr)",PlotRange->All]];
Print[ListPlot[{"l",Sqrt[("v")^2-("L")^2/(("m\[Chi]")^2 ("r")^2)]}/.temptrajdict["traj"]/.temptrajdict,PlotLabel->"vr(L)",PlotRange->All]];
Print[ListPlot[{"l",Sqrt[("v")^2]}/.temptrajdict["traj"]/.temptrajdict,PlotLabel->"v(L)"]];
Print[ListPlot[{"l",Sqrt[("L")^2/(("m\[Chi]")^2 ("r")^2)]}/.temptrajdict["traj"]/.temptrajdict,PlotLabel->"\!\(\*SubscriptBox[\(v\), \(\[Perpendicular]\)]\)(L)",PlotRange->All]];
Print[ListPlot[{"l",Log10@(("v")/(("L")/("m\[Chi]""r")))}/.temptrajdict["traj"]/.temptrajdict,PlotLabel->"\!\(\*FractionBox[\(v\), SubscriptBox[\(v\), \(\[Perpendicular]\)]]\)",PlotRange->All]];
Print[ListPlot[{"l",("\[CapitalDelta]E")/(("m\[Chi]")/2 ("v")^2+"V")}/.temptrajdict["traj"]/.temptrajdict,PlotLabel->"\!\(\*FractionBox[\(\[CapitalDelta]E \((l)\)\), \(H \((l)\)\)]\)"]];
Print[Plot[temptrajdict["v(l)"][l],{l,0,temptrajdict["v(l)"]["Domain"][[1,2]]},PlotLabel->"v(l) - interpolated"]];
Print[Plot[temptrajdict["r(l)"][r],{r,0,temptrajdict["r(l)"]["Domain"][[1,2]]},PlotLabel->"r(l) - interpolated"]];
Print[ListPlot[{"l",("m\[Chi]")/2 ("v")^2+"V"(*-"\[CapitalDelta]E"*)}/.temptrajdict["traj"]/.temptrajdict,PlotRange->All,PlotLabel->"H(l)"]];
Print["check if FE is conserving \!\(\*SubscriptBox[\(L\), \(0\)]\) and \!\(\*SubscriptBox[\(E\), \(0\)]\)"];
Print[ListPlot[{"l",(("m\[Chi]")/2 ("v")^2+"V"+"\[CapitalDelta]E")/("E0")}/.temptrajdict["traj"]/.temptrajdict,PlotLabel->"\!\(\*FractionBox[\(\[CapitalDelta]E\\\  + \\\ H\), SubscriptBox[\(E\), \(0\)]]\)"]];
Print[ListPlot[{"l",(("m\[Chi]")/2 ("v")^2+"V")}/.temptrajdict["traj"]/.temptrajdict,PlotLabel->"H"]];
Print[ListPlot[{"l",("\[CapitalDelta]E")}/.temptrajdict["traj"]/.temptrajdict,PlotLabel->"\[CapitalDelta]E"]];
Print[ListPlot[{"l",("L"Exp["Log\[CapitalDelta]L"])/("L0")}/.temptrajdict["traj"]/.temptrajdict,PlotLabel->"\!\(\*FractionBox[\(L\\\ \*SuperscriptBox[\(e\), \(-\(\[Integral]\*FractionBox[\(dl\), \(2  T\)] \*FractionBox[\(dE\), \(dl\)]\)\)]\), SubscriptBox[\(L\), \(0\)]]\)"]];
Print[ListPlot[{"l","L"}/.temptrajdict["traj"]/.temptrajdict,PlotLabel->"L"]];
Print[ListPlot[{"l",("m\[Chi]""r" Sqrt[("v")^2-("vr")^2])/("L0")}/.temptrajdict["traj"]/.temptrajdict,PlotLabel->("L")/("L0")]];
Print[ListPlot[{"l",Exp@"Log\[CapitalDelta]L"}/.temptrajdict["traj"]/.temptrajdict,PlotLabel->"\!\(\*SuperscriptBox[\(e\), \(-\(\[Integral]\*FractionBox[\(dl\), \(2  T\)] \*FractionBox[\(dE\), \(dl\)]\)\)]\)"]];
Print["investigate conditions for capture"];
Print[ListPlot[{"l",(("m\[Chi]")/2 ("v")^2+"V")}/.temptrajdict["traj"]/.temptrajdict,PlotLabel->"H(l)"]];
Print[ListPlot[{"l",("E0"-"V"-"\[CapitalDelta]E")}/.temptrajdict["traj"]/.temptrajdict,PlotLabel->"\!\(\*FractionBox[\(m\[Chi]\), \(2\)]\)\!\(\*SuperscriptBox[\(v\), \(2\)]\)(l)"]];
Print["\!\(\*SubscriptBox[\(r\), \(final\)]\) and \!\(\*SubscriptBox[\(r\), \(E\)]\)"];
Print["r"/.temptrajdict["traj"][[-1]]];
Print["rE"/.EarthRepl];
]


(* ::Subsection:: *)
(*MC Processes*)


(* ::Subsubsection::Closed:: *)
(*Gravitational Potential*)


(*Vgrav[r_]:=Piecewise[{{-(("G" "ME")/r),r>"rE"/. Constants`EarthRepl},{-(("G" "ME")/(2("rE")^3))( 3("rE")^2-r^2),r<"rE"/. Constants`EarthRepl}},-(("G" "ME")/("rE"/. Constants`EarthRepl))]/.Constants`SIConstRepl/.Constants`EarthRepl*)
Vgrav[r_]:=Piecewise[{{-(("rE" ("vesc")^2/2)/r),r>"rE"/. Constants`EarthRepl},{-(("rE" ("vesc")^2/2)/(2("rE")^3))( 3("rE")^2-r^2),r<"rE"/. Constants`EarthRepl}},-(("rE" ("vesc")^2/2)/("rE"/. Constants`EarthRepl))]/.Constants`SIConstRepl/.Constants`EarthRepl


(* ::Input:: *)
(*"we use \!\(\*FractionBox[\(1\), \(2\)]\) \!\(\*SubscriptBox[\(r\), \(E\)]\) \!\(\*SuperscriptBox[SubscriptBox[\(v\), \(esc\)], \(2\)]\) instead of G \!\(\*SubscriptBox[\(M\), \(E\)]\) in the gravitational potential so it's exactly consistent with \!\(\*SubscriptBox[\(v\), \(esc\)]\) assumed elsewhere."*)
(*"rE" ("vesc")^2/2/.Constants`EarthRepl*)
(*("G" "ME")/.Constants`EarthRepl/.SIConstRepl*)


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
(*Weak Capture Regime Trajectory - w Potential - Runge Kutta*)


(* ::Input::Initialization:: *)
Clear[GetWeakCaptureTrajectoryinvRK2wPotential]
GetWeakCaptureTrajectoryinvRK2wPotential[m\[Chi]_,v\[Chi]E_,bE_,\[Kappa]_,ELfunc_,N_:30,Vofr_:Vgrav]:=Module[{rE,\[CapitalDelta]l,E0,L0,steps,traj,ti,l,r,v,vr,\[CapitalDelta]E,L,Log\[CapitalDelta]L,Lcomped,vsol,ldom,lsol,rdom,voft,rfinal(*,vfinal,vesc*),lMC1,lMC2,lMC},
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

\[CapitalDelta]E = ti["\[CapitalDelta]E"] + \[CapitalDelta]l \[Kappa]^2 ELfunc[r,v];


vr = ti["vr"] + \[CapitalDelta]l ((v^2-ti["vr"]^2)/(r v) - (m\[Chi] Vofr'[r])/(m\[Chi] v)-ti["vr"]/(m\[Chi] v^2) \[Kappa]^2 ELfunc[r,v]);(*assumes ELfunc > 0, so we include - for dissipation*)

L=m\[Chi] r Sqrt[v^2-vr^2];(*this choice of incrementation, preserves energy but wildly violates L conservation / behaviour*)

Log\[CapitalDelta]L = ti["Log\[CapitalDelta]L"]+ (\[CapitalDelta]l \[Kappa]^2 ELfunc[r,v])/(m\[Chi] v^2);


AppendTo[traj,<|"l"->l,"r"->r,"v"->v,"vr"->vr,"\[CapitalDelta]E"->\[CapitalDelta]E,"V"->m\[Chi] Vofr[r],"L"->L,"Log\[CapitalDelta]L"->Log\[CapitalDelta]L(*,"Lcomped"->Lcomped*)|>];

If[m\[Chi]/2 v^2+ m\[Chi] Vofr[r]< 0,Break[]];(*particle is captured, ie. bound to the earth, placed after the append statement so that at least one additional point is added. This is needed for the captured condition used further down the pipeline in the case where initial conditions already imply the particle is bound.*)

If[r>rE,Break[]];
,{s,3N}];(*will break once reaches surface (<~ 2 N) this is just a buffer*)

If[Max[Abs[(("m\[Chi]")/2 ("v")^2+"V"+"\[CapitalDelta]E")/E0-1]/.traj]>0.20,Print["Energy conservation has been violated by more than 20% in Forward Euler: ",Max@Abs[(("m\[Chi]")/2 ("v")^2+"V"+"\[CapitalDelta]E")/E0-1]/.traj]];

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
lMC1 = l/.Quiet@FindRoot[lsol[l]-"rcore"/.Constants`EarthRepl,{l,rdom[[2]]/4,rdom[[1]],rdom[[2]]/2}];
lMC2 = l/.Quiet@FindRoot[lsol[l]-"rcore"/.Constants`EarthRepl,{l,(3rdom[[2]])/4,rdom[[2]]/2,rdom[[2]]}];
lMC={lMC1,lMC2};,
(*doesn't pass through the core (No Core)*)
lMC={"NC"}
],
(*It's captured somewhere, so we won't need to distinguish between crust and core for since we won't compute hard capture.*)
lMC = {"Cap"}
];

<|"v(l)"->vsol,"domain"->ldom,"r(l)"->lsol,"rdom"->rdom,"traj"->traj,"\[Kappa]"->\[Kappa],"bE"->bE,"m\[Chi]"->m\[Chi],"v\[Chi]E"->v\[Chi]E,"E0"->E0,"L0"->L0,"steps"->steps,"rfinal"->rfinal,"lMCs"->lMC,"V(r)"->Function[{r},m\[Chi] Vofr[r]]|>
]


(* ::Subsection:: *)
(*Weak Capture Probabilities*)


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

\[Lambda]invof\[Chi]=If[#>v\[Chi]min && \[Omega]maxofv\[Chi]t[#]>\[Omega]capture[#]&&(\[Xi]of\[Omega]andv\[Chi]t[\[Omega]capture[#],#])>10^\[Sigma]dom[[2,1]],("ne"/.params)\[Sigma]ofv\[Chi]and\[Omega]min[#,\[Omega]capture[#]],0]&; (*must have v > minimum velocity allowed, must have maximum possible energy transfer greater than energy needed to capture, and must have that energy transfer be in the domain of the interpolation function*)

\[Lambda]invof\[Chi]
]


(* ::Input::Initialization:: *)
Clear[GetTotaleWeakCaptureInteractionLength]
GetTotaleWeakCaptureInteractionLength[\[Sigma]Dicts_]:=Module[{},
Total@Table[GetWeakCaptureInteractionLength[\[Sigma]Dicts[[i]]][#],{i,Length[\[Sigma]Dicts]}]&
]


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
(*Weak Capture Regime Capture Probability - w Potential - position dependent Subscript[v, esc] - e*)


(* ::Input::Initialization:: *)
Clear[GetWeakCaptureProbabilitywPot]
GetWeakCaptureProbabilitywPot[vDict_,\[Sigma]Dict_]:= Module[{v\[Chi]min,v\[Chi]raw,v\[Chi],\[Chi]traj,bE,m\[Chi],\[Omega]min,params,ldom,\[Kappa],regionindex,\[Sigma]dom,vesc,\[Omega]capture,\[Delta]\[Xi],\[Xi]of\[Omega]andl,\[Omega]maxofl,d\[Sigma]dERoftand\[Omega],\[Sigma]oftand\[Omega]min,rCore,rE,lcore,lcrust,\[Lambda]invoft,vbound1,vbound2,\[Lambda]integrand,integralof\[Lambda]invwrtr},
\[Omega]min="\[Omega]min"/.\[Sigma]Dict;
params=\[Sigma]Dict["params"];
\[Delta]\[Xi]=(10^"\[Xi]"/.\[Sigma]Dict)[[1]];
m\[Chi]=("m\[Chi]"/.vDict);
ldom = "domain"/.vDict;
\[Kappa]=("\[Kappa]"/.vDict);
regionindex = "regionindex"/.\[Sigma]Dict; (*1 in crust, 2 in core*)
\[Sigma]dom=("\[Sigma]of\[Xi]f"/.\[Sigma]Dict)["Domain"];

vesc=Function[{l},Sqrt[-((2("V(r)"/.vDict)[("r(l)"/.vDict)[l]])/m\[Chi])]];

v\[Chi]min=Max[3Sqrt[("\[HBar]" \[Omega]min)/(2 m\[Chi])]/.params,vesc[#]]&;
v\[Chi]raw=("v(l)"/.vDict);
v\[Chi]=If[v\[Chi]raw[#]>v\[Chi]min[#],("v(l)"/.vDict)[#],0]&;
bE=("bE"/.vDict);(*impact parameter*)

\[Xi]of\[Omega]andl=EnergyLoss`\[Xi]of\[Omega]andv\[Chi][#1,\[Omega]min,("v(l)"/.vDict)[#2],m\[Chi],params]&;
\[Omega]maxofl=EnergyLoss`\[Omega]maxofv\[Chi][("v(l)"/.vDict)[#1],m\[Chi],params]&;

\[Sigma]oftand\[Omega]min=10^(Re[("\[Sigma]of\[Xi]f"/.\[Sigma]Dict)[Log10@v\[Chi][#1],Log10@\[Xi]of\[Omega]andl[#2,#1]]])&;

\[Omega]capture = 1/(2 ("\[HBar]"/.params)) m\[Chi] (v\[Chi][#]^2 - vesc[#]^2)&;

\[Lambda]invoft=If[v\[Chi]raw[#]>v\[Chi]min[#] && \[Omega]maxofl[#]>\[Omega]capture[#]&&(\[Xi]of\[Omega]andl[\[Omega]capture[#],#])>10^\[Sigma]dom[[2,1]],\[Kappa]^2 ("ne"/.params)\[Sigma]oftand\[Omega]min[#,\[Omega]capture[#]],0]&;

(*find the points where capture goes to 0 due to large velocity (a result of acceleration), if the region is too close to the boundary, we need to integrate over them more carefully since mathematica can miss it. *)
vbound1 =Quiet@Check[ l/.FindRoot[\[Omega]maxofl[#]-\[Omega]capture[#]&[l],{l,ldom[[2]]/4,ldom[[1]],ldom[[2]]/2}],ldom[[2]]/2];
vbound2= Quiet@Check[l/.FindRoot[\[Omega]maxofl[#]-\[Omega]capture[#]&[l],{l,(3ldom[[2]])/4,ldom[[2]]/2,ldom[[2]]}],ldom[[2]]/2];

rCore="rcore"/.Constants`EarthRepl;
rE = "rE"/.Constants`EarthRepl;

\[Lambda]integrand[l_?NumericQ]:=\[Lambda]invoft[l];


If[vDict["lMCs"]=={ "NC"},
(*Doesn't pass through the core, just integrate*)
If[regionindex==2,
(*material is in the core, return 0*)
Return[0,Module],
(*material is in the mantle, just integrate*)
integralof\[Lambda]invwrtr=Check[#[12],Print[{m\[Chi],"v\[Chi]E"/.vDict,bE,\[Kappa],"materialname"/.\[Sigma]Dict}];#[15]]&[Function[{MR},If[vbound1/(ldom[[2]]-ldom[[1]])<10^-2||1-vbound2/(ldom[[2]]-ldom[[1]])<10^-2,NIntegrate[\[Lambda]integrand[l],{l,ldom[[1]],vbound1},PrecisionGoal->3,MaxRecursion->MR,AccuracyGoal->3]+NIntegrate[\[Lambda]integrand[l],{l,vbound2,ldom[[2]]},PrecisionGoal->3,MaxRecursion->MR,AccuracyGoal->3],NIntegrate[\[Lambda]integrand[l],{l,ldom[[1]],ldom[[2]]},PrecisionGoal->3,MaxRecursion->MR,AccuracyGoal->3]]]];
],
(*Passes through the core*)
If[regionindex==2,
(*material is in the core, integrate between the two Mantle-Core boundary points*)
integralof\[Lambda]invwrtr=Check[#,Print[{m\[Chi],"v\[Chi]E"/.vDict,bE,\[Kappa],"materialname"/.\[Sigma]Dict}];#]&[NIntegrate[\[Lambda]integrand[l],{l,vDict["lMCs"][[1]],vDict["lMCs"][[2]]},PrecisionGoal->3,MaxRecursion->12,AccuracyGoal->3]];,
(*material is in the mantle, integrate over the regions before and after the core, then add*)
integralof\[Lambda]invwrtr=Check[#,Print[{m\[Chi],"v\[Chi]E"/.vDict,bE,\[Kappa],"materialname"/.\[Sigma]Dict}];#]&[NIntegrate[\[Lambda]integrand[l],{l,ldom[[1]],vDict["lMCs"][[1]]},PrecisionGoal->3,MaxRecursion->12,AccuracyGoal->3]+NIntegrate[\[Lambda]integrand[l],{l,vDict["lMCs"][[2]],ldom[[2]]},PrecisionGoal->3,MaxRecursion->12,AccuracyGoal->3]];
]
];

integralof\[Lambda]invwrtr
]


(* ::Subsubsection::Closed:: *)
(*Weak Capture Regime Capture Probability - w Potential - position dependent Subscript[v, esc] - Nuc*)


(* ::Input::Initialization:: *)
Clear[GetWeakCaptureProbabilityNucwPot]
GetWeakCaptureProbabilityNucwPot[vDict_,\[Sigma]Dict_]:= Module[{v\[Chi]min,v\[Chi]raw,v\[Chi],\[Chi]traj,bE,m\[Chi],\[Omega]min,params,NucleusParams,ldom,\[Kappa],regionindex,\[Sigma]dom,mN,\[Omega]capture,vesc,\[Delta]\[Xi],\[Xi]of\[Omega]andl,\[Omega]maxofl,d\[Sigma]dERoftand\[Omega],\[Sigma]oftand\[Omega]min,\[Lambda]invoft,\[Lambda]integrand,vbound1,vbound2,rCore,rE,lcore,lcrust,integralof\[Lambda]invwrtr},
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

vesc=Function[{l},Sqrt[-((2("V(r)"/.vDict)[("r(l)"/.vDict)[l]])/m\[Chi])]]; 

v\[Chi]min=Max[3Sqrt[("\[HBar]" \[Omega]min)/(2 m\[Chi])]/.params,vesc[#]]&;(*minimum allowable velocity given the lower bound on energy transfer, times an O(1) fudge factor for numerical stability (3) *)
v\[Chi]raw=("v(l)"/.vDict);
v\[Chi]=If[v\[Chi]raw[#]>v\[Chi]min[#],("v(l)"/.vDict)[#],0]&;

bE=("bE"/.vDict);(*impact parameter*)

\[Xi]of\[Omega]andl=EnergyLoss`\[Xi]of\[Omega]Nuc[#1,\[Omega]min,("v(l)"/.vDict)[#2],m\[Chi],mN,params]&;(*[\[Omega],l]*)
\[Omega]maxofl=EnergyLoss`\[Omega]maxNuc[("v(l)"/.vDict)[#1],m\[Chi],mN,params]&;

\[Sigma]oftand\[Omega]min=10^(Re[("\[Sigma]of\[Xi]f"/.\[Sigma]Dict)[Log10@v\[Chi][#1],Log10@\[Xi]of\[Omega]andl[#2,#1]]])&;(*[l,\[Omega]]*)

\[Omega]capture = 1/(2 ("\[HBar]"/.params)) m\[Chi] (v\[Chi][#]^2 - vesc[#]^2)&; 

\[Lambda]invoft=If[v\[Chi]raw[#]>v\[Chi]min[#] && \[Omega]maxofl[#]>\[Omega]capture[#]&&(\[Xi]of\[Omega]andl[\[Omega]capture[#],#])>10^\[Sigma]dom[[2,1]],\[Kappa]^2 ("nI"/.NucleusParams)\[Sigma]oftand\[Omega]min[#,\[Omega]capture[#]],0]&;

vbound1 =Quiet@Check[ l/.FindRoot[\[Omega]maxofl[#]-\[Omega]capture[#]&[l],{l,ldom[[2]]/4,ldom[[1]],ldom[[2]]/2}],ldom[[2]]/2];
vbound2= Quiet@Check[l/.FindRoot[\[Omega]maxofl[#]-\[Omega]capture[#]&[l],{l,(3ldom[[2]])/4,ldom[[2]]/2,ldom[[2]]}],ldom[[2]]/2];

rCore="rcore"/.Constants`EarthRepl;
rE = "rE"/.Constants`EarthRepl;

\[Lambda]integrand[l_?NumericQ]:=\[Lambda]invoft[l];

If[vDict["lMCs"]=={ "NC"},(*could use -1 and -2 for No core and Captured*)
(*Doesn't pass through the core, just integrate*)
If[regionindex==2,
(*material is in the core, return 0*)
Return[0,Module],
(*material is in the mantle, just integrate*)
integralof\[Lambda]invwrtr=If[vbound1/(ldom[[2]]-ldom[[1]])<10^-2||1-vbound2/(ldom[[2]]-ldom[[1]])<10^-2,NIntegrate[\[Lambda]integrand[l],{l,ldom[[1]],vbound1},PrecisionGoal->3,MaxRecursion->12,AccuracyGoal->3]+NIntegrate[\[Lambda]integrand[l],{l,vbound2,ldom[[2]]},PrecisionGoal->3,MaxRecursion->12,AccuracyGoal->3],NIntegrate[\[Lambda]integrand[l],{l,ldom[[1]],ldom[[2]]},PrecisionGoal->3,MaxRecursion->12,AccuracyGoal->3]];
],
(*Passes through the core*)
If[regionindex==2,
(*material is in the core, integrate between the two Mantle-Core boundary points*)
integralof\[Lambda]invwrtr=Check[#,Print[{m\[Chi],"v\[Chi]E"/.vDict,bE,\[Kappa],"materialname"/.\[Sigma]Dict}];#]&[NIntegrate[\[Lambda]integrand[l],{l,vDict["lMCs"][[1]],vDict["lMCs"][[2]]},PrecisionGoal->3,MaxRecursion->12,AccuracyGoal->3]];,
(*material is in the mantle, integrate over the regions before and after the core, then add*)
integralof\[Lambda]invwrtr=Check[#,Print[{m\[Chi],"v\[Chi]E"/.vDict,bE,\[Kappa],"materialname"/.\[Sigma]Dict}];#]&[NIntegrate[\[Lambda]integrand[l],{l,ldom[[1]],vDict["lMCs"][[1]]},PrecisionGoal->3,MaxRecursion->12,AccuracyGoal->3]+NIntegrate[\[Lambda]integrand[l],{l,vDict["lMCs"][[2]],ldom[[2]]},PrecisionGoal->3,MaxRecursion->12,AccuracyGoal->3]];
]
];

integralof\[Lambda]invwrtr
]


(* ::Subsection:: *)
(*Initial Conditions*)


(* ::Subsubsection::Closed:: *)
(*vmax*)


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
Clear[Getv\[Chi]atEarth]
Getv\[Chi]atEarth[v\[Chi]Max_,nbs_:10,nv\[Chi]s_:40,lowerbrat_:(*10^-6*)10^-5,lowervrat_:10^-8]:=Module[{rE,vmin,v\[Chi]s,bs,bEarth,vrEarth,v\[Chi],\[Chi]speedEarth},

    rE = "rE"/.Constants`EarthRepl;
    vmin = "vesc"/.Constants`EarthRepl;
   
v\[Chi]s = vmin+ 10^Subdivide[ Log10[lowervrat(v\[Chi]Max-vmin)],Log10@(v\[Chi]Max-vmin),nv\[Chi]s];

bs = Subdivide[lowerbrat rE,rE,nbs];

Table[<|"\[Chi]speedEarth"->v\[Chi]s[[i]],"bEarth"->bs[[j]]|>,{i,nv\[Chi]s},{j,nbs}]
]


(* ::Subsection:: *)
(*Get Capture Probability*)


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
(*Interpolate Capture Probability distributions*)


(* ::Input::Initialization:: *)
Clear[InterpolatePc]
InterpolatePc[PDictList_]:=Module[{Intertable,Interf,\[Kappa]GatheredDicts,mpD,meD,mpDcap,m\[Chi],v\[Chi]Max,\[Beta]D,v0},

{mpD,meD,m\[Chi],mpDcap,v\[Chi]Max,\[Beta]D,v0}=Union[{"mpD","meD","m\[Chi]","mpDcap","v\[Chi]Max","\[Beta]D","v0"}/.PDictList][[1]];
\[Kappa]GatheredDicts=Gather[PDictList,(#1["\[Kappa]"]==#2["\[Kappa]"])&];
Interf=Table[<|"Pf"->Interpolation[Log10@{{"v\[Chi]E","bEarth"},"Pcap"}/.\[Kappa]GatheredDicts[[i]],InterpolationOrder->1],"\[Kappa]"->(Union["\[Kappa]"/.\[Kappa]GatheredDicts[[i]]][[1]]),"mpD"->mpD,"meD"->meD,"m\[Chi]"->m\[Chi],"mpDcap"->mpDcap,"v\[Chi]Max"->v\[Chi]Max,"\[Beta]D"->\[Beta]D,"v0"->v0,"PfE"->Interpolation[Log10@{{"v\[Chi]E","bEarth"},"PcapE"}/.\[Kappa]GatheredDicts[[i]],InterpolationOrder->1],"PfNuc"->Interpolation[Log10@{{"v\[Chi]E","bEarth"},"PcapNuc"}/.\[Kappa]GatheredDicts[[i]],InterpolationOrder->1],"y"->Interpolation[Log10@{{"v\[Chi]E","bEarth"},"y"}/.\[Kappa]GatheredDicts[[i]],InterpolationOrder->1],"CSDcap"->Interpolation[{Log10@{"v\[Chi]E","bEarth"},1-2"CSDcap"}/.\[Kappa]GatheredDicts[[i]],InterpolationOrder->1],"int\[Lambda]invE"->Interpolation[Log10@{{"v\[Chi]E","bEarth"},Max["int\[Lambda]invE",10^-100]}/.\[Kappa]GatheredDicts[[i]],InterpolationOrder->1],"int\[Lambda]invNuc"->Interpolation[Log10@{{"v\[Chi]E","bEarth"},Max["int\[Lambda]invNuc",10^-100]}/.\[Kappa]GatheredDicts[[i]],InterpolationOrder->1]|>,{i,Length[\[Kappa]GatheredDicts]}];
Interf
]


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

\[Rho]DM0= "\[CapitalOmega]c"/"\[CapitalOmega]b" n0 mSMtot/.\[CapitalLambda]CDMrepl;
\[Rho]DM0/mDMtot fD
]


(* ::Input::Initialization:: *)
naDM[fD_,mDMtot_] := Module[{\[Rho]DM},
(*fD - [ ] ionized mass fraction of local DM density in the galaxy
mDMtot - [kg] sum of masses of DM species *)

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
m\[Chi]="m\[Chi]"/.InterPDicts[[i]]; 
\[Beta]D="\[Beta]D"/.InterPDicts[[i]];
Pfdom=Union@Table[Pfs[[i,k]]["Domain"],{k,Length[Pfs[[i]]]}];

If[(Dimensions@Pfdom)[[1]]!=1,Print["The interpolation functions corresponding to these keys have inequivalent domains, so will not be integrated over."];Return[0,Module]];

Pfdom = Pfdom[[1]];
fintegrand[bE_?NumericQ,v\[Chi]_?NumericQ]:=bE^2 E^(-(1/2) m\[Chi] (v\[Chi]^2-vesc^2) \[Beta]D)  v\[Chi]^3/Sqrt[1 + bE^2/rE^2] Sum[10^(Pfs[[i,k]][Log10@v\[Chi],Log10@bE]),{k,Length[Pfs[[i]]]}];(*[(s^-1)] goes as <Subscript[v, \[Chi]]> Subscript[n, aDM] Subscript[A, Earth] so this rate over Subscript[V, Earth] is <Subscript[v, \[Chi]]> Subscript[n, aDM] Subscript[R, Earth]^-1 which is what was used in the appendix*)

intregboundaries=If[mpDcap,Table[(10^Pfdom[[1,2]]-10^Pfdom[[1,1]])10^c,{c,-6,-1}],Table[(10^Pfdom[[1,2]]-10^Pfdom[[1,1]])10^c,{c,-8,-1}]]; 
dNcMaxdt=AEarth (3 "nD" Sqrt[2/\[Pi]]  (m\[Chi] \[Beta]D)^(3/2))/rE^3 Quiet@Check[#,Print[{m\[Chi],\[Beta]D,"\[Kappa]"/.InterPDicts[[i]]}];#]&[(NIntegrate[fintegrand[bE,v\[Chi]],{v\[Chi],10^Pfdom[[1,1]],10^Pfdom[[1,1]]+intregboundaries[[1]]},
{bE,10^Pfdom[[2,1]],10^Pfdom[[2,2]]},PrecisionGoal->3,MaxRecursion->8,AccuracyGoal->3]+Sum[With[{ct=c},NIntegrate[fintegrand[bE,v\[Chi]],{v\[Chi],10^Pfdom[[1,1]]+intregboundaries[[ct-1]],10^Pfdom[[1,1]]+intregboundaries[[ct]]},
{bE,10^Pfdom[[2,1]],10^Pfdom[[2,2]]},PrecisionGoal->3,MaxRecursion->8,AccuracyGoal->3]],{c,2,Length[intregboundaries]}]+NIntegrate[fintegrand[bE,v\[Chi]],{v\[Chi],10^Pfdom[[1,1]]+intregboundaries[[-1]],10^Pfdom[[1,2]]},
{bE,10^Pfdom[[2,1]],10^Pfdom[[2,2]]},PrecisionGoal->3,MaxRecursion->8,AccuracyGoal->3])];   (*we break up the integration region because Nuclear capture probability is EXTREMELY sharply peaked near v\[Chi]=vesc. Mathematica will ignore this contribution for most of the parameter space if we don't divide up the integration region like this. *)

(* May want to define a function of the integration bounds so this expression isn't so disgusting. *)

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


(* ::Subsubsection::Closed:: *)
(*Evaporation Rate - w Evap Volume *)


(* ::Input:: *)
(*StringQ[Function[{x},x]]*)


(* ::Input::Initialization:: *)
Clear[GetEvaporationRateVevap]
GetEvaporationRateVevap[\[Sigma]Dict_,\[Kappa]_,Vevapfunc_:"NA",mpDcap_:True]:= Module[{v\[Chi]min,v\[Chi]max,m\[Chi],\[Omega]min,params,regionindex,NucleusParams,mN,nT,\[Beta]E,\[Omega]capture,\[Delta]\[Xi],\[Xi]min,\[Xi]of\[Omega]andv\[Chi],\[Omega]maxofv\[Chi],d\[Sigma]dERofv\[Chi]and\[Omega],\[Sigma]ofv\[Chi]and\[Omega]min,\[CapitalOmega]evapperparticle,\[CapitalGamma]evappercapturedparticle,v\[Chi]st\[Omega]capis\[Omega]max,RegionVolume,\[CapitalOmega]evapperparticlewvT,EarthVolume,EvapRescale,fintegrand,intregboundaries,Nv\[Chi]=100},
\[Omega]min="\[Omega]min"/.\[Sigma]Dict;
params=\[Sigma]Dict["params"];
NucleusParams=\[Sigma]Dict["NucleusParams"];
\[Delta]\[Xi]=(10^"\[Xi]"/.\[Sigma]Dict)[[1]];
m\[Chi]=("m\[Chi]"/.\[Sigma]Dict);
mN = "mN"/.\[Sigma]Dict;
regionindex="regionindex"/.\[Sigma]Dict;
\[Beta]E = {"\[Beta]crust","\[Beta]core"}[[regionindex]]/.Constants`EarthRepl;(*Earth Temeperature*)
nT = "nI"/.NucleusParams;
EarthVolume =  (4\[Pi])/3 ("rE"/.Constants`EarthRepl)^3;

\[Xi]of\[Omega]andv\[Chi]=EnergyLoss`\[Xi]of\[Omega]Nuc[#1,\[Omega]min,#2,m\[Chi],mN,params]&;

\[Sigma]ofv\[Chi]and\[Omega]min=10^(Re[("\[Sigma]of\[Xi]f"/.\[Sigma]Dict)[Log10@#1,Log10@\[Xi]of\[Omega]andv\[Chi][#2,#1]]])&;

v\[Chi]min = 10^(("\[Sigma]of\[Xi]f"/.\[Sigma]Dict)["Domain"][[1,1]]);
\[Xi]min = 10^(("\[Sigma]of\[Xi]f"/.\[Sigma]Dict)["Domain"][[2,1]]);

(*Post \[Omega]capture*)
\[Omega]maxofv\[Chi]=EnergyLoss`\[Omega]maxNuc[#1,m\[Chi],mN,params]&; (*maximum energy transferable in a collision*)

\[Omega]capture = 1/(2 ("\[HBar]"/.params)) m\[Chi] (("vesc"/.Constants`EarthRepl)^2-#^2 )&;(*energy transfered to aDM such that it has v > Subscript[v, esc]*)

v\[Chi]st\[Omega]capis\[Omega]max=v\[Chi]/.NSolve[\[Omega]maxofv\[Chi][v\[Chi]]==\[Omega]capture[v\[Chi]],v\[Chi]][[2]];(*find the lowest velocity such that aDM can be evaporated through a single scatter*)

EvapRescale=1/(-If[x^2<500,E^(-(1/2) x^2) Sqrt[2/\[Pi]] x,0]+Erf[x/Sqrt[2]])/.x->Sqrt[m\[Chi]]"vesc" Sqrt[\[Beta]E]/.Constants`EarthRepl; (*account for the fact that the MB distribution is truncated, renormalize*)

fintegrand[vT_?NumericQ,v\[Chi]_?NumericQ]:=If[StringQ[Vevapfunc],vT^2 E^(-\[Beta]E 1/2 mN vT^2) ((vT+v\[Chi])^3-((vT-v\[Chi])^2)^(3/2))/(3 vT v\[Chi]) v\[Chi]^2 E^(-\[Beta]E 1/2 m\[Chi] v\[Chi]^2) If[\[Omega]maxofv\[Chi][v\[Chi]]>\[Omega]capture[v\[Chi]]&&(\[Xi]of\[Omega]andv\[Chi][\[Omega]capture[v\[Chi]],v\[Chi]])>\[Xi]min,1,0]\[Sigma]ofv\[Chi]and\[Omega]min[v\[Chi],\[Omega]capture[v\[Chi]]],Vevapfunc[v\[Chi],\[Kappa],regionindex]/EarthVolume vT^2 E^(-\[Beta]E 1/2 mN vT^2) ((vT+v\[Chi])^3-((vT-v\[Chi])^2)^(3/2))/(3 vT v\[Chi]) v\[Chi]^2 E^(-\[Beta]E 1/2 m\[Chi] v\[Chi]^2) If[\[Omega]maxofv\[Chi][v\[Chi]]>\[Omega]capture[v\[Chi]]&&(\[Xi]of\[Omega]andv\[Chi][\[Omega]capture[v\[Chi]],v\[Chi]])>\[Xi]min,1,0]\[Sigma]ofv\[Chi]and\[Omega]min[v\[Chi],\[Omega]capture[v\[Chi]]]];

intregboundaries=If[mpDcap,Table[(("vesc"/.Constants`EarthRepl)-v\[Chi]st\[Omega]capis\[Omega]max)10^c,{c,-5,-1}],Table[(("vesc"/.Constants`EarthRepl)-v\[Chi]st\[Omega]capis\[Omega]max)10^c,{c,-7,-1}]];(* DIFF !!! *)

\[CapitalOmega]evapperparticlewvT=EvapRescale  nT  ((Sqrt[m\[Chi] mN ]\[Beta]E)/(2 \[Pi]))^3 4 \[Pi] 2 \[Pi] Quiet@(NIntegrate[fintegrand[vT,v\[Chi]],{v\[Chi],v\[Chi]st\[Omega]capis\[Omega]max,v\[Chi]st\[Omega]capis\[Omega]max + intregboundaries[[1]]},{vT,0,\[Infinity]},PrecisionGoal->3,MaxRecursion->5,AccuracyGoal->3]+Sum[NIntegrate[fintegrand[vT,v\[Chi]],{v\[Chi],v\[Chi]st\[Omega]capis\[Omega]max+ intregboundaries[[c-1]],v\[Chi]st\[Omega]capis\[Omega]max + intregboundaries[[c]]},{vT,0,\[Infinity]},PrecisionGoal->3,MaxRecursion->5,AccuracyGoal->3],{c,2,Length[intregboundaries]}]+NIntegrate[fintegrand[vT,v\[Chi]],{v\[Chi],v\[Chi]st\[Omega]capis\[Omega]max+ intregboundaries[[-1]],"vesc"/.Constants`EarthRepl},{vT,0,\[Infinity]},PrecisionGoal->3,MaxRecursion->5,AccuracyGoal->3]);

\[CapitalGamma]evappercapturedparticle=If[StringQ[Vevapfunc],RegionVolume =  4\[Pi] Integrate[r^2 Region\[CapitalTheta][r,regionindex],{r,0,"rE"/.Constants`EarthRepl}]; \[CapitalOmega]evapperparticlewvT RegionVolume/EarthVolume,\[CapitalOmega]evapperparticlewvT];  (*\[CapitalGamma]\[Kappa]^2 corresponds to Revp in appendix*)(*Assumes homogeneously distributed aDM in the Earth*)

<|"\[CapitalGamma]"->\[CapitalGamma]evappercapturedparticle,"regionindex"->regionindex|>
]


(* ::Input::Initialization:: *)
Clear[IterateGetEvaporationRateVevap]
IterateGetEvaporationRateVevap[Nuclear\[Sigma]Dicts_,\[Kappa]_,Vevapfunc_:"NA"]:=Module[{\[CapitalGamma]list={},\[CapitalGamma]Total},
Monitor[
Do[
AppendTo[\[CapitalGamma]list,GetEvaporationRateVevap[Nuclear\[Sigma]Dicts[[n]],\[Kappa],Vevapfunc]],
{n,Length[Nuclear\[Sigma]Dicts]}
];
,n];
\[CapitalGamma]Total=Total[\[CapitalGamma]list[[;;,1]]];

<|"\[CapitalGamma]Total"->\[CapitalGamma]Total,"\[CapitalGamma]s"->\[CapitalGamma]list|>
]


(* ::Subsubsection::Closed:: *)
(*Evaporation Rate Electronic - w Evap Volume*)


(* ::Input::Initialization:: *)
Clear[GetEvaporationRateendVevap]
GetEvaporationRateendVevap[\[Sigma]Dict_,\[Kappa]_,Vevapfunc_:"NA",mpDcap_:True]:= Module[{v\[Chi]min,v\[Chi]max,m\[Chi],\[Omega]min,params,regionindex,mT,nT,vF,vesc,vTm,vTp,\[Beta]E,\[Omega]capture,\[Delta]\[Xi],\[Xi]min,\[Xi]of\[Omega]andv\[Chi]t,\[Omega]maxofv\[Chi]t,d\[Sigma]dERofv\[Chi]and\[Omega],\[Sigma]ofv\[Chi]and\[Omega]min,\[CapitalOmega]evapperparticle,\[CapitalGamma]evappercapturedparticle,v\[Chi]st\[Omega]capis\[Omega]max,RegionVolume,\[CapitalOmega]evapperparticlewvT,EarthVolume,EvapRescale,fintegrand,intregboundaries,PLOT=False},

params=\[Sigma]Dict["params"];
\[Delta]\[Xi]=(10^"\[Xi]"/.\[Sigma]Dict)[[1]];
m\[Chi]=("m\[Chi]"/.\[Sigma]Dict);
mT = "m"/.\[Sigma]Dict;
regionindex="regionindex"/.\[Sigma]Dict;
\[Beta]E = {"\[Beta]crust","\[Beta]core"}[[regionindex]]/.Constants`EarthRepl;(*Earth Temeperature*)
nT = "ne"/.params;
vF="vF"/.params;
vesc = ("vesc"/.Constants`EarthRepl);
EarthVolume =  (4\[Pi])/3 ("rE"/.Constants`EarthRepl)^3;

v\[Chi]min = Max[10^(("\[Sigma]of\[Xi]f"/.\[Sigma]Dict)["Domain"][[1,1]]),Re[Sqrt[vesc^2 - 8/(m\[Chi] \[Beta]E)]]];
\[Xi]min = 10^(("\[Sigma]of\[Xi]f"/.\[Sigma]Dict)["Domain"][[2,1]]);


\[Omega]min="\[Omega]min"/.\[Sigma]Dict;
\[Xi]of\[Omega]andv\[Chi]t=EnergyLoss`\[Xi]of\[Omega]andv\[Chi]nd[#1,\[Omega]min,#2,m\[Chi],params]&;

\[Sigma]ofv\[Chi]and\[Omega]min=10^(("\[Sigma]of\[Xi]f"/.\[Sigma]Dict)[Log10@#1,Log10@\[Xi]of\[Omega]andv\[Chi]t[#2,#1]])&;

\[Omega]maxofv\[Chi]t=EnergyLoss`\[Omega]maxofv\[Chi]nd[#1,m\[Chi],params]&; (*maximum energy transferable in a collision*)

\[Omega]capture = 1/(2 ("\[HBar]"/.params)) m\[Chi] (vesc^2-#^2 )&;(*energy transfered to aDM such that it has v > Subscript[v, esc]*)

v\[Chi]st\[Omega]capis\[Omega]max=v\[Chi]min;

EvapRescale=1/(-If[x^2<500,E^(-(1/2) x^2) Sqrt[2/\[Pi]] x,0]+Erf[x/Sqrt[2]])/.x->Sqrt[m\[Chi]]"vesc" Sqrt[\[Beta]E]/.Constants`EarthRepl; (*account for the fact that the aDM MB distribution is truncated, normalize to 1*)


fintegrand[vT_?NumericQ,v\[Chi]_?NumericQ]:=If[StringQ[Vevapfunc],vT^2 (fnd[vT/vF]/.params)(((vT+v\[Chi])^3-((vT-v\[Chi])^2)^(3/2))/(3 vT v\[Chi]))v\[Chi]^2 E^(-\[Beta]E 1/2 m\[Chi] v\[Chi]^2) If[\[Omega]maxofv\[Chi]t[v\[Chi]]>\[Omega]capture[v\[Chi]]&&(\[Xi]of\[Omega]andv\[Chi]t[\[Omega]capture[v\[Chi]],v\[Chi]])>\[Xi]min,1,0]\[Sigma]ofv\[Chi]and\[Omega]min[v\[Chi],\[Omega]capture[v\[Chi]]],Vevapfunc[v\[Chi],\[Kappa],regionindex]/EarthVolume vT^2 (fnd[vT/vF]/.params)(((vT+v\[Chi])^3-((vT-v\[Chi])^2)^(3/2))/(3 vT v\[Chi]))v\[Chi]^2 E^(-\[Beta]E 1/2 m\[Chi] v\[Chi]^2) If[\[Omega]maxofv\[Chi]t[v\[Chi]]>\[Omega]capture[v\[Chi]]&&(\[Xi]of\[Omega]andv\[Chi]t[\[Omega]capture[v\[Chi]],v\[Chi]])>\[Xi]min,1,0]\[Sigma]ofv\[Chi]and\[Omega]min[v\[Chi],\[Omega]capture[v\[Chi]]]]; (* DIFF !!! *)

intregboundaries=If[mpDcap,Table[(("vesc"/.Constants`EarthRepl)-v\[Chi]st\[Omega]capis\[Omega]max)10^c,{c,-5,-1}],Table[(("vesc"/.Constants`EarthRepl)-v\[Chi]st\[Omega]capis\[Omega]max)10^c,{c,-7,-1}]];

{vTm,vTp}={vF Dielectrics`\[Zeta]m,vF Dielectrics`\[Zeta]p}/.params;

\[CapitalOmega]evapperparticlewvT=EvapRescale  nT (("D"/.params)/Sqrt[2 \[Pi]]) (Sqrt[m\[Chi] \[Beta]E]/vF)^3 (*Quiet@*)(NIntegrate[fintegrand[vT,v\[Chi]],{v\[Chi],v\[Chi]st\[Omega]capis\[Omega]max,v\[Chi]st\[Omega]capis\[Omega]max + intregboundaries[[1]]},{vT,vTm,vTp},PrecisionGoal->3,MaxRecursion->5,AccuracyGoal->3]+Sum[With[{ct=c},NIntegrate[fintegrand[vT,v\[Chi]],{v\[Chi],v\[Chi]st\[Omega]capis\[Omega]max+ intregboundaries[[ct-1]],v\[Chi]st\[Omega]capis\[Omega]max + intregboundaries[[ct]]},{vT,vTm,vTp},PrecisionGoal->3,MaxRecursion->5,AccuracyGoal->3]],{c,2,Length[intregboundaries]}]+NIntegrate[fintegrand[vT,v\[Chi]],{v\[Chi],v\[Chi]st\[Omega]capis\[Omega]max+ intregboundaries[[-1]],"vesc"/.Constants`EarthRepl},{vT,vTm,vTp},PrecisionGoal->3,MaxRecursion->5,AccuracyGoal->3]);

\[CapitalGamma]evappercapturedparticle=If[StringQ[Vevapfunc],RegionVolume =  4\[Pi] Integrate[r^2 Region\[CapitalTheta][r,regionindex],{r,0,"rE"/.Constants`EarthRepl}]; \[CapitalOmega]evapperparticlewvT RegionVolume/EarthVolume,\[CapitalOmega]evapperparticlewvT]; 

<|"\[CapitalGamma]"->\[CapitalGamma]evappercapturedparticle,"regionindex"->regionindex|>
]


(* ::Input::Initialization:: *)
Clear[IterateGetEvaporationRateendVevap]
IterateGetEvaporationRateendVevap[Electronic\[Sigma]Dicts_,\[Kappa]_,Vevapfunc_:"NA"]:=Module[{\[CapitalGamma]list={},\[CapitalGamma]Total},
(*Monitor[*)
Do[
If[Electronic\[Sigma]Dicts[[n]]!=<|"NA"->Null|>,
AppendTo[\[CapitalGamma]list,GetEvaporationRateendVevap[Electronic\[Sigma]Dicts[[n]],\[Kappa],Vevapfunc]]],
{n,Length[Electronic\[Sigma]Dicts]}
];
(*,n];*)
\[CapitalGamma]Total=Total[\[CapitalGamma]list[[;;,1]]];
<|"\[CapitalGamma]Total"->\[CapitalGamma]Total,"\[CapitalGamma]s"->\[CapitalGamma]list|>
]


(* ::Subsubsection::Closed:: *)
(*Get Penetration Depth*)


(* ::Input::Initialization:: *)
Clear[GetPenetrationDepthfnofv]
GetPenetrationDepthfnofv[v\[Chi]_,Electronic\[Sigma]Dicts_,Nuclear\[Sigma]Dicts_,vesc_:"vesc"/.Constants`EarthRepl]:=Module[{m\[Chi],\[Beta]C,\[Beta]M,\[Lambda]inve,\[Lambda]invNuc,eregiontable,eregiondict,vdom,nT,Nucregiontable,Nucregiondict},

	m\[Chi]=Union@Join["m\[Chi]"/.Electronic\[Sigma]Dicts,"m\[Chi]"/.Nuclear\[Sigma]Dicts];
	If[Length[m\[Chi]]>1,Print["Masses in all \[Sigma]Dicts must be equal."];Return[<|"N/A"->Null|>,Module](*,m\[Chi]=m\[Chi][[1]]*)];
	
	{\[Beta]C,\[Beta]M}={"\[Beta]core","\[Beta]crust"}/.Constants`EarthRepl;
	
	eregiontable = Gather[Electronic\[Sigma]Dicts,#1["regionindex"]==#2["regionindex"]&];
	eregiondict=Association@Table[Union["regionindex"/.eregiontable[[i]]][[1]]->eregiontable[[i]],{i,Length[eregiontable]}];
	
	\[Lambda]inve =Association@Table[r->Sum[vdom=("\[Sigma]f"/.eregiondict[#2][[o]])["Domain"][[1]];
	nT=("ne"/.eregiondict[#2][[o]]["params"]);If[#1>vdom[[1]]&&#1<vdom[[2]],nT 10^("\[Sigma]f"/.eregiondict[#2][[o]])[#1],0],{o,Length[eregiondict[#2]]}]&[Log10@v\[Chi],r],{r,2}];(*\[Lambda]inve[v\[Chi]][r] - r is the regionindex*)
	
	Nucregiontable = Gather[Nuclear\[Sigma]Dicts,#1["regionindex"]==#2["regionindex"]&];
	Nucregiondict=Association@Table[Union["regionindex"/.Nucregiontable[[i]]][[1]]->Nucregiontable[[i]],{i,Length[eregiontable]}];
	
	\[Lambda]invNuc=Association@Table[r->Sum[vdom=("\[Sigma]f"/.Nucregiondict[#2][[o]])["Domain"][[1]];
	nT=("nI"/.Nucregiondict[#2][[o]]["NucleusParams"]);If[#1>vdom[[1]]&&#1<vdom[[2]],nT 10^("\[Sigma]f"/.Nucregiondict[#2][[o]])[#1],0],{o,Length[Nucregiondict[#2]]}]&[Log10@v\[Chi],r],{r,2}];(*\[Lambda]invNuc[v\[Chi]][r] - r is the regionindex*)

<|"m\[Chi]"->m\[Chi],"\[Lambda]C"->Quiet@Check[(\[Lambda]inve[2]+\[Lambda]invNuc[2])^-1,10^4 "rE"/.EarthRepl],"\[Lambda]M"->Quiet@Check[(\[Lambda]inve[1]+\[Lambda]invNuc[1])^-1,10^4 "rE"/.EarthRepl],"v\[Chi]"->v\[Chi],"\[Lambda]inve"->\[Lambda]inve,"\[Lambda]invNuc"->\[Lambda]invNuc|>
]


(* ::Input::Initialization:: *)
Clear[GetPenetrationDepth]
GetPenetrationDepth[Electronic\[Sigma]Dicts_,Nuclear\[Sigma]Dicts_,vesc_:"vesc"/.Constants`EarthRepl]:=Module[{m\[Chi],\[Beta]C,\[Beta]M,vbar,\[Lambda]inve,\[Lambda]invNuc,eregiontable,eregiondict,vdom,nT,Nucregiontable,Nucregiondict},

(* *\[Sigma]Dicts should be a single list format (ie. for 1 mass) *)

m\[Chi]=Union@Join["m\[Chi]"/.Electronic\[Sigma]Dicts,"m\[Chi]"/.Nuclear\[Sigma]Dicts];
If[Length[m\[Chi]]>1,Print["Masses in all \[Sigma]Dicts must be equal."];Return[<|"N/A"->Null|>,Module](*,m\[Chi]=m\[Chi][[1]]*)];

{\[Beta]C,\[Beta]M}={"\[Beta]core","\[Beta]crust"}/.Constants`EarthRepl;

vbar = <|2->Log10@Max[Sqrt[2/(m\[Chi] \[Beta]C)],vesc],1->Log10@Max[Sqrt[2/(m\[Chi] \[Beta]M)],vesc]|>;

eregiontable = Gather[Electronic\[Sigma]Dicts,#1["regionindex"]==#2["regionindex"]&];
eregiondict=Association@Table[Union["regionindex"/.eregiontable[[i]]][[1]]->eregiontable[[i]],{i,Length[eregiontable]}];

\[Lambda]inve =Association@Table[r->Sum[vdom=("\[Sigma]f"/.eregiondict[#2][[o]])["Domain"][[1]];
nT=("ne"/.eregiondict[#2][[o]]["params"]);If[#1>vdom[[1]]&&#1<vdom[[2]],nT 10^("\[Sigma]f"/.eregiondict[#2][[o]])[#1],0],{o,Length[eregiondict[#2]]}]&[vbar[r],r],{r,2}];(*\[Lambda]inve[v\[Chi]][r] - r is the regionindex*)

Nucregiontable = Gather[Nuclear\[Sigma]Dicts,#1["regionindex"]==#2["regionindex"]&];
Nucregiondict=Association@Table[Union["regionindex"/.Nucregiontable[[i]]][[1]]->Nucregiontable[[i]],{i,Length[eregiontable]}];

\[Lambda]invNuc=Association@Table[r->Sum[vdom=("\[Sigma]f"/.Nucregiondict[#2][[o]])["Domain"][[1]];
nT=("nI"/.Nucregiondict[#2][[o]]["NucleusParams"]);If[#1>vdom[[1]]&&#1<vdom[[2]],nT 10^("\[Sigma]f"/.Nucregiondict[#2][[o]])[#1],0],{o,Length[Nucregiondict[#2]]}]&[vbar[r],r],{r,2}];(*\[Lambda]invNuc[v\[Chi]][r] - r is the regionindex*)

<|"m\[Chi]"->m\[Chi],"\[Lambda]C"->Quiet@Check[(\[Lambda]inve[2]+\[Lambda]invNuc[2])^-1,10^4 "rE"/.EarthRepl],"\[Lambda]M"->Quiet@Check[(\[Lambda]inve[1]+\[Lambda]invNuc[1])^-1,10^4 "rE"/.EarthRepl],"vbar"->vbar,"\[Lambda]inve"->\[Lambda]inve,"\[Lambda]invNuc"->\[Lambda]invNuc|>
]


(* ::Input::Initialization:: *)
Clear[EvaporationwMFP]
EvaporationwMFP[\[CapitalGamma]total_,\[Lambda]MFPs_,\[Kappa]_,regionindex_]:=Module[{rC,rE,\[Lambda]eff,EvaporationVolume,RegionVolume},

rC = "rcore"/.Constants`EarthRepl;
rE="rE"/.Constants`EarthRepl;

\[Lambda]eff = If[\[Lambda]MFPs["\[Lambda]M"]/\[Kappa]^2<rE-rC,\[Lambda]MFPs["\[Lambda]M"]/\[Kappa]^2,rE-rC + \[Lambda]MFPs["\[Lambda]C"]/\[Kappa]^2];

EvaporationVolume =  4\[Pi] Integrate[r^2 Region\[CapitalTheta][r,regionindex],{r,Max[rE-\[Lambda]eff,0],rE}]; 

RegionVolume =  4\[Pi] Integrate[r^2 Region\[CapitalTheta][r,regionindex],{r,0,"rE"/.Constants`EarthRepl}]; 

\[CapitalGamma]total EvaporationVolume/RegionVolume

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
if \[Lambda]eff > rE-rC and region index = 2, total - Mantle volume,
if \[Lambda]eff > rE, region volume. *)

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

v\[Chi]min=\[Sigma]f["Domain"][[1,1]];

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

regionindex="regionindex"/.\[Sigma]Dict;

RegionVolume =  4\[Pi] Integrate[r^2 Region\[CapitalTheta][r,regionindex],{r,0,"rE"/.Constants`EarthRepl}];
EarthVolume =  (4\[Pi])/3 ("rE"/.Constants`EarthRepl)^3;

GetinvMFT[\[Sigma]Dict]GetPejection[\[Sigma]Dict] RegionVolume/EarthVolume (*function of v\[Chi]*)
]


(* ::Input::Initialization:: *)
GetTotalAppendixEvaporationRate[\[Sigma]Dicts_]:=Module[{},
Sum[GetAppendixEvaporationRate[\[Sigma]Dicts[[i]]],{i,Length[\[Sigma]Dicts]}]
]


(* ::Subsubsection::Closed:: *)
(*Get Subscript[N, c] in the presence of evaporation w Integrated Evap Volume*)


(* ::Input::Initialization:: *)
Clear[GetNc]
GetNc[NcMaxDictList_,Nuc\[Sigma]Dict_,Electronic\[Sigma]Dicts_,ElectronicEvapDict_,keys_:{"Pf"},UseMFP_:True]:=Module[{Vevapfunc,\[CapitalGamma]DictNuc,\[CapitalGamma]Dicte,\[CapitalGamma]total,\[CapitalGamma]totalNuc,\[CapitalGamma]totale,indexforcapturethreshold,\[Kappa],(*UseMFP=True,*)m\[Chi],eind,Nucind,\[Lambda]MFPs,\[CapitalGamma]totalMFP,regionindex,NcMax,dNcMaxdt,namesepchar,keysuffix,tE,timeenoughforequil,Nc,timeenoughforequilMFP, NcMFP,tempNcDict={},NcDicts={}},

Vevapfunc = Function[{v\[Chi]t,\[Kappa],rI},EvaporationVolume[GetPenetrationDepthfnofv[v\[Chi]t,Electronic\[Sigma]Dicts,Nuc\[Sigma]Dict],\[Kappa],rI]["Vevap"]];

keysuffix =If[keys=={"Pf"},"",StringJoin[Table["_"<>keys[[k]],{k,Length[keys]}]]];

Do[
\[Kappa]="\[Kappa]"/.NcMaxDictList[[i]];

If[UseMFP,
\[CapitalGamma]DictNuc=IterateGetEvaporationRateVevap[Nuc\[Sigma]Dict,\[Kappa],Vevapfunc];
\[CapitalGamma]Dicte=IterateGetEvaporationRateendVevap[ElectronicEvapDict,\[Kappa],Vevapfunc]; ,\[CapitalGamma]DictNuc=IterateGetEvaporationRateVevap[Nuc\[Sigma]Dict,\[Kappa]];
\[CapitalGamma]Dicte=IterateGetEvaporationRateendVevap[ElectronicEvapDict,\[Kappa]];
];

\[CapitalGamma]totalNuc="\[CapitalGamma]Total"/.\[CapitalGamma]DictNuc;
\[CapitalGamma]totale="\[CapitalGamma]Total"/.\[CapitalGamma]Dicte;(*evaporation rate per captured particle assuming homogeneous distribution in the Earth*)

NcMax="NcMax"<>keysuffix/.NcMaxDictList[[i]];
dNcMaxdt="dNcMaxdt"<>keysuffix/.NcMaxDictList[[i]];

tE=NcMax/dNcMaxdt;(*[s] age of the Earth*)

\[CapitalGamma]total = \[CapitalGamma]totalNuc+\[CapitalGamma]totale;

timeenoughforequil=1/tE<(\[Kappa]^2 \[CapitalGamma]total);
Nc=If[timeenoughforequil,dNcMaxdt/(\[Kappa]^2 \[CapitalGamma]total),NcMax ];

tempNcDict = Append[NcMaxDictList[[i]] ,{"timeenoughforequil"<>keysuffix->timeenoughforequil,"Nc"<>keysuffix->Nc,"\[CapitalGamma]total"->\[CapitalGamma]total,"\[CapitalGamma]totalNuc"->\[CapitalGamma]totalNuc,"\[CapitalGamma]totale"->\[CapitalGamma]totale}];
AppendTo[NcDicts,tempNcDict];

,{i,Length[NcMaxDictList]}];

NcDicts
]


(* ::Subsubsection::Closed:: *)
(*Scan parameter space*)


(* ::Input::Initialization:: *)
Clear[ScanGetCapture]
ScanGetCapture[meDratios_,\[Kappa]s_,v0s_,Electronic\[Sigma]Dicts_,Nuclear\[Sigma]Dicts_,mpDcap_:True,ElectronicEvapDict_:<|"NA"->Null|>]:=Module[{PDictList,InterPDicts,PDictLists={},ScanDictList={},meDratio,v0,mpD,meD,m\[Chi],\[Beta]D,NcMaxDict,NcDict},
m\[Chi]= "m\[Chi]"/.Electronic\[Sigma]Dicts[[1]];

Do[
meDratio=("meD/mpD"/.meDratios[[me]]);
If[mpDcap,
mpD = m\[Chi];meD = meDratio m\[Chi];,
meD = m\[Chi];mpD =  m\[Chi]/meDratio;
];
v0 = v0s[[vs]];
\[Beta]D ="\[Beta]D"/.v0to\[Beta]D[v0,meD,mpD];

PDictList=IterateGetCaptureProbability[meDratio,\[Kappa]s,v0,Electronic\[Sigma]Dicts,Nuclear\[Sigma]Dicts,mpDcap];(*uses Forward Euler*)

PlotPDist[PDictList];

InterPDicts=InterpolatePc[PDictList];

NcMaxDict=GetNcMax[InterPDicts];

NcDict=GetNc[NcMaxDict,Nuclear\[Sigma]Dicts,Electronic\[Sigma]Dicts,ElectronicEvapDict];

AppendTo[ScanDictList, NcDict];

,{me,Length[meDratios]},{vs,Length[v0s]}];

(*ScanDictList=IterateKEchange[Flatten@ScanDictList,{0.05},"\[Alpha]"/.Constants`SIConstRepl];*)

ScanDictList
]


(* ::Subsection:: *)
(*Get Significant Capture Regimes*)


(* ::Subsubsection::Closed:: *)
(*KE Change *)


(* ::Input::Initialization:: *)
VCoulombE[N_,\[Alpha]D_]:=\[Alpha]D/("\[Alpha]") ("e")^2/(4 \[Pi] "\[Epsilon]0") N/("rE")/.Constants`EarthRepl/.Constants`SIConstRepl


(* ::Input::Initialization:: *)
Clear[KEchange]
KEchange[ScanDict_,fD_,\[Alpha]D_]:=Module[{mpD,meD,nD0,Nc,Ntot,pDvescchange,eDvescchange,pDKEfracchange,eDKEfracchange},
{mpD,meD}={"mpD","meD"}/.ScanDict;

nD0=naDM[fD,mpD +meD];

Nc="Nc"/.ScanDict;
Ntot = Nc/."nD"->nD0;

pDvescchange=Sqrt[(2 VCoulombE[Ntot,\[Alpha]D])/mpD]; 
eDvescchange=Sqrt[(2 VCoulombE[Ntot,\[Alpha]D])/meD];
pDKEfracchange=pDvescchange/Sqrt[("vesc")^2+3/("\[Beta]D" mpD)]/.Constants`EarthRepl/.ScanDict;
eDKEfracchange=eDvescchange/Sqrt[("vesc")^2+3/("\[Beta]D" meD)]/.Constants`EarthRepl/.ScanDict;

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
