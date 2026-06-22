(* ::Package:: *)

BeginPackage["DarkAtomicFormFactor`"];

DarkReducedMass::usage =
  "DarkReducedMass[m1,m2] gives m1 m2/(m1+m2). The two masses must use the same units.";

DarkHydrogenElectronFF::usage =
  "DarkHydrogenElectronFF[q,alphaD,meD,mpD] gives the rescaled hydrogen electron-cloud form factor f_D(q). q, meD, and mpD are in eV.";

DarkHydrogenChargeFF::usage =
  "DarkHydrogenChargeFF[q,alphaD,meD,mpD] gives the neutral dark-hydrogen charge amplitude Zfit-f_D(q). q, meD, and mpD are in eV.";

DarkHydrogenAFF2::usage =
  "DarkHydrogenAFF2[q,alphaD,meD,mpD] gives |DarkHydrogenChargeFF|^2, the multiplicative atomic screening factor for a rate or cross section.";

DarkHydrogenScreeningMomentum::usage =
  "DarkHydrogenScreeningMomentum[alphaD,meD,mpD] gives alphaD times the dark reduced mass, in the same mass units as meD and mpD.";

HydrogenFitCoefficients::usage =
  "HydrogenFitCoefficients is the association containing the four-Gaussian SM hydrogen fit.";

Begin["`Private`"];

(* ---------- constants ---------- *)

alphaEM = 1/137.035999084;
meSM = 510998.95069;          (* eV *)
mpSM = 938272088.16;          (* eV *)
hbarcEVAngstrom = 1973.269804; (* eV Angstrom *)

HydrogenFitCoefficients = <|
  "a" -> {0.489918, 0.262003, 0.196767, 0.049879},
  "b" -> {20.6593, 7.74039, 49.5519, 2.20159},
  "c" -> 0.001305
|>;

(* Use the fitted value at q=0 rather than setting Z=1 by hand, so that
   the neutral charge form factor vanishes exactly at q=0. *)
zFit = Total[HydrogenFitCoefficients["a"]] + HydrogenFitCoefficients["c"];
muH = meSM mpSM/(meSM + mpSM);

(* ---------- public functions ---------- *)

DarkReducedMass[m1_?NumericQ, m2_?NumericQ] := m1 m2/(m1 + m2);

DarkHydrogenScreeningMomentum[
    alphaD_?NumericQ,
    meD_?NumericQ,
    mpD_?NumericQ
] := alphaD DarkReducedMass[meD, mpD];

DarkHydrogenElectronFF[
    qEV_?NumericQ,
    alphaD_?NumericQ,
    meD_?NumericQ,
    mpD_?NumericQ
] := Module[
  {muD, scale, qAInv, x, a, b, c},

  If[alphaD <= 0 || meD <= 0 || mpD <= 0,
    Return[$Failed]
  ];

  muD = DarkReducedMass[meD, mpD];

  (* Rescale the physical momentum into the argument of the SM hydrogen fit:
       q_SM,effective = q (alphaEM muH)/(alphaD muD). *)
  scale = (alphaEM muH)/(alphaD muD);

  (* Convert q from eV to inverse Angstrom. *)
  qAInv = qEV/hbarcEVAngstrom;
  x = scale qAInv/(4 Pi);

  a = HydrogenFitCoefficients["a"];
  b = HydrogenFitCoefficients["b"];
  c = HydrogenFitCoefficients["c"];

  Total[a Exp[-b x^2]] + c
];

DarkHydrogenChargeFF[
    qEV_?NumericQ,
    alphaD_?NumericQ,
    meD_?NumericQ,
    mpD_?NumericQ
] := zFit - DarkHydrogenElectronFF[qEV, alphaD, meD, mpD];

DarkHydrogenAFF2[
    qEV_?NumericQ,
    alphaD_?NumericQ,
    meD_?NumericQ,
    mpD_?NumericQ
] := Abs[DarkHydrogenChargeFF[qEV, alphaD, meD, mpD]]^2;

End[];
EndPackage[];
