(* ::Package:: *)

(*
 * $Rev: 18 $
 * Author: Jesse Perla (c) 2015
 * Use, modification and distribution are subject to the 
 * Boost Software License, Version 1.0. (See accompanying file 
 * LICENSE_ 1_ 0.txt or copy at http://www.boost.org/LICENSE_1_ 0.txt)
*)
(* Deriving approximations of moments from a MGF and taylor series.*)

If[Max[StringCount[$Path, "etk_mathematica"]] == 0, AppendTo[$Path, "c:\\working\\libraries\\etk_mathematica"]];
<< ETK`;
<< ToMatlab`;

$Assumptions = \[Lambda] > 0 && \[Alpha] > 0;
(* Example with Poisson *)
poisson$mgf= ConvertToPureFunction[MomentGeneratingFunction[PoissonDistribution[\[Lambda]],z],{z}];
f = #^\[Alpha]&
EFromMGF[ Exp[\[Lambda](Exp[#]-1)]&, f, 5]
EFromMGF[ poisson$mgf, #^\[Alpha]&, 5]

(* Example with Truncated Poisson *)
truncated$poisson$mgf = (Exp[\[Lambda] Exp[#]]-1)/(Exp[\[Lambda]]-1)&;
EFromMGF[ truncated$poisson$mgf, #^\[Alpha]&, 1]

(* Checking algebraic solutions *)
EFromMGF[ truncated$poisson$mgf, #&, 3] == (\[Lambda] Exp[\[Lambda]])/(Exp[\[Lambda]]-1)//FullSimplify
(EFromMGF[ truncated$poisson$mgf, #^2&, 3] - EFromMGF[ truncated$poisson$mgf, #&, 3]^2)== (\[Lambda] Exp[\[Lambda]])/(Exp[\[Lambda]]-1) (1 - \[Lambda]/(Exp[\[Lambda]]-1))//FullSimplify
(* PMF for the zero-truncated poisson*)
ztp$pmf = \[Lambda]^#/(#!(Exp[\[Lambda]]-1))&;

E$ztp$MGF = EFromMGF[ truncated$poisson$mgf, #^\[Alpha]&, 3]
E$ztp$PMF = EFromPMF[ztp$pmf, f, 8,1]
{E$ztp$MGF, E$ztp$PMF } /. {\[Lambda]-> .1, \[Alpha]->.2}



