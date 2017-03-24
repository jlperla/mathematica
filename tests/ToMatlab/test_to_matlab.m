(* ::Package:: *)

(*
% $Rev: 1 $
% $Date: 2015-01-29 15:42:05 -0800 (Thu, 29 Jan 2015) $
% $LastChangedBy: jlperla $
% Author: Jesse Perla (c) 2012
% Use, modification and distribution are subject to the 
% Boost Software License, Version 1.0. (See accompanying file 
% LICENSE_ 1_ 0.txt or copy at http://www.boost.org/LICENSE_1_ 0.txt)
*)

<<ETK`;
<< ToMatlab`;

CaptionPrint["Test for a few string rules"];
val = "f(l)(2)"
val2= "z(l)"
val3= "f(l,d)"
val5= "f(l)(1,2)"
string$rule = {f_ ~~ "("~~ sub__ ~~ ")(" ~~ arg__ ~~ ")" -> f~~"_" ~~sub ~~ "(" ~~ arg ~~")"};

StringMatchQ[val, f_ ~~ "("~~ sub__ ~~ ")(" ~~ arg__ ~~ ")" ]
StringMatchQ[val2, f_ ~~ "("~~ sub__ ~~ ")(" ~~ arg__ ~~ ")" ]

StringReplace[val, string$rule]
StringReplace[val2, string$rule]
StringReplace[val3, string$rule]
StringReplace[val5, string$rule]

constant$subs = {"l", "d"};
InConstantSubs[s_]:= StringMatchQ[s, constant$subs];
subscript$string$rule = {f_ ~~ "("~~ sub_?InConstantSubs ~~ ")" -> f~~"_" ~~ sub};
val = "z(l)"
val2 = "z(2)"
StringReplace[val, subscript$string$rule]
StringReplace[val2, subscript$string$rule]


<<ETK`;
<< ToMatlab`;
CaptionPrint["Some matlab testing"];
$constant$subscripts = {"l"};	
exp = x + \[Gamma] + z[bar] + f[l][z] + f[2][1] + z[l] + 2 * (f[bar][z] + q$temp + f[z]) + Sqrt[\[Alpha][l] + \[Beta][2]]
str = ToMatlab@exp
rule = {h -> x, h[z]-> z, h[z,x]-> z, h[x,z,y]-> y, h[l]-> z, h[bar]-> z, h[bar][z] -> y,
h[x,y]-> z[bar], h[z]->x^2 , h[l]-> z, h[z] :> l, h[l][z_] -> 4}
RulesToMatlab@rule



