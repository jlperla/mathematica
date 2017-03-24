(* ::Package:: *)

(*
 * $Rev: 1 $
 * Author: Jesse Perla (c) 2013
 * Use, modification and distribution are subject to the 
 * Boost Software License, Version 1.0. (See accompanying file 
 * LICENSE_ 1_ 0.txt or copy at http://www.boost.org/LICENSE_1_ 0.txt)
 
 The display scripted was modified from http://mathematica.stackexchange.com/questions/30884/displaying-index-as-subscript-on-output-e-g-ci-c-i-with-notation-or
 
*)


(* The following will convert from dollar signs to scripted, for display.  take a$b -> a[b], etc.*)
ToScriptedNotation[exp_] := 
exp /. s_Symbol /; StringMatchQ[SymbolName[s], "*$*"] :> ToExpression[
          StringJoin @@ Riffle[StringSplit[SymbolName[s], "$" -> "["], "]", {4, -1, 3}]];

  
symbols := Sequence[bar, hat, tilde, vec, dot, underbar, plus, minus, star];
(*
(* Probably should protect the symbols, but errors seem to popup*)
  
With[{symbols := Sequence[bar, hat, tilde, vec, underbar, plus, minus, star,dot]},
  Unprotect @ symbols; ClearAll @ symbols; Protect @ symbols;
]
*)

(*
makeDef[pat_, body_] := (
  MakeBoxes[a : pat, fmt_] := ToBoxes[Interpretation[body, a], fmt] /;
    StringMatchQ[ToString @ Unevaluated @ h, Union[$scriptedconstants, $scriptedfunctions]];
      MakeBoxes[a : pat[sub_], fmt_] := ToBoxes[Interpretation[body[sub], a], fmt] /;
        StringMatchQ[ToString @ Unevaluated @ h, $scriptedfunctions]
)
*)


NotScriptedVarQ[z_] := ! MemberQ[$scriptedfunctionsvars, HoldPattern @ z, Infinity];


makeDef[pat_, body_] := (
  MakeBoxes[a : pat, fmt_] := ToBoxes[Interpretation[body, a], fmt] /;
    MemberQ[Union[$scriptedconstants, $scriptedfunctions], Unevaluated @ h];
      MakeBoxes[a : pat[sub_], fmt_] := ToBoxes[Interpretation[body[sub], a], fmt] /;
        MemberQ[$scriptedfunctions, Unevaluated @ h]
)

set1 = {
  bar      -> OverBar,
  dot      -> OverDot,
  hat      -> OverHat,
  tilde    -> OverTilde,
  vec      -> OverVector,
  underbar -> UnderBar,(* After here subsuper is ncessary to have both annotations*)
  plus     -> SuperPlus,
  minus    -> SuperMinus,
  star     -> SuperStar
 };

set2 = {
  plus  -> "+",
  minus -> "-",
  star  -> "*"
 };

makeDef[h_Symbol[#], #2[h]] & @@@ set1;
makeDef[h_Symbol[argssub__][#], Subscript[#2[h], argssub]] & @@@ Take[set1, 5];
makeDef[h_Symbol[argsub_][#], Subsuperscript[h, argsub, #2]] & @@@ set2;
makeDef[h_Symbol[argssub__?NotScriptedVarQ], Subscript[h, argssub]]; (* This feature having trouble with Mma10*)



