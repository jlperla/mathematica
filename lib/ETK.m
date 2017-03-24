(* ::Package:: *)

(*
 * $Rev: 25 $
 * Author: Jesse Perla (c) 2009, 2010, 2013
 * Use, modification and distribution are subject to the 
 * Boost Software License, Version 1.0. (See accompanying file 
 * LICENSE_ 1_ 0.txt or copy at http://www.boost.org/LICENSE_1_ 0.txt)
 
 The Equality threading was modified from  © Copyright 1997, Roman E. Maeder.
 The display scripted was modified from http://mathematica.stackexchange.com/questions/30884/displaying-index-as-subscript-on-output-e-g-ci-c-i-with-notation-or
 
*)

(*
(* Usage *)
TestFunction::"usage"="Identity function as test for packages";
CaptionPrint::"usage"="Prints the argument as large caption with formatting";
BigCaptionPrint::"usage"="Even larger and bold";
DumpGlobals::"usage"="Outputs all global data in a dump file";
Second::"usage"="Gets the second element in a list";
FlipInequalitySign::"usage"="Changes the signs of all inequalities.  e.g. multiplying by < 0";
ReverseSubstitutions::"usage"="Reverses a list of substitutions.  Useful for constant groupings and verification";
NegateSubstitutions::"usage"="Reverses sign of a list of substitutions.  Useful for constant groupings and verification";
FlipSubstitutions::"usage"="Flips the  substitutions";
SubstitutionVariations::"usage"="Adds a bunch of variations of the substitutions to help in simplification";
DisplayDefinitions::"usage"="Takes a list of rules and displays as definitions"
*)

	Print["Clearing and Loading etk libraries"]
	ClearAll["Global`*"]
	
	
	(* A simple function for testing packages *)
	TestFunction[x_] := 2 * x;

	CaptionPrint[x_] := Print[Style[x, Bold, FontSize-> 16 ]];
	BigCaptionPrint[x_] := Print[Style[x, Bold, Underlined, FontSize-> 18 ]];
	
	ReverseSubstitutions[subs_]:= Reverse[subs, 2]; (* Reverses the direction.  e.g. {a->b, c->d} --> {b->a, d->c}  Use to get rid of subs to check results *)	
	NegateSubstitutions[subs_]:= Map[-#[[1]] -> -#[[2]] &, subs] (* Reverses the sign.  e.g. {a->b} --> {-a -> -b};*)
	FlipSubstitutions[subs_]:= Map[1/#[[1]] -> 1/#[[2]] &, subs] (* {a->b}	--> {1/a -> 1/b} *)
	SubstitutionVariations[subs_]:= Module[{all$subs=Union[subs, NegateSubstitutions@subs, FlipSubstitutions@subs, NegateSubstitutions@FlipSubstitutions@subs]},
		Union[all$subs, FullSimplify@all$subs]]; (* Applies all sorts of variations to the substitutions which may show up in patterns *)
	DisplayDefinitions[undo$subs_]:= #[[1]] \[Congruent]  #[[2]] & /@ undo$subs;
		

	DumpGlobals[save$name_, sub$directory_: "", postfix_ : ""] := Module[{notebookDirectoryExists,outputDirectory},
		notebookDirectoryExists = Check[NotebookDirectory[],0];
		outputDirectory = If[notebookDirectoryExists == 0, Directory[], notebookDirectoryExists];
		If[Not[DirectoryQ[outputDirectory <> sub$directory]],CreateDirectory[outputDirectory <> sub$directory],0];
		SetDirectory[outputDirectory <> sub$directory];
		CaptionPrint["Saving Output to " <> outputDirectory <> sub$directory <>"/"<> save$name <> postfix <>".mx"];
		DumpSave[save$name <> postfix <>".mx", "Global`"];
		SetDirectory[outputDirectory];
	];	
	
	SetWorkingDirectory[] := Module[{outputDirectory,notebookDirectoryExists},
		notebookDirectoryExists = Check[NotebookDirectory[], 0];
		outputDirectory = If[notebookDirectoryExists == 0, Directory[], 
		notebookDirectoryExists];
		SetDirectory[outputDirectory];
		]
	
	PowerExpandSimplify := Composition[FullSimplify, PowerExpand];
	
	
	Second[seq_] := seq[[2]]; %Simple utility to get the second element


	(* The following is modified from the Equation Threading package *)
	listableQ[f_] := MemberQ[Attributes[f], Listable]

	protected = Unprotect[Equal]

	Equal/: lhs:f_Symbol?listableQ[___, _Equal, ___] :=
		Thread[ Unevaluated[lhs], Equal ]

	Protect[Evaluate[protected]]

	(* VERY IMPORTANT!  The transformation of inequalities doesn't take into account multiplications by negative numbers (it can't except to generate multiple cases which can be refined).
	Do  manually by the FlipSign
	
	Also, equations with mismatched equalities and inequalities cannot be added directly, etc.  To do this, multiply by -1 and then flip the sign so they match.  e.g.
	ineq = a  x > b; ineq2 = x < c; ineq2 + (-ineq // FlipInequalitySign)
	*)

	protected = Unprotect[Less]
	Less/: lhs:f_Symbol?listableQ[___, _Less, ___] :=
		Thread[ Unevaluated[lhs], Less ]
	Protect[Evaluate[protected]]

	protected = Unprotect[LessEqual]
	LessEqual/: lhs:f_Symbol?listableQ[___, _LessEqual, ___] :=
		Thread[ Unevaluated[lhs], LessEqual ]
	Protect[Evaluate[protected]]

	protected = Unprotect[Greater]
	Greater/: lhs:f_Symbol?listableQ[___, _Greater, ___] :=
		Thread[ Unevaluated[lhs], Greater ]
	Protect[Evaluate[protected]]

	protected = Unprotect[GreaterEqual]
	GreaterEqual/: lhs:f_Symbol?listableQ[___, _GreaterEqual, ___] :=
		Thread[ Unevaluated[lhs], GreaterEqual ]
	Protect[Evaluate[protected]]

	(* Use if multiplying by a negative number, etc.*)
	FlipInequalitySign[exp_]:= exp /. {Less -> Greater, LessEqual -> GreaterEqual, Greater -> Less, GreaterEqual -> LessEqual}

 CollectFullSimplify[exp_, pat_] :=
  Replace[
    Collect[exp, pat, FullSimplify],
    x : pat :> RuleCondition @ FullSimplify @ x,
    {0, -1}
  ]
(* Why Fixed point? Want to allow the recollection after simplification of exponents, but having trouble finding a minimal example. Might be unncessary... *)
(* CollectFullSimplify[exp_, pat_]:= FixedPoint[CollectFullSimplifyImpl[#, pat]&, exp, 10] (* 10 arbitrary just in case non-monotone *) *)

(* The following will convert from dollar signs to scripted, for display.  take a$b -> a[b], etc.*)

ToScriptedNotation[exp_] := 
exp /. s_Symbol /; StringMatchQ[SymbolName[s], "*$*"] :> ToExpression[
          StringJoin @@ Riffle[StringSplit[SymbolName[s], "$" -> "["], "]", {4, -1, 3}]];


(* See http://mathematica.stackexchange.com/questions/9570/how-do-i-reassign-canonical-ordering-of-symbols *)
$canonicalorder = {_Integer};
CanonicalDisplay[expr_] :=
  Module[{h, rls},
    rls = MapIndexed[x : # :> h[#2, x] &,  $canonicalorder];
    HoldForm @@ {ToScriptedNotation[expr] /. rls} /. h[_, x_] :> x
  ]
 
  (* Utility function.  Necessary if using DisplayScripted as a workaround for a bug with Interpretation.  Copy/past as plain text to ignore the escape characters*)
  ToLatex[exp_] := Module[{},
  HoldForm[exp];
  str = ToString[exp, TeXForm];
   ReleaseHold[exp];
  str] 

  (* From http://mathematica.stackexchange.com/questions/8142/how-can-i-separate-a-separable-function/8146#8146 *)
  ClearAll[getGX];
getGX[expr_, xvar_, yvar_] :=
  With[{dlogg = D[Log[expr], xvar] // FullSimplify},
     Exp[Integrate[dlogg, xvar]] /; FreeQ[dlogg, yvar]];

Clear[getHY];
getHY[expr_, xvar_, yvar_] := FullSimplify[(#/getGX[#, xvar, yvar]) &[expr]]
SeparateFunction[expr_, xvar_, yvar_] := {getGX[expr, xvar, yvar], getHY[expr, xvar, yvar]};
  
assert$threshold = 1.*^-9;
AssertSmall[exp_] := Assert[Abs[exp] < assert$threshold]
AssertTrue[exp_] := If[TrueQ[exp //FullSimplify], {}, BigCaptionPrint["ERROR!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"]]
On[Assert]

MapToFunctions[f_, vals_] := MapThread[#1 -> #2 &, {f, vals}] (*e.g. MapToFunctions[{v$l[z], v$h[z]}, {0, 1}] -> {v$l[z] -> 0, v$h[z] -> 1} *)
ChopSmallTerms[exp_, threshold_: Exp[-10]] := Expand@exp //. {(Times[z_, __] /; Abs[z] < threshold) -> 0}

Integer$moment$from$M[M_, n_]:= D[M[z], {z, n}]/. {z-> 0};
(* Find an approximation of a moment form a taylor series of the expectation*)
EFromMGF[mgf_, f_, n$max_]:= ((CoefficientList[ Series[h[x], {x, E$x[1], n$max}] //Normal, x].Union[{1}, E$x /@ Range[1,n$max]]//FullSimplify) //.Table[ E$x[n] -> (Integer$moment$from$M[mgf, n]), {n, 1, n$max}]) /. {h -> f} //FullSimplify

EFromPMF[pfm_, f_, n_, start_: 0] := 
 (ztp$pmf /@ Range[start, n]).(f /@ Range[start, n]) // FullSimplify


SetAttributes[ConvertToPureFunction, HoldAll];
ConvertToPureFunction[expr_, vars_List] := 
 With[{variables = Unevaluated@vars}, 
  Block[variables, 
   Evaluate@(Hold[expr] /. Thread[vars -> Slot /@ Range@Length@vars]) & // ReleaseHold
   ]]

CalculateElasticity[f_, var_] := 
 D[Log[f //. {var -> Exp[log$var]}] // PowerExpandSimplify, 
   log$var] /. {log$var -> Log[var]}
ToElasticityRule[elast_, f_, var_] := {f'[var] -> elast f[var] / var}
