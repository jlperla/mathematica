(* :Title: Ito's Lemma *)

(* :Author: Mark Fisher *)

(* :Context: ItosLemma` *)

(* :Mathematica Version: 3.0 *)

(* :History: Version 1.0, June 1999.
			 Version 1.1, June 2000.
			 	Simplified ItoD.
			 	Added option Scalarize to Diffusion.
			 	Added new versions of ItoMake for simplified entry that
			 		use new global symbols DriftSymbol and DiffusionSymbol.
			 Version 1.1.1, April 2003. Minor changes.
*)

(* :Sources:
A preliminary version was inspired by the package "Diffusion" written
by Steele and Stine (described in "Mathematica and Diffusions",
chapter 9 of "Economic and Financial Modeling with Mathematica", 1993,
Springer-Verlag, edited by Hal R. Varian).

For an introduction to stochastic differential equations and Ito's lemma,
see Bernt Øksendal, Stochastic Differential Equations: An Introduction with
Applications, Springer-Verlag (currently in its 5th edition).
*)

(* :Keywords: Ito's lemma, stochastic calculus, stochastic differential
equations, Brownian motion *)

(* :Summary:
This package implements Ito's lemma for any number of Ito processes
with any number of arbitrarily-correlated Brownian motions.

There are two main functions in this package: ItoMake and ItoD. ItoMake
should be used to declare an Ito process prior to using ItoD to compute the
stochastic derivative. For example (assuming the global symbols have been
set as in the housekeeping example below), ItoMake[x[t], \[Mu], \[Sigma]]
creates the global rule x[t + dt] -> x[t] + \[Mu] dt + \[Sigma] Subscript[dB, 1].
ItoD[f[x[t], t]] constructs a Taylor series for f[x[t + dt], t + dt] around
dt = 0 and Subscript[dB, 1] = 0. The series is first-order in dt and
second-order in Subscript[dB, 1], where "Ito multiplication rules" are
applied to the second-order term.

The functions Drift and Diffusion can be used to extract the drift and
diffusion from the output of ItoD.

(***** Housekeeping Notice *****)

The package relies on six global symbols: TimeSymbol, TimeIncrement,
BrownianIncrement, CorrelationSymbol, DriftSymbol, and DiffusionSymbol.
They can be set at the beginning of each session (after loading ItosLemma)
to more convenient symbols. Here is an example:

Clear[t, dt, dB, \[Rho], \[Mu], \[Sigma]]
{TimeSymbol, TimeIncrement, BrownianIncrement,
	CorrelationSymbol, DriftSymbol, DiffusionSymbol} =
	{t, dt, dB, \[Rho], \[Mu], \[Sigma]}

One way to make these assignments automatic is to copy the preceding
code to the bottom of this file, after the EndPackage[] statement.

There is a "brief explanation to the code" in the Notebook ItosLemma.nb.

*)

(* :Examples: Supplied in the Notebook ItosLemma.nb. *)

BeginPackage["ItosLemma`"]

ItosLemma::usage = "ItosLemma.m is a package that implements Ito's Lemma.
The package uses six global symbols: TimeSymbol, TimeIncrement,
BrownianIncrement, CorrelationSymbol, DriftSymbol, and DiffusionSymbol.
They can be defined in terms of more convenient symbols.\n
Example:\n
Clear[t, dt, dB, \[Rho], \[Mu], \[Sigma]];\n
{TimeSymbol, TimeIncrement, BrownianIncrement, CorrelationSymbol,
	DriftSymbol, DiffusionSymbol} = {t, dt, dB, \[Rho], \[Mu], \[Sigma]}"

TimeSymbol::usage = "TimeSymbol is the symbol that represents time."

TimeIncrement::usage = "TimeIncrement is the symbol that represents an
infinitesimal change in time."

BrownianIncrement::usage = "BrownianIncrement is the symbol that
represents an infinitesimal change in a Brownian motion."

CorrelationSymbol::usage = "CorrelationSymbol is the symbol that
represents the correlation between Brownian motions."

DriftSymbol::usage = "DriftSymbol is the symbol that represents the
drift."

DiffusionSymbol::usage = "DiffusionSymbol is the symbol that represents
the diffusion."

ItoD::usage = "ItoD[expr] applies Ito's lemma to expr. ItoD takes the
option OrthogonalBrownians."

RelativeItoD::usage = "RelativeItoD[expr] computes ItoD[expr]/expr."

OrthogonalBrownians::usage = "OrthogonalBrownians is an option for ItoD. The
default setting is OrthogonalBrownians -> True."

SuppressTime::usage = "SuppressTime is an option for ItoD that specifies
whether to suppress the dependence on time in the output. The default
setting is SuppressTime -> True."

ItoMake::usage = "ItoMake[x[args], mu, sigma] associates a rule with the
(ultimate) head of the x[args]. The TimeSymbol should be one of the
arguments; for example, ItoMake[x[t], \[Mu], \[Sigma]]. For multiple
Brownians, sigma may be a list."

ItoMakeFromItoD::usage = "ItoMakeFromItoD[name, expr] calls
ItoMake[x[args], Drift[dx], Diffusion[dx]], where dx is the output from
ItoD. ItoMakeFromItoD passes the option BrownianList to Diffusion."

IncludeArguments::usage = "IncludeArguments is an option for ItoMake (and
related functions) that specifies whether to include the arguments in the
drift and diffusion specifications. IncludeArguments is only used in
vector versions of ItoMake-related functions."

RelativeItoMake::usage = "RelativeItoMake[name, mu, sigma] calls
ItoMake[name, name mu, name sigma]."

RelativeItoMakeFromItoD::usage = "RelativeItoMakeFromItoD[name, expr] calls
ItoMake[name, name Drift[expr], name Diffusion[expr]]. RelativeItoMakeFromItoD
passes the option BrownianList to Diffusion."

ExponentialItoMake::usage = "ExponentialItoMake[name, mu, sigma] calls
ItoMake[name, name (mu + sigma^2/2), name sigma]."

ExponentialItoMakeFromItoD::usage = "ExponentialItoMakeFromItoD[name, expr]
calls ItoMake[name (Drift[expr] + (1/2)Diffusion[expr].Diffusion[expr]),
Diffusion[expr]]. ExponentialItoMakeFromItoD passes the option BrownianList
to Diffusion."

VectorItoMake::usage = "VectorItoMake[name, n] makes a vector of n Ito
processes."

Drift::usage = "Drift[expr] returns the drift of an Ito process."

Diffusion::usage = "Diffusion[expr] returns the diffusion of an Ito
process. Diffusion takes the options BrownianList and Scalarize."

BrownianList::usage = "BrownianList is an option for Diffusion. The default
setting is BrownianList -> Automatic."

Scalarize::usage = "Scalarize is an option for Diffusion. The default
setting is Scalarize -> True, which specifies that a length-one Diffusion
should be returned as a scalar."

Begin["`Private`"]
(*****************  code starts here **********************)

(***** main engine *****)

Options[ItoD] = {OrthogonalBrownians -> True, SuppressTime -> True}

SetAttributes[ItoD, Listable]

ItoD[expr_, opts___?OptionQ] :=
	Module[{t, dt, dB, rho, suppress, ortho, zerorule,
		dtexpr, brownians, dtime, dbrown, d2brown, timedrift,
		diffusion, jensen},
	(* get the global values for these symbols *)
	{t, dt, dB, rho} = {TimeSymbol, TimeIncrement, BrownianIncrement,
		CorrelationSymbol};
	{suppress, ortho} = {SuppressTime, OrthogonalBrownians} /.
		{opts} /. Options[ItoD];
	zerorule = Subscript[dB, _] | dt -> 0;
	(* identify the Ito processes *)
	dtexpr = expr /. t -> t + dt;
	(* find the brownians *)
	(* if there are no brownians, make one up *)
	brownians = Union[Cases[dtexpr, Subscript[dB, _], {0, Infinity}]];
	If[brownians === {}, brownians = {Subscript[dB, 1]}];
	(* compute the derivatives *)
	dtime = D[dtexpr, dt];
    dbrown = D[dtexpr, #]& /@ brownians;
	d2brown = D[dbrown, #]& /@ brownians;
	(* compute the "ito-taylor series" *)
	timedrift = (dtime /. zerorule) dt;
	diffusion = (dbrown /. zerorule) . brownians;
	jensen = (1/2) *
		Expand[brownians . (d2brown /. zerorule) . brownians] /.
			{Subscript[dB, _]^2 -> dt,
			 Subscript[dB, i_] Subscript[dB, j_] ->
				If[TrueQ[ortho],
					0,
					Subscript[rho, i, j] dt]};
	(* assemble the parts *)
	(diffusion + Collect[timedrift + jensen, dt]) /.
		If[TrueQ[suppress], f_[t] -> f, {}]
	]

(* related function to simplify output in special cases *)

RelativeItoD[expr_, opts___?OptionQ] :=
	Collect[ItoD[expr, opts]/expr /.
		If[TrueQ[SuppressTime /. {opts} /. Options[ItoD]],
			(* then *)
			f_[TimeSymbol] -> f,
			(* else *)
			{}],
		TimeIncrement, Simplify]

(***** extractors *****)

Drift[expr_] := Coefficient[expr, TimeIncrement]

Diffusion::lost = "BrownianList omitted ``."

Options[Diffusion] = {BrownianList -> Automatic, Scalarize -> True}

Diffusion[expr_, opts___?OptionQ] :=
	Module[{dB, blist, scalar, found, brownians, lost},
	(* get the global value for this symbol *)
	dB = BrownianIncrement;
	{blist, scalar} = {BrownianList, Scalarize} /. {opts} /.
		Options[Diffusion];
	found = Union[Cases[expr, Subscript[dB, _], {0, Infinity}]];
	If[TrueQ[blist === Automatic],
		(* then *)
		brownians = found,
		(* else *)
		brownians = blist;
		lost = Complement[found, blist];
		If[lost =!= {}, Message[Diffusion::lost, lost]]];
	If[TrueQ[scalar] && Length[brownians] === 1,
		(* then *)
		Coefficient[expr, First @ brownians],
		(* else *)
		Coefficient[expr, #]& /@ brownians]
	]

(***** constructor *****)

ItoMake::notime = "TimeSymbol `1` does not appear as an argument in `2`."
ItoMake::badhead = "The Head of `1` matches one of `2`, `3`, or `4`."
ItoMake::badarg = "One or more arguments of `1` matches either `2` or `3`."

(* only used in vector functions (see below) *)
Options[ItoMake] = {IncludeArguments -> False}

ItoMake[x_[args__], mu_, sig_] :=
	Module[{t, dt, dB, diffusion, head, lhs},
	(* get the global values for these symbols *)
	{t, dt, dB} = {TimeSymbol, TimeIncrement, BrownianIncrement};
	If[FreeQ[{args}, t],
		Message[ItoMake::notime, t, x[args]]; Return[$Failed]];
	If[MatchQ[x, t | dt | dB],
		Message[ItoMake::badhead, x[args], t, dt, dB]; Return[$Failed]];
	If[MemberQ[{args}, dt | dB, Infinity],
		Message[ItoMake::badarg, x[args], dt, dB]; Return[$Failed]];
	diffusion = Switch[sig,
		_List, 	sig . Array[Subscript[dB, #]&, Length[sig]],
		_, 		sig Subscript[dB, 1]];
	(* find the ultimate Head *)
	head = FixedPointList[Head, x][[-3]];
	lhs = Block[Evaluate[{head}],
		(* construct the pattern and Hold it *)
		Hold @ Evaluate[
			If[# === t, # + dt, Pattern[#, Blank[]]]& /@ x[args]
			]];
	(* define the rule and return the SDE *)
	(Set @@ Append[lhs, x[args] + mu dt + diffusion]) - x[args]
	]

(*****************************************************************)
(* ItoMake-related functions to simplify input for special cases *)
(*****************************************************************)

RelativeItoMake[x_[args__], mu_, sig_] :=
	ItoMake[x[args], x[args] mu, x[args] sig]

ExponentialItoMake[x_[args__], mu_, sig_]:=
	ItoMake[x[args], x[args] (mu + sig.sig/2) /. Dot -> Times, x[args] sig]

(* make Ito process from output of ItoD *)

ItoMakeFromItoD[x_[args__], expr_, opts___?OptionQ] :=
	With[{diff = Diffusion[expr, opts]},
	ItoMake[x[args], Drift[expr], If[diff === {}, {0}, diff]]
	]

RelativeItoMakeFromItoD[x_[args__], expr_, opts___?OptionQ] :=
	ItoMake[x[args], x[args] Drift[expr],
		x[args] Diffusion[expr, opts]]

ExponentialItoMakeFromItoD[x_[args__], expr_, opts___?OptionQ] :=
	With[{diff = Diffusion[expr, Scalarize -> False, opts]},
		ItoMake[x[args], x[args] (Drift[expr] + (1/2) diff.diff),
			x[args] diff]
		]

(* scalar process with vector browians *)

ItoMake[x_[args__], mu_, sig_, n_Integer, opts___?OptionQ] :=
	ItoMake[x[args], ##]& @@ makeargs[x[args], mu, sig, n, opts]

RelativeItoMake[x_[args__], mu_, sig_, n_Integer, opts___?OptionQ] :=
	RelativeItoMake[x[args], ##]& @@ makeargs[x[args], mu, sig, n, opts]

ExponentialItoMake[x_[args__], mu_, sig_, n_Integer, opts___?OptionQ] :=
	ExponentialItoMake[x[args], ##]& @@ makeargs[x[args], mu, sig, n, opts]

(* auxiliary function *)
makeargs[x_[args__], mu_, sig_, n_, opts___?OptionQ] :=
	{mu[x][args], Array[Subscript[sig, #][x][args]&, n]} /. 
		If[TrueQ[IncludeArguments /. {opts} /. Options[ItoMake]], {},
			f_[x][args] :> f[x]]

(* vector Ito processes *)

(* n Brownian motion processes *)
VectorItoMake[x_[args__], 0, 1, n_Integer] :=
	Table[ItoMake[Subscript[x, i][args], 0,
		Table[If[i == j, 1, 0], {j, n}]], {i, n}]

(* n Ito processes, n Brownians *)
VectorItoMake[x_[args__], mu_, sig_, n_Integer, opts___?OptionQ] :=
	Table[ItoMake[Subscript[x, i][args], mu, sig, n, opts], {i, n}]

(* m Ito processes, n Brownians *)
VectorItoMake[x_[args__], mu_, sig_, {m_Integer, n_Integer},
		opts___?OptionQ] :=
	Table[ItoMake[Subscript[x, i][args], mu, sig, n, opts], {i, m}]

(* use default symbols for drift and diffusion *)

(#[x_[args__]] := #[x[args],
	Subscript[DriftSymbol, x], Subscript[DiffusionSymbol, x]])& /@
		{ItoMake, RelativeItoMake}

Options[ExponentialItoMake] = {Global`OverTilde -> True}

ExponentialItoMake[x_[args__], opts___?OptionQ] :=
	Module[{tilde, fun},
	tilde = Global`OverTilde /. {opts} /. Options[ExponentialItoMake];
	fun = If[TrueQ[tilde], Global`OverTilde, Identity];
	ExponentialItoMake[x[args],
			Subscript[fun[DriftSymbol], x],
			Subscript[DiffusionSymbol, x]]
	]

(#[x_[args__], n_Integer, opts___?OptionQ] := #[x[args],
	DriftSymbol, DiffusionSymbol, n, opts])& /@
		{ItoMake, RelativeItoMake, ExponentialItoMake}

VectorItoMake[x_[args__], n_Integer, opts___?OptionQ] :=
	VectorItoMake[x[args], DriftSymbol, DiffusionSymbol, n, opts]

VectorItoMake[x_[args__], {m_Integer, n_Integer}, opts___?OptionQ] :=
	VectorItoMake[x[args], DriftSymbol, DiffusionSymbol, {m, n}, opts]

End[]
EndPackage[]
