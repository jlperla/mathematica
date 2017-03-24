(* ::Package:: *)

(*
 * $Rev: 1 $
 * Original code: http://www.ericswanson.us/perturbation.html
 * Modified by: Jesse Perla
*)

Print["This is Perturbation AIM Full Version 2.7"] ;
(*                                                                                                           *)
(* Credits: original algorithm designed by Gary Anderson, 2002-3;                                            *)
(*   corrected, optimized, and recoded by Eric Swanson, 2004.                                                *)
(*   Additional features and optimizations by Eric Swanson, 2008, 2009.                                      *)
(* Code is currently maintained by Eric Swanson, see http://www.ericswanson.org or                           *)
(*   http://www.ericswanson.us for most recent version and for background, instructions, capabilities,       *)
(*   examples, and tips.                                                                                     *)
(* This version of the code has only been tested with Mathematica version 6.  Perturbation AIM version 2.4   *)
(*   will work in both Mathematica 5 and Mathematica 6, but lacks some of the functionality and numerical    *)
(*   robustness of the current perturbation AIM version.  While it is not too hard to modify the code to     *)
(*   work in Mathematica 4 and earlier, those earlier versions are *much* slower at solving linear systems   *)
(*   of equations, and a few of the function calls are also a bit more awkward.                              *)
(* This is the "advanced" or "Full" version of the perturbationAIM code.  It is faster and more numerically  *)
(*   robust than older versions of the code, but the underlying algorithm is also a little bit more opaque   *)
(*   as a result of these optimizations.  If you are trying to understand how the code works, there is a     *)
(*   "simple" version of perturbation AIM 2.4 available at the web site above, which does not feature many   *)
(*   of these optimizations.  To understand the code, it is recommended that you go through the "simple"     *)
(*   version first, and then go through the code below.                                                      *)
(* No matter which version you are using, you should compute your final set of results using the arbitrary   *)
(*   precision capabilities of the code, as these are *much* more likely to be free of numerical problems.   *)
(*   As a general rule, I have found that 50 digits of precision (AIMPrecision = 50) and a zero tolerance    *)
(*   of 10^-10 (AIMZeroTol = 10^-10) provide adequate robustness for most models.  However, very highly      *)
(*   curved models, such as those that arise in finance, may require even higher levels of precision to      *)
(*   generate reliable answers.                                                                              *)
(*                                                                                                           *)

SetOptions["stdout",PageWidth->120] ; (* most people don't use 80 columns anymore; set the default to 120 *)
<<Combinatorica`                      (* we'll use the Compositions function *)

AIMZeroTol = 10^-10 ;
AIMPrecision = $MachinePrecision ;

(* Useful utilities that we will repeatedly call.                                                            *)
(* The primary improvement here with respect to the "simple" version is the introduction of a new routine,   *)
(*   AIMBackLookingVars, that searches through the model for any variables that are clearly "backward-       *)
(*   looking " (i.e., obviously a function of just a few state variables, typically one or two, rather than   *)
(*   the whole set of state variables of the model).  You could call these clearly "predetermined"           *)
(*   variables if you like, but I've avoided using that term.  This helps speed up the code by reducing      *)
(*   the amount of differentiation that must be done and imposing a number of zero restrictions on the       *)
(*   solution.                                                                                               *)

AIMGenericArgs[eqns_]:= AIMGenericArgs[eqns] = Table[Unique["arg"],{AIMNArgs[eqns]}] ;

AIMLagVars[eqns_]:= AIMLagVars[eqns] =
Flatten[Map[Table[#[t+i],{i,Min[Cases[eqns,#[t+j_.]->j,Infinity]],-1}]&, AIMVarNames[eqns]]] ;

AIMMaxLag[eqns_]:= AIMMaxLag[eqns] = -Min[Cases[eqns,x_Symbol[t+i_.]->i,Infinity]] ;

AIMMaxLead[eqns_]:= AIMMaxLead[eqns] = Max[Cases[eqns,x_Symbol[t+i_.]->i,Infinity]] ;

AIMNArgs[eqns_]:= AIMNArgs[eqns] = Length[AIMStateVars[eqns]] ;

AIMNEqns[eqns_]:= AIMNEqns[eqns] = Length[Flatten[{eqns}]] ;

AIMNVars[eqns_]:= AIMNVars[eqns] = Length[AIMVarNames[eqns]] ;

AIMShocks[eqns_]:= AIMShocks[eqns] = Union[Cases[eqns,eps[_][t],Infinity]] ;

AIMSSSubs={Sigma->0, eps[_][_]->0, x_[t+_.]:>Symbol[SymbolName[x]<>"AIMSS"]} ;

AIMSSVars[eqns_]:= AIMSSVars[eqns] = Map[Symbol[SymbolName[#]<>"AIMSS"]&, AIMVarNames[eqns]] ;

AIMStateVars[eqns_]:= AIMStateVars[eqns] = Join[AIMLagVars[eqns], AIMShocks[eqns], {Sigma}] ;

AIMVarNames[eqns_]:= AIMVarNames[eqns] = Union[Cases[eqns,x_Symbol[t+i_.]->x,Infinity]] ;

(* Here is the AIMBackLookingVars routine.  It first finds all equations that have no leads.  For those      *)
(*   equations that have only one variable dated t, it is then obvious that the solution to this variable    *)
(*   is entirely a function of the lagged variables and shocks in that particular equation.  This            *)
(*   particular variable can then be solved in terms of a reduced state vector, so that its particular       *)
(*   bFunc has fewer arguments, hence fewer derivatives and fewer undetermined coefficients that must be     *)
(*   solved.                                                                                                 *)
(* The routine then iterates, so that a variable is also counted as clearly backward-looking if the only     *)
(*   other variables in its equation are lagged variables, date t shocks, and other variables dated t that   *)
(*   have already been shown to be clearly backward-looking.                                                 *)

AIMBackLookingVars[eqns_]:= AIMBackLookingVars[eqns] =
Module[{eqnsnoleads,symbols,nsv,nsvargs,possiblensv,nextpos},
eqnsnoleads = Select[Flatten[{eqns}],Cases[#,_[t+i_/;i>0],Infinity] =={}&] ;
symbols = Map[Join[AIMLagVars[#],Union[Cases[#,_[t],Infinity]]]&, eqnsnoleads] ;
{nsv,nsvargs} = FixedPoint[
	Function[nsv,
	  symbols = symbols /.Thread[nsv[[1]]->nsv[[2]]] ;
	  possiblensv = Map[Cases[#,_Symbol[t]]&, symbols] ;
	  nextpos = Flatten[Position[Map[Length,possiblensv], 1]] ;
	 {Join[nsv[[1]], Flatten[possiblensv[[nextpos]]]],
	    Join[nsv[[2]], DeleteCases[symbols[[nextpos]],_Symbol[t],Infinity]]}
	], {{},{}}] ;
{nsv /.x_[t]:>x, Map[Union[Flatten[#]]&,nsvargs]}
] ;

AIMBLVPos[eqns_]:= AIMBLVPos[eqns] = Flatten[Map[Position[AIMVarNames[eqns],#]&, AIMBackLookingVars[eqns][[1]]]] ;


(* Calculate the steady state *)

AIMSS[eqns_,parameters_List,opts___Rule]:= AIMSS[eqns,parameters,AIMZeroTol,AIMPrecision] =
Module[{ssguess,sseqns,symbols,ss,zeropos,nonzeropos,sseqns2,ss2},
 ssguess = AIMSSGuess /.{opts} /.AIMSSGuess->AIMSSGuess[eqns] /. x_Real:>SetPrecision[x,AIMPrecision] ;
 sseqns = eqns /.Equal->Subtract //.parameters /.AIMSSSubs /. x_Real:>SetPrecision[x,AIMPrecision] ;
 AIMModelDiagnostics[eqns] ;
 If[Length[ssguess]=!=AIMNVars[eqns], Print["Model has ",AIMNVars[eqns]," variables but AIMSSGuess has ",Length[ssguess]," elements"]; $Failed,
  Print["Finding steady state, time is ", Date[],", starting from initial guess: ", Thread[AIMVarNames[eqns]->N[ssguess]]] ;
  symbols = Complement[Union[Cases[sseqns,_Symbol,Infinity]], Join[AIMSSVars[eqns],{E}]] ;
  If[Length[symbols]>0, Print["Warning: found symbols ",symbols," in equations"]] ;
  ss = Chop[If[AIMPrecision==$MachinePrecision, FindRoot[sseqns,Transpose[{AIMSSVars[eqns],ssguess}],MaxIterations->5000],
	FindRoot[sseqns,Transpose[{AIMSSVars[eqns],ssguess}],MaxIterations->5000,WorkingPrecision->.96*AIMPrecision]],
	AIMZeroTol] ; 
  If[Head[ss]===FindRoot, $Failed, (* return $Failed if FindRoot cannot solve, else continue: *)
  zeropos = Flatten[Position[Map[Last,ss], 0]] ;
  AIMSSGuess[eqns] = Chop[Map[Last,ss]] ; (* update default SS guess *)
  AIMRefineSS[eqns,sseqns,zeropos]
]]] ;

AIMSS[eqns_,opts__Rule]:= AIMSS[eqns,{},opts] ; (* if AIMSS called without separate parameters, insert {} *)
AIMSS[eqns_]:= AIMSS[eqns,{}] ;
AIMSS[eqns_,parameters_List,zerotol_,precision_]:= AIMSS[eqns,parameters] ;

AIMSSGuess[eqns_]:= Table[0,{AIMNVars[eqns]}] ; (* sets the default steady-state initial guess *)

AIMRefineSS[eqns_,sseqns_,zeropos_]:=
Module[{nonzeropos,sseqns2,ss2,zeropos2},
 nonzeropos = Complement[Range[AIMNVars[eqns]], zeropos] ;
 If[Length[zeropos]>0,
   Print["Found steady-state values < AIMZeroTol for ", Length[zeropos]," variables: ", AIMVarNames[eqns][[zeropos]]] ;
   Print["Imposing hard zeros for these variables and re-solving, time is ", Date[]]] ;
 sseqns2 = DeleteCases[sseqns /.Thread[AIMSSVars[eqns][[zeropos]]->0], True] ;
 ss2 = Chop[If[AIMPrecision==$MachinePrecision, FindMinimum[sseqns2.sseqns2, Transpose[{AIMSSVars[eqns][[nonzeropos]],
     						AIMSSGuess[eqns][[nonzeropos]]+AIMZeroTol}],MaxIterations->5000],
					FindMinimum[sseqns2.sseqns2, Transpose[{AIMSSVars[eqns][[nonzeropos]],
					  SetPrecision[AIMSSGuess[eqns][[nonzeropos]]+AIMZeroTol,AIMPrecision]}],
						MaxIterations->5000, WorkingPrecision->.96*AIMPrecision]], AIMZeroTol] ;
 If[NumberQ[First[ss2]]==False, $Failed,
 Print["Precision of steady state is estimated to be ",Precision[Last[ss2]]," digits, sum of squared errors is ", First[ss2]] ;
 zeropos2 = Flatten[Position[Map[Last,Last[ss2]],0]] ;
 If[Length[zeropos2]==0, ReplacePart[Thread[AIMSSVars[eqns]->0],Thread[nonzeropos->Last[ss2]]],
   AIMRefineSS[eqns,sseqns2, Union[zeropos, Flatten[Map[Position[AIMSSVars[eqns],#]&,Map[First,Last[ss2][[zeropos2]]]]]]]
]]] ;


(* Front end that does error-checking and formats the output *)

AIMSeries[eqns_,parameters_List,deg_Integer]:=
Module[{soln,const,argsubs},
 If[AIMSS[eqns,parameters,AIMZeroTol,AIMPrecision] === $Failed ||
	(soln=AIMSoln[eqns,parameters,deg,AIMZeroTol,AIMPrecision]) === $Failed, $Failed,
 Print["Formatting output, time is ", Date[]] ;
 const = AIMSSVars[eqns] /.AIMSS[eqns,parameters,AIMZeroTol,AIMPrecision] ;
 argsubs = Thread[AIMGenericArgs[eqns]->AIMStateVars[eqns]] ;
 Thread[Through[AIMVarNames[eqns][t]] == const + (soln /.argsubs)]
]] ;

AIMSeries[eqns_,deg_Integer]:= AIMSeries[eqns,{},deg] ; (* if AIMSeries called without separate parameters, insert {} *)



(* Front end for Linear AIM *)

AIMSoln[eqns_,parameters_,1,zerotol_,precision_]:= AIMSoln[eqns,parameters,1,zerotol,precision] =
Module[{eqnseq0,allvars,hmat,epsmat,cofb,s0inv,alllagvars},
 eqnseq0 = Flatten[{eqns /.x_Real:>SetPrecision[x,precision]}] //.(parameters /.x_Real:>SetPrecision[x,precision]) /.Equal->Subtract ;
 allvars = Flatten[Map[Through[AIMVarNames[eqns] [t+#]]&, Range[-Max[AIMMaxLag[eqns],1],AIMMaxLead[eqns]]]] ;
 If[!ValueQ[AIMSS[eqns,parameters,zerotol,precision]], AIMSS[eqns,parameters]] ;
 Print["Linearizing model around steady state, time is ", Date[]] ;
 hmat = Chop[Outer[D,eqnseq0,allvars] /.AIMSSSubs /.AIMSS[eqns,parameters,zerotol,precision], zerotol] ;
 epsmat = Chop[Outer[D,eqnseq0,AIMShocks[eqns]] /.AIMSSSubs /.AIMSS[eqns,parameters,zerotol,precision], zerotol] ;
 If[({cofb,s0inv} = AIMLinearSoln[hmat,AIMMaxLead[eqns]]) === {{},{}}, $Failed,
 Print["Precision of linear solution is estimated to be ",Precision[{cofb,s0inv}]," digits, time is ", Date[]] ;
 alllagvars = allvars[[Range[Max[AIMMaxLag[eqns],1]*AIMNVars[eqns]]]] ;
 Chop[cofb .alllagvars - s0inv .epsmat .AIMShocks[eqns] /.Thread[AIMStateVars[eqns]->AIMGenericArgs[eqns]], zerotol]
]] ;


AIMSoln[eqns_,parameters_,0,zerotol_,precision_]:=
Table[0,{AIMNVars[eqns]}] ; (* this is unnecessary but allows the user to call AIMSeries[eqns,0] or AIMSoln[eqns,0] *)


(* This is the heart of the program.  Derivatives are evaluated at steady state, the expectation is taken,   *)
(*   and coefficients are solved using the method of undetermined coefficients.                              *)
(* As in the simple version, we solve a certainty-equivalent version of the problem first.  In addition,     *)
(*   for the clearly backward-looking variables, we delete terms and coefficients that are clearly zero,     *)
(*   and thus solve only those coefficients that are in general nonzero.                                     *)

AIMSoln[eqns_,parameters_,deg_Integer,zerotol_,precision_]:= AIMSoln[eqns,parameters,deg,zerotol,precision] =
Module[{args,blvpos,bfuncargs,drvindxs,cedrvindxs,stdrvindxs,cecoeffs,stcoeffs,nextceTerms,nextstTerms,
	cesubs,cesystem,cesoln,cematrix,dum,stsubs,stsystem,stsoln},
 args = AIMGenericArgs[eqns] ;
 blvpos = AIMBLVPos[eqns] ;
 bfuncargs = Table[args,{AIMNVars[eqns]}] ;
 bfuncargs[[blvpos]] = AIMBackLookingVars[eqns][[2]] /.Thread[AIMStateVars[eqns]->AIMGenericArgs[eqns]] ;
 drvindxs = Compositions[deg,AIMNArgs[eqns]] ;
  cedrvindxs = Select[drvindxs, #[[-1]]==0 &] ;
  stdrvindxs = Select[drvindxs, #[[-1]] >1 &] ;
 cecoeffs = Table[Unique[coeff,{Temporary}],{AIMNVars[eqns]},{Length[cedrvindxs]}] ;
  stcoeffs = Table[Unique[coeff,{Temporary}],{AIMNVars[eqns]},{Length[stdrvindxs]}] ;
 nextceTerms = cecoeffs .Map[Apply[Times,Power[args,#]]&, cedrvindxs] ;
  nextstTerms = stcoeffs .Map[Apply[Times,Power[args,#]]&, stdrvindxs] ;
 nextceTerms[[blvpos]] = Map[ nextceTerms[[#]] /.Thread[Complement[args,bfuncargs[[#]]]->0]&, blvpos] ;
  nextstTerms[[blvpos]] = 0 ;
 cecoeffs[[blvpos]] = Map[Complement[Cases[nextceTerms[[#]],_Symbol,Infinity], bfuncargs[[#]]]&, blvpos] ;
  stcoeffs[[blvpos]] = {} ;
 cecoeffs = Flatten[cecoeffs] ;
  stcoeffs = Flatten[stcoeffs] ;
 cesubs = Thread[Map[bFunc,AIMVarNames[eqns]] -> MapThread[Apply[Function,{#1,#2}]&,
				{bfuncargs, AIMSoln[eqns,parameters,deg-1,zerotol,precision] + nextceTerms}]] ;
 Print["Calculating CE Derivatives, time is ", Date[]] ;
 cesoln = If[AIMNArgs[eqns]===1, PrependTo[args,dum]; nextceTerms,
   cesystem = Chop[Flatten[CoefficientArrays[AIMCEDerivatives[eqns,deg][0]
	//.parameters /.AIMSS[eqns,parameters,zerotol,precision] /.cesubs, Drop[args,-2]]], zerotol] ;
   Print["Undetermined CE coefficients to solve: ",Length[cecoeffs]] ;
   Print["Calculating CE solution, time is ", Date[]] ;
   cematrix = CoefficientArrays[DeleteCases[cesystem,0], cecoeffs] ;
   nextceTerms /.Thread[cecoeffs -> Chop[LinearSolve[cematrix[[2]],-cematrix[[1]]],zerotol]]] ;
 Print["Precision of solution is estimated to be ",Precision[cesoln]," digits, time is ", Date[]] ;
 If[stcoeffs==={}, stsoln=0,
  stsubs = Thread[Map[bFunc,AIMVarNames[eqns]] -> MapThread[Apply[Function,{#1,#2}]&,
				{bfuncargs, AIMSoln[eqns,parameters,deg-1,zerotol,precision] + cesoln + nextstTerms}]] ;
  Print["Calculating Stoch Derivatives, time is ", Date[]] ;
  stsystem = Chop[Expand[Flatten[CoefficientArrays[Chop[
	AIMStDerivatives[eqns,deg][0] //.parameters /.AIMSS[eqns,parameters,zerotol,precision]
		/.stsubs, zerotol],Drop[args,-1]]]], zerotol] /.eps[x_][_]^n_->mom[x,n] /.eps[_][_]->0 ;
  Print["Undetermined Stoch coefficients to solve: ",Length[stcoeffs]] ;
  Print["Calculating Stoch solution, time is ", Date[]] ;
  stsoln = nextstTerms /.(Flatten[NSolve[stsystem, stcoeffs]] /.x_/;Abs[N[x]]<zerotol ->0) ;
 ] ;
 Apply[Remove,Join[cecoeffs,stcoeffs]] ;
 Print["Precision of solution is estimated to be ",Precision[stsoln]," digits, time is ", Date[]] ;
 AIMSoln[eqns,parameters,deg-1,zerotol,precision] + cesoln + stsoln
] ;

(* That's essentially it.  The following routine calculates derivatives of the equations composed with the   *)
(*   (unknown) solution functions b.	                                                                     *)
(* Like the simple version, we use univariate differentiation to compute the multivariate derivatives of     *)
(*   the function Fob.  Unlike the simple version, we compute the certainty-equivalent derivatives (those    *)
(*   not involving Sigma) first, and the stochastic derivatives (those involving Sigma) second.  This        *)
(*   allows us to avoid ever computing the first derivatives of Fob with respect to Sigma, which we know     *)
(*   from theory must be zero.                                                                               *)

AIMCEDerivatives[eqns_,2]:= AIMCEDerivatives[eqns,2] =
Derivative[2][Function[tAIM, Evaluate[AIMSubBFuncsIntoEqns[eqns] /.
 Thread[AIMStateVars[eqns]->tAIM*Join[Drop[AIMGenericArgs[eqns],-2],{1,0}]]]]] ;

AIMCEDerivatives[eqns_,deg_Integer]:= AIMCEDerivatives[eqns,deg] = AIMCEDerivatives[eqns,deg-1]' ;


AIMStDerivatives[eqns_,2]:= AIMStDerivatives[eqns,2] =
Function[tAIM, Evaluate[Derivative[0,2] [Function[{tAIM1,tAIM2},
	Evaluate[AIMSubBFuncsIntoEqns[eqns] /.Thread[AIMStateVars[eqns]->
	Append[tAIM1*Drop[AIMGenericArgs[eqns],-1],tAIM2]]]]] [tAIM,tAIM]]] ;

AIMStDerivatives[eqns_,deg_Integer]:= AIMStDerivatives[eqns,deg] = AIMStDerivatives[eqns,deg-1]' ;


AIMSubBFuncsIntoEqns[eqns_]:= AIMSubBFuncsIntoEqns[eqns] =
With[{deveqns = eqns /.x_Symbol[t+i_.]:>Symbol[SymbolName[x]<>"AIMSS"] +x[t+i] /.Equal->Subtract},
 AIMSubOutLeadVars[eqns,deveqns,AIMMaxLead[eqns]]
] ;


AIMBFuncs[eqns_]:= AIMBFuncs[eqns] =
Module[{bfuncs},
 bfuncs = Map[Apply[bFunc[#],AIMStateVars[eqns]]&, AIMVarNames[eqns]] ;
 bfuncs[[AIMBLVPos[eqns]]] = MapThread[Apply[bFunc[#1],#2]&, AIMBackLookingVars[eqns]] ;
 bfuncs
] ;


AIMSubOutLeadVars[origeqns_,eqnssofar_,0]:= 
eqnssofar /.Thread[Through[AIMVarNames[origeqns][t]]->AIMBFuncs[origeqns]] /.eps[x_][t+i_]->Sigma*eps[x][t+i] ;

AIMSubOutLeadVars[origeqns_,eqnssofar_,lead_Integer]:=
With[{reducedeqns=eqnssofar /.Thread[Through[AIMVarNames[origeqns][t+lead]]-> (AIMBFuncs[origeqns]/.t->t+lead)]},
 AIMSubOutLeadVars[origeqns, reducedeqns, lead-1]
] ;


(* Print out model diagnostics (obviously) *)

AIMModelDiagnostics[eqns_] := (
Print["Model Diagnostics:"] ;
Print["Number of equations: ", AIMNEqns[eqns]] ;
Print["Number of variables: ", AIMNVars[eqns]] ;
Print["Number of shocks:    ", Length[AIMShocks[eqns]]] ;
Print["Maximum lag:         ", AIMMaxLag[eqns]] ;
Print["Maximum lead:        ", AIMMaxLead[eqns]] ;
Print["Lagged variables: ", AIMLagVars[eqns]] ;
Print["Shocks: ", AIMShocks[eqns]] ;
Print[" together with Sigma, these yield ", AIMNArgs[eqns], " state variables"] ;
Print["Variables ", AIMBackLookingVars[eqns][[1]], " are clearly backward-looking"] ;
Print[" using the following reduced state variable sets for these variables:"] ;
Print[AIMBackLookingVars[eqns][[2]]] ;
Print["List of all variables: ", AIMVarNames[eqns]] ;
Print["Using internal working precison of ", AIMPrecision, " digits (AIMPrecision)"] ;
Print["Treating numbers < ", N[AIMZeroTol], " as zero (AIMZeroTol)"] ;
) ;


(* Everything that follows is Linear AIM.  See Anderson and Moore (1985) for a description of the algorithm.  *)
(* In contrast to the simple version of the code, here we perform "exact" ShiftRights first (i.e., shift      *)
(*   equations forward that have a complete row of zeros in the lead matrix, which reduces the number of      *)
(*   calls to the QR decomposition and thus reduces potential numerical problems.)                            *)
(* Also, we compute the stability conditions of the model using the Schur decomposition if the eigenvectors   *)
(*   corresponding to the large eigenvalues of the model are linearly dependent.  The Schur decomposition     *)
(*   is guaranteed to exist for all models, while linearly independent eigenvectors are not.  The             *)
(*   disadvantage of using Schur vectors is that computing them is slow: first, the computation is            *)
(*   inherently slower, and second, the fact that the code must swap Schur blocks using a routine that I      *)
(*   (ets) wrote, rather than internal optimized machine code, slows it down quite a bit.                     *)

AIMLinearSoln[hmat_,nleads_Integer]:= AIMLinearSoln[hmat,nleads] =
Module[{(*hrows,hcols,*)shiftedh(*,qmat,stabconds,bmat,smat*)},
 {hrows,hcols} = Dimensions[hmat] ;
 {shiftedh,qmat} = FixedPoint[AIMExactShiftRight, {hmat,{}}, nleads] ;
 {shiftedh,qmat} = FixedPoint[AIMShiftRight, {shiftedh,qmat}, nleads*hrows] ;
 Print["Lead matrix is full rank; computing stability conditions, time is ", Date[]] ;
 stabconds = AIMLinearStabConds[Transpose[AIMAR1Form[shiftedh]]] ;
 (* CAN TURN ON FOR MORE INFO ON EIGENVECTORS Print["Stability conditions, somehow related to the transpose of the AIMAR1Form's eigenvectors with eigenvectors within a numeric threshold rounded to 0", stabconds]; *)
 Print["Model has ",Length[stabconds]," unstable roots, time is ", Date[]] ;
 qmat = Join[qmat,stabconds] ;
(* CAN TURN ON FOR MORE INFO ON EIGENVECTORS 
Print["qmat:", qmat];
Print["Eigenvalues of AIMAR1Form[shiftedh]'", Eigenvalues[Transpose[AIMAR1Form[shiftedh]]]];
Print["Eigenvectors of AIMAR1Form[shiftedh]'", Eigenvectors[Transpose[AIMAR1Form[shiftedh]]]];
*)
 If[Length[qmat]<hrows*nleads, Print["Multiple Linear Solutions"]; {{},{}},
 If[Length[qmat]>hrows*nleads, Print["No Linear Solutions"]; {{},{}},
 bmat = LinearSolve[AIMLastBlock[qmat],-AIMFirstBlocks[qmat]] ;
 smat = hmat .Join[IdentityMatrix[hcols-hrows*nleads], MapThread[Join,{Table[0,{hrows*nleads},{hrows}],bmat}]] ;
 Chop[{Take[bmat,hrows], Inverse[AIMLastBlock[smat]]}, AIMZeroTol]
]]] ;


AIMLinearSoln[hmat_,0] :=
With[{cofb=LinearSolve[AIMLastBlock[hmat],-AIMFirstBlocks[hmat]]}, Chop[{cofb,Inverse[AIMLastBlock[hmat]]}, AIMZeroTol]] ;


AIMExactShiftRight[{hmatold_,qmatsofar_}]:=
Module[{hrows=Length[hmatold],zerorows,firstblocks},
If[ Length[zerorows = Position[Map[Max,Abs[AIMLastBlock[hmatold]]],0]] ==0, {hmatold,qmatsofar},
 firstblocks = AIMFirstBlocks[hmatold] ;
 Print["Shifting ",Length[zerorows]," equations ",zerorows," forward"] ;
 {ReplacePart[hmatold, ArrayFlatten[{{Table[0,{hrows},{hrows}],firstblocks}}], zerorows,zerorows],
   Join[qmatsofar,firstblocks[[Flatten[zerorows]]]]}
]] ;


AIMShiftRight[{hmatold_,qmatsofar_}]:=
Module[{hrows=Length[hmatold],q,r,zerorows,hmatnew,firstblocks},
 {q,r} = Take[QRDecomposition[AIMLastBlock[hmatold], Pivoting->True], 2] ;
If[ Length[zerorows = Position[Abs[Tr[r,List]], x_/;x<AIMZeroTol]] ==0, {hmatold,qmatsofar},
 Print["Shifting ",Length[zerorows]," linear combinations of equations forward"] ;
 firstblocks = AIMFirstBlocks[hmatnew = q .hmatold] ;
 Chop[{ReplacePart[hmatnew,ArrayFlatten[{{Table[0,{hrows},{hrows}],firstblocks}}], zerorows,zerorows],
   Join[qmatsofar,firstblocks[[Flatten[zerorows]]]]}, AIMZeroTol]
]] ;


AIMAR1Form[shiftedh_]:=
With[{hrows=Length[shiftedh],hcols=Length[First[shiftedh]]},
 Chop[Join[MapThread[Join,{Table[0,{hcols-2*hrows},{hrows}],IdentityMatrix[hcols-2*hrows]}],
	LinearSolve[AIMLastBlock[shiftedh],-AIMFirstBlocks[shiftedh]]], AIMZeroTol]
] ;


AIMLastBlock[matrix_]:= With[{dims=Dimensions[matrix]}, Take[matrix,{1,dims[[1]]},{dims[[2]]-dims[[1]]+1,dims[[2]]}]
] ;


AIMFirstBlocks[matrix_]:= With[{dims=Dimensions[matrix]}, Take[matrix,{1,dims[[1]]},{1,dims[[2]]-dims[[1]]}]
] ;


(* Mathematica's implementation of the Schur decomposition does not allow the user to specify any particular *)
(*   order for the eigenvalues to appear along the diagonal.  To compute the stability conditions of the     *)
(*   model, we require all eigenvalues greater than one in absolute value to appear in the top part of the   *)
(*   diagonal.  The following two routines swap adjacent Schur blocks until this ordering is achieved.       *)
(* The routine is based on a 2002 paper by Jan Brandts entitled "Matlab Code for Sorted Real Schur Forms,"   *)
(*   available at http://staff.science.uva.nl/~brandts.  One advantage of this particular algorithm is that  *)
(*   it allows for swapping 1x1 and 2x2 Schur blocks, so that we can avoid complex numbers entirely.  This   *)
(*   keeps annoying (albeit tiny) imaginary numbers from cropping up in our stability conditions and hence   *)
(*   our solution.                                                                                           *)
(* One can also use a Generalized Schur Decomposition instead of the plain vanilla Schur Decomposition, but  *)
(*   in Linear AIM the gains from this are small because of the efficient way that AIM handles leads (using  *)
(*   ShiftRights to render the lead matrix invertible, and avoiding stacking equations into the more         *)
(*   numerically unstable AR1 form until the very last minute).                                              *)

AIMLinearStabConds[matrix_]:=
With[{n=Length[matrix]}, Module[{stabconds,q,t,subdiagelts,twoblockpos,blockpos,eigvals,bigrtpos,swaplist},
 stabconds = Chop[Eigenvectors[matrix, Length[Cases[Abs[Eigenvalues[matrix]],i_/;i>1+AIMZeroTol]]], AIMZeroTol] ;
 If[Min[Map[Max,Abs[stabconds]]]>0, stabconds,
  Print["Found degenerate Blanchard-Kahn eigenvector, using Schur decomposition instead, time is ", Date[]] ;
  {q,t} = Chop[SchurDecomposition[matrix], AIMZeroTol] ;
  subdiagelts = Map[t[[#+1,#]]&, Range[n-1]] ;
  twoblockpos = Flatten[Position[subdiagelts, i_/;i!=0]] ;
  blockpos = Complement[Range[n+1], twoblockpos+1] ;
  eigvals = Map[t[[#,#]]&, Drop[blockpos,-1]] ;
  eigvals[[Flatten[Map[Position[blockpos,#]&, twoblockpos]]]] += Map[Sqrt[t[[#,#+1]] *t[[#+1,#]]]&, twoblockpos] ;
  bigrtpos = Flatten[Position[Abs[eigvals], i_/;i>1+AIMZeroTol]] ;
  swaplist = Flatten[Map[Range[bigrtpos[[#]]-1,#,-1]&, Range[Length[bigrtpos]]]];
  {q,t,blockpos} = Fold[AIMSwapSchurBlocks, {Transpose[q],t,blockpos}, swaplist];
  stabconds = Take[q, blockpos[[Length[bigrtpos]+1]] -1]
]]] ;


AIMSwapSchurBlocks[{q_,t_,blockpos_},swappos_]:=
Module[{index1,index2,n1,n2,a11,a12,a22,xmat,ans,xq,index12,qswpd,tprime,tswpd,blockposswpd},
 index1 = Range[blockpos[[swappos]], blockpos[[swappos+1]]-1] ;
 index2 = Range[blockpos[[swappos+1]], blockpos[[swappos+2]]-1] ;
 {n1,n2} = {Length[index1],Length[index2]} ;
 {a11,a12,a22} = {Take[t,index1,index1], Take[t,index1,index2], Take[t,index2,index2]} ;
 xmat = Table[Unique[x,{Temporary}],{n1},{n2}] ;
 ans = Flatten[NSolve[Flatten[a11.xmat-xmat.a22]-Flatten[a12], Flatten[xmat]]];
 xq = First[QRDecomposition[ ArrayFlatten[{{-xmat /.ans,IdentityMatrix[n1]},{IdentityMatrix[n2],Table[0,{n2},{n1}]}}]]] ;
 index12 = Join[index1,index2] ;
 qswpd = q ; (* premultiply q by xq *)
  qswpd[[index12]] = xq .qswpd[[index12]] ;
 tprime = Transpose[t] ; (* transform t to xq .t .xq' *)
  tprime[[index12]] = xq .tprime[[index12]] ;
  tswpd = Transpose[tprime] ;
  tswpd[[index12]] = xq .tswpd[[index12]] ;
 blockposswpd = ReplacePart[blockpos, blockpos[[swappos+1]] +n2-n1, swappos+1] ;
 Apply[Remove,Flatten[xmat]] ;
 Chop[{qswpd,tswpd,blockposswpd}, AIMZeroTol]
] ;

(* JP ADDITION: This evaluates the solution at a set of substitutions, returning the new state.*)
EvaluteSolution[soln_, subs_]:= Table[Part[soln /. subs, i, 2], {i,  1, AIMNVars[eqns]}];
