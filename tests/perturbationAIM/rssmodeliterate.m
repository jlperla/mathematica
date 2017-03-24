

(*
 * $Rev: 1 $
 * Original code: http://www.ericswanson.us/perturbation.html
 * Modified by: Jesse Perla
*)

(* Load perturbationAIM *)
<<perturbationAIM.m


(* builds on the nkmodellogint.m from the RSS St. Louis Fed conference paper, with the
   following modifications:
     1. put interest rate rule in terms of output gap rather than output growth
     2. modify the consol to give it a duration of 10 yrs instead of 25 yrs
     3. adjust parameter values to get steady-state ratios roughly calibrated to US data
     4. introduce an adjustment cost parameter psi that we'll set to 0 here but play
	around with later.

Eric Swanson, 8/2007.
*)

(* objective = (C[t] - habitsize *(C[t-1]))^(1-phi) /(1-phi) - chi0 *L[t]^(1+chi) /(1+chi)
 *)

eqns = {
(* Marginal Utility of Consumption and Consumer's Euler Equation *)
MUc[t] == (C[t] - habitsize *C[t-1])^-phi,
MUc[t] == beta *(Exp[Int[t]]/pi[t+1]) *MUc[t+1],

(* Price-setting equations *)
zn[t] == (1+theta) *MUc[t] *MC[t] *Y[t] + beta *xi *pi[t+1]^((1+theta)/theta/(1-alpha)) *zn[t+1],
zd[t] == Y[t] *MUc[t] + beta *xi *pi[t+1]^(1/theta) *zd[t+1],
p0[t]^(1+(1+theta)/theta *alpha/(1-alpha)) == zn[t] /zd[t],
pi[t]^(-1/theta) == (1-xi) *(p0[t]*pi[t])^(-1/theta) + xi,

(* Marginal cost and quadratic adjustment costs to labor, psi/2 *Log[L[t]/L[t-1]]^2 *)
MC[t] == wreal[t] /(1-alpha) *Y[t]^(alpha/(1-alpha)) /A[t]^(1/(1-alpha)) /KBar^(alpha/(1-alpha)),
chi0 *L[t]^chi /MUc[t] == wreal[t] - psi /L[t] *Log[L[t]/L[t-1]]
						+ beta *psi /L[t] *Log[L[t+1]/L[t]] *MUc[t+1]/MUc[t],

(* Output equations *)
Y[t] == Disp[t]^(-1) *A[t] *KBar^alpha *L[t]^(1-alpha),
Disp[t]^(1/(1-alpha)) == (1-xi) *p0[t]^(-(1+theta)/theta/(1-alpha))
	+ xi *pi[t]^((1+theta)/theta/(1-alpha)) *Disp[t-1]^(1/(1-alpha)),
C[t] == Y[t] - G[t] - IBar - psi/2 *Log[L[t]/L[t-1]]^2, (* aggregate resource constraint w/ adj costs *)

(* Monetary Policy Rule *)
piavg[t] == rhoinflavg *piavg[t-1] + (1-rhoinflavg) *pi[t],
(*piavg[t] == (pi[t] + pi[t-1] + pi[t-2] + pi[t-3]) /4,*)
4*Int[t] == (1-taylrho) * ( 4*Log[1/beta] + 4*Log[piavg[t]]
	+ taylpi * (4*Log[piavg[t]] - piBar)
	+ tayly*(Y[t]-YBar)/YBar )
	+ taylrho * 4*Int[t-1] + eps[Int][t],  (* multiply Int, infl by 4 to put at annual rate *)
Intr[t] == Log[Exp[Int[t-1]]/pi[t]], (* ex post real short rate *)

(* Exogenous Shocks *)
Log[A[t]/ABar] == rhoa * Log[A[t-1]/ABar] + eps[A][t],
Log[G[t]/GBar] == rhog * Log[G[t-1]/GBar] + eps[G][t],
pistar[t] == (1-rhopistar) *piBar + rhopistar *pistar[t-1] + 0*eps[pistar][t],

(* Long-Term Bond Price, Yield, and Term Premium *)
pricebond[t] == 1/100 + beta *consoldelta *MUc[t+1] /MUc[t] /pi[t+1] *pricebond[t+1],
pricebondrn[t] == 1/100 + consoldelta *pricebondrn[t+1] /Exp[Int[t]], (* risk-neutral bond price *)
ytm[t] == Log[consoldelta*pricebond[t]/(pricebond[t]-1/100)] *400,  (* yield in annualized percent *)
ytmrn[t] == Log[consoldelta*pricebondrn[t]/(pricebondrn[t]-1/100)] *400,  (* yield in annualized percent *)
termprem[t] == 100 * (ytm[t] - ytmrn[t]) (* term premium in annualized basis points *)
}


loglinearizevars = {A, C, Disp, G, L, MUc, p0, pi, piavg, wreal, Y, zd, zn} ;
logRules = Map[#[x_]->E^(#[x])&, loglinearizevars] ;

(* List of all variables:
{A, C, Disp, G, Int, Intr, L, MC, MUc, p0,
  pi, piavg, pistar, pricebond, pricebondrn, termprem, wreal, Y, ytm, ytmrn, zd, zn}
*)

SSGuess = {0, .5, 0, -.8, .01, .01, 0, .8, 1.1, 0,
  0, 0, 0, .4, .4, 0, .5, 1, 4, 4, 3.5, 3.5}

AIMPrecision = 50 ;
AIMZeroTol = 10^-10 (* 10^-50 *) ;  (* for psi = 1000, AIMZeroTol = 10^-8 and Precision 80 seems to work *)

parametersfixed = SetPrecision[{alpha->.3, beta->.99, delta->.02, chi0->Exp[wrealAIMSS +MUcAIMSS], psi->0,
  theta->.2, ABar->1, GBar ->.17 *YBar, KBar->10 *YBar, IBar->delta *KBar, YBar->Exp[YAIMSS], piBar->0,
  taylrho->.73, taylpi->.53, tayly->.93, rhog->.9, rhopistar->0, rhoinflavg->.7,
  consoldelta -> Exp[ytmAIMSS/400] *(1-1/40)}, AIMPrecision] ;

rssmodel = eqns /.logRules //. parametersfixed ;

parameterlist = SetPrecision[Flatten[Outer[{phi->#1, habitsize->#2, chi->#3, rhoa->#4, xi->#5}&,
  {5,5.5}, {.75,.8}, {2.25,2.5}, {.9,.95}, {.5,.55}], 4],  (* note: this set has 32 parameter vectors *)
(* {5}, {.75}, {2.25,2.5}, {.9}, {.5}], 4], here's a much smaller set of parameter vectors (only 2) *)
  AIMPrecision] ; (* list of all parameter value combinations to search over *)

(* Note: it takes several minutes to solve the model for each parameter vector above, so if you want to just
   test the file, you should reduce the number of parameter values considered in the parameterlist, e.g., uncomment
   the smaller list of 2 vectors instead of using the full 32. *)
   
momSubs = {Sigma->1, mom[_,3]->0, mom[G,2]->.004^2, mom[Int,2]->.004^2} ;
momalist = {.005,.0075,.01,.015,.02,.025}^2 ; (* all possible values for mom[A,2] to search over; note that
	   				      	 it's much faster to search over mom[_,_] values than other parameters. *)

(* print out the steady-state term premium estimate: *)

tppos = Position[AIMVarNames[rssmodel], termprem][[1,1]] ;

parameterlistorig = parameterlist ;

AIMSSGuess[rssmodel] = SSGuess ; (* initialize steady state guess *)

(* This routine loops through all of the parameter value combinations in the Outer[_] function above *)
momentlist = Flatten[Map[Function[param, moments = Map[Function[moma,
  (Print[param]; Print["mom[A,2] -> ", moma] ;
  ss = AIMSS[rssmodel,param,AIMZeroTol,AIMPrecision] ;
  If[ss===$Failed, ss = AIMSS[rssmodel,param,AIMSSGuess->SSGuess]] ; (* try a second SSguess if needed *)
  If[ss===$Failed, 10^6*{1,1,1,1,1,1,1,1,1}, (* if still fail, then record bad fit, else continue *)
  eqnseq0 = rssmodel //.SetPrecision[param,AIMPrecision] /.Equal->Subtract ;
  allvars = Flatten[Map[Through[AIMVarNames[rssmodel] [t+#]]&,
		      Range[-Max[AIMMaxLag[rssmodel],1], AIMMaxLead[rssmodel]]]] ;
  hmat = Outer[D,eqnseq0,allvars] /.AIMSSSubs /.ss ;
  epsmat = -Outer[D,eqnseq0,AIMShocks[rssmodel]] /.AIMSSSubs /.ss ;
  {cofb,s0inv} = AIMLinearSoln[hmat,AIMMaxLead[rssmodel]] ;
  If[cofb==={}, 10^6*{1,1,1,1,1,1,1,1,1}, (* else continue *)
  s0invsmall = s0inv .epsmat ;
  Print[tp = AIMSeries[rssmodel,param,2][[tppos]] /.momSubs /.mom[A,2]->moma /.x_Real:> N[x]] ;
  sigmaeps = s0invsmall .
           DiagonalMatrix[AIMShocks[rssmodel] /.eps[x_][_]->mom[x,2] /.momSubs /.mom[A,2]->moma] .
           Transpose[s0invsmall] ;
  bigb = Chop[N[cofb]] ;
  sigmax = FixedPoint[(s = # + bigb.#.Transpose[bigb]; bigb = bigb.bigb; s)&, sigmaeps, 8] ;
  moments = Sqrt[Map[sigmax[[#,#]]&,Range[AIMNEqns[rssmodel]]]] [[
	Flatten[Map[Position[AIMVarNames[rssmodel], #] &, {C, Y, L, wreal, pi, Int, Intr, ytm}]] ]
	] * {100, 100, 100, 100, 400, 400, 400, 1} ;
  Join[moments, {Chop[Last[tp],10^-4]}]]]
) ], momalist] ;
  Unset[AIMSoln[rssmodel,param,1,AIMZeroTol,AIMPrecision]] ; (* free up memory *)
  Unset[AIMSoln[rssmodel,param,2,AIMZeroTol,AIMPrecision]] ; (* free up memory *)
  Unset[AIMLinearSoln[hmat,1]] ; (* free up memory *)
  Unset[AIMSS[rssmodel,param,AIMZeroTol,AIMPrecision]] ; (* free up memory *)
  moments
], parameterlist] ,1] ;

momenttargets = {1.19, 1.50, 1.71, .82, 2.52, 2.71, 2.30, 2.37, 106} ; (* the moments to match *)
momentweights = {1, 0, 1, 1, 1, 1, 1, 1, .01^2} ; (* diagonal weighting matrix with equal weights, account for units *)
distancelist = Map[((#-momenttargets) *momentweights) .(#-momenttargets) &, momentlist] ; (* distance from target *)

mindist = Min[distancelist] ;
bestfitpos = Position[distancelist,mindist][[1,1]] ;
bestfitpos1 = Mod[Flatten[Position[distancelist,mindist]] -1, Length[momalist]] +1 ;
bestfitpos2 = Floor[Flatten[Position[distancelist,mindist]]/Length[momalist]] +1 ;
bestparameters = N[Flatten[parameterlist[[bestfitpos2]]]] ; (* best-fitting set of parameters *)
bestmoma = momalist[[bestfitpos1]] ; (* best-fitting set of moments *)

Print["\n"] ;
Print["best parameters: ", bestparameters] ;
Print["best mom[A,2]: ", bestmoma] ;
Print["distance to target moments: ", mindist] ;
Print["memory used: ", MaxMemoryUsed[]] ;
Print["\n"] ;
Print["best-fitting moments: ", momentlist[[bestfitpos]]] ;

