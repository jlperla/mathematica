(* ::Package:: *)

(*
 * $Rev: 1 $
 * Original code: http://www.ericswanson.us/perturbation.html
 * Modified by: Jesse Perla
*)

(* Load perturbationAIM *)
<<perturbationAIM.m


(* Builds on the model from Rudebusch and Swanson (2008, FRB St. Louis Economic Review) and FRBSL conference, with the
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
};


parametervalues = {
  alpha -> .3,
  beta -> .99,
  delta -> .02,
  phi -> 2,
  habitsize -> .66,
  chi0 -> (* 4.73921522897550247856008201926876616780052304558537869572838202224763 *)
    Exp[wrealAIMSS +MUcAIMSS], (* normalize L in the nonstochastic steady state to be 1; note the use of *AIMSS
				variables here to accomplish this.  The AIMSS suffix denotes the nonstochastic
				steady-state value of a variable, as defined by the AIMSS[_] function call. *)
  chi -> 1.5,

  psi -> 0 *YBar, (* quadratic adjustment costs to labor *)
  theta -> .2,
  xi -> .75,

  ABar -> 1,
  GBar -> .17 *YBar,
  KBar -> 10 *YBar,
  IBar -> delta *KBar,
  YBar -> Exp[YAIMSS], (* nonstochastic steady-state level of Y *)
  piBar -> 0,

  taylrho -> .73,
  taylpi -> .53,
  tayly -> .93,

  rhoa -> .9,
  rhog -> .9,
  rhopistar -> 0,
  rhoinflavg -> .7,

  consoldelta -> Exp[ytmAIMSS/400] *(1-1/40)  (* consol depreciation, set to make duration of consol =10 yrs *) 
}


loglinearizevars = {A, C, Disp, G, L, MUc, p0, pi, piavg, pricebond, pricebondrn, wreal, Y, zd, zn} ;
logRules = Map[#[x_]->E^(#[x])&, loglinearizevars] ;

(* List of all variables:
{A, C, Disp, G, Int, Intr, L, MC, MUc, p0,
  pi, piavg, pistar, pricebond, pricebondrn, termprem, wreal, Y, ytm, ytmrn, zd, zn}
*)

SSGuess = {0, .5, 0, -.8, .01, .01, 0, .8, 1.1, 0,
  0, 0, 0, -1, -1, 0, .5, 1, 4, 4, 3.5, 3.5};

AIMPrecision = 50 ; (* 250 *)
AIMZeroTol = 10^-10 (* 10^-50 *) ;  (* for psi = 1000, AIMZeroTol = 10^-8 and Precision 80 seems to work *)

rssmodel = eqns /.logRules //.SetPrecision[parametervalues,AIMPrecision] ;
(* Just to simplify output
Print[rssmodel];
Print[SSGuess];
*)
ss = AIMSS[rssmodel,{},AIMSSGuess->SSGuess] ;

momSubs = {Sigma->1, mom[_,3]->0, mom[A,2]->.01^2, mom[G,2]->.004^2, mom[Int,2]->.004^2, mom[pistar,2]->0^2}


(* print out the steady-state term premium estimate: *)

tppos = Position[AIMVarNames[rssmodel], termprem][[1,1]] ;

Print[AIMSeries[rssmodel,2][[tppos]] /.momSubs /.x_Real:> N[x]] ;



(* Simulate solution forward nperiods periods *)

Print["\nThe Part[] function will issue some warnings here; fear not, everything is ok..."]

nperiods = 1000 ; (* number of periods *)
simdegree = 1 ;  (* degree of approximation *)

reportvars = {C, Y, L, wreal, pi, Int, Intr, ytm, termprem} ; (* variables of interest *)
simvars = Join[ Map[Head, AIMLagVars[rssmodel]], reportvars, {pricebond, ytmrn}] ; (* variables to simulate *)

lagvarpos = Range[Length[AIMLagVars[rssmodel]]] ;
reportvarpos = Length[AIMLagVars[rssmodel]] + Range[Length[reportvars]] ;

shocks = RandomReal[NormalDistribution[0,1],{nperiods,Length[AIMShocks[rssmodel]]}] .
					DiagonalMatrix[AIMShocks[rssmodel] /.eps[x_][_]->Sqrt[mom[x,2]] /.momSubs] ;

soln = AIMSoln[rssmodel,{},simdegree,AIMZeroTol,AIMPrecision][[Flatten[Map[Position[AIMVarNames[rssmodel], #] &,
													simvars]]]] ;

soln = Chop[N[soln /.Last[AIMGenericArgs[rssmodel]]->1 /. momSubs]] /.
	Thread[Drop[AIMGenericArgs[rssmodel],-Length[AIMShocks[rssmodel]]-1] -> Map[inputvar[[#]]&,lagvarpos]] /.
	Thread[Take[AIMGenericArgs[rssmodel],{-Length[AIMShocks[rssmodel]]-1,-2}] -> 
							Map[inputshock[[#]] &, Range[Length[AIMShocks[rssmodel]]]]] ;

iterator = Apply[Function,{ {inputvar,inputshock}, soln} ] ;


syndata = Transpose[FoldList[iterator, Table[0,{Length[simvars]}], shocks]] ;

synpricebond = Flatten[syndata[[Flatten[Position[simvars, pricebond]]]]] + (pricebondAIMSS /.ss) ;
synytm = Flatten[syndata[[Flatten[Position[simvars, ytm]]]]] + (ytmAIMSS /.ss) ;
synytmrn = Flatten[syndata[[Flatten[Position[simvars, ytmrn]]]]] + (ytmAIMSS /.ss) ;
synInt = Flatten[syndata[[Last[Position[simvars,Int]]]]] + (IntAIMSS /.ss) ;
syntp = Flatten[syndata[[Flatten[Position[simvars, termprem]]]]] + (termpremAIMSS /.ss) ;

meantp = Mean[syntp] ; (* mean term prem in annualized basis points *)

stddevs = Map[StandardDeviation,syndata[[reportvarpos]]] * {100, 100, 100, 100, 400, 400, 400, 1, 1}
		(* real variable std devs in percentage points, not logs; inflation std devs in
			annualized percentage points; interest rate std devs in annualized basis points *)

meanslope = (Mean[synytm] - 400*Mean[synInt]) * 100 ; (* slope in annualized basis points *)
meanslopern = (Mean[synytmrn] - 400*Mean[synInt]) * 100 ; (* slope in annualized basis points *)
stdslope = StandardDeviation[synytm - 400*synInt] * 100 ;

synehpr = (((consoldelta /.parametervalues/.ss) *Drop[Exp[synpricebond],1] + Drop[Exp[synInt],-1] *1/100) /
		Drop[Exp[synpricebond],-1] - Drop[Exp[synInt],-1]) * 40000 ; (* ehpr in annualized basis points *)
meanehpr = Mean[synehpr] ;
stdehpr = StandardDeviation[synehpr] ;

synDpricebond = Drop[synpricebond,1] - Drop[synpricebond,-1] ;
cscoeff = Covariance[-synDpricebond,Drop[synytm/400 - synInt,-1]] / Variance[synytm/400 - synInt] ;



Print["simulated mean tp = ", meantp]
Print["simulated mean slope = ", meanslope]
Print["simulated mean ehpr = ", meanehpr]
Print["simulated std dev tp = ", stdtp]
Print["simulated std dev slope = ", stdslope]
Print["simulated std dev ehpr = ", stdehpr]
Print["Solution Complete"]
