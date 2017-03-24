(* ::Package:: *)

(*
 * $Rev: 1 $
 * Original code: http://www.ericswanson.us/perturbation.html
 * Modified by: Jesse Perla
*)

(* Load perturbationAIM *)
<<perturbationAIM.m


eqns={
  Y[t] == A[t] *K[t-1]^alpha,
  Log[A[t]] == rho *Log[A[t-1]] + eps[a][t],
  K[t] == (1-delta) *K[t-1] + Inv[t], (* Mathematica treats "I" as Sqrt[-1] *)
  Y[t] == C[t] + Inv[t],
  C[t]^-gamma == beta *(1+r[t+1]) *C[t+1]^-gamma,
  r[t] == alpha *A[t] *K[t-1]^(alpha-1) - delta,
  Welf[t] == C[t]^(1-gamma) /(1-gamma) + beta *Welf[t+1]
}

(* variables listed here are transformed to logs and then approximated in
  logs rather than in levels (analogous to log-linearizing a linear model):
*)
logvars = {A, C, K, Y}
logrules = Map[#[x_]->E^(#[x])&, logvars]

(* parameter values *)
parametervals={
  alpha->0.3,
  beta->0.99,
  gamma->1.1,
  delta->0.1,
  rho->0.8
}

(* substitute log transformation rules and parameter values into equations: *)
sgmodel = eqns /.logrules //.parametervals


(* complete variable list (Mathematica places in alphabetical order):
{A, C, Inv, K, r, Welf, Y}
*)

(* find the steady state *)
AIMSS[sgmodel]

(* Find the third-order approximation: *)
AIMSeries[sgmodel,3]
