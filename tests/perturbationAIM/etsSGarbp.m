
(*
 * $Rev: 1 $
 * Original code: http://www.ericswanson.us/perturbation.html
 * Modified by: Jesse Perla
*)

(* Load perturbationAIM *)
<<perturbationAIM.m


(* Equations of the model *)
eqns={
  Y[t] == A[t] *K[t-1]^alpha,
  Log[A[t]] == rho *Log[A[t-1]] + eps[a][t],
  K[t] == (1-delta) *K[t-1] + Inv[t], (* Mathematica treats "I" as Sqrt[-1] *)
  Y[t] == C[t] + Inv[t],
  C[t]^-gamma == beta *(1+r[t+1]) *C[t+1]^-gamma,
  r[t] == alpha *A[t] *K[t-1]^(alpha-1) - delta,
  Welf[t] == C[t]^(1-gamma) /(1-gamma) + beta *Welf[t+1]
}

(* variables to transform to logs and then approximate *)
logvars = {A, C, K, Y}
logrules = Map[#[x_]->E^(#[x])&, logvars]

(* parameter values *)
parametervals={
  alpha->0.3`70,   (* 70 digits of precision: in general, you need to       *)
  beta->0.99`70,   (*  input more precision than you desire for the output. *)
  gamma->1.1`70,   (* Also, the larger the model and the higher the order   *)
  delta->0.1`70,   (*  of approximation, the more precision you must input  *)
  rho->0.8`70      (*  to achieve a given desired precision of output.      *)
}

(* substitute log transformation rules and parameter values into equations: *)
sgmodel = eqns /.logrules //.parametervals


(* complete variable list (Mathematica places in alphabetical order):
{A, C, Inv, K, r, Welf, Y}
*)

(* Now find the steady state to about 50 significant digits. *)

(* To truly get 50 significant digits of precision, one should set the
  AIMZeroTol parameter to 10^-50.  I have not done this here because I
  want to avoid giving users the misimpression that this is what they
  should do for their models in general.  As discussed in the comments to
  perturbationAIMsimple.m, the AIMZeroTol parameter should be thought of
  essentially as an "economic zero" rather than a "numerical zero" (the
  latter is left to the internals of Mathematica to handle appropriately).
  Since we're using arbitrary-precision numbers here, we have more freedom
  to set AIMZeroTol conservatively, but values much smaller than 10^-20 are
  probably not very meaningful economically.
*)

AIMZeroTol=10^-20 ;


(* The AIMPrecision option essentially determines the precision of the
  output of AIMSS; it must be less than the precision of the inputs or you
  will get a warning.  You must calculate the steady state to greater
  precision than you desire for the output of AIMSeries.
*)

AIMSS[sgmodel, AIMPrecision->65] ;


(* Find first- and second-order coefficients to 50 significant digits *)

Print[AIMSeries[sgmodel,1] //TableForm]

AIMSeries[sgmodel,2] //TableForm


(* Truncate the number of digits reported in the answer, if you prefer
  (the following reports only the first 6 significant digits, even though
  the answer was computed to about 50):

NumberForm[TableForm[AIMSeries[sgmodel,2]],6]
*)


(* Check the accuracy of the output (reports number of significant digits
  to the right of decimal point in the answer):

Accuracy[AIMSeries[sgmodel,2]]
*)
