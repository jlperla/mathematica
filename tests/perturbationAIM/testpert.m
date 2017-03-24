(*
 * $Rev: 1 $
 * Original code: http://www.ericswanson.us/perturbation.html
 * Modified by: Jesse Perla
*)

(* Load perturbationAIM *)
<<perturbationAIM.m


(* Note: these models are somewhat contrived (not much economic content)
 and are meant to test out various aspects of the code.  Nonetheless,
 they may serve as useful examples of what the code is capable of.  The
 arbitrary-precision calculations at the end may be of particular
 interest.
*)


Print["\n\n"] ;
Print["Test 1: A simple, contrived New Keynesian model"] ;
eqns={
  Y[t] == alpha *Y[t+1] - gamma *(E^(r[t+1]-pi[t+1]) -1) + eps[y][t],
  pi[t] == beta *pi[t+1] + delta *Y[t] + eps[pi][t],
  r[t] == tayly *Y[t-1] + taylpi *pi[t-1]
}

logvars = {}
logrules = Map[#[x_]->E^(#[x])&, logvars]

parametervalues={
  alpha->1,
  beta->0.99,
  gamma->0.15,
  delta->0.25,
  tayly->0.5,
  taylpi->1.5
}

test1 = eqns /.logrules //.parametervalues
Print[AIMSeries[test1,4] //TableForm]





Print["\n\n"] ;
Print["Test 2: A version of the model with no shocks\n"] ;
Print["Note: some Windows versions of Mathematica will crash on this test"] ;
Print["  due to improper handling of multiplication by an empty matrix"] ;
Print["  that is caused by the lack of shocks in the model."] ;
Print["If you are using Windows, just comment out this test (Test 2) and re-run"] ;
Print["  testpert.m.  You will also have to make sure that your models always have"] ;
Print["  at least one shock in them to avoid a crash."] ;
eqns={
  Y[t] == alpha *Y[t+1] - gamma *(E^(r[t+1]-pi[t+1]) -1),
  pi[t] == beta *pi[t+1] + delta *Y[t],
  r[t] == tayly *Y[t-1] + taylpi *pi[t-1]
}

logvars = {}
logrules = Map[#[x_]->E^(#[x])&, logvars]

parametervalues={
  alpha->1,
  beta->0.99,
  gamma->0.15,
  delta->0.25,
  tayly->0.5,
  taylpi->1.5
}

test2 = eqns /.logrules //.parametervalues
Print[AIMSeries[test2,4] //TableForm]





Print["\n\n"] ;
Print["Test 3: A version of the model with no lags"] ;
eqns={
  Y[t] == alpha *Y[t+1] - gamma *(E^(r[t+1]-pi[t+1]) -1) + eps[y][t],
  pi[t] == beta *pi[t+1] + delta *Y[t] + eps[pi][t],
  r[t] == tayly *Y[t] + taylpi *pi[t]
}

logvars = {}
logrules = Map[#[x_]->E^(#[x])&, logvars]

parametervalues={
  alpha->1,
  beta->0.99,
  gamma->0.15,
  delta->0.25,
  tayly->0.5,
  taylpi->1.5
}

test3 = eqns /.logrules //.parametervalues
Print[AIMSeries[test3,4] //TableForm]





Print["\n\n"] ;
Print["Test 4: A version of the model with no leads"] ;
eqns={
  Y[t] == alpha *Y[t-1] - gamma *(E^(r[t]-pi[t]) -1) + eps[y][t],
  pi[t] == beta *pi[t-1] + delta *Y[t-1] + eps[pi][t],
  r[t] == tayly *Y[t-1] + taylpi *pi[t-1]
}

logvars = {}
logrules = Map[#[x_]->E^(#[x])&, logvars]

parametervalues={
  alpha->1,
  beta->0.99,
  gamma->0.15,
  delta->0.25,
  tayly->0.5,
  taylpi->1.5
}

test4 = eqns /.logrules //.parametervalues
Print[AIMSeries[test4,4] //TableForm]






Print["\n\n"] ;
Print["Test 5: A version of the model with multiple lags"] ;
eqns={
  Y[t] == alpha *Y[t+1] - gamma *(E^(r[t+1]-pi[t+1]) -1) + eps[y][t],
  pi[t] == beta *pi[t-4] + delta *Y[t] + eps[pi][t],
  r[t] == tayly * Y[t-1] + taylpi *pi[t-1]
}

logvars = {}
logrules = Map[#[x_]->E^(#[x])&, logvars]

parametervalues={
  alpha->1,
  beta->0.99,
  gamma->0.15,
  delta->0.25,
  tayly->0.5,
  taylpi->1.5
}

test5 = eqns /.logrules //.parametervalues
Print[AIMSeries[test5,3] //TableForm]






Print["\n\n"] ;
Print["Test 6: A version of the model with multiple leads"] ;
eqns={
  Y[t] == alpha *Y[t+1] - gamma *(E^(r[t+1]-pi[t+1]) -1) + eps[y][t],
  pi[t] == beta *pi[t+4] + delta *Y[t] + eps[pi][t],
  r[t] == tayly *Y[t-1] + taylpi *pi[t-1]
}

parametervalues={
  alpha->1,
  beta->0.99,
  gamma->0.15,
  delta->0.25,
  tayly->0.5,
  taylpi->1.5
}

logvars = {}
logrules = Map[#[x_]->E^(#[x])&, logvars]

test6 = eqns /.logrules //.parametervalues
Print[AIMSeries[test6,3] //TableForm]






Print["\n\n"] ;
Print["Test 7: A version of the model with one equation"] ;
eqns = Y[t] == alpha *Y[t+1] + beta *Y[t-1] + E^eps[y][t]

test7 = eqns /. {alpha->.9, beta->.05}
Print[AIMSeries[test7,4] //TableForm]







Print["\n\n"] ;
Print["Test 8: A version of the model with eps[_][t+1]"] ;
eqns={
  Y[t] == alpha *Y[t+1] - gamma *(E^(r[t+1]-pi[t+1]) -1) + eps[y][t],
  pi[t] == beta *pi[t+1] + delta *Y[t] + eps[pi][t+1],
  r[t] == tayly * Y[t-1] + taylpi *pi[t-1]
}

parametervalues={
  alpha->1,
  beta->0.99,
  gamma->0.15,
  delta->0.25,
  tayly->0.5,
  taylpi->1.5
}

logvars = {}
logrules = Map[#[x_]->E^(#[x])&, logvars]

test8 = eqns /.logrules /.parametervalues
Print[AIMSeries[test8,4] //TableForm]






Print["\n\n"] ;
Print["Test 9: A version of the model with arbitrary precision coefficients"];
eqns={
  Y[t] == alpha *Y[t+1] - gamma *(E^(r[t+1]-pi[t+1]) -1) + eps[y][t],
  pi[t] == beta *pi[t+1] + delta *Y[t] + eps[pi][t],
  r[t] == tayly * Y[t-1] + taylpi *pi[t-1]
}

parametervalues={
  alpha->1`60, (* 60 digits of precision *)
  beta->0.99`60,
  gamma->0.15`60,
  delta->0.25`60,
  tayly->0.5`60,
  taylpi->1.5`60
}

logvars = {}
logrules = Map[#[x_]->E^(#[x])&, logvars]

test9 = eqns /.logrules /.parametervalues
AIMZeroTol=10^-15
AIMSS[test9,AIMPrecision->50]
Print[AIMSeries[test9,4] //TableForm]

