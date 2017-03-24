(* ::Package:: *)

<<ETK`;
<<CSE`;
(* I think I want to do a depthfirst MapAll so simplify the expressions from the bottom up if they are recursive?*)
CollectFullSimplifyImpl[exp_, pat_]:= If[MatchQ[#,pat], FullSimplify@#, #]& //@ Collect[exp, pat, FullSimplify];
(* Why Fixed point? Want to allow the recollection//further simplication if the exponents have simplied, but having trouble finding a minimal example. Might be unncessary... *)
CollectFullSimplify[exp_, pat_]:= FixedPoint[CollectFullSimplifyImpl[#, pat]&, exp, 10] (* 10 arbitrary just in case non-monotone *)

exp = ((b - a)/a - b/a )x^(((b - a)/a) - (b/a)  )+(1 + a)x + (b - 1)x;
exp2 = ((b - a)/a - b/a )x^(((b - a)/a) - (b/a)  )+ x^x^((b - a)/a - b/a-1/x+x^((b - a)/a - b/a)); (* Test out the recursion *)
Print["Collect, then Collect with simplification, then CollectFullSimplify"]
Collect[exp, x^_]
Collect[exp, x^_, FullSimplify]
CollectFullSimplify[exp, x^_]

Print["Collect, then Collect with simplification, then CollectFullSimplify"]
Collect[exp2, x^_]
Collect[exp2, x^_, FullSimplify]
CollectFullSimplify[exp, x^_]


   $Assumptions =  
 g > 0 && r > 0 && g < r && b[l] > 0 && b[h] > 0 && 
  g + b[l] + b[h]
exp3 = (z^(((g - r)^2 - g b[h] + (-g + r) b[l])/(
     g (g - r - b[h]))) (g z^(-((r b[h])/(
         g (-g + r + b[h])))) (-g + r + b[h]) + (g - 
          r) z^(((g - r) (r + b[l]))/(
        g (g - r - b[h]))) (g - 
          r - b[h] - b[l])))/((g - 
       r) (-r b[h] + (g - r) (r + b[l])))
Collect[exp3, z^_]
Collect[exp3, z^_, Simplify]
CollectFullSimplifyImpl[exp3, z^_]
CollectFullSimplify[exp3, z^_]
ExtractCommonSubexpressions@%




