(*
 * $Rev: 1 $
 * Original code: http://www.internationalmathematicasymposium.org/IMS99/paper46/FractionalCalculus.nb
 * Modified by: Jesse Perla
*)

Unprotect[D];
D[x_^m_.,{x_,q_}]:=Gamma[m+1]/Gamma[m-q+1] x^(m-q)/;Head[q]==Real || Head[q]==Rational||Head[q]==Complex
Protect[D];


Remove[RiemannLiouville];
(*--- main function --- *)
RiemannLiouville[f_, {x_, order_, a_: 0}] := Block[{n, int, y},
  		If[NumericQ[order] && Simplify[order > 0], n = Floor[order]; 
   q = order - n];
  		int = 
   Integrate[(x - y)^(-q - 1) (f /. x -> y), {y, a, x}, 
    GenerateConditions -> True];
  		D[int / Gamma[-q], {x, n}] /; FreeQ[int, y]
  		]
		
RiemannLiouville[c_ f_, {x_, order_, a_: 0}] := 
  c RiemannLiouville[f, {x, order, a}] /; FreeQ[c, x];
  
RiemannLiouville[f_ + g_, {x_, order_, a_: 0}] := 
 RiemannLiouville[f, {x, order, a}] + 
  RiemannLiouville[g, {x, order, a}]
  
RiemannLiouville[ 
   RiemannLiouville[f_, {x_, order1_, a_: 0}], {x_, order2_, a_: 0}] :=
   RiemannLiouville[f, {x, order1 + order2, a}] /; order1 < 0;
  
(*
Options[RiemannLiouville] 
SetOptions[RiemannLiouville,ShowConditions->False]
*)

WeylPlus[f_,{t_,q_}]:=RiemannLiouville[(-1)^q f,{t,q,\[Infinity]}]
WeylMinus[f_,{t_,q_}]:=RiemannLiouville[f,{t,q,-\[Infinity]}]
WeylPlus[c_ f_,{x_,order_}]:=c WeylPlus[f,{x,order}] /; FreeQ[c,x];
WeylPlus[f_ + g_,{x_,order_}]:=WeylPlus[f,{x,order}] + WeylPlus[g,{x,order}]
WeylMinus[c_ f_,{x_,order_}]:=c WeylMinus[f,{x,order}] /; FreeQ[c,x];
WeylMinus[f_ + g_,{x_,order_}]:=WeylMinus[f,{x,order}] + WeylMinus[g,{x,order}]
WeylPlus[ WeylPlus[f_,{x_,order1_}],{x_,order2_}]:=WeylPlus[f,{x,order1+order2,a}] /;(order1<0&&order2<0)	
WeylMinus[ WeylMinus[f_,{x_,order1_}],{x_,order2_}]:=WeylMinus[f,{x,order1+order2,a}] /;(order1<0&&order2<0)	
WeylPlus[f_,{x_,0}]:=f 
WeylPlus[ WeylPlus[f_,{x_,order1_}],{x_,order2_}]:=WeylPlus[f,{x,order1+order2,a}] /;(order1>0&&order2>0)	
WeylMinus[ WeylMinus[f_,{x_,order1_}],{x_,order2_}]:=WeylMinus[f,{x,order1+order2,a}] /;(order1>0&&order2>0)
(*
SetOptions[WeylMinus,ShowConditions->False]
SetOptions[WeylPlus,ShowConditions->False];
*)