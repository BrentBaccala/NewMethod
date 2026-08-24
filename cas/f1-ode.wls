#!/usr/bin/env wolframscript

(* NewMethod paper: the xi ("first") separated equation

     d/dxi ( xi df1/dxi ) + [ (1/2) E xi - (1/4) m^2/xi + beta1 ] f1 = 0

   Everything except xi and f1 is a constant:
     EE     -- the energy E   (E is Protected in WL: Euler's number)
     m      -- the magnetic quantum number
     beta1  -- the first separation constant

   NOTE on formatting: this is a .wls script, so every top-level line that
   parses as a COMPLETE expression ends the statement.  Multi-line
   expressions must therefore be wrapped in ( ... ).
*)

ClearAll[f1, xi, EE, m, beta1];

(* Constant is the attribute D/Dt consult to know dEE/dxi == 0.  DSolve
   infers it anyway (only f1 is the dependent function), but this states it. *)
SetAttributes[{EE, m, beta1}, Constant];

(* xi > 0 holds identically: xi = r + z is a parabolic coordinate
   (NewMethod.tex, eq. above (\ref{f1ode})) and r >= |z|.  It is NOT needed
   to get the solution -- DSolve returns the identical result with
   $Assumptions = {} -- it only pins the branch of xi^(m/2) = E^(m Log[xi]/2)
   and of the 1/Sqrt[xi] in the Whittaker form.  Under xi < 0 those pick up
   the I Pi branch shift. *)
$Assumptions = {EE \[Element] Reals, beta1 \[Element] Reals,
                m \[Element] Integers, xi > 0};

(* Print does NOT apply a wrapper like InputForm/TeXForm to its argument --
   it prints the wrapper.  ToString[expr, form] is what actually formats. *)
show[label_, expr_] := (
  Print["\n== ", label, " =="];
  Print["  ", ToString[expr, InputForm]];
  Print["  TeX: ", ToString[expr, TeXForm]];
);

Print["Dt[EE, xi] = ", Dt[EE, xi], "   (0 => EE is treated as a constant)"];

(* ---- the equation, in the divergence form printed in the paper ---- *)
ode = ( D[xi D[f1[xi], xi], xi]
        + ((1/2) EE xi - (1/4) m^2/xi + beta1) f1[xi] == 0 );

show["ODE, divergence form", ode];
show["ODE, expanded and multiplied through by xi",
     Collect[Expand[xi First[ode]], f1[xi] | f1'[xi] | f1''[xi]] == 0];

(* ---- general solution ----
   Solve for f1, NOT f1[xi]: DSolve[ode, f1[xi], xi] returns a rule for the
   symbol f1[xi] alone, so f1'[xi] and f1''[xi] survive the replacement and
   the verification cannot close.  DSolve[ode, f1, xi] returns
   f1 -> Function[{xi}, ...], which substitutes into the derivatives too. *)
sol = DSolve[ode, f1, xi];

(* Print the ENTIRE DSolve result, exactly as returned: a list (one element
   per solution branch) of lists of rules f1 -> Function[{xi}, ...]. *)
show["DSolve raw result (all branches)", sol];
Print["\n   number of solution branches: ", Length[sol]];
Do[show["  branch " <> ToString[i] <> ": f1[xi]", f1[xi] /. sol[[i]]],
   {i, Length[sol]}];

(* Verification.  FullSimplify on the symbolic HypergeometricU/LaguerreL
   solution with symbolic m segfaults the 14.0 kernel, so verify numerically:
   substitute the closed form back and evaluate the residual to 40 digits at
   assorted (E, beta1, m, xi).  Residual ~ 0 to working precision == solved. *)
resid[f_] := D[xi D[f, xi], xi] + ((1/2) EE xi - (1/4) m^2/xi + beta1) f;

$MaxExtraPrecision = 400;   (* else N::meprec floods the log with the raw expr *)

check[label_, f_] := (
  Print["\n== residual of ", label, " (want 0) =="];
  Quiet[Do[Print["   E=", ev, " beta1=", bv, " m=", mv, " xi=", xv, "  ->  ",
       Chop[N[resid[f] /. {C[1] -> 1, C[2] -> 1} /. {EE -> ev, beta1 -> bv,
                           m -> mv, xi -> xv}, 40], 10^-25]],
     {ev, {-1, 3}}, {bv, {0, 7/5}}, {mv, {0, 2}}, {xv, {13/10}}], N::meprec];
);

Do[check["general solution, branch " <> ToString[i], f1[xi] /. sol[[i]]],
   {i, Length[sol]}];

(* ---- Bessel limit: beta1 -> 0 ---- *)
ode0 = ode /. beta1 -> 0;
sol0 = DSolve[ode0, f1, xi];
show["beta1 = 0, raw result (all branches)", sol0];
Do[show["  beta1=0 branch " <> ToString[i] <> ": f1[xi]", f1[xi] /. sol0[[i]]],
   {i, Length[sol0]}];
Print["\n== verification of the beta1=0 Bessel solution (want True) =="];
Do[Print["   branch ", i, ": ", FullSimplify[ode0 /. sol0[[i]], xi > 0]],
   {i, Length[sol0]}];

(* ---- the textbook Whittaker form ----
   xi^2 f1'' + xi f1' + (E xi^2/2 + beta1 xi - m^2/4) f1 = 0  becomes the
   Whittaker equation under f1 = xi^(-1/2) w(2 k xi), k = Sqrt[-E/2], with
   kappa = beta1/(2k) and mu = m/2. *)
kk    = Sqrt[-EE/2];
kappa = beta1/(2 kk);
whitW = xi^(-1/2) WhittakerW[kappa, m/2, 2 kk xi];
whitM = xi^(-1/2) WhittakerM[kappa, m/2, 2 kk xi];
show["Whittaker W form", whitW];
show["Whittaker M form", whitM];

check["Whittaker W form", whitW];
check["Whittaker M form", whitM];

(* ---- general case rewritten ---- *)
Do[show["general solution under FunctionExpand, branch " <> ToString[i],
        FunctionExpand[f1[xi] /. sol[[i]]]], {i, Length[sol]}];
