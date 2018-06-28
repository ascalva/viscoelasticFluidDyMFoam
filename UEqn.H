// Momentum predictor

tmp<fvVectorMatrix> UEqn
(
    fvm::ddt(U)
  + fvm::div(phi, U)
  - visco.divTau(U)
);

if( momentumPredictor ) {
    solve(UEqn() == -fvc::grad(p));
}
