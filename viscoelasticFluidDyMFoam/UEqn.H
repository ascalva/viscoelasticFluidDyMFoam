tmp<fvVectorMatrix> UEqn
(
    fvm::ddt(U)
  + fvm::div(phi, U)
  - visco.divTau(U)
);

UEqn().relax();

solve(UEqn() == -fvc::grad(p));
