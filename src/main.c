#include <petscdm.h>
#include <petscdmda.h>
#include <petscdmdatypes.h>
#include <petscdmtypes.h>
#include <petscerror.h>
#include <petscksp.h>
#include <petsclog.h>
#include <petscmat.h>
#include <petscoptions.h>
#include <petscpc.h>
#include <petscsys.h>
#include <petscsystypes.h>
#include <petsctime.h>
#include <petscvec.h>
#include <petscviewer.h>
#include <petscviewerhdf5.h>
#include <slepceps.h>

#include "system.h"

int main(int argc, char **argv) {
  PetscCall(SlepcInitialize(&argc, &argv, (char *)0, "MMG\n"));
  PCCtx test;
  Mat A;
  Vec rhs, sol;
  KSP ksp;
  PC pc;

  PetscCall(PC_init(&test));
  PetscCall(DMCreateMatrix(test.dm, &A));
  PetscCall(DMCreateGlobalVector(test.dm, &rhs));
  PetscCall(DMCreateGlobalVector(test.dm, &sol));

  PetscCall(formBoundary(&test));
  PetscCall(formMatrix(&test, A));
  PetscCall(formRHS(&test, rhs));
  PetscCall(KSPCreate(PETSC_COMM_WORLD, &ksp));
  PetscCall(KSPSetOperators(ksp, A, A));
  PetscCall(KSPSetFromOptions(ksp));
  PetscCall(KSPSetUp(ksp));
  PetscCall(KSPSolve(ksp, rhs, sol));
  PetscCall(VecView(sol, PETSC_VIEWER_STDOUT_WORLD));

  PetscCall(SlepcFinalize());
}
