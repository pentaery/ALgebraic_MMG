#include <petscao.h>
#include <petscdm.h>
#include <petscdmda.h>
#include <petscdmdatypes.h>
#include <petscdmtypes.h>
#include <petscerror.h>
#include <petscis.h>
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
#include <time.h>

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
  PetscCall(MatView(A, PETSC_VIEWER_STDOUT_WORLD));


  MatPartitioning part;
  IS is, isg;
  AO ao;
  PetscCall(MatPartitioningCreate(PETSC_COMM_WORLD, &part));
  PetscCall(MatPartitioningSetAdjacency(part, A));
  PetscCall(MatPartitioningSetType(part, MATPARTITIONINGPARMETIS));
  PetscCall(MatPartitioningSetFromOptions(part));
  PetscCall(MatPartitioningApply(part, &is));
  PetscCall(MatPartitioningDestroy(&part));
  PetscCall(ISView(is, PETSC_VIEWER_STDOUT_WORLD));
  PetscCall(ISPartitioningToNumbering(is, &isg));
  PetscCall(ISView(isg, PETSC_VIEWER_STDOUT_WORLD));
  PetscCall(AOCreateBasicIS(isg, NULL, &ao));
  PetscCall(AOView(ao, PETSC_VIEWER_STDOUT_WORLD));
  PetscInt size;
  PetscCall(ISGetSize(is, &size));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Size of IS: %d\n", size));
  PetscCall(ISGetSize(isg, &size));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Size of ISG: %d\n", size));

  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Finished partitioning!\n"));


  PetscCall(KSPCreate(PETSC_COMM_WORLD, &ksp));
  PetscCall(KSPSetOperators(ksp, A, A));
  PetscCall(KSPSetFromOptions(ksp));
  PetscCall(KSPSetUp(ksp));
  PetscCall(KSPSolve(ksp, rhs, sol));

  PetscCall(SlepcFinalize());
}
