#pragma once
#include <petscdm.h>
#include <petscdmda.h>
#include <petscdmdatypes.h>
#include <petscdmtypes.h>
#include <petscerror.h>
#include <petscksp.h>
#include <petscmat.h>
#include <petscsys.h>
#include <petscsystypes.h>
#include <petscvec.h>
#include <petscviewer.h>

#define tD 100.0
#define xCont 1e-6
#define volfrac 0.1
#define f0 1
#define DIM 3

#define kH 1
#define kL 1e-3
#define xlow 0

PetscErrorCode formkappa(PCCtx *s_ctx, Vec x, PetscInt penal);
PetscErrorCode formMatrix(PCCtx *s_ctx, Mat A);
PetscErrorCode formRHS(PCCtx *s_ctx, Vec rhs, Vec x, PetscInt penal);

