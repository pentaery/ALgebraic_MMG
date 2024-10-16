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
#define f0 1

#define kH 1
#define kL 1e-6

#define DIM 3

typedef struct preconditioner_context {
  DM dm;
  Vec kappa[DIM], *ms_bases_c, *ms_bases_cc;
  Vec cost;
  Vec boundary;
  PetscInt *coarse_startx, *coarse_lenx, *coarse_starty, *coarse_leny,
      *coarse_startz, *coarse_lenz;
  PetscInt sub_domains, lv2_eigen_op;
  PetscInt max_eigen_num_lv1, *eigen_num_lv1;
  PetscInt max_eigen_num_lv2, eigen_num_lv2;
  PetscScalar H_x, H_y, H_z, L, W, H;
  // 总的长宽高，网格的长宽高
  PetscScalar *eigen_max_lv1, *eigen_min_lv1, eigen_bd_lv1, eigen_max_lv2,
      eigen_min_lv2, eigen_bd_lv2;
  PetscInt M, N, P;
  // 网格数
  Mat Rc, Rcc;
  PetscScalar meas_elem, meas_face_xy, meas_face_yz, meas_face_zx;
  PetscInt coarse_elem_num;
  Vec *ms_bases_c_tmp;
} PCCtx;

PetscErrorCode formBoundary(PCCtx *s_ctx);
PetscErrorCode formkappa(PCCtx *s_ctx);
PetscErrorCode formMatrix(PCCtx *s_ctx, Mat A);
PetscErrorCode formRHS(PCCtx *s_ctx, Vec rhs);
PetscErrorCode PC_init(PCCtx *s_ctx);