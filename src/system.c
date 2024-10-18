#include "system.h"

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

PetscErrorCode formBoundary(PCCtx *s_ctx) {
  PetscFunctionBeginUser;
  PetscInt startx, starty, startz, nx, ny, nz, ex, ey, ez;
  PetscScalar ***array;
  PetscCall(DMCreateGlobalVector(s_ctx->dm, &s_ctx->boundary));
  PetscCall(
      DMDAGetCorners(s_ctx->dm, &startx, &starty, &startz, &nx, &ny, &nz));
  PetscCall(DMDAVecGetArray(s_ctx->dm, s_ctx->boundary, &array));
  for (ez = startz; ez < startz + nz; ++ez) {
    for (ey = starty; ey < starty + ny; ++ey) {
      for (ex = startx; ex < startx + nx; ++ex) {
        if (ex >= PetscFloorReal(0.45 * s_ctx->M) &&
            ex <= PetscCeilReal(0.55 * s_ctx->M) - 1 &&
            ey >= PetscFloorReal(0.45 * s_ctx->N) &&
            ey <= PetscCeilReal(0.55 * s_ctx->N) - 1 && ez == 0) {
          array[ez][ey][ex] = 1;
        } else {
          array[ez][ey][ex] = 0;
        }
      }
    }
  }

  PetscCall(DMDAVecRestoreArray(s_ctx->dm, s_ctx->boundary, &array));

  PetscFunctionReturn(0);
}

PetscErrorCode formkappa(PCCtx *s_ctx) {
  PetscFunctionBeginUser;
  PetscInt startx, starty, startz, nx, ny, nz, ex, ey, ez, i;
  PetscScalar ***arraykappa[DIM];
  PetscCall(
      DMDAGetCorners(s_ctx->dm, &startx, &starty, &startz, &nx, &ny, &nz));
  for (i = 0; i < DIM; ++i) {
    PetscCall(DMDAVecGetArray(s_ctx->dm, s_ctx->kappa[i], &arraykappa[i]));
  }
  for (ez = startz; ez < startz + nz; ++ez) {
    for (ey = starty; ey < starty + ny; ++ey) {
      for (ex = startx; ex < startx + nx; ++ex) {
        for (i = 0; i < DIM; ++i) {
          arraykappa[i][ez][ey][ex] = 1.0;
        }
      }
    }
  }
  for (i = 0; i < DIM; ++i) {
    PetscCall(DMDAVecRestoreArray(s_ctx->dm, s_ctx->kappa[i], &arraykappa[i]));
  }

  PetscFunctionReturn(0);
}

PetscErrorCode formMatrix(PCCtx *s_ctx, Mat A) {
  PetscFunctionBeginUser;
  PetscCall(formkappa(s_ctx));
  Vec kappa_loc[DIM];  // Destroy later.
  PetscInt startx, starty, startz, nx, ny, nz, ex, ey, ez, i;
  PetscScalar ***arr_kappa_3d[DIM], ***arrayBoundary, val_A[2][2], avg_kappa_e;
  MatStencil row[2], col[2];

  for (i = 0; i < DIM; ++i) {
    PetscCall(DMGetLocalVector(s_ctx->dm, &kappa_loc[i]));
    PetscCall(DMGlobalToLocal(s_ctx->dm, s_ctx->kappa[i], INSERT_VALUES,
                              kappa_loc[i]));
    PetscCall(DMDAVecGetArrayRead(s_ctx->dm, kappa_loc[i], &arr_kappa_3d[i]));
  }
  PetscCall(DMDAVecGetArrayRead(s_ctx->dm, s_ctx->boundary, &arrayBoundary));

  PetscCall(
      DMDAGetCorners(s_ctx->dm, &startx, &starty, &startz, &nx, &ny, &nz));
  PetscCall(MatZeroEntries(A));
  for (ez = startz; ez < startz + nz; ++ez)
    for (ey = starty; ey < starty + ny; ++ey)
      for (ex = startx; ex < startx + nx; ++ex) {
        if (ex >= 1) {
          row[0] = (MatStencil){.i = ex - 1, .j = ey, .k = ez};
          row[1] = (MatStencil){.i = ex, .j = ey, .k = ez};
          col[0] = (MatStencil){.i = ex - 1, .j = ey, .k = ez};
          col[1] = (MatStencil){.i = ex, .j = ey, .k = ez};
          avg_kappa_e = 2.0 / (1.0 / arr_kappa_3d[0][ez][ey][ex - 1] +
                               1.0 / arr_kappa_3d[0][ez][ey][ex]);
          val_A[0][0] = s_ctx->H_y * s_ctx->H_z / s_ctx->H_x * avg_kappa_e;
          val_A[0][1] = -s_ctx->H_y * s_ctx->H_z / s_ctx->H_x * avg_kappa_e;
          val_A[1][0] = -s_ctx->H_y * s_ctx->H_z / s_ctx->H_x * avg_kappa_e;
          val_A[1][1] = s_ctx->H_y * s_ctx->H_z / s_ctx->H_x * avg_kappa_e;
          PetscCall(MatSetValuesStencil(A, 2, &row[0], 2, &col[0], &val_A[0][0],
                                        ADD_VALUES));
        }
        if (ey >= 1) {
          row[0] = (MatStencil){.i = ex, .j = ey - 1, .k = ez};
          row[1] = (MatStencil){.i = ex, .j = ey, .k = ez};
          col[0] = (MatStencil){.i = ex, .j = ey - 1, .k = ez};
          col[1] = (MatStencil){.i = ex, .j = ey, .k = ez};
          avg_kappa_e = 2.0 / (1.0 / arr_kappa_3d[1][ez][ey - 1][ex] +
                               1.0 / arr_kappa_3d[1][ez][ey][ex]);
          val_A[0][0] = s_ctx->H_x * s_ctx->H_z / s_ctx->H_y * avg_kappa_e;
          val_A[0][1] = -s_ctx->H_x * s_ctx->H_z / s_ctx->H_y * avg_kappa_e;
          val_A[1][0] = -s_ctx->H_x * s_ctx->H_z / s_ctx->H_y * avg_kappa_e;
          val_A[1][1] = s_ctx->H_x * s_ctx->H_z / s_ctx->H_y * avg_kappa_e;
          PetscCall(MatSetValuesStencil(A, 2, &row[0], 2, &col[0], &val_A[0][0],
                                        ADD_VALUES));
        }
        if (ez >= 1) {
          row[0] = (MatStencil){.i = ex, .j = ey, .k = ez - 1};
          row[1] = (MatStencil){.i = ex, .j = ey, .k = ez};
          col[0] = (MatStencil){.i = ex, .j = ey, .k = ez - 1};
          col[1] = (MatStencil){.i = ex, .j = ey, .k = ez};
          avg_kappa_e = 2.0 / (1.0 / arr_kappa_3d[2][ez - 1][ey][ex] +
                               1.0 / arr_kappa_3d[2][ez][ey][ex]);
          val_A[0][0] = s_ctx->H_x * s_ctx->H_y / s_ctx->H_z * avg_kappa_e;
          val_A[0][1] = -s_ctx->H_x * s_ctx->H_y / s_ctx->H_z * avg_kappa_e;
          val_A[1][0] = -s_ctx->H_x * s_ctx->H_y / s_ctx->H_z * avg_kappa_e;
          val_A[1][1] = s_ctx->H_x * s_ctx->H_y / s_ctx->H_z * avg_kappa_e;
          PetscCall(MatSetValuesStencil(A, 2, &row[0], 2, &col[0], &val_A[0][0],
                                        ADD_VALUES));
        }
        if (arrayBoundary[ez][ey][ex] == 1) {
          col[0] = (MatStencil){.i = ex, .j = ey, .k = ez};
          row[0] = (MatStencil){.i = ex, .j = ey, .k = ez};
          val_A[0][0] = 2.0 * s_ctx->H_x * s_ctx->H_y / s_ctx->H_z *
                        arr_kappa_3d[2][ez][ey][ex];
          PetscCall(MatSetValuesStencil(A, 1, &col[0], 1, &row[0], &val_A[0][0],
                                        ADD_VALUES));
        }
      }
  // A的赋值
  PetscCall(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY));
  for (i = 0; i < DIM; ++i) {
    PetscCall(
        DMDAVecRestoreArrayRead(s_ctx->dm, kappa_loc[i], &arr_kappa_3d[i]));
    PetscCall(DMRestoreLocalVector(s_ctx->dm, &kappa_loc[i]));
    PetscCall(VecDestroy(&kappa_loc[i]));
  }
  PetscCall(
      DMDAVecRestoreArrayRead(s_ctx->dm, s_ctx->boundary, &arrayBoundary));
  PetscFunctionReturn(0);
}

PetscErrorCode formRHS(PCCtx *s_ctx, Vec rhs) {
  PetscFunctionBeginUser;
  PetscScalar ***array, ***arraykappa, ***arrayBoundary;
  PetscInt startx, starty, startz, nx, ny, nz, ex, ey, ez;
  PetscCall(
      DMDAGetCorners(s_ctx->dm, &startx, &starty, &startz, &nx, &ny, &nz));

  // Set RHS
  PetscCall(DMDAVecGetArrayRead(s_ctx->dm, s_ctx->kappa[2], &arraykappa));
  PetscCall(DMDAVecGetArrayRead(s_ctx->dm, s_ctx->boundary, &arrayBoundary));
  PetscCall(VecSet(rhs, 0));
  PetscCall(DMDAVecGetArray(s_ctx->dm, rhs, &array));
  for (ez = startz; ez < startz + nz; ++ez) {
    for (ey = starty; ey < starty + ny; ey++) {
      for (ex = startx; ex < startx + nx; ex++) {
        array[ez][ey][ex] += s_ctx->H_x * s_ctx->H_y * s_ctx->H_z * f0;
        if (arrayBoundary[ez][ey][ex] > 0.5) {
          // array[ez][ey][ex] += 2 * arraykappa[ez][ey][ex] * tD * s_ctx->H_x *
          //                      s_ctx->H_y / s_ctx->H_z;
          array[ez][ey][ex] = 0;
        }
      }
    }
  }

  PetscCall(DMDAVecRestoreArrayRead(s_ctx->dm, s_ctx->kappa[2], &arraykappa));
  PetscCall(
      DMDAVecRestoreArrayRead(s_ctx->dm, s_ctx->boundary, &arrayBoundary));
  PetscCall(DMDAVecRestoreArray(s_ctx->dm, rhs, &array));

  PetscFunctionReturn(0);
}

PetscErrorCode PC_init(PCCtx *s_ctx) {
  PetscFunctionBeginUser;
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Initializing the PC context...\n"));
  PetscInt grid = 50;
  PetscCall(PetscOptionsGetInt(NULL, NULL, "-grid", &grid, NULL));
  PetscInt mesh[DIM] = {grid, grid, grid};
  PetscScalar dom[DIM] = {1.0, 1.0, 1.0};
  PetscCheck((dom[0] > 0.0 && dom[1] > 0.0 && dom[2] > 0.0), PETSC_COMM_WORLD,
             PETSC_ERR_ARG_WRONG,
             "Errors in dom(L, W, H)=[%.5f, %.5f, %.5f].\n", dom[0], dom[1],
             dom[2]);
  PetscCheck((mesh[0] > 0 && mesh[1] > 0 && mesh[2] > 0), PETSC_COMM_WORLD,
             PETSC_ERR_ARG_WRONG, "Errors in mesh(M, N, P)=[%d, %d, %d].\n",
             mesh[0], mesh[1], mesh[2]);
  s_ctx->L = dom[0];
  s_ctx->W = dom[1];
  s_ctx->H = dom[2];
  s_ctx->M = mesh[0];
  s_ctx->N = mesh[1];
  s_ctx->P = mesh[2];

  s_ctx->sub_domains = 2;
  PetscCall(PetscOptionsGetInt(NULL, NULL, "-sd", &s_ctx->sub_domains, NULL));
  PetscCheck(s_ctx->sub_domains >= 1, PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,
             "Error in sub_domains=%d.\n", s_ctx->sub_domains);
  PetscCall(
      PetscPrintf(PETSC_COMM_WORLD, "sub_domains=%d\n", s_ctx->sub_domains));

  s_ctx->max_eigen_num_lv1 = 4;
  PetscCall(PetscOptionsGetInt(NULL, NULL, "-en_lv1", &s_ctx->max_eigen_num_lv1,
                               NULL));
  PetscCheck(s_ctx->max_eigen_num_lv1 >= 1, PETSC_COMM_WORLD,
             PETSC_ERR_ARG_WRONG,
             "Error in max_eigen_num_lv1=%d for the level-1 problem.\n",
             s_ctx->max_eigen_num_lv1);
  PetscCall(
      PetscPrintf(PETSC_COMM_WORLD, "en_lv1=%d\n", s_ctx->max_eigen_num_lv1));

  s_ctx->max_eigen_num_lv2 = 4;

  PetscCall(PetscOptionsGetInt(NULL, NULL, "-en_lv2", &s_ctx->max_eigen_num_lv2,
                               NULL));
  PetscCheck(s_ctx->max_eigen_num_lv2 >= 1, PETSC_COMM_WORLD,
             PETSC_ERR_ARG_WRONG,
             "Error in max_eigen_num_lv2=%d for the level-2 problem.\n",
             s_ctx->max_eigen_num_lv2);
  PetscCall(
      PetscPrintf(PETSC_COMM_WORLD, "en_lv2=%d\n", s_ctx->max_eigen_num_lv2));

  s_ctx->H_x = dom[0] / (double)mesh[0];
  s_ctx->H_y = dom[1] / (double)mesh[1];
  s_ctx->H_z = dom[2] / (double)mesh[2];

  PetscCall(DMDACreate3d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,
                         DM_BOUNDARY_NONE, DMDA_STENCIL_STAR, s_ctx->M,
                         s_ctx->N, s_ctx->P, PETSC_DECIDE, PETSC_DECIDE,
                         PETSC_DECIDE, 1, 1, NULL, NULL, NULL, &(s_ctx->dm)));
  // If oversampling=1, DMDA has a ghost point width=1 now, and this will
  // change the construction of A_i in level-1.
  PetscCall(DMSetUp(s_ctx->dm));

  for (PetscInt i = 0; i < DIM; ++i)
    PetscCall(DMCreateGlobalVector(s_ctx->dm, &s_ctx->kappa[i]));

  if (s_ctx->sub_domains == 1) {
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,
                          "The subdomain may not be proper (too "
                          "long/wide/high), reset sub_domains from %d to 2.\n",
                          s_ctx->sub_domains));
    s_ctx->sub_domains = 2;
  }
  s_ctx->coarse_elem_num =
      s_ctx->sub_domains * s_ctx->sub_domains * s_ctx->sub_domains;

  PetscCall(PetscMalloc1(s_ctx->sub_domains, &s_ctx->coarse_startx));
  PetscCall(PetscMalloc1(s_ctx->sub_domains, &s_ctx->coarse_starty));
  PetscCall(PetscMalloc1(s_ctx->sub_domains, &s_ctx->coarse_startz));
  PetscCall(PetscMalloc1(s_ctx->sub_domains, &s_ctx->coarse_lenx));
  PetscCall(PetscMalloc1(s_ctx->sub_domains, &s_ctx->coarse_leny));
  PetscCall(PetscMalloc1(s_ctx->sub_domains, &s_ctx->coarse_lenz));
  PetscCall(PetscMalloc1(s_ctx->coarse_elem_num, &s_ctx->eigen_max_lv1));
  PetscCall(PetscMalloc1(s_ctx->coarse_elem_num, &s_ctx->eigen_min_lv1));

  s_ctx->meas_elem = s_ctx->H_x * s_ctx->H_y * s_ctx->H_z;
  s_ctx->meas_face_yz = s_ctx->H_y * s_ctx->H_z;
  s_ctx->meas_face_zx = s_ctx->H_z * s_ctx->H_x;
  s_ctx->meas_face_xy = s_ctx->H_x * s_ctx->H_y;
  s_ctx->ms_bases_c_tmp = (Vec *)malloc(s_ctx->max_eigen_num_lv1 * sizeof(Vec));

  PetscInt proc_nx, proc_ny, proc_nz;
  PetscCall(DMDAGetCorners(s_ctx->dm, NULL, NULL, NULL, &proc_nx, &proc_ny,
                           &proc_nz));

  PetscInt eigen_len_lv1, eigen_len_lv1_max;
  eigen_len_lv1 =
      s_ctx->coarse_lenx[0] * s_ctx->coarse_leny[0] * s_ctx->coarse_lenz[0];
  PetscCallMPI(MPI_Allreduce(&eigen_len_lv1, &eigen_len_lv1_max, 1, MPI_INT,
                             MPI_MAX, PETSC_COMM_WORLD));
  PetscCall(MatCreateAIJ(PETSC_COMM_WORLD, proc_nx * proc_ny * proc_nz,
                         s_ctx->coarse_elem_num * s_ctx->max_eigen_num_lv1,
                         PETSC_DEFAULT, PETSC_DEFAULT, eigen_len_lv1_max, NULL,
                         0, NULL, &s_ctx->Rc));

  PetscCall(MatCreateAIJ(PETSC_COMM_WORLD,
                         s_ctx->coarse_elem_num * s_ctx->max_eigen_num_lv1,
                         s_ctx->max_eigen_num_lv2, PETSC_DEFAULT, PETSC_DEFAULT,
                         s_ctx->coarse_elem_num * s_ctx->max_eigen_num_lv1,
                         NULL, 0, NULL, &s_ctx->Rcc));

  PetscInt m, n;
  PetscCall(MatGetSize(s_ctx->Rc, &m, &n));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "m1=%d, n1=%d\n", m, n));
  PetscCall(MatGetSize(s_ctx->Rcc, &m, &n));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "m1=%d, n1=%d\n", m, n));

  PetscFunctionReturn(0);
}