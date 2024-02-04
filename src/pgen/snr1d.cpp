//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file snr1d.cpp
//! \brief Problem generator for spherical blast wave problem.  Works in Cartesian,
//!        cylindrical, and spherical coordinates.  Contains post-processing code
//!        to check whether blast is spherical for regression tests
//!
//! REFERENCE: P. Londrillo & L. Del Zanna, "High-order upwind schemes for
//!   multidimensional MHD", ApJ, 530, 508 (2000), and references therein.

// C headers

// C++ headers
#include <algorithm>
#include <cmath>
#include <cstdio>     // fopen(), fprintf(), freopen()
#include <cstring>    // strcmp()
#include <sstream>
#include <stdexcept>
#include <string>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../globals.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"

namespace {
  Real rout, rin, pa, da, prat, drat;
  Real mdot, edot, t_m, t_e;
  int mdot_flag, edot_flag;
  void SNSource(MeshBlock *pmb, const Real time, const Real dt,
                const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
                const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
                AthenaArray<Real> &cons_scalar);
} // namespace

void MyBoundary_ix1(MeshBlock *pmb, Coordinates *pco,
                    AthenaArray<Real> &prim, FaceField &b,
                    Real time, Real dt,
                    int il, int iu, int jl, int ju, int kl, int ku, int ngh);

void MySource(MeshBlock *pmb, const Real time, const Real dt,
              const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
              const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
              AthenaArray<Real> &cons_scalar);

void Mesh::InitUserMeshData(ParameterInput *pin) {
  printf("Initializing user mesh data...\n");
  rout = pin->GetReal("problem", "radius");
  rin  = rout - pin->GetOrAddReal("problem", "ramp", 0.0);
  pa   = pin->GetOrAddReal("problem", "pamb", 1.0);
  da   = pin->GetOrAddReal("problem", "damb", 1.0);
  prat = pin->GetReal("problem", "prat");
  drat = pin->GetOrAddReal("problem", "drat", 1.0);
  mdot = pin->GetOrAddReal("problem", "mdot", 0.0);
  edot = pin->GetOrAddReal("problem", "edot", 0.0);
  t_m  = pin->GetOrAddReal("problem", "t_m",  0.0);
  t_e  = pin->GetOrAddReal("problem", "t_e",  0.0);
  mdot_flag = pin->GetOrAddInteger("problem", "mdot_flag", 0);
  edot_flag = pin->GetOrAddInteger("problem", "edot_flag", 0);
  if (pin->GetString("mesh", "ix1_bc") == "user") {
    EnrollUserBoundaryFunction(BoundaryFace::inner_x1, MyBoundary_ix1);
  }
  EnrollUserExplicitSourceFunction(MySource);
  //AllocateUserHistoryOutput(1);
  //EnrollUserHistoryOutput(0, radius, "radius");
  return;
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//! \brief Spherical blast wave test problem generator
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  Real rout = pin->GetReal("problem", "radius");
  Real rin  = rout - pin->GetOrAddReal("problem", "ramp", 0.0);
  Real pa   = pin->GetOrAddReal("problem", "pamb", 1.0);
  Real da   = pin->GetOrAddReal("problem", "damb", 1.0);
  Real prat = pin->GetReal("problem", "prat");
  Real drat = pin->GetOrAddReal("problem", "drat", 1.0);
  Real gamma = peos->GetGamma();
  Real gm1 = gamma - 1.0;

  // get coordinates of center of blast, and convert to Cartesian if necessary
  Real x1_0   = pin->GetOrAddReal("problem", "x1_0", 0.0);
  Real x2_0   = pin->GetOrAddReal("problem", "x2_0", 0.0);
  Real x3_0   = pin->GetOrAddReal("problem", "x3_0", 0.0);
  Real x0, y0, z0;
  if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
    x0 = x1_0;
    y0 = x2_0;
    z0 = x3_0;
  } else if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
    x0 = x1_0*std::cos(x2_0);
    y0 = x1_0*std::sin(x2_0);
    z0 = x3_0;
  } else if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
    x0 = x1_0*std::sin(x2_0)*std::cos(x3_0);
    y0 = x1_0*std::sin(x2_0)*std::sin(x3_0);
    z0 = x1_0*std::cos(x2_0);
  } else {
    // Only check legality of COORDINATE_SYSTEM once in this function
    std::stringstream msg;
    msg << "### FATAL ERROR in blast.cpp ProblemGenerator" << std::endl
        << "Unrecognized COORDINATE_SYSTEM=" << COORDINATE_SYSTEM << std::endl;
    ATHENA_ERROR(msg);
  }

  // setup uniform ambient medium with spherical over-pressured region
  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {
        Real rad;
        if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
          Real x = pcoord->x1v(i);
          Real y = pcoord->x2v(j);
          Real z = pcoord->x3v(k);
          rad = std::sqrt(SQR(x - x0) + SQR(y - y0) + SQR(z - z0));
        } else if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
          Real x = pcoord->x1v(i)*std::cos(pcoord->x2v(j));
          Real y = pcoord->x1v(i)*std::sin(pcoord->x2v(j));
          Real z = pcoord->x3v(k);
          rad = std::sqrt(SQR(x - x0) + SQR(y - y0) + SQR(z - z0));
        } else { // if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0)
          Real x = pcoord->x1v(i)*std::sin(pcoord->x2v(j))*std::cos(pcoord->x3v(k));
          Real y = pcoord->x1v(i)*std::sin(pcoord->x2v(j))*std::sin(pcoord->x3v(k));
          Real z = pcoord->x1v(i)*std::cos(pcoord->x2v(j));
          rad = std::sqrt(SQR(x - x0) + SQR(y - y0) + SQR(z - z0));
        }

        Real den = da;
        if (rad < rout) {
          if (rad < rin) {
            den = drat*da;
          } else {   // add smooth ramp in density
            Real f = (rad-rin) / (rout-rin);
            Real log_den = (1.0-f) * std::log(drat*da) + f * std::log(da);
            den = std::exp(log_den);
          }
        }

        phydro->u(IDN,k,j,i) = den;
        phydro->u(IM1,k,j,i) = 0.0;
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = 0.0;
        if (NON_BAROTROPIC_EOS) {
          Real pres = pa;
          if (rad < rout) {
            if (rad < rin) {
              pres = prat*pa;
            } else {  // add smooth ramp in pressure
              Real f = (rad-rin) / (rout-rin);
              Real log_pres = (1.0-f) * std::log(prat*pa) + f * std::log(pa);
              pres = std::exp(log_pres);
            }
          }
          phydro->u(IEN,k,j,i) = pres/gm1;
          if (RELATIVISTIC_DYNAMICS)  // this should only ever be SR with this file
            phydro->u(IEN,k,j,i) += den;
        }
      }
    }
  }
  return;
}

void MyBoundary_ix1(MeshBlock *pmb, Coordinates *pco,
                    AthenaArray<Real> &prim, FaceField &b,
                    Real time, Real dt,
                    int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  // zero velocity boundary conditions
  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=1; i<=ngh; ++i) {
        prim(IDN,k,j,il-i) = prim(IDN,k,j,il);
        prim(IVX,k,j,il-i) = 0.0;
        prim(IVY,k,j,il-i) = 0.0;
        prim(IVZ,k,j,il-i) = 0.0;
        prim(IPR,k,j,il-i) = prim(IPR,k,j,il);
      }
    }
  }
  return;
}

void MySource(MeshBlock *pmb, const Real time, const Real dt,
              const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
              const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
              AthenaArray<Real> &cons_scalar) {
  SNSource(pmb, time, dt, prim, prim_scalar, bcc, cons, cons_scalar);
  return;
}

namespace {
  void SNSource(MeshBlock *pmb, const Real time, const Real dt,
                const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
                const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
                AthenaArray<Real> &cons_scalar) {
    Real g = pmb->peos->GetGamma();
    Real epsilon = 1.15167;
    Real K = 1.527096;
    Real energy = 4.0/3.0*M_PI*rout*rout*rout*pa*prat/(g-1.0);; 
    Real r_st = epsilon*std::pow(energy,0.2)*std::pow(da,-0.2)*std::pow(time,0.4);
    Real v_st = 0.4*epsilon*std::pow(energy,0.2)*std::pow(da,-0.2)*std::pow(time,-0.6);
    Real p_st = K*energy/std::pow(r_st,3.0)/2.0/M_PI;
    //std::cout << "r_st=" << r_st << " v_st=" << v_st << " p_st=" << p_st << std::endl;
    Real mdot_st = 3.0*da*v_st/r_st;
    Real edot_st = 3.0*da*v_st*v_st*v_st/r_st;
    auto pcoord = pmb->pcoord;
    for (int k = pmb->ks; k <= pmb->ke; ++k) {
      for (int j = pmb->js; j <= pmb->je; ++j) {
        for (int i = pmb->is; i <= pmb->ie; ++i) {
          Real temp = prim(IPR,k,j,i) / prim(IDN,k,j,i);
          if (std::abs(prim(IVX,k,j,i))>1e-4) {
            Real rad;
            if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
              Real x = pcoord->x1v(i);
              Real y = pcoord->x2v(j);
              Real z = pcoord->x3v(k);
              rad = std::sqrt(SQR(x) + SQR(y) + SQR(z));
            } else if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
              Real x = pcoord->x1v(i)*std::cos(pcoord->x2v(j));
              Real y = pcoord->x1v(i)*std::sin(pcoord->x2v(j));
              Real z = pcoord->x3v(k);
              rad = std::sqrt(SQR(x) + SQR(y) + SQR(z));
            } else { // if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0)
              Real x = pcoord->x1v(i)*std::sin(pcoord->x2v(j))*std::cos(pcoord->x3v(k));
              Real y = pcoord->x1v(i)*std::sin(pcoord->x2v(j))*std::sin(pcoord->x3v(k));
              Real z = pcoord->x1v(i)*std::cos(pcoord->x2v(j));
              rad = std::sqrt(SQR(x) + SQR(y) + SQR(z));
            }

            //std::cout << i << std::endl;
            if (time > t_m) {
              if (mdot_flag==0) {
                cons(IDN,k,j,i) += dt * mdot * pow(prim(IPR,k,j,i),5.0/6.0) * 1.8616894573743115;
              } else {
                //cons(IDN,k,j,i) += dt * mdot * mdot_st;
                //cons(IDN,k,j,i) += dt * mdot * mdot_st * rad/r_st; // (rad/r_st-0.5);// * prim(IDN,k,j,i) / 4.0 / da ;;
                cons(IDN,k,j,i) += dt * mdot * mdot_st  * pow(prim(IPR,k,j,i) / p_st, 5.0/6.0);
              }
            }
            if (time > t_e) {
              if (edot_flag==0) {
                cons(IEN,k,j,i) += dt * edot * pow(time,-11.0/5.0) * 1.7787463061055937;
              } else {
                cons(IEN,k,j,i) += dt * edot * edot_st; // * prim(IPR,k,j,i) / p_st;
              }
              
            }
          }
        }
      }
    }
    return;
  }
}

//TODO(@mhguo): add history output: Ehot, Mhot, radius, etc.
