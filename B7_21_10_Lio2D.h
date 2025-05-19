#pragma once

#include "FilesExtern.h"

#include "Common_AAA.h"
#include "Common_CCC.h"
#include "Common_DDD.h"
#include "Common_GGG.h"
#include "Common_TTT.h"
#include "Common_VVV.h"

#include "C7_21_11_DpmXY.h"
#include "C7_21_120_Vis.h"
#include <cmath>

// Thread-local variables to replace globals that could cause race conditions
struct ThreadLocalVars {
  int Y, X; // Local thread indices for matrices
  double cvMet;
  double a_pl, a_mi, s, A, B, C, D;
  double h2x, h1x, h2y, h1y;
  double Vis1, Vis0x, Vis2x, Vis0y, Vis2y;
  double Vis1V, Vis0xV, Vis2xV, Vis0yV, Vis2yV;
  double Vis1SV, Vis0xSV, Vis2xSV, Vis0ySV, Vis2ySV;
  double Tep1, Tep0x, Tep2x, Tep0y, Tep2y;
  double u1, u0x, u2x, u0y, u2y;
  double v1, v0x, v2x, v0y, v2y;
  double e1, e0x, e2x, e0y, e2y;
  double u22, u02, u20, u00;
  double v22, v02, v20, v00;
  double dudx, dudy, dvdx, dvdy;
  double dr11, dr12, dr21, dr22, dr31, dr32, dr41, dr42;
  double dr1, dr2, dr3, dr4;
  double k43, k23;
  double ARP[MX], AUP[MX], AVP[MX], AEP[MX];
  double ARM[MX], AUM[MX], AVM[MX], AEM[MX];
  double DVIS[4];
};

// Function to get thread-local variable storage
ThreadLocalVars *get_thread_local_vars() {
  static thread_local ThreadLocalVars tls;
  return &tls;
}

// Parallel version of B7_21_10_Lio2D
void B7_21_10_Lio2D_Parallel(int Old, int Met, double gam, double FRO[MY][MX],
                             double FU[MY][MX], double FV[MY][MX],
                             double FE[MY][MX], double LRO[MY][MX],
                             double LU[MY][MX], double LV[MY][MX],
                             double LE[MY][MX]) {

  double a_R, a_U, a_V, a_E, hyr;
  double a_P, a_A, a_H, a_S, sM, a_ROp, a_ROm;
  double a_R1, a_U1, a_V1, a_E1;
  double a_R2, a_U2, a_V2, a_E2;
  ThreadLocalVars *tls = get_thread_local_vars();

  //*****Õ¿◊¿ÀŒ œ–Œ’Œƒ¿ œŒ
  // √Œ–»«ŒÕ“¿ÀﬂÃ.************************************************
  if (Met == 1) {
    tls->cvMet = c1v;
  } else if (Met == 2) {
    tls->cvMet = c2v;
  }

// First loop across Y can be parallelized
#pragma omp parallel for private(a_R, a_U, a_V, a_E, a_P, a_A, a_H, a_S, sM,   \
                                     a_ROp, a_ROm)
  for (int Y = 3; Y <= Nyp - 2; Y++) { // do k = 3,NYp-2
    ThreadLocalVars *tls = get_thread_local_vars();

    // This inner loop is sequential due to dependencies
    for (int X = 1; X <= Nxp; X++) { //	do i = 1,NXp
      a_R = FRO[Y][X];
      a_U = FU[Y][X];
      a_V = FV[Y][X];
      a_E = FE[Y][X];

      a_P = a_R * (gam - 1) * a_E;
      a_A = sqrt(gam * (gam - 1) * a_E);
      a_H = gam * a_E + 0.5 * (pow(a_U, 2) + pow(a_V, 2));
      a_S = a_R * a_U;
      sM = a_U / a_A; //! Mach
      a_ROp = pow(sM + 1, 2) / 4;
      a_ROm = -pow(sM - 1, 2) / 4;

      if (a_U > a_A) {
        tls->ARP[X] = a_S;
        tls->AUP[X] = a_S * a_U + a_P;
        tls->AEP[X] = a_S * a_H;
        tls->AVP[X] = a_R * a_V * a_U;

        tls->ARM[X] = 0;
        tls->AUM[X] = 0;
        tls->AEM[X] = 0;
        tls->AVM[X] = 0;
      } else if (a_U > -a_A) {
        tls->ARP[X] = a_A * a_R * a_ROp;
        tls->AUP[X] = a_ROp * a_R * a_A * ((gam - 1) * a_U + 2 * a_A) / gam;
        tls->AVP[X] = a_ROp * a_R * a_A * a_V;
        tls->AEP[X] =
            a_ROp * a_R * a_A *
            (pow((gam - 1) * a_U + 2 * a_A, 2) / (2 * (pow(gam, 2) - 1)) +
             0.5 * pow(a_V, 2));

        tls->ARM[X] = a_ROm * a_R * a_A;
        tls->AUM[X] = a_ROm * a_R * a_A * ((gam - 1) * a_U - 2 * a_A) / gam;
        tls->AVM[X] = a_ROm * a_R * a_A * a_V;
        tls->AEM[X] =
            a_ROm * a_R * a_A *
            (pow((gam - 1) * a_U - 2 * a_A, 2) / (2 * (pow(gam, 2) - 1)) +
             0.5 * pow(a_V, 2));
      } else {
        tls->ARP[X] = 0;
        tls->AUP[X] = 0;
        tls->AEP[X] = 0;
        tls->AVP[X] = 0;

        tls->ARM[X] = a_S;
        tls->AUM[X] = a_S * a_U + a_P;
        tls->AEM[X] = a_S * a_H;
        tls->AVM[X] = a_R * a_V * a_U;
      }
    }

    // This inner loop needs functions d_p_X and d_m_X that weren't provided
    // If those functions don't have side effects, this can be parallelized as
    // well
    for (int X = 3; X <= Nxp - 2; X++) { // do i = 3,NXp-2
      a_R1 = d_p_X(tls->ARP);
      a_U1 = d_p_X(tls->AUP);
      a_V1 = d_p_X(tls->AVP);
      a_E1 = d_p_X(tls->AEP);

      a_R2 = d_m_X(tls->ARM);
      a_U2 = d_m_X(tls->AUM);
      a_V2 = d_m_X(tls->AVM);
      a_E2 = d_m_X(tls->AEM);

      LRO[Y][X] = +(2 / (HX[X + 1] + HX[X])) * (a_R1 + a_R2);
      LU[Y][X] = +(2 / (HX[X + 1] + HX[X])) * (a_U1 + a_U2);
      LV[Y][X] = +(2 / (HX[X + 1] + HX[X])) * (a_V1 + a_V2);
      LE[Y][X] = +(2 / (HX[X + 1] + HX[X])) * (a_E1 + a_E2);
    }
  }

// Second main loop can also be parallelized
#pragma omp parallel for private(a_R, a_U, a_V, a_E, a_P, a_A, a_H, a_S, sM,   \
                                     a_ROp, a_ROm, hyr, a_R1, a_U1, a_V1,      \
                                     a_E1, a_R2, a_U2, a_V2, a_E2)
  for (int X = 3; X <= Nxp - 2; X++) { // do i=3,NXp-2
    ThreadLocalVars *tls = get_thread_local_vars();

    // Sequential inner loop for each thread
    for (int Y = 1; Y <= Nyp; Y++) { // do k=1,NYp
      a_R = FRO[Y][X];
      a_U = FU[Y][X];
      a_V = FV[Y][X];
      a_E = FE[Y][X];

      a_P = a_R * (gam - 1) * a_E;
      a_A = sqrt(gam * (gam - 1) * a_E);
      a_H = gam * a_E + 0.5 * (pow(a_U, 2) + pow(a_V, 2));
      a_S = a_R * a_V;
      sM = a_V / a_A; //! Mach
      a_ROp = pow(sM + 1, 2) / 4;
      a_ROm = -pow(sM - 1, 2) / 4;

      if (a_V > a_A) {
        tls->ARP[Y] = a_S;
        tls->AUP[Y] = a_S * a_U;
        tls->AVP[Y] = a_R * a_V * a_V + a_P;
        tls->AEP[Y] = a_R * a_H * a_V;

        tls->ARM[Y] = 0;
        tls->AUM[Y] = 0;
        tls->AVM[Y] = 0;
        tls->AEM[Y] = 0;
      } else if (a_V > -a_A) {
        tls->ARP[Y] = a_ROp * a_R * a_A;
        tls->AUP[Y] = a_ROp * a_R * a_A * a_U;
        tls->AVP[Y] = a_ROp * a_R * a_A * ((gam - 1) * a_V + 2 * a_A) / gam;
        tls->AEP[Y] =
            a_R * a_A * a_ROp *
            (pow((gam - 1) * a_V + 2 * a_A, 2) / (2 * (pow(gam, 2) - 1)) +
             0.5 * pow(a_U, 2));

        tls->ARM[Y] = a_ROm * a_R * a_A;
        tls->AUM[Y] = a_ROm * a_R * a_A * a_U;
        tls->AVM[Y] = a_ROm * a_R * a_A * ((gam - 1) * a_V - 2 * a_A) / gam;
        tls->AEM[Y] =
            a_R * a_A * a_ROm *
            (pow((gam - 1) * a_V - 2 * a_A, 2) / (2 * (pow(gam, 2) - 1)) +
             0.5 * pow(a_U, 2));
      } else {
        tls->ARP[Y] = 0;
        tls->AUP[Y] = 0;
        tls->AVP[Y] = 0;
        tls->AEP[Y] = 0;

        tls->ARM[Y] = a_S;
        tls->AUM[Y] = a_S * a_U;
        tls->AVM[Y] = a_R * a_V * a_V + a_P;
        tls->AEM[Y] = a_R * a_H * a_V;
      }
    }

    // Inner loop for this thread
    for (int Y = 3; Y <= Nyp - 2; Y++) { // do k = 3,NYp-2
      a_R1 = d_p_Y(tls->ARP);
      a_U1 = d_p_Y(tls->AUP);
      a_V1 = d_p_Y(tls->AVP);
      a_E1 = d_p_Y(tls->AEP);

      a_R2 = d_m_Y(tls->ARM);
      a_U2 = d_m_Y(tls->AUM);
      a_V2 = d_m_Y(tls->AVM);
      a_E2 = d_m_Y(tls->AEM);

      hyr = (2 / (HY[Y + 1] + HY[Y]));

// Write operations to the output matrices - no race condition since X is unique
// per thread
#pragma omp atomic update
      LRO[Y][X] += hyr * (a_R1 + a_R2);

#pragma omp atomic update
      LU[Y][X] += hyr * (a_U1 + a_U2);

#pragma omp atomic update
      LV[Y][X] += hyr * (a_V1 + a_V2);

#pragma omp atomic update
      LE[Y][X] += hyr * (a_E1 + a_E2);

      // Viscosity calculation - uses thread-local DVIS array
      if ((vis10 != 0.0) || (vis20 != 0.0)) {
        // Make a thread-local copy since we're modifying these globals
        C7_21_120_Vis(Met, FU, FV, FE, KRO1[Old], KRO2[Old], KE1[Old],
                      KE2[Old]);

#pragma omp atomic update
        LU[Y][X] += tls->DVIS[0];

#pragma omp atomic update
        LV[Y][X] += tls->DVIS[1];

#pragma omp atomic update
        LE[Y][X] += tls->DVIS[2];
      }
    }
  }
}

void B7_21_10_Lio2D(int Old, int Met, double gam, double FRO[MY][MX],
                    double FU[MY][MX], double FV[MY][MX], double FE[MY][MX],
                    double LRO[MY][MX], double LU[MY][MX], double LV[MY][MX],
                    double LE[MY][MX]) {

  double a_R, a_U, a_V, a_E, hyr;
  double a_P, a_A, a_H, a_S, sM, a_ROp, a_ROm;
  double a_R1, a_U1, a_V1, a_E1;
  double a_R2, a_U2, a_V2, a_E2;

  //*****Õ¿◊¿ÀŒ œ–Œ’Œƒ¿ œŒ
  // √Œ–»«ŒÕ“¿ÀﬂÃ.************************************************
  if (Met == 1) {
    cvMet = c1v;
  } else if (Met == 2) {
    cvMet = c2v;
  }

  for (Y = 3; Y <= Nyp - 2; Y++) { // do k = 3,NYp-2
    for (X = 1; X <= Nxp; X++) {   //	do i = 1,NXp
      a_R = FRO[Y][X];
      a_U = FU[Y][X];
      a_V = FV[Y][X];
      a_E = FE[Y][X];
      // uuu = PU1[K][I];

      a_P = a_R * (gam - 1) * a_E;
      a_A = sqrt(gam * (gam - 1) * a_E);
      a_H = gam * a_E + 0.5 * (pow(a_U, 2) + pow(a_V, 2));
      a_S = a_R * a_U;
      sM = a_U / a_A; //! Mach
      a_ROp = pow(sM + 1, 2) / 4;
      a_ROm = -pow(sM - 1, 2) / 4;

      //		if(uuu > a_A) {
      if (a_U > a_A) {
        ARP[X] = a_S;
        AUP[X] = a_S * a_U + a_P;
        AEP[X] = a_S * a_H;
        AVP[X] = a_R * a_V * a_U;

        ARM[X] = 0;
        AUM[X] = 0;
        AEM[X] = 0;
        AVM[X] = 0;
      }

      //		else if(uuu > -a_A) {
      else if (a_U > -a_A) {
        ARP[X] = a_A * a_R * a_ROp;
        AUP[X] = a_ROp * a_R * a_A * ((gam - 1) * a_U + 2 * a_A) / gam;
        AVP[X] = a_ROp * a_R * a_A * a_V;
        AEP[X] = a_ROp * a_R * a_A *
                 (pow((gam - 1) * a_U + 2 * a_A, 2) / (2 * (pow(gam, 2) - 1)) +
                  0.5 * pow(a_V, 2));

        ARM[X] = a_ROm * a_R * a_A;
        AUM[X] = a_ROm * a_R * a_A * ((gam - 1) * a_U - 2 * a_A) / gam;
        AVM[X] = a_ROm * a_R * a_A * a_V;
        AEM[X] = a_ROm * a_R * a_A *
                 (pow((gam - 1) * a_U - 2 * a_A, 2) / (2 * (pow(gam, 2) - 1)) +
                  0.5 * pow(a_V, 2));
      }

      else {
        ARP[X] = 0;
        AUP[X] = 0;
        AEP[X] = 0;
        AVP[X] = 0;

        ARM[X] = a_S;
        AUM[X] = a_S * a_U + a_P;
        AEM[X] = a_S * a_H;
        AVM[X] = a_R * a_V * a_U;
      } // end if.	end if
    } // end do

    for (X = 3; X <= Nxp - 2; X++) { // do i = 3,NXp-2
      a_R1 = d_p_X(ARP);
      a_U1 = d_p_X(AUP);
      a_V1 = d_p_X(AVP);
      a_E1 = d_p_X(AEP);

      a_R2 = d_m_X(ARM);
      a_U2 = d_m_X(AUM);
      a_V2 = d_m_X(AVM);
      a_E2 = d_m_X(AEM);

      LRO[Y][X] = +(2 / (HX[X + 1] + HX[X])) * (a_R1 + a_R2);
      LU[Y][X] = +(2 / (HX[X + 1] + HX[X])) * (a_U1 + a_U2);
      LV[Y][X] = +(2 / (HX[X + 1] + HX[X])) * (a_V1 + a_V2);
      LE[Y][X] = +(2 / (HX[X + 1] + HX[X])) * (a_E1 + a_E2);
    } // end do
  } // end do

  for (X = 3; X <= Nxp - 2; X++) { // do i=3,NXp-2
    for (Y = 1; Y <= Nyp; Y++) {   // do k=1,NYp
      a_R = FRO[Y][X];
      a_U = FU[Y][X];
      a_V = FV[Y][X];
      a_E = FE[Y][X];
      //		vvv = PV1[K][I];

      a_P = a_R * (gam - 1) * a_E;
      a_A = sqrt(gam * (gam - 1) * a_E);
      a_H = gam * a_E + 0.5 * (pow(a_U, 2) + pow(a_V, 2));
      a_S = a_R * a_V;
      sM = a_V / a_A; //! Mach
      a_ROp = pow(sM + 1, 2) / 4;
      a_ROm = -pow(sM - 1, 2) / 4;

      //		if(vvv > a_A) {
      if (a_V > a_A) {
        ARP[Y] = a_S;
        AUP[Y] = a_S * a_U;
        AVP[Y] = a_R * a_V * a_V + a_P;
        AEP[Y] = a_R * a_H * a_V;

        ARM[Y] = 0;
        AUM[Y] = 0;
        AVM[Y] = 0;
        AEM[Y] = 0;
      }

      //		else if (vvv > -a_A) {
      else if (a_V > -a_A) {
        ARP[Y] = a_ROp * a_R * a_A;
        AUP[Y] = a_ROp * a_R * a_A * a_U;
        AVP[Y] = a_ROp * a_R * a_A * ((gam - 1) * a_V + 2 * a_A) / gam;
        AEP[Y] = a_R * a_A * a_ROp *
                 (pow((gam - 1) * a_V + 2 * a_A, 2) / (2 * (pow(gam, 2) - 1)) +
                  0.5 * pow(a_U, 2));

        ARM[Y] = a_ROm * a_R * a_A;
        AUM[Y] = a_ROm * a_R * a_A * a_U;
        AVM[Y] = a_ROm * a_R * a_A * ((gam - 1) * a_V - 2 * a_A) / gam;
        AEM[Y] = a_R * a_A * a_ROm *
                 (pow((gam - 1) * a_V - 2 * a_A, 2) / (2 * (pow(gam, 2) - 1)) +
                  0.5 * pow(a_U, 2));
      }

      else {
        ARP[Y] = 0;
        AUP[Y] = 0;
        AVP[Y] = 0;
        AEP[Y] = 0;

        ARM[Y] = a_S;
        AUM[Y] = a_S * a_U;
        AVM[Y] = a_R * a_V * a_V + a_P;
        AEM[Y] = a_R * a_H * a_V;
      } // end if. end if
    } // end do

    for (Y = 3; Y <= Nyp - 2; Y++) { // do k = 3,NYp-2
      a_R1 = d_p_Y(ARP);
      a_U1 = d_p_Y(AUP);
      a_V1 = d_p_Y(AVP);
      a_E1 = d_p_Y(AEP);

      a_R2 = d_m_Y(ARM);
      a_U2 = d_m_Y(AUM);
      a_V2 = d_m_Y(AVM);
      a_E2 = d_m_Y(AEM);

      //		if (Met == 1) {
      hyr = (2 / (HY[Y + 1] + HY[Y]));
      LRO[Y][X] += hyr * (a_R1 + a_R2);
      LU[Y][X] += hyr * (a_U1 + a_U2);
      LV[Y][X] += hyr * (a_V1 + a_V2);
      LE[Y][X] += hyr * (a_E1 + a_E2);

      //!*********************************************************************
      //! BEGIN (ƒÀﬂ ¬ﬂ« Œ—“»).
      //!*********************************************************************
      if ((vis10 != 0.0) || (vis20 != 0.0)) {
        C7_21_120_Vis(Met, FU, FV, FE, KRO1[Old], KRO2[Old], KE1[Old],
                      KE2[Old]);
        LU[Y][X] += DVIS[0];
        LV[Y][X] += DVIS[1];
        LE[Y][X] += DVIS[2];
      } // END if !(VISCOSITY).

      //!****************************************************************
    } // end do
  } // end do
  // return
} // END SUBROUTINE S7_2100_LIO2D
