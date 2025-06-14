#pragma once

#include "C11_DpmXY.h"
#include "C120_Vis.h"

void C10_Lio2D(int Met, double GRO[MY][MX], double GU[MY][MX],
               double GV[MY][MX], double GE[MY][MX], double LRO[MY][MX],
               double LU[MY][MX], double LV[MY][MX], double LE[MY][MX]) {

  double a_R, a_U, a_V, a_E, hyr, gam;
  double a_P, a_A, a_H, a_S, sM, a_ROp, a_ROm;
  double a_R1, a_U1, a_V1, a_E1;
  double a_R2, a_U2, a_V2, a_E2;

  //*****Õ¿◊¿ÀŒ œ–Œ’Œƒ¿ œŒ
  //√Œ–»«ŒÕ“¿ÀﬂÃ.************************************************
  if (Met == 1) {
    cvMet = c1v;
    gam = gam1;
  } else if (Met == 2) {
    cvMet = c2v;
    gam = gam2;
  }

  for (Y = 3; Y <= Nyp - 2; Y++) { // do k = 3,NYp-2
    for (X = 1; X <= Nxp; X++) {   //	do i = 1,NXp
      a_R = GRO[Y][X];
      a_U = GU[Y][X];
      a_V = GV[Y][X];
      a_E = GE[Y][X];
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
      a_R = GRO[Y][X];
      a_U = GU[Y][X];
      a_V = GV[Y][X];
      a_E = GE[Y][X];
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
        C120_Vis(Met, GU, GV, GE, GRO1, GRO2, GE1, GE2);
        LU[Y][X] += DVIS[0];
        LV[Y][X] += DVIS[1];
        LE[Y][X] += DVIS[2];
      } // END if !(VISCOSITY).

      //!****************************************************************
    } // end do
  } // end do
  // return
} // END SUBROUTINE S7_2100_LIO2D
