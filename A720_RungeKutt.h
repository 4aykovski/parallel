#pragma once
#include "B0_A721_RK_Step.h"
// #include "B22_RK_LayerUp.h"
//
//

void A722_RK_LayerUp_Parallel();

// Параллельная версия A720_RungeKutt
void A720_RungeKutt_Parallel(double TAU) {
  int x, y;
  if (TypeTau == '1') {
    B0_A721_RK_Step_Parallel(1, TAU);
  } else if (TypeTau == '2') {
    B0_A721_RK_Step_Parallel(1, TAU / 2);
    B0_A721_RK_Step_Parallel(2, TAU);
  } else if (TypeTau == '3') {
    B0_A721_RK_Step_Parallel(1, TAU / 2);
    B0_A721_RK_Step_Parallel(2, TAU / 2);
    B0_A721_RK_Step_Parallel(3, TAU / 3);
  } else {
    cout << "The schema type is incorrectly specified!!!" << endl;
  }
  A722_RK_LayerUp_Parallel();
}

// Параллельная версия A722_RK_LayerUp
void A722_RK_LayerUp_Parallel() {
  int Nxt, x, y;

#pragma omp parallel for private(x)
  for (y = 1; y <= Nyp; y++) {
    for (x = 1; x <= Nxp; x++) {
      RO1[y][x] = GRO1[y][x];
      U1[y][x] = GU1[y][x];
      V1[y][x] = GV1[y][x];
      E1[y][x] = GE1[y][x];
      RO2[y][x] = GRO2[y][x];
      U2[y][x] = GU2[y][x];
      V2[y][x] = GV2[y][x];
      E2[y][x] = GE2[y][x];
    }
  }
}

void A722_RK_LayerUp();

void A720_RungeKutt(double TAU) {
  int x, y;

  if (TypeTau == '1') {
    B0_A721_RK_Step(1, TAU);
  }

  else if (TypeTau == '2') {
    B0_A721_RK_Step(1, TAU / 2);
    B0_A721_RK_Step(2, TAU);
  }

  else if (TypeTau == '3') {
    B0_A721_RK_Step(1, TAU / 2);
    B0_A721_RK_Step(2, TAU / 2);
    B0_A721_RK_Step(3, TAU / 3);
  }

  else {
    cout << "The schema type is incorrectly specified!!!" << endl;
  }

  A722_RK_LayerUp();
}

void A722_RK_LayerUp() {
  int Nxt, x, y;
  // if (TypeTau == '2') Nxt = 0;
  // else                Nxt = 1;
  for (y = 1; y <= Nyp; y++) {   // DO K = 1,NYp
    for (x = 1; x <= Nxp; x++) { // DO I = 1,NXp
      RO1[y][x] = GRO1[y][x];
      U1[y][x] = GU1[y][x];
      V1[y][x] = GV1[y][x];
      E1[y][x] = GE1[y][x];
      RO2[y][x] = GRO2[y][x];
      U2[y][x] = GU2[y][x];
      V2[y][x] = GV2[y][x];
      E2[y][x] = GE2[y][x];
    } // ENDDO
  } // ENDDO
}
