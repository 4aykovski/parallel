#pragma once

#include "B3_BndCnd.h"
#include "C10_Lio2D.h"
#include "D20_PKtau.h"

void B0_A721_RK_Step(int Step, double Tau) {
  //	int Km = K - 1,
  int x, y;
  double Save;
  int Old, Nxt;

  C10_Lio2D(1, GRO1, GU1, GV1, GE1, LRO1, LU1, LV1,
            LE1); // Расчёт нового слоя "Первого" газа.
  C10_Lio2D(2, GRO2, GU2, GV2, GE2, LRO2, LU2, LV2,
            LE2); // Расчёт нового слоя "Второго" газа.

  if ((Step == 1) or (TypeTau == '2' and Step == 2)) {
    for (y = 1; y <= Nyp; y++) {   // do k = 1,nyp
      for (x = 1; x <= Nxp; x++) { // do i = 1,nxp
        KRO1[y][x] = LRO1[y][x];
        KRO2[y][x] = LRO2[y][x];
        KU1[y][x] = LU1[y][x];
        KU2[y][x] = LU2[y][x];
        KV1[y][x] = LV1[y][x];
        KV2[y][x] = LV2[y][x];
        KE1[y][x] = LE1[y][x];
        KE2[y][x] = LE2[y][x];
      } // enddo
    } // enddo
  } // endif

  else {
    for (y = 1; y <= Nyp; y++) {   // do k = 1,nyp
      for (x = 1; x <= Nxp; x++) { // do i = 1,nxp
        KRO1[y][x] += LRO1[y][x];
        KRO2[y][x] += LRO2[y][x];
        KU1[y][x] += LU1[y][x];
        KU2[y][x] += LU2[y][x];
        KV1[y][x] += LV1[y][x];
        KV2[y][x] += LV2[y][x];
        KE1[y][x] += LE1[y][x];
        KE2[y][x] += LE2[y][x];
      } // enddo
    } // enddo
  } // endelse

  D20_PKtau(Tau, RO1, U1, V1, E1, RO2, U2, V2, E2, GRO1, GU1, GV1, GE1, GRO2,
            GU2, GV2, GE2); // Первый (из двух) шаг Рунге-Кутта.

  B3_BndCnd();
}
