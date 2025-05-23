#pragma once

#include "FilesExtern.h"

#include "Common_AAA.h"
#include "Common_CCC.h"
#include "Common_DDD.h"
#include "Common_GGG.h"
#include "Common_TTT.h"
#include "Common_VVV.h"

#include "A7_21_0_RK_Step.h"
#include "A7_22_RK_LayerUp.h"

void A7_20_RungeKutt(double Tau) {
  if (TypeTau == '5') {
    A7_21_0_RK_Step(1, Tau / 4);
    A7_21_0_RK_Step(2, Tau / 6);
    A7_21_0_RK_Step(3, 3.0 * Tau / 8);
    A7_21_0_RK_Step(4, Tau / 2);
    A7_21_0_RK_Step(5, Tau);
  }

  else if (TypeTau == '3') {
    A7_21_0_RK_Step(1, Tau / 2);
    A7_21_0_RK_Step(2, Tau / 2);
    A7_21_0_RK_Step(3, Tau);
  }

  else if (TypeTau == '2') {
    A7_21_0_RK_Step(1, Tau / 2);
    A7_21_0_RK_Step(2, Tau);
  }

  else if (TypeTau == '1') {
    A7_21_0_RK_Step(1, Tau);
  }

  else {
    cout << "Nepravilno ukazan tip shemy!!!" << endl;
  }

  A7_2_2_RK_LayerUp();
}

void A7_20_RungeKutt_Parallel(double Tau) {
  if (TypeTau == '5') {
    A7_21_0_RK_Step_Parallel(1, Tau / 4);
    A7_21_0_RK_Step_Parallel(2, Tau / 6);
    A7_21_0_RK_Step_Parallel(3, 3.0 * Tau / 8);
    A7_21_0_RK_Step_Parallel(4, Tau / 2);
    A7_21_0_RK_Step_Parallel(5, Tau);
  } else if (TypeTau == '3') {
    A7_21_0_RK_Step_Parallel(1, Tau / 2);
    A7_21_0_RK_Step_Parallel(2, Tau / 2);
    A7_21_0_RK_Step_Parallel(3, Tau);
  } else if (TypeTau == '2') {
    A7_21_0_RK_Step_Parallel(1, Tau / 2);
    A7_21_0_RK_Step_Parallel(2, Tau);
  } else if (TypeTau == '1') {
    A7_21_0_RK_Step_Parallel(1, Tau);
  } else {
    cout << "Nepravilno ukazan tip shemy!!!" << endl;
  }
  A7_2_2_RK_LayerUp();
}
