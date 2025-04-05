#pragma once

#include "B7_21_10_Lio2D.h"
#include "B7_21_20_PKtau.h"

void A7_21_0_RK_Step(int Step, double Tau) {
//	int Km = K - 1, 
int x, y;
double Save;
int Old, Nxt;

if (Step == 1) {
	B7_21_10_Lio2D(0,1,gam1, RO1,U1,V1,E1,  LRO1,LU1,LV1,LE1); //Расчёт нового слоя "Первого" газа.
	B7_21_10_Lio2D(0,2,gam2, RO2,U2,V2,E2,  LRO2,LU2,LV2,LE2); //Расчёт нового слоя "Второго" газа.
	
	B7_21_20_PKtau(Tau,      RO1,U1,V1,E1,  RO2,U2,V2,E2,
							 KRO1[1],KU1[1],KV1[1],KE1[1],  
							 KRO2[1],KU2[1],KV2[1],KE2[1]); //Первый (из двух) шаг Рунге-Кутта.
}
	
else {
	Old = (Step+1) % 2;
	Nxt = Step % 2;
	
	B7_21_10_Lio2D(Old,1,gam1, KRO1[Old],KU1[Old],KV1[Old],KE1[Old],  LRO1,LU1,LV1,LE1); //Расчёт нового слоя "Первого" газа.
	B7_21_10_Lio2D(Old,2,gam2, KRO2[Old],KU2[Old],KV2[Old],KE2[Old],  LRO2,LU2,LV2,LE2); //Расчёт нового слоя "Второго" газа.
	
	B7_21_20_PKtau(Tau,      
		KRO1[Old],KU1[Old],KV1[Old],KE1[Old],  KRO2[Old],KU2[Old],KV2[Old],KE2[Old],
		KRO1[Nxt],KU1[Nxt],KV1[Nxt],KE1[Nxt],  KRO2[Nxt],KU2[Nxt],KV2[Nxt],KE2[Nxt]); //Первый (из двух) шаг Рунге-Кутта.
}

}

