#pragma once

#include "Ret7_21_210_AK.h"
#include "B7_21_22_BndCnd.h"

//void M7_4_0_Modeling(double Tau, double om, double c1v, double c2v) {
void B7_21_20_PKtau(double Tau,
	double OldRO1[MY][MX],double OldU1[MY][MX],double OldV1[MY][MX],double OldE1[MY][MX],
	double OldRO2[MY][MX],double OldU2[MY][MX],double OldV2[MY][MX],double OldE2[MY][MX],
	double NewRO1[MY][MX],double NewU1[MY][MX],double NewV1[MY][MX],double NewE1[MY][MX],
	double NewRO2[MY][MX],double NewU2[MY][MX],double NewV2[MY][MX],double NewE2[MY][MX]) {

//M7_200_Lio2D(1,K1RO1,K1RO2,K1U1,K1U2,K1V1,K1V2,K1E1,K1E2); //Расчёт нового слоя "Первого" газа.
//M7_200_Lio2D(2,K1RO1,K1RO2,K1U1,K1U2,K1V1,K1V2,K1E1,K1E2); //Расчёт нового слоя "Второго" газа.

//M7_3100_Lio2D1(1); //Расчёт нового слоя "Первого" газа.
//M7_3100_Lio2D1(2); //Расчёт нового слоя "Второго" газа.

double a_U, a_V, a_E, vu, vun;
double c1UN, c1VN, e1NN;
double c2UN, c2VN, e2NN;
double f1, f2, q1, q2, b1, b2;
double delt, delt1, delt2;
double ak1, beta, p1N, p2N;

for(Y=3; Y<=Nyp-2; Y++) { //DO K = 3,Nyp-2
	for(X=3; X<=Nxp-2; X++) { //DO I = 3,Nxp-2

		NewRO1[Y][X] = RO1[Y][X] - Tau * LRO1[Y][X];
		NewRO2[Y][X] = RO2[Y][X] - Tau * LRO2[Y][X];

ak1 = F7_21_210_AK(OldRO1[Y][X],OldRO2[Y][X],OldU1[Y][X],OldU2[Y][X],
				   OldE1[Y][X]/c1v,OldE2[Y][X]/c2v);

		beta = OldE1[Y][X]/c1v/(OldE1[Y][X]/c1v+1/om*OldE2[Y][X]/c2v);

		p1N = NewRO1[Y][X];
		p2N = NewRO2[Y][X];
//!			write(*,*)p1N,p2N

		if(p1N/(p1N+p2N/om) < concmin) {
			a_U = (RO2[Y][X]*U2[Y][X]-Tau*LU2[Y][X]) / p2N;
			a_V = (RO2[Y][X]*V2[Y][X]-Tau*LV2[Y][X]) / p2N;
			vu  = pow(U2[Y][X],2) + pow(V2[Y][X],2);
			vun = pow(a_U,2)+pow(a_V,2);
			a_E = (RO2[Y][X]*(E2[Y][X]+vu/2)-Tau*LE2[Y][X]) / p2N - vun/2.0;

			c2UN = a_U;
			c2VN = a_V;
			e2NN = a_E;

			c1UN = c2UN;
			c1VN = c2VN;
			e1NN = c1v/c2v*e2NN;
//c				write(*,*)'1',e1NN,e2NN
		}

		else if(p2N/(p1N+p2N) < concmin) {
			a_U = (RO1[Y][X]*U1[Y][X]-Tau*LU1[Y][X]) / p1N;
			a_V = (RO1[Y][X]*V1[Y][X]-Tau*LV1[Y][X]) / p1N;
			vu  = pow(U1[Y][X],2) + pow(V1[Y][X],2);
			vun = pow(a_U,2) + pow(a_V,2);
a_E = (RO1[Y][X]*(E1[Y][X]+vu/2)-Tau*LE1[Y][X]) / p1N - vun/2.0;

			c1UN = a_U;
			c1VN = a_V;
			e1NN = a_E;

			c2UN = c1UN;
			c2VN = c1VN;
			e2NN = c2v/c1v*e1NN;
		}

		else {
			f1 = -LU1[Y][X];
			f2 = -LU2[Y][X];

delt = (p1N+ak1*Tau) * (p2N+ak1*Tau) - pow(ak1*Tau,2);

delt1 = (RO1[Y][X]*U1[Y][X]+Tau*f1) * (p2N+ak1*Tau) +
	 ak1*Tau * (RO2[Y][X]*U2[Y][X]+f2*Tau);

delt2 = (RO2[Y][X]*U2[Y][X]+Tau*f2) * (p1N+ak1*Tau) +
	 ak1*Tau * (RO1[Y][X]*U1[Y][X]+f1*Tau);

			c1UN = delt1/delt;
			c2UN = delt2/delt;

			f1 = -LV1[Y][X];
			f2 = -LV2[Y][X];

delt  = (p1N+ak1*Tau) * (p2N+ak1*Tau) - pow(ak1*Tau,2);

delt1 = (RO1[Y][X]*V1[Y][X]+Tau*f1) * (p2N+ak1*Tau) +
	 ak1*Tau * (RO2[Y][X]*V2[Y][X]+f2*Tau);

delt2 = (RO2[Y][X]*V2[Y][X]+Tau*f2) * (p1N+ak1*Tau) +
	 ak1*Tau * (RO1[Y][X]*V1[Y][X]+f1*Tau);

			c1VN = delt1/delt;
			c2VN = delt2/delt;

q1 = -LE1[Y][X] + ak1*c1UN*(c2UN-c1UN) + ak1*beta*pow(c1UN-c2UN,2) +
	 ak1*c1VN*(c2VN-c1VN) + ak1*beta*pow(c1VN-c2VN,2);

q2 = -LE2[Y][X] + ak1*c2UN*(c1UN-c2UN) + ak1*(1-beta)*pow(c1UN-c2UN,2) +
	 ak1*c2VN*(c1VN-c2VN) + ak1*(1-beta)*pow(c1VN-c2VN,2);

delt = (p1N+qk*ak1*Tau/c1v) * (p2N+qk*ak1*Tau/c2v) -
	 pow(qk*ak1*Tau,2) / (c1v*c2v);

b1 = Tau*q1 + RO1[Y][X] * (E1[Y][X]+(pow(U1[Y][X],2)+pow(V1[Y][X],2))/2) -
	 p1N * (pow(c1UN,2)+pow(c1VN,2)) / 2;

b2 = Tau*q2 + RO2[Y][X] * (E2[Y][X]+(pow(U2[Y][X],2)+pow(V2[Y][X],2))/2) -
	 p2N * (pow(c2UN,2)+pow(c2VN,2)) / 2;

delt1 = b1*(p2N+qk*ak1*Tau/c2v) + qk*ak1*Tau/c2v*b2;
delt2 = b2*(p1N+qk*ak1*Tau/c1v) + qk*ak1*Tau/c1v*b1;

			e1NN = delt1/delt;
			e2NN = delt2/delt;
		} //end if


		NewU1[Y][X] = c1UN;
		NewU2[Y][X] = c2UN;
		NewV2[Y][X] = c2VN;
		NewV1[Y][X] = c1VN;
		NewE1[Y][X] = e1NN;
		NewE2[Y][X] = e2NN;
	} //ENDDO
} //enddo

B7_21_22_BndCnd(NewRO1,NewU1,NewV1,NewE1,  NewRO2,NewU2,NewV2,NewE2); //Первый (из двух) шаг Рунге-Кутта.


} //END SUBROUTINE S7_3000_MODELING



