#pragma once

double F7_1_Tau() {
double am, am1, am2, am3, am4, aa;
//int k,i;

	am  = 0.0;
	am1 = 0.0;
	am2 = 0.0;
	am3 = 0.0;
	am4 = 0.0;

for(Y=1; Y<=Nyp; Y++) { //DO K = 1,NYp
	for(X=1; X<=Nxp; X++) { //do  i = 1,nXp
//c		aa=	abs(U1[K][I])+sqrt(gam1*(gam1-1)*E1[K][I])
//!		aa=abs(U1[K][I])/(2*HX[I])+abs(V1[K][I])/(2*HY[K])+
//!     *sqrt(gam1*(gam1-1)*E1[K][I])*sqrt(4/HX[I]**2+4/HY[K]**2)

aa = abs(U1[Y][X])/HX[X] + abs(V1[Y][X])/HY[Y] +
	 sqrt(gam1*(gam1-1)*E1[Y][X]) * sqrt(1.0/pow(HX[X],2)+1.0/pow(HY[Y],2));
	
		if(aa > am1) {
			am1 = aa;
		} //end if

aa = abs(U2[Y][X])/HX[X] + abs(V2[Y][X])/HY[Y] +
	 sqrt(gam2*(gam2-1)*E2[Y][X]) * sqrt(1.0/pow(HX[X],2)+1.0/pow(HY[Y],2));

			if(aa > am3) {
				am3 = aa;
			} //end if

		} //ENDDO
	} //ENDDO

	am = max(am1,am3);

	return Akur/am;


} //END F7_1_TaU1().


