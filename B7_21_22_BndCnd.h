#pragma once

//void M7_4_0_Modeling(double Tau, double om, double c1v, double c2v) {
void B7_21_22_BndCnd(
	double CRO1[MY][MX],double CU1[MY][MX],double CV1[MY][MX],double CE1[MY][MX],
	double CRO2[MY][MX],double CU2[MY][MX],double CV2[MY][MX],double CE2[MY][MX]) {

//!	DO K=3,Nyp-2
for(Y=1; Y<=Nyp; Y++) { //DO K = 1,Nyp
	if(BoundCondLeft == 1) {
		CRO1[Y][1]= RO1[Y][1]; CRO1[Y][2]= RO1[Y][2];
		CU1[Y][1] = U1[Y][1];  CU1[Y][2] = U1[Y][2];
		CV1[Y][1] = V1[Y][1];  CV1[Y][2] = V1[Y][2];
		CE1[Y][1] = E1[Y][1];  CE1[Y][2] = E1[Y][2];
		
		CRO2[Y][1]= RO2[Y][1]; CRO2[Y][2]= RO2[Y][2];
		CU2[Y][1] = U2[Y][1];  CU2[Y][2] = U2[Y][2];
		CV2[Y][1] = V2[Y][1];  CV2[Y][2] = V2[Y][2];
		CE2[Y][1] = E2[Y][1];  CE2[Y][2] = E2[Y][2];
	} //end if
	else if(BoundCondLeft == 0) {
		CRO1[Y][1]= CRO1[Y][5]; CRO1[Y][2]= CRO1[Y][4];
		CU1[Y][1] = -CU1[Y][5]; CU1[Y][2] =-CU1[Y][4];
		CV1[Y][1] = CV1[Y][5];  CV1[Y][2] = CV1[Y][4];
		CE1[Y][1] = CE1[Y][5];  CE1[Y][2] = CE1[Y][4];
//c		U1N[K][2)=0
		CRO2[Y][1]= CRO2[Y][5]; CRO2[Y][2]= CRO2[Y][4];
		CU2[Y][1] = -CU2[Y][5]; CU2[Y][2] =-CU2[Y][4];
		CV2[Y][1] = CV2[Y][5];  CV2[Y][2] = CV2[Y][4];
		CE2[Y][1] = CE2[Y][5];  CE2[Y][2] = CE2[Y][4];
//c		U2N[K][2)=0
	}


	if(BoundCondRight == 1) { // ÂÍÈÌÀÍÈÅ - ÐÀÍÜØÅ ÁÛË Npr
		CRO1[Y][Nxp]= CRO1[Y][Nxp-1]= CRO1[Y][Nxp-2];
		CU1[Y][Nxp] = CU1[Y][Nxp-1] = CU1[Y][Nxp-2];
		CV1[Y][Nxp] = CV1[Y][Nxp-1] = CV1[Y][Nxp-2];
		CE1[Y][Nxp] = CE1[Y][Nxp-1] = CE1[Y][Nxp-2];
		
		CRO2[Y][Nxp]= CRO2[Y][Nxp-1]= CRO2[Y][Nxp-2];
		CU2[Y][Nxp] = CU2[Y][Nxp-1] = CU2[Y][Nxp-2];
		CV2[Y][Nxp] = CV2[Y][Nxp-1] = CV2[Y][Nxp-2];
		CE2[Y][Nxp] = CE2[Y][Nxp-1] = CE2[Y][Nxp-2];
	}
	else if(BoundCondRight == 2) {
		CRO1[Y][Nxp]= CRO1[Y][Nxp-4]; CRO1[Y][Nxp-1]= CRO1[Y][Nxp-3];  
		CU1[Y][Nxp] =-CU1[Y][Nxp-4];  CU1[Y][Nxp-1] =-CU1[Y][Nxp-3];   
		CV1[Y][Nxp] = CV1[Y][Nxp-4];  CV1[Y][Nxp-1] = CV1[Y][Nxp-3];   
		CE1[Y][Nxp] = CE1[Y][Nxp-4];  CE1[Y][Nxp-1] = CE1[Y][Nxp-3];   
		CRO2[Y][Nxp]= CRO2[Y][Nxp-4]; CRO2[Y][Nxp-1]= CRO2[Y][Nxp-3];  
		CU2[Y][Nxp] =-CU2[Y][Nxp-4];  CU2[Y][Nxp-1] =-CU2[Y][Nxp-3];   
		CV2[Y][Nxp] = CV2[Y][Nxp-4];  CV2[Y][Nxp-1] = CV2[Y][Nxp-3];   
		CE2[Y][Nxp] = CE2[Y][Nxp-4];  CE2[Y][Nxp-1] = CE2[Y][Nxp-3];   
	} //end if
} //ENDDO


//c******************************************************
//c**************************************************
for(X=1; X<=Nxp; X++) { //DO I = 1,Nxp
	CRO1[1][X]= CRO1[5][X]; CRO1[2][X]= CRO1[4][X];
	CU1[1][X] = CU1[5][X];  CU1[2][X] = CU1[4][X];
	CV1[1][X] =-CV1[5][X];  CV1[2][X] =-CV1[4][X]; CV1[3][X] = 0;
	CE1[1][X] = CE1[5][X];  CE1[2][X] = CE1[4][X];

	CRO2[1][X]= CRO2[5][X]; CRO2[2][X]= CRO2[4][X];
	CU2[1][X] = CU2[5][X];  CU2[2][X] = CU2[4][X];
	CV2[1][X] =-CV2[5][X];  CV2[2][X] =-CV2[4][X]; CV2[3][X] = 0;
	CE2[1][X] = CE2[5][X];  CE2[2][X] = CE2[4][X];
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	CRO1[Nyp][X]= CRO1[Nyp-1][X]= CRO1[Nyp-2][X];
	CU1[Nyp][X] = CU1[Nyp-1][X] = CU1[Nyp-2][X];
	CV1[Nyp][X] = CV1[Nyp-1][X] = CV1[Nyp-2][X];
	CE1[Nyp][X] = CE1[Nyp-1][X] = CE1[Nyp-2][X];
//c		V1N(NY-1][I]=0
	CRO2[Nyp][X]= CRO2[Nyp-1][X]= CRO2[Nyp-2][X];
	CU2[Nyp][X] = CU2[Nyp-1][X] = CU2[Nyp-2][X];
	CV2[Nyp][X] = CV2[Nyp-1][X] = CV2[Nyp-2][X];
	CE2[Nyp][X] = CE2[Nyp-1][X] = CE2[Nyp-2][X];
//c		V2N(NY-1][I]=0
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
} //ENDDO

} //END SUBROUTINE S7_3000_MODELING



