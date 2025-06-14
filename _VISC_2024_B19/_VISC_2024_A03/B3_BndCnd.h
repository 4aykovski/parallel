#pragma once

//void M7_4_0_Modeling(double Tau, double om, double c1v, double c2v) {
void B3_BndCnd()
//	double CRO1[MY][MX],double CU1[MY][MX],double CV1[MY][MX],double CE1[MY][MX],
//	double CRO2[MY][MX],double CU2[MY][MX],double CV2[MY][MX],double CE2[MY][MX]) 
{

//!	DO K=3,Nyp-2
for(Y=1; Y<=Nyp; Y++) { //DO K = 1,Nyp
	if(BoundCondLeft == 1) {
		GRO1[Y][1]= RO1[Y][1]; GRO1[Y][2]= RO1[Y][2];
		GU1[Y][1] = U1[Y][1];  GU1[Y][2] = U1[Y][2];
		GV1[Y][1] = V1[Y][1];  GV1[Y][2] = V1[Y][2];
		GE1[Y][1] = E1[Y][1];  GE1[Y][2] = E1[Y][2];
		
		GRO2[Y][1]= RO2[Y][1]; GRO2[Y][2]= RO2[Y][2];
		GU2[Y][1] = U2[Y][1];  GU2[Y][2] = U2[Y][2];
		GV2[Y][1] = V2[Y][1];  GV2[Y][2] = V2[Y][2];
		GE2[Y][1] = E2[Y][1];  GE2[Y][2] = E2[Y][2];
	} //end if
	else if(BoundCondLeft == 0) {
		GRO1[Y][1]= GRO1[Y][5]; GRO1[Y][2]= GRO1[Y][4];
		GU1[Y][1] = -GU1[Y][5]; GU1[Y][2] =-GU1[Y][4];
		GV1[Y][1] = GV1[Y][5];  GV1[Y][2] = GV1[Y][4];
		GE1[Y][1] = GE1[Y][5];  GE1[Y][2] = GE1[Y][4];
//c		U1N[K][2)=0
		GRO2[Y][1]= GRO2[Y][5]; GRO2[Y][2]= GRO2[Y][4];
		GU2[Y][1] = -GU2[Y][5]; GU2[Y][2] =-GU2[Y][4];
		GV2[Y][1] = GV2[Y][5];  GV2[Y][2] = GV2[Y][4];
		GE2[Y][1] = GE2[Y][5];  GE2[Y][2] = GE2[Y][4];
//c		U2N[K][2)=0
	}


	if(BoundCondRight == 1) { // ÂÍÈÌÀÍÈÅ - ÐÀÍÜØÅ ÁÛË Npr
		GRO1[Y][Nxp]= GRO1[Y][Nxp-1]= GRO1[Y][Nxp-2];
		GU1[Y][Nxp] = GU1[Y][Nxp-1] = GU1[Y][Nxp-2];
		GV1[Y][Nxp] = GV1[Y][Nxp-1] = GV1[Y][Nxp-2];
		GE1[Y][Nxp] = GE1[Y][Nxp-1] = GE1[Y][Nxp-2];
		
		GRO2[Y][Nxp]= GRO2[Y][Nxp-1]= GRO2[Y][Nxp-2];
		GU2[Y][Nxp] = GU2[Y][Nxp-1] = GU2[Y][Nxp-2];
		GV2[Y][Nxp] = GV2[Y][Nxp-1] = GV2[Y][Nxp-2];
		GE2[Y][Nxp] = GE2[Y][Nxp-1] = GE2[Y][Nxp-2];
	}
	else if(BoundCondRight == 2) {
		GRO1[Y][Nxp]= GRO1[Y][Nxp-4]; GRO1[Y][Nxp-1]= GRO1[Y][Nxp-3];  
		GU1[Y][Nxp] =-GU1[Y][Nxp-4];  GU1[Y][Nxp-1] =-GU1[Y][Nxp-3];   
		GV1[Y][Nxp] = GV1[Y][Nxp-4];  GV1[Y][Nxp-1] = GV1[Y][Nxp-3];   
		GE1[Y][Nxp] = GE1[Y][Nxp-4];  GE1[Y][Nxp-1] = GE1[Y][Nxp-3];   
		GRO2[Y][Nxp]= GRO2[Y][Nxp-4]; GRO2[Y][Nxp-1]= GRO2[Y][Nxp-3];  
		GU2[Y][Nxp] =-GU2[Y][Nxp-4];  GU2[Y][Nxp-1] =-GU2[Y][Nxp-3];   
		GV2[Y][Nxp] = GV2[Y][Nxp-4];  GV2[Y][Nxp-1] = GV2[Y][Nxp-3];   
		GE2[Y][Nxp] = GE2[Y][Nxp-4];  GE2[Y][Nxp-1] = GE2[Y][Nxp-3];   
	} //end if
} //ENDDO


//c******************************************************
//c**************************************************
for(X=1; X<=Nxp; X++) { //DO I = 1,Nxp
	GRO1[1][X]= GRO1[5][X]; GRO1[2][X]= GRO1[4][X];
	GU1[1][X] = GU1[5][X];  GU1[2][X] = GU1[4][X];
	GV1[1][X] =-GV1[5][X];  GV1[2][X] =-GV1[4][X]; GV1[3][X] = 0;
	GE1[1][X] = GE1[5][X];  GE1[2][X] = GE1[4][X];

	GRO2[1][X]= GRO2[5][X]; GRO2[2][X]= GRO2[4][X];
	GU2[1][X] = GU2[5][X];  GU2[2][X] = GU2[4][X];
	GV2[1][X] =-GV2[5][X];  GV2[2][X] =-GV2[4][X]; GV2[3][X] = 0;
	GE2[1][X] = GE2[5][X];  GE2[2][X] = GE2[4][X];
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	GRO1[Nyp][X]= GRO1[Nyp-1][X]= GRO1[Nyp-2][X];
	GU1[Nyp][X] = GU1[Nyp-1][X] = GU1[Nyp-2][X];
	GV1[Nyp][X] = GV1[Nyp-1][X] = GV1[Nyp-2][X];
	GE1[Nyp][X] = GE1[Nyp-1][X] = GE1[Nyp-2][X];
//c		V1N(NY-1][I]=0
	GRO2[Nyp][X]= GRO2[Nyp-1][X]= GRO2[Nyp-2][X];
	GU2[Nyp][X] = GU2[Nyp-1][X] = GU2[Nyp-2][X];
	GV2[Nyp][X] = GV2[Nyp-1][X] = GV2[Nyp-2][X];
	GE2[Nyp][X] = GE2[Nyp-1][X] = GE2[Nyp-2][X];
//c		V2N(NY-1][I]=0
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
} //ENDDO

} //END SUBROUTINE S7_3000_MODELING



