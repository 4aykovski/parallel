#pragma once

void A7_2_2_RK_LayerUp() {

int Nxt,x,y;

if (TypeTau == '2') Nxt = 0;
else                Nxt = 1;
	
for (y = 1; y <= Nyp; y++) { //DO K = 1,NYp
	for(x=1; x<=Nxp; x++) { //DO I = 1,NXp
		RO1[y][x] = KRO1[Nxt][y][x];
		U1[y][x] = KU1[Nxt][y][x];
		V1[y][x] = KV1[Nxt][y][x];
		E1[y][x] = KE1[Nxt][y][x];
		RO2[y][x] = KRO2[Nxt][y][x];
		U2[y][x] = KU2[Nxt][y][x];
		V2[y][x] = KV2[Nxt][y][x];
		E2[y][x] = KE2[Nxt][y][x];
	} //ENDDO
} //ENDDO



}

