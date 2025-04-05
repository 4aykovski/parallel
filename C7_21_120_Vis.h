#pragma once

#include "C7_21_121_VTdxy.h"

void C7_21_120_Vis(int Met, double Umet[MY][MX], double Vmet[MY][MX], double Emet[MY][MX],
	double TrRO1[MY][MX], double TrRO2[MY][MX], double TrE1[MY][MX], double TrE2[MY][MX]) {

//double cvMet;
//if (K == 8 && I == 8) cout << PRO1 - K1RO1 << endl;

if(Met==1) {
	Vis1  = vis1(TrE1[Y][X]  /c1v, TrE2[Y][X]  /c2v, TrRO1[Y][X],   TrRO2[Y][X]);
	Tep1  = tep1(TrE1[Y][X]  /c1v, TrE2[Y][X]  /c2v, TrRO1[Y][X],   TrRO2[Y][X]);
	Vis0x = vis1(TrE1[Y][X-1]/c1v, TrE2[Y][X-1]/c2v, TrRO1[Y][X-1], TrRO2[Y][X-1]);
	Vis2x = vis1(TrE1[Y][X+1]/c1v, TrE2[Y][X+1]/c2v, TrRO1[Y][X+1], TrRO2[Y][X+1]);
	Tep2x = tep1(TrE1[Y][X+1]/c1v, TrE2[Y][X+1]/c2v, TrRO1[Y][X+1], TrRO2[Y][X+1]);
	Vis0y = vis1(TrE1[Y-1][X]/c1v, TrE2[Y-1][X]/c2v, TrRO1[Y-1][X], TrRO2[Y-1][X]);
	Vis2y = vis1(TrE1[Y+1][X]/c1v, TrE2[Y+1][X]/c2v, TrRO1[Y+1][X], TrRO2[Y+1][X]);
	Tep2y = tep1(TrE1[Y+1][X]/c1v, TrE2[Y+1][X]/c2v, TrRO1[Y+1][X], TrRO2[Y+1][X]);
	visVol = visVol1;
}
else if(Met==2) {
	Vis1  = vis2(TrE1[Y][X]  /c1v, TrE2[Y][X]  /c2v, TrRO1[Y][X],   TrRO2[Y][X]);
	Tep1 =  tep2(TrE1[Y][X]  /c1v, TrE2[Y][X]  /c2v, TrRO1[Y][X],   TrRO2[Y][X]);
	Vis0x = vis2(TrE1[Y][X-1]/c1v, TrE2[Y][X-1]/c2v, TrRO1[Y][X-1], TrRO2[Y][X-1]);
	Vis2x = vis2(TrE1[Y][X+1]/c1v, TrE2[Y][X+1]/c2v, TrRO1[Y][X+1], TrRO2[Y][X+1]);
	Tep2x = tep2(TrE1[Y][X+1]/c1v, TrE2[Y][X+1]/c2v, TrRO1[Y][X+1], TrRO2[Y][X+1]);
	Vis0y = vis2(TrE1[Y-1][X]/c1v, TrE2[Y-1][X]/c2v, TrRO1[Y-1][X], TrRO2[Y-1][X]);
	Vis2y = vis2(TrE1[Y+1][X]/c1v, TrE2[Y+1][X]/c2v, TrRO1[Y+1][X], TrRO2[Y+1][X]);
	Tep2y = tep2(TrE1[Y+1][X]/c1v, TrE2[Y+1][X]/c2v, TrRO1[Y+1][X], TrRO2[Y+1][X]);
	visVol = visVol2;
}
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	
	h1x = HX[X];   h1y = HY[Y];
	h2x = HX[X+1]; h2y = HY[Y+1];
	
	u1  = Umet[Y][X];     v1  = Vmet[Y][X];
	u0x = Umet[Y][X-1];   v0x = Vmet[Y][X-1];
	u2x = Umet[Y][X+1];   v2x = Vmet[Y][X+1];
	u0y = Umet[Y-1][X];   v0y = Vmet[Y-1][X];
	u2y = Umet[Y+1][X];   v2y = Vmet[Y+1][X];

	u00 = Umet[Y-1][X-1]; v00 = Vmet[Y-1][X-1];
	u02 = Umet[Y-1][X+1]; v02 = Vmet[Y-1][X+1];
	u20 = Umet[Y+1][X-1]; v20 = Vmet[Y+1][X-1];
	u22 = Umet[Y+1][X+1]; v22 = Vmet[Y+1][X+1];

	Vis1V  = Vis1  * visVol;
	Vis0xV = Vis0x * visVol;
	Vis2xV = Vis2x * visVol;
	Vis0yV = Vis0y * visVol;
	Vis2yV = Vis2y * visVol;

	Vis1SV = Vis1V   + k43*Vis1;
	Vis0xSV = Vis0xV - k23*Vis0x;
	Vis2xSV = Vis2xV + k43*Vis2x;
	Vis0ySV = Vis0yV - k23*Vis0y;
	Vis2ySV = Vis2yV + k43*Vis2y;

dr1 = dHdH(Vis2xSV,Vis1SV,  u2x,u1,u0x,      h2x,h1x);
dr3 = dHdH(Vis2y,  Vis1,    u2y,u1,u0y,      h2y,h1y);
	Vis2xSV = Vis2xV - k23*Vis2x;
dr2 = dXdY(Vis2xSV,Vis0xSV, v22,v02,v20,v00, h2x,h1x,h2y,h1y);
dr4 = dYdX(Vis2y,  Vis0y,   v22,v20,v02,v00, h2x,h1x,h2y,h1y);
DVIS[0] = -(dr1+dr2+dr3+dr4);
//UUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUU

dr1 = dHdH(Vis2ySV,Vis1SV,  v2y,v1,v0y,      h2y,h1y);
dr4 = dHdH(Vis2x,  Vis1,    v2x,v1,v0x,      h2x,h1x);
	Vis2ySV = Vis2yV - k23 * Vis2y;
dr2 = dYdX(Vis2ySV,Vis0ySV, u22,u20,u02,u00, h2x,h1x,h2y,h1y);
dr3 = dXdY(Vis2x,  Vis0x,   u22,u02,u20,u00, h2x,h1x,h2y,h1y);
DVIS[1] = -(dr1+dr2+dr3+dr4);
//VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV

Vis1SV = Vis1V - k23 * Vis1;

dudx = dH(Umet[Y][X+1],Umet[Y][X-1], HX[X+1],HX[X-1]);
dvdx = dH(Vmet[Y][X+1],Vmet[Y][X-1], HX[X+1],HX[X-1]);
dudy = dH(Umet[Y+1][X],Umet[Y-1][X], HY[Y+1],HY[Y-1]);
dvdy = dH(Vmet[Y+1][X],Vmet[Y-1][X], HY[Y+1],HY[Y-1]);

dr1 = 2.0 * (pow(dudx,2) + pow(dvdy,2));
dr2 = pow(dudy+dvdx, 2);
dr3 = pow(dudx+dvdy, 2);
DVIS[2] = -(Vis1*(dr1+dr2) + Vis1SV*dr3);

e1  = Emet[Y][X]   /cvMet;
e2x = Emet[Y][X+1] /cvMet;
e0x = Emet[Y][X-1] /cvMet;
e2y = Emet[Y+1][X] /cvMet;
e0y = Emet[Y-1][X] /cvMet;
	
dr1 = dHdH(Tep2x,Tep1, e2x,e1,e0x, h2x,h1x); 
dr2 = dHdH(Tep2y,Tep1, e2y,e1,e0y, h2y,h1y);

DVIS[2] -= (dr1 + dr2);
//EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE

}

