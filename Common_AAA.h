#pragma once
//C     parameter(L = 801, M = 5001, MN = 5001)
//const int n1001 = 1001;

#define MR 2
#define MX 1001
#define MY 1001

double RO1[MY][MX], LRO1[MY][MX], KRO1[MR][MY][MX]; //BRO1[MY][MX];
double RO2[MY][MX], LRO2[MY][MX], KRO2[MR][MY][MX]; //BRO2[MY][MX];
double U1[MY][MX],  LU1[MY][MX],  KU1[MR][MY][MX];  //BU1[MY][MX];
double U2[MY][MX],  LU2[MY][MX],  KU2[MR][MY][MX];  //BU2[MY][MX];
double V1[MY][MX],  LV1[MY][MX],  KV1[MR][MY][MX];  //BV1[MY][MX];
double V2[MY][MX],  LV2[MY][MX],  KV2[MR][MY][MX];  //BV2[MY][MX];
double E1[MY][MX],  LE1[MY][MX],  KE1[MR][MY][MX];  //BE1[MY][MX];
double E2[MY][MX],  LE2[MY][MX],  KE2[MR][MY][MX];  //BE2[MY][MX];

double P[MX];

double ARP[MX], AUP[MX], AVP[MX], AEP[MX];
double ARM[MX], AUM[MX], AVM[MX], AEM[MX];

double HX[MX], HY[MY], XVV[MX], YVV[MY];



/*COMMON / AA / U1(1001, 3001), RO1(1001, 3001), E1(1001, 3001),
* V1(1001, 3001), U2(1001, 3001), RO2(1001, 3001), E2(1001, 3001),
* V2(1001, 3001), U1N(1001, 3001), RO1N(1001, 3001), E1N(1001, 3001),
* V1N(1001, 3001), U2N(1001, 3001), RO2N(1001, 3001),
* E2N(1001, 3001), V2N(1001, 3001), U1N1(1001, 3001),
* RO1N1(1001, 3001), E1N1(1001, 3001), V1N1(1001, 3001),
* U2N1(1001, 3001), RO2N1(1001, 3001),
* E2N1(1001, 3001), V2N1(1001, 3001),
* P(3001)
*/