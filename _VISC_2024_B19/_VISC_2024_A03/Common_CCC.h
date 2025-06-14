#pragma once
//!common / ccc / dl, an0, am10m, am20m, om, tpaz, r1, r2

int Y, X; //ИНДЕКСЫ ДЛЯ ПРОХОДА ПО МАТРИЦАМ!!!!! (БЫЛО K, I).
double an0, am10m, am20m, al1m, al2m, al1p, al2p, p1p, p2p, om;
double d0, d0Nusl;

double r1, r2, c1v, c2v, qk, cvMet;
double Gap; //common/a/GAP
double a_pl, a_mi, s, A, B, C, D;

double DVIS[4]; //На всякий случай, увеличил размер массива на 1.

double h2x, h1x, h2y, h1y;
double Vis1,  Vis0x,  Vis2x,  Vis0y,  Vis2y; 
double Vis1V, Vis0xV, Vis2xV, Vis0yV, Vis2yV;
double Vis1SV,Vis0xSV,Vis2xSV,Vis0ySV,Vis2ySV;

double Tep1, Tep0x, Tep2x, Tep0y, Tep2y;

double u1, u0x, u2x, u0y, u2y;
double v1, v0x, v2x, v0y, v2y;
double e1, e0x, e2x, e0y, e2y;
double u22, u02, u20, u00;
double v22, v02, v20, v00;

double dudx, dudy, dvdx, dvdy;
double dr11, dr12, dr21, dr22, dr31, dr32, dr41, dr42;
double dr1, dr2, dr3, dr4;

double k43, k23;