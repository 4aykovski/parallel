#pragma once
void Razmaz_SW();

void A5_InitalCond() {
//cout << "M5_0_InitalCond()" << endl;

//int i, k;
int x, y;
double gam, pp, ro;
double d2d1, dav, tem;
double pm, rop, rom, up, um, vv0, ucd;

//cout << "cmax = " << cmax << endl;
//cout << "p1p = " << p1p << endl;
//cout << "p2p = " << p2p << endl;
//cout << "c1v = " << c1v << endl;
//cout << "c2v = " << c1v << endl;
//cout << "pp = " << p1p + p2p << endl;
//	subroutine raspad(pp,pm,up,um,rop,rom,gam1,gam2,ucd,al1m,c1v,al2m,c2v)

for(y=3; y<=Nyp; y++) { //do k = 3,nyp
	for(x=1; x<=Nxp; x++) { //do i = 1,nxp
		U1[y][x] = 0.0;
		U2[y][x] = 0.0;
		V1[y][x] = 0.0;
		V2[y][x] = 0.0;
		E1[y][x] = c1v;
		E2[y][x] = c2v;
			
		if((XVV[x] < x3-0.01*h0) && (abs(cmax-1) > 0.0001)) {
gam = (al1p*gam1*c1v+al2p*gam2*c2v)/(al1p*c1v+al2p*c2v);
			pp = p1p+p2p;
			ro = pp*(gam+1)/(gam-1+2/pow(cmax,2));
			RO2[y][x] = al2p*ro;
			RO1[y][x] = al1p*ro;
	
d2d1 = (2*gam)/(gam+1)*pow(cmax,2)-(gam-1)/(gam+1);
			dav = d2d1*(al1p*r1+al2p*r2)*pp;
			tem = dav/((al1p*r1+al2p*r2)*ro);

			E1[y][x] = c1v*tem;
			E2[y][x] = c2v*tem;
U2[y][x] = sqrt(2*(al1p*r1+al2p*r2)/((gam+1)*d2d1+gam-1))*(d2d1-1);
			U1[y][x] = U2[y][x];
//!				vv0=U1[k][i]
		} //end if
//!			xx=xx+hx[i]
	} //enddo
} //enddo



//Размазывание Фронта Ударной Волны.
if(dSW > 0) Razmaz_SW();
//Размазывание Фронта Ударной Волны.



for(x=1; x<=Nxp; x++) { //do i = 1,nxp
	RO1[1][x]= RO1[5][x]; RO1[2][x]= RO1[4][x];
	U1[1][x] = U1[5][x];  U1[2][x] = U1[4][x];
	V1[1][x] =-V1[5][x];  V1[2][x] =-V1[4][x];
	E1[1][x] = E1[5][x];  E1[2][x] = E1[4][x];

	RO2[1][x]= RO2[5][x]; RO2[2][x]= RO2[4][x];
	U2[1][x] = U2[5][x];  U2[2][x] = U2[4][x];
	V2[1][x] =-V2[5][x];  V2[2][x] =-V2[4][x];
	E2[1][x] = E2[5][x];  E2[2][x] = E2[4][x];
} //enddo
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	!

pp  = 1;
pm  = (gam1-1)*RO1[1][1]*E1[1][1];
rop = om;
rom = RO1[1][1];
up  = 0;
um  = U1[1][1];

//vv0 = 0.57*tpaz/(dL*1000);
//vv0 = 0.79*tpaz/(dL*1000);
vv0 = 0;
//	call raspad(pp,pm,up,um,rop,rom,gam1,gam2,ucd,al1m,c1v,al2m,c2v)
//vv0 = ucd;
vv0 = (U1[1][1]+vv0)/2;
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
cout << "U_shift = " << (vv0-U1[1][1])/(tpaz*1E-4) << "; ";
vv0 = U1[1][1] + U_shift;
cout << U_shift / (tpaz * 1E-4) << endl;
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


for(y=1; y<=Nyp; y++) { //do k = 1,nyp
	for(x=1; x<=Nxp; x++) { //do i = 1,nxp
		U1[y][x] = U1[y][x]-vv0;
		U2[y][x] = U2[y][x]-vv0;
//c		U1[k][i]=U1[k][i]+cmax*sqrt(gam1*r1)/6
//c		U2[k][i]=U1[k][i]
	} //enddo
} //enddo


} //END SUBROUTINE S5_0_INITAL_COND


void Razmaz_SW() {
//int i,k;
int x,y;
double xPr; //, dSW; Координата профиля (продольного).
double aR1, bR1, dR1;
double aR2, bR2, dR2;
double aE1, bE1, dE1;
double aE2, bE2, dE2;
double aU, bU, dU;

//dSW = 1.0 / 100; //dSW = 0.0;
	aR1 = RO1[3][1]; bR1 = RO1[3][Nxp];
	aR2 = RO2[3][1]; bR2 = RO2[3][Nxp];
	aE1 = E1[3][1];  bE1 = E1[3][Nxp];
	aE2 = E2[3][1];  bE2 = E2[3][Nxp];
	aU  = U1[3][1];  bU  = U1[3][Nxp];
	//cout << "aR = " << aR << endl;
	//cout << "bR = " << bR << endl;
	dR1 = (aR1 - bR1) / dSW;
	dR2 = (aR2 - bR2) / dSW;
	dE1 = (aE1 - bE1) / dSW;
	dE2 = (aE2 - bE2) / dSW;
	dU  = (aU  - bU)  / dSW;
	for(x=1; x<=Nxp; x++) {
		xPr = XVV[x];
//		if ((x>x3) && (x<x3+dSW)) {
		if((xPr<x3) && (xPr>x3-dSW)) {
			for(y=3; y<=Nyp; y++) {
				RO1[y][x] = bR1 + (x3 - xPr) * dR1;
				RO2[y][x] = bR2 + (x3 - xPr) * dR2;
				E1[y][x]  = bE1 + (x3 - xPr) * dE1;
				E2[y][x]  = bE2 + (x3 - xPr) * dE2;
				U1[y][x] =  bU + (x3 - xPr) * dU;
				U2[y][x]  = U1[y][x];
			}
		}
	}
}



/*
//subroutine raspad(pp, pm, up, um, rop, rom, gam1, gam2, ucd, al1m, c1v, al2m, c2v)
void raspad() {

double eps, pcd;
int i;
	eps = 1.e-6;
	i = 0;
	pcd = (pp + pm) / 2;
	ucd = (up + um) / 2;
	gam = (al1m * gam1 * c1v + al2m * gam2 * c2v) / (al1m * c1v + al2m * c2v);
//222	continue
Label222: //continue
//	if (Pcd.GT.Pm) then
	if (Pcd > Pm) {
		//A = SQRT(((Gam1 + 1.) / 2. * Pcd + (Gam1 - 1) / 2. * Pm) * ROm)
		A = SQRT(((Gam1 + 1.) / 2. * Pcd + (Gam1 - 1) / 2. * Pm) * ROm);
	}
	else {
		//A = (Gam1 - 1) / 2. / Gam1 * SQRT(Gam1 * Pm * ROm)
		A = (Gam1 - 1) / 2. / Gam1 * SQRT(Gam1 * Pm * ROm);
		//A = A*(1.-Pcd/Pm) / (1.-(Pcd/Pm)**((Gam1-1)/2./Gam1))
		A = A * (1. - Pcd / Pm) / (1. - (Pcd / Pm) * *((Gam1 - 1) / 2. / Gam1));
	//end if
	}

	//if (Pcd .GT. Pp) then
	if (Pcd > Pp) {
//		B = SQRT(((Gam + 1.) / 2. * Pcd + (Gam - 1) / 2. * Pp) * ROp)
		B = SQRT(((Gam + 1.) / 2. * Pcd + (Gam - 1) / 2. * Pp) * ROp);
	}
	else {
//		B = (Gam-1)/2./Gam * SQRT(Gam*Pp*ROp)
		B = (Gam - 1) / 2. / Gam * SQRT(Gam * Pp * ROp);
//		B = B*(1.-Pcd/Pp) / (1.-(Pcd/Pp)**((Gam-1)/2./Gam))
		B = B * (1. - Pcd / Pp) / (1. - (Pcd / Pp) * *((Gam - 1) / 2. / Gam));
	//end if
	}

//	Pcd1 = (B*Pm+A*Pp+A*B*(Um-Up)) / (A+B)
	Pcd1 = (B * Pm + A * Pp + A * B * (Um - Up)) / (A + B);
//	Ucd  = (A*Um + B*Up + Pm-Pp) / (A+B)
	Ucd  = (A*Um + B*Up + Pm-Pp) / (A+B);

//	if (ABS(Pcd-Pcd1) < EPS) goto 212
	if (ABS(Pcd - Pcd1) < EPS) goto Label212

//	Pcd = Pcd1
	Pcd = Pcd1;
//	I = I+1
	I = I + 1;
//	goto 222
	goto Label222;

//212	continue
Label212: //continue
//return
//end subroutine raspad
}
	*/
