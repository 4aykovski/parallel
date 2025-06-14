#pragma once

void A2_Param() {
//cout << "M4_Param()" << endl;
//!write(*, *) 'S4_Param.for'
double p0, am20p, p1m, p2m;
	
	om = om2 / om1;	
	dL = 0.1;
	qk = 3.0 / (1 + om);
	Gap = 1.0 / 3.0;

	tpaz = dL/(sqrt(1.380662*3*60.22045/om1)*100)*1e6;
	
	U_shift = U_shift*tpaz*1E-4;

	d0Nusl = d0;

	d0 = d0*1.e-10*tpaz/pow(dL,2);
	p0  = 0.8;
	an0 = p0*1.01325*10 / (3*1.380662);

	vis10 *= visS;
	vis20 *= visS;
	vis10 = vis10*tpaz / (om1/6.022045*an0*pow(dL,2))*1.0E-11;
	vis20 = vis20*tpaz / (om2/6.022045*an0*pow(dL,2))*1.0E-11;

//!am10m=alfa10
//c	am10m=0
	alfa10 = alfa101;
	am10m = alfa10*om/(1+(om-1)*alfa10);
	am20m = 1-am10m;

	al1m  = am10m/(am10m+am20m*om);
	al2m  = 1-al1m;
//c	am10p=1
	am20p = 1-am10p;
	al1p  = am10p/(am10p+am20p*om);
	al2p  = 1-al1p;
	p1m   = am10m;
	p2m   = (1-p1m)*om;
	p1p   = am10p;
	p2p   = om*am20p;
	r1 = 1;
	r2 = r1/om;

//c	gam1=5./3.
//c	gam2=5./3.
	c1v = r1/(gam1-1);
	c2v = r2/(gam2-1);

//!AKUR=0.95


//!x3=-dd/2-4.0*sigma
//!x3=-2*dd
//!x3 =  xk/2
//!x3 =-dd/2-3.0*sigma

//!x3=am

//!TIME_MAX  = 50

	Time = 0.0;

	TIME_MAX = TIME_MAX / tpaz;

//!Обезразмеривание Частоты записей для Restart.
	DtRestart = DtRestart / tpaz;
//!Новое время для записи Restart.
	TmRestart = DtRestart;

//!Обезразмеривание Частоты вывода на Консоль(Экран).
	DtConsol  = DtConsol / (tpaz/1);
//!Новое время для  вывода на Консоль(Экран).
	TmConsol = DtConsol;

//!Обезразмеривание Частоты выводов Результатов.
	DtWriteStep  = DtWriteStep / tpaz;
//!Новое время для вывода Результатов.
	TmWriteStep = DtWriteStep;

	DtEmptyString = DtEmptyString / tpaz;
	TmEmptyString = DtEmptyString;
//!*************************************************************

	k43 = 4.0 / 3.0;
	k23 = 2.0 / 3.0;
//	cout << "tpaz=" << tpaz << endl;

} //END SUBROUTINE S4_Param
