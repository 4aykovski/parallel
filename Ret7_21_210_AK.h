#pragma once

#include "Ret7_21_211_Diffu.h"

double F7_21_210_AK(double pl1, double pl2, double ck1, double ck2, double t1, double t2)
{
//c	ak_1=dl*an0*sig1**2*sqrt(1/om)
//c     **4/3.0*pl1*pl2/(1+om)*sqrt(3.14/2)*sqrt(t1*om+t2)*(1+sig2/sig1)**2
double x1, dav, d1;
	x1  = pl1 / (pl1+pl2/om);
	dav = r1 * pl1 * t1 + r2 * pl2 * t2;
//!	d1 = diffu(t1,t2,pl1,pl2,anu12,d0)
	d1 = F_21_211_Diffu(t1, t2, pl1, pl2);

//	S7_220_AK = x1*(1-x1)*dav / d1
//!	ak=0
	return x1 * (1 - x1) * dav / d1;
} //END FUNCTION S7_220_AK






