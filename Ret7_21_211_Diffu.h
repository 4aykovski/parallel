#pragma once

double F_21_211_Diffu(double t1, double t2, double pl1, double pl2) {

//!	d0=3.0/(4*(dl*an0*sig1**2*1.e5))
//!     **sqrt((1+om)/(om*3.14/2))/(1+sig2/sig1)**2
//!   d0=0.097*1.e-8*tpaz
//!	x1=pl1/(pl1+pl2/om)

double tem;
	tem = t1*om+t2;

//	S7_221_DIFFU = d0 * (tem/(om+1))**anu12
	return 	d0 * pow(tem/(om+1),anu12);

} //END FUNCTION S7_221_DIFFU



