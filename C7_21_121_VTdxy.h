#pragma once

#include "Ret7_21_211_Diffu.h"

double vis1(double t1, double t2, double ro1, double ro2) {
	//!w=1.0/(dl*an0*sig1**2*1.e5)
	//!vis1=w*5.0/32.0*sqrt(t1/3.14)/(1+ro2/ro1/om*sqrt(2.0)/4.0
	//!**sqrt(1+t2/t1/om)*(1+sig2/sig1)**2)
double vis, x2, d1;		
	vis = vis10 * pow(t1,anu1);
	x2  = ro2 / (ro1*om+ro2);
    //d1 = S7_221_DIFFU(t1,t2,r01,ro2);   //ÁÛËÎ r01, גלוסעמ r01!!!!!!!!!
	d1 = F_21_211_Diffu(t1, t2, ro1, ro2);  //ÁÛËÎ r01, גלוסעמ r01!!!!!!!!!
return vis*ro1 / (ro1+3.0/5*vis/d1*x2);
} //end function vis1

double vis2(double t1, double t2, double ro1, double ro2) {
	//!w=1.0/(dl*an0*sig1**2*1.e5)
	//!vis2=w*5.0/32.0*sqrt(t2/3.14*om)*(sig1/sig2)**2/
	//!*(1+5.0/12.0*ro1/ro2*sqrt(1+t1/t2*om)*(1+sig1/sig2)**2)
double vis, x2, d1;	
	vis = vis20 * pow(t2,anu2);
	x2 = ro2 / (ro1 * om + ro2);
	d1 = F_21_211_Diffu(t1, t2, ro1, ro2);
return vis*ro2 / (ro2+6*vis/d1*(1-x2));
} //end function vis2

double tep1(double t1, double t2, double ro1, double ro2) {
	//!w=1.0/(dl*an0*sig1**2*1.e5)
	//!tepl1=w*75.0/128.0*sqrt(t1/3.14)/(1+5*ro2/ro1/om*sqrt(2.0)/8.0
	//!**sqrt(1+t2/t1/om)*(1+sig2/sig1)**2)
double vis, tep, x2, d1;
	vis = vis10 * pow(t1,anu1);
	tep = vis * 5.0 / 2.0 / (4.0 / 3.0);
	x2 = ro2 / (ro1 * om + ro2);
//d1 = S7_221_DIFFU(t1, t2, r01, ro2);     //ÁÛËÎ r01, גלוסעמ r01!!!!!!!!!
	d1 = F_21_211_Diffu(t1, t2, ro1, ro2);   //ÁÛËÎ r01, גלוסעמ r01!!!!!!!!!
return tep*ro1 / (ro1+3.0/2.0*vis/d1*x2);
} //end function tepl1

double tep2(double t1, double t2, double ro1, double ro2) {
	//!w=1.0/(dl*an0*sig1**2*1.e5)
	//!tepl2=w*(sig1/sig2)**2*sqrt(t2/3.14/om)/(1+15.0*ro1/ro2
	//!**sqrt(2.0)/16.0*sqrt(1+t1/t2*om)*(1+sig1/sig2)**2)
double vis, tep, x2, d1; 
	vis = vis20 * pow(t2,anu2);
	tep = vis * 5.0 / 2.0 / (4.0 / 3.0) / om;
	x2 = ro2 / (ro1 * om + ro2);
	//!d1 = diffu(t1,t2,ro1,ro2,anu12,d0)
	d1 = F_21_211_Diffu(t1, t2, ro1, ro2);
return tep*ro2 / (ro2+9.0/2.0*vis/d1*(1-x2));
	
} //end function tepl2


double dH(double uv2, double uv0,	double h2, double h1) {
	return (uv2 - uv0) / (h1 + h2);
} //end function dX


double dHdH(double vis2xy, double vis1xy, 
	        double uve2, double uve1, double uve0, 
	        double h2, double h1) {
return (vis2xy*(uve2-uve1)/h2 - vis1xy*(uve1-uve0)/h1) / h2;
} //end function dxdx


/*double dYdY(double fi2, double fi1,
	        double psi2, double psi1, double psi0, 
	        double h2, double h1) {
return (fi2*(psi2-psi1)/h2 - fi1*(psi1-psi0)/h1) / h2;
} //end function dydy
*/

double dXdY(double vis2x, double vis0x, 
	        double uv22, double uv02, double uv20, double uv00, 
	        double h2x, double h1x, double h2y, double h1y) {
return (vis2x*(uv22-uv02)/(h2y+h1y) - vis0x*(uv20-uv00)/(h2y+h1y)) / (h2x + h1x);
	
} //end function dxdy


double dYdX(double vis2y, double vis0y, 
	        double uv22, double uv20, double uv02, double uv00, 
	        double h2x, double h1x, double h2y, double h1y) {
return (vis2y*(uv22-uv20)/(h2x+h1x) - vis0y*(uv02-uv00)/(h2x+h1x)) / (h2y + h1y);
} //end function dydx
	


/*
double vis1t(double t1, double t2, double ro1, double ro2) {
	//!w = 1.0 / (dl * an0 * sig1 * *2 * 1.e5)
	//!vis1 = w * 5.0 / 32.0 * sqrt(t1 / 3.14) / (1 + ro2 / ro1 / om * sqrt(2.0) / 4.0
	//!**sqrt(1 + t2 / t1 / om) * (1 + sig2 / sig1) * *2)
	vis = vis10 * t1 * *anu1;
	x2 = ro2 / (ro1 * om + ro2);
	//!d1 = diffu(t1, t2, r01, ro2, anu12, d0)
	d1 = S7_221_DIFFU(t1, t2, r01, ro2);

	vis1t = vis / (ro1 + 3.0 / 5 * vis / d1 * x2);
	return
} //end function vis1t
	

function vis2t(t1, t2, ro1, ro2)
//!w = 1.0 / (dl * an0 * sig1 * *2 * 1.e5)
//!vis2 = w * 5.0 / 32.0 * sqrt(t2 / 3.14 * om) * (sig1 / sig2) * *2 /
//!*(1 + 5.0 / 12.0 * ro1 / ro2 * sqrt(1 + t1 / t2 * om) * (1 + sig1 / sig2) * *2)
	vis = vis20 * t2 * *anu2;
	x2 = ro2 / (ro1 * om + ro2);
//!d1 = diffu(t1, t2, ro1, ro2, anu12, d0)
	d1 = S7_221_DIFFU(t1, t2, ro1, ro2);

	vis2t = vis / (ro2 + 6 * vis / d1 * (1 - x2));
	return
	end function vis2t

	END MODULE M7_2120_VIS_TEPL_DXY

*/
