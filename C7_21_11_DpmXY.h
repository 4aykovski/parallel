#pragma once
//#include "M7_2121_PLUS_DINS.h"
 
double d_p_X(double f[]) {
//double a_pl, a_mi, s, A, B;
    a_pl = f[X+1] - f[X];   //a_pl_i = plus_i(f,I);
    a_mi = f[X]   - f[X-1]; //a_mi_i = dins_i(f,I);

    s = (2*a_pl*a_mi+1e-06) / (pow(a_pl,2)+pow(a_mi,2)+1e-06);
    A = f[X] + s*(1.-Gap*s)/4*a_mi + s*(1.+Gap*s)/4*a_pl;
    
    a_pl = f[X]   - f[X-1]; //a_pl_i = plus_i(f,ii);
    a_mi = f[X-1] - f[X-2]; //a_mi_i = dins_i(f,ii);

    s = (2*a_pl*a_mi+1e-06) / (pow(a_pl,2)+pow(a_mi,2)+1e-06);
    B = f[X-1] + s*(1.-Gap*s)/4*a_mi + s*(1.+Gap*s)/4*a_pl;

    return A - B;
} //end function d_p_i

 

double d_m_X(double f[]) {
//double a_pl, a_mi, s, C, D;
    a_pl = f[X+2] - f[X+1]; //a_pl_i = plus_i(f,ii);
    a_mi = f[X+1] - f[X];   //a_mi_i = dins_i(f,ii);

    s = (2*a_pl*a_mi+1e-06) / (pow(a_pl,2)+pow(a_mi,2)+1e-06);
    C = f[X+1] - s*(1.-Gap*s)/4.0*a_pl - s*(1.+Gap*s)/4.0*a_mi;
    
    a_pl = f[X+1] - f[X];   //a_pl_i = plus_i(f,I);
    a_mi = f[X]   - f[X-1]; //a_mi_i = dins_i(f,I);

    s = (2*a_pl*a_mi+1e-06) / (pow(a_pl,2)+pow(a_mi,2)+1e-06);
    D = f[X] - s*(1.-Gap*s)/4.0*a_pl - s*(1.+Gap*s)/4.0*a_mi;

    return C - D;
} //end function d_m_i



double d_p_Y(double f[]) {
//double a_pl, a_mi, s, A, B;

    a_pl = f[Y+1] - f[Y];   //a_pl_k = plus_k(f,K);
    a_mi = f[Y]   - f[Y-1]; //a_mi_k = dins_k(f,K);

    s = (2*a_pl*a_mi+1e-06) / (pow(a_pl,2)+pow(a_mi,2)+1e-06);
    A = f[Y] + s*(1.-Gap*s)/4*a_mi + s*(1.+Gap*s)/4*a_pl;

    a_pl = f[Y]   - f[Y-1]; //a_pl_k = plus_k(f,kk);
    a_mi = f[Y-1] - f[Y-2]; //a_mi_k = dins_k(f,kk);

    s = (2*a_pl*a_mi+1e-06) / (pow(a_pl,2)+pow(a_mi,2)+1e-06);
    B = f[Y-1] + s*(1.-Gap*s)/4*a_mi + s*(1.+Gap*s)/4*a_pl;

    return A - B;
} //end function d_p_k



double d_m_Y(double f[]) {
//double a_pl, a_mi, s, C, D;

    a_pl = f[Y+2]-f[Y+1]; //a_pl_k = plus_k(f,kk);
    a_mi = f[Y+1]  -f[Y]; //a_mi_k = dins_k(f,kk);

    s = (2*a_pl*a_mi+1e-06) / (pow(a_pl,2)+pow(a_mi,2)+1e-06);
    C = f[Y+1] - s*(1.-Gap*s)/4.0*a_pl - s*(1.+Gap*s)/4.0*a_mi;

    a_pl = f[Y+1]-f[Y];   //a_pl_k = plus_k(f,K);
    a_mi = f[Y]  -f[Y-1]; //a_mi_k = dins_k(f,K);

    s = (2*a_pl*a_mi+1e-06) / (pow(a_pl,2)+pow(a_mi,2)+1e-06);
    D = f[Y] - s*(1.-Gap*s)/4.0*a_pl - s*(1.+Gap*s)/4.0*a_mi;

    return C - D;
} //end function d_m_k


//!**************************************************************
//    END MODULE M7_2110_D_PKMP_PIMI
