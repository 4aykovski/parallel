#pragma once
//COMMON / GGG / 

double xN, xK;					//ЛЕВАЯ И ПРАВАЯ ГРАНИЦЫ ОСНОВНОЙ ОБЛАСТИ.
double yN, yK;					//НИЖНЯЯ И ВЕРХНЯЯ ГРАНИЦЫ ОСНОВНОЙ ОБЛАСТИ.
double Xmin, Xmax, Ymin, Ymax;	//ГРАНИЦЫ ВСЕЙ ОБЛАСТИ.
//Но Xmax должна ещё зависеть от времени расчётов и числа Маха.

double h0;						//ШАГ ОСНОВНОЙ СЕТКИ.
double x3;						//ПОЛОЖЕНИЕ УДАРНОЙ ВОЛНЫ.
double dSW;                     //Ширина размазки ударной волны.

double Gxy1, Gxy2, Gxy3;	    //ШИРИНЫ ДОПОЛНИТЕЛЬНЫХ СЛОЁВ.
double qm1, qm2, qm3;           //!Коэффициенты для шагов в дополнительных областях.

double dL;						//Обезразмеривание геометрии.

int MN;				     //max(Nxp, Nyp).
int Nx, Ny;			     //КОЛИЧЕСТВО ШАГОВ ОСНОВНОЙ СЕТКИ.
int Nxp, Nyp, Nyp0;		 //КОЛИЧЕСТВО ШАГОВ ВСЕЙ СЕТКИ.
int Nxl, Nxr, Nyn, Nyv;

int IndexXN, IndexXK, IndexYN, IndexYK; //Возможно, нужно для диффузии!!

int BoundCondLeft;  //В фортране это была переменная Nl.
int BoundCondRight; //В фортране это была переменная Nr.
int BoundCondLower; //В фортране это была переменная Nd.


string HeadFiles;


