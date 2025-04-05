#pragma once
#include "Timer.h"
#include "Ret7_1_Tau.h"
#include "A7_20_RungeKutt.h"
#include "A7_39_WritRestart.h"



int IndexSW();

void A7_0_Evolution() {

int TmRestartOld     = int(Time / DtRestart);
int TmWriteStepOld   = int(Time / DtWriteStep);
int TmConsolOld      = int(Time / DtConsol);
int TmEmptyStringOld = int(Time / DtEmptyString);

//int k, i;
int y, x;
double Tau;
NT = 1;
Timer t;
//Timer t1;

//cout << "t(mkc) = " << TimePaz << "   " << NT << "   " << t.get() << endl;
cout << "t(mkc)=" << TimePaz << "    " << t.get() << "\n\n";
FTime << t.get() << endl;
FTime.precision(6);

Label6:

Tau = F7_1_Tau();

A7_20_RungeKutt(Tau);

Time = Time + Tau;
TimePaz = Time * tpaz;
NT = NT + 1;

int TmRestartNew = int(Time / DtRestart);
int TmWriteStepNew = int(Time / DtWriteStep);
int TmConsolNew = int(Time / DtConsol);
int TmEmptyStringNew = int(Time / DtEmptyString);


//!*********************************************************
if(TmRestartNew > TmRestartOld) {
	A7_39_WritRestart();
	TmRestartOld = TmRestartNew;
} //end if
//CALL S7_5_AMAX(metr)

//!*********************************************************

if(TmWriteStepNew > TmWriteStepOld) { //if(TmWriteStepNew > TmWriteStepOld) then
//	CALL S6_1_WRITE_STEP(Time,tpaz,om,c1v,c2v,gam1,gam2)
	A6_2_WriteStep();
//	FTime << Time * tpaz << "   " << Tau * tpaz << "   " << t.get() << endl;
	FTime << TimePaz << "\t  " << Tau * tpaz << "   " << t.get() << endl;
	TmWriteStepOld = TmWriteStepNew;
//	write(24, *)	Time* tpaz, Tau* TPAZ
} //end if

if(TmConsolNew > TmConsolOld) { //if(TmConsolNew > TmConsolOld) then
	TmConsolOld = TmConsolNew;
//	double time = t.stop();
//	std::cout << std::fixed << time << std::endl;
//	cout << "t(mkc) = " << Time * tpaz << "   " << NT << "   " << t.get() << endl;
//	cout << "t(mkc) = " << TimePaz << "   " << NT << "   " << t.get() << endl;
	cout << "t(mkc)=" << TimePaz << "  "
		<< t.get() << "  Nxp-iSW=" << Nxp-IndexSW() << endl;
} //end if


if(TmEmptyStringNew > TmEmptyStringOld) { //if(TmEmptyStringNew > TmEmptyStringOld) then
	TmEmptyStringOld = TmEmptyStringNew;
	cout << int(Time * tpaz) << endl;
} //end if


if(Time < TIME_MAX) { //if(Time < TIME_MAX) then
//	goto 6
	goto Label6;
} //end if

//7	format(8e14.4)
//1	FORMAT(A5,500E14.4)
//2	FORMAT (500E18.6)

//STOP
return;
}

//cout << "M7_0000_Evolution()" << endl;
//	for (Time = 1; Time < 3; Time++) M6_2_WriteStep();
//	cout << "int(8.8) =" << int(8.8) << endl;


int IndexSW() {
	int i = 1;
	double RoRight = RO1[1][Nxp];
	double RoSW = RO1[1][i];
	while(RO1[1][i] > RoRight * 1.01) {
		i++;
		if(i <= Nxp) RoSW = RO1[1][i];
	}
	return i;
}
