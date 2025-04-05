#pragma once

#include "FilesExtern.h"

#include "Common_AAA.h"
#include "Common_CCC.h"
#include "Common_DDD.h"
#include "Common_GGG.h"
#include "Common_TTT.h"
#include "Common_VVV.h"

//#include "Timer.h"
#include "FilesCreate.h"

#include "A1_Vvod.h"
#include "A2_Param.h"
#include "A3_0_Grid.h"
#include "A4_Diff_T0.h"
#include "A5_InitalCond.h"
#include "A6_0_WriteData0.h"
#include "A7_0_Evolution.h"

int main() {

	cout << "vvod from files: ";  //<< endl;
	cin >> nf;
	
	//cout << "RO1=" << RO1<<endl;
	A1_Vvod();
	A2_Param();
	
	A3_0_Grid();
	A4_Diff_T0();
	//M4_Param();
	A5_InitalCond();
	A6_0_WriteData0();

//	Timer t;
	A7_0_Evolution();

	cout << endl;
	return 0;
}




//////////!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

/*	cout << endl;
	cout << sd << " -- sd" << endl;	   //!Расст.м - ду.центрами(безразмерная).
	cout << dd << " -- dd" << endl;	   //!Диаметр капель(мм.)..
	cout << sigma << " -- sigma" << endl;	   //!Размазка.
	cout << KONF << " -- KONF" << endl;     //!Тип конфиг. - 1(1С), 2(2С), 3(3С), 4(1С - UP), 5(2С - UP).
	cout << X3C << " -- X3C" << endl;      //!Отступ удалённых от ударной волны капель(1).
	cout << am10p << " -- am10p" << endl;    //!Конценр.лёг.газа перед и за волной(0 1).
	cout << alfa101 << ' ' << alfa102 << ' ' << alfa103 << " -- alfa101 alfa102 alfa103" << endl;
	//!конценр.лег.газа перед и за волной(0.0 0.0 0.0).

	cout << xN << ' ' << xK << " -- xN xK" << endl;	   //!Границы.осн.обл.:лев., пр., верх.
	cout << yN << ' ' << yK << " -- yN yK" << endl;	   //!Границы.осн.обл. : лев., пр., верх.
	cout << x3 << " -- x3" << endl;	   //!Границы.осн.обл. : лев., пр., верх.
	cout << h0 << " -- h0" << endl;	   //!Шаг сетки, мм.

	cout << Nl1 << ' ' << Nl2 << ' ' << Nl3 << " -- Nl1 Nl2 Nl3" << endl;  //!Кол - во.доп.слоёв Слева.
	cout << Nr1 << ' ' << Nr2 << ' ' << Nr3 << " -- Nr1 Nr2 Nr3" << endl;  //!Кол - во.доп.слоёв Справа.
	cout << Nu1 << ' ' << Nu2 << ' ' << Nu3 << " -- Nu1 Nu2 Nu3" << endl;  //!Кол - во.доп.слоёв Сверху.
	cout << qm1 << ' ' << qm2 << ' ' << qm3 << " -- qm1 qm2 qm3" << endl;  //!Коэффициент для шага.
	cout << Nl << ' ' << Nr << " -- Nl Nr" << endl;	   //!Усл.на гр.по х(лев., пр)(1 - своб., 0 - тв.стенка).
	cout << Nd << " -- Nd" << endl;		   //!УСЛОВИЕ НА НИЖНЕЙ ГРАНИЦЕ(1 - СВ. 2 - СИММЕТР.УСЛ.).

	//!*****************************************************************
	//!Хранятся в COMMON / TTT /
	cout << DtRestart << " -- DtRestart" << endl;      //!Частота записей для Restart(мкс).
	cout << DtConsol << " -- DtConsol" << endl;       //!Частота сообщений на Консоли(Экране)(dt).
	cout << DtWriteStep << " -- DtWriteStep" << endl;    //!Частота вывода результатов(мкс)(dxp).
	cout << DtEmptyString << " -- DtEmptyString" << endl;  //!Частота пропусков строки.
	cout << TimeDif << " -- TimeDif" << endl;	     //!Время формир.диффузионной области в миллисек.
	cout << TIME_MAX << " -- TIME_MAX" << endl;	     //!Время окончания расчётов.

	//!*****************************************************************
	//!Хранятся в COMMON / CCC /

	cout << cmax << " -- cmax" << endl;	       //!Число Маха(1.095).
	cout << sig1 << ' ' << sig2 << " -- sig1 sig2" << endl;	   //!Эффектив.сечения молекул.(3.66 4.94)(3.19 5.66).
	cout << gam1 << ' ' << gam2 << " -- gam1 gam2" << endl;	   //!Показатели адиабаты газов.(1.667 1.667)(1.4 1.09).
	cout << om1 << ' ' << om2 << " -- om1 om2" << endl;	   //!Масса.молек.1 - 2 - го.газа.(40 131.2)(29.04 146.07).

	cout << d0 << ' ' << anu12 << " -- d0 anu12" << endl;	   //!Диффузия, коэффициент.
	cout << vis10 << ' ' << anu1 << " -- vis10 anu1" << endl;   //!Вязкость 1, степень 1.
	cout << vis20 << ' ' << anu2 << " -- vis10 anu1" << endl;   //!Вязкость 2, степень 2.
	cout << Akur << " -- AKUR" << endl;		   //!Курант.
	cout << concmin << " -- concmin" << endl;	   //!Минимальня концентрация.

	//!*******************************************************************************
	cout << U_shift << " -- U_shift" << endl;      //!Скорость сдвига[м / с](-2.0).
	
*/
