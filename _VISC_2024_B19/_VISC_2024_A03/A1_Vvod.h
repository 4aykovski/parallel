#pragma once

void next(ifstream& vvod) {
	vvod.clear(); vvod.ignore(INT_MAX, '\n');
}

void A1_Vvod() {
//	HX[5] = 2.5;
ifstream vvod("__Vvod_Cpp.dat");
/*if (!vvod) {
	cout << "No file"; return -1;
}*/
vvod >> TypeTau; next(vvod);	
	//Тип реализации по времени (p-ПредикторКорректор, 1-Эйлер, 2,3,4,5 Рунге-Кутт).
vvod >> Akur; next(vvod);  	    //(0.9)- Курант.
vvod >> h0; next(vvod);			//Шаг сетки, мм (0.25).

vvod >> xN >> xK; next(vvod);	//Лев.и прав.гр.осн.области, мм (-7.0 4.0).
vvod >> yN >> yK; next(vvod);	//Ниж.и верх.гр.осн.области, мм (-7.0 7.0).
vvod >> dd; next(vvod);			//Диаметр капель до диффузии, мм (3.0).
vvod >> x3 >> dSW; next(vvod);			//Положение ударной волны, мм (-6.0).

vvod >> sd; next(vvod);			//Расст.м-ду.центрами (безразмерная).
vvod >> KONF; next(vvod);		//Тип конфигурации - 1-(1С) 2-(2С) 3-(3С) 4-(1С-UP) 5-(2С-UP).
vvod >> X3C; next(vvod);		//Отступ удалённых от ударной волны капель (относительный).
vvod >> am10p; next(vvod);		//Конценр.лёг.газа перед и за волной  (0 1).
vvod >> alfa101 >> alfa102 >> alfa103; next(vvod);
	
vvod >> Gxy1 >> Gxy2 >> Gxy3; next(vvod);	//Ширина доп.слоёв, мм (5 10 200).
vvod >> qm1 >> qm2 >> qm3; next(vvod);		//(1.05 1.1 1.2) - коэффициент для шага.
vvod >> BoundCondLeft >> BoundCondRight; next(vvod);				//(1 1) - Усл.на гр.по х(лев., пр)(1 - своб., 0 - тв.стенка).
vvod >> BoundCondLower; next(vvod);						//(2)- УСЛОВИЕ НА НИЖНЕЙ ГРАНИЦЕ        1 - СВ. 2 - СИММЕТР.УСЛ.

vvod >> DtRestart; next(vvod);          //(10) - Частота записей для Restart(мкс).
vvod >> DtConsol; next(vvod);	        //(10) - Частота сообщений на Консоли(Экране)(dt).
vvod >> DtWriteStep; next(vvod);	    //(10) - Частота вывода результатов(мкс)(dxp).
vvod >> DtEmptyString; next(vvod);	    //(100) - Частота пропусков строки.
vvod >> TimeDif; next(vvod);	        //(35) - Время формир.диффузионной области в миллисек.
vvod >> TIME_MAX; next(vvod);	        //(300) - Время окончания расчётов.

vvod >> cmax; next(vvod);               //(1.2) - Число Маха(1.095).
vvod >> sig1 >> sig2; next(vvod);       //(3.01 5.05) - Эффективные сечения молекул.(3.66 4.94)(3.19 5.66).
vvod >> gam1 >> gam2; next(vvod);       //(1.4 1.09) - Показатели адиабаты газов.(1.667 1.667)(1.4 1.09).
vvod >> om1 >> om2; next(vvod);         //(29.04 146.07) - Масса молекул 1 - го, 2 - го газа(40 131.2)(29.04 146.07).
vvod >> d0 >> anu12; next(vvod);        //(0.097  0.786) - Диффузия, коэффициент.

vvod >> vis10 >> vis20; next(vvod);     //(18.26 14.75) - Вязкость 1, Вязкость 2.
vvod >> anu1 >> anu2; next(vvod);       //(0.786 0.8)   - степень 1, степень 2.
vvod >> visS >> visVol1 >> visVol2; next(vvod); //(0.0 0.0) - Вторая (объёмная) Вязкость первого и второго газов.

vvod >> concmin; next(vvod);	        //(0.0001) - Минимальня концентрация.
vvod >> U_shift; next(vvod);            //(1.788391) - скорость сдвига[м / с].


cout << "TypeTau = " << TypeTau << endl;
cout << "Akur = " << Akur << endl;
cout << "h0 = " << h0 << " mm." << endl;
cout << "xSW = " << x3 << " mm." << endl;

}


