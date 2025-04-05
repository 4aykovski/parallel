#pragma once

//#define W 50
void A6_1_SaverGrid() {
	int x, y;

	ofstream hx("hx.dat");
	hx.precision(8);
	hx << Nxp << endl;

	for(x=1; x<=Nxp; x++) { 
//		hx << scientific << XVV[i] * 100 << endl;
		hx << fixed << XVV[x] * 100 << endl;
	}
	hx.close();  // Закрываем файл

	ofstream hy("hy.dat");
	hy.precision(8);
	hy << Nyp << endl;

	for(y=3; y<=Nyp; y++) {
		hy << scientific << YVV[y] * 100 << endl;
	}
	hy.close();  // Закрываем файл

	ofstream hxX("hxX.dat");
	hxX.precision(8);
	hxX << Nxp << endl;

	for(x=1; x<=Nxp; x++) {
		hxX << scientific << HX[x] * 100 << endl;
	}
	hxX.close();  // Закрываем файл

	ofstream hyY("hyY.dat");
	hyY.precision(8);
	hyY << Nyp << endl;

	for(y=3; y<=Nyp; y++) {
		hyY << scientific << HY[y] * 100 << endl;
	}
	hyY.close();  // Закрываем файл


	/*
	open(12, file = 'hx.dat')
	open(13, file = 'hy.dat')
	open(120, file = 'hxX.dat')
	open(130, file = 'hyY.dat')

	write(12, *) nxp
	DO I = 1, NXP
		write(12, 2) xvv(i)* dl * 1000
	END DO
	write(13, *) nyp
	DO K = 3, NYP
		write(13, 2) yvv(K) * dl * 1000
	END DO

	write(120, *) nxp
	DO I = 1, NXP
		write(120, 2) HX(i) * dl * 1000
	END DO
	write(130, *) nyp
	DO K = 3, NYP
		write(130, 2) HY(K) * dl * 1000
	END DO
	//!  		write(12,*) nxp
	//!  		write(12,2) (xvv(i)*dl*1000,i=1,nxp)
	//!  		write(13,*) nyp
	//!  		write(13,2) (yvv(K)*dl*1000,K=3,NYp)
	*/
}
