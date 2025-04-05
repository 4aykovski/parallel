#pragma once

#include "A6_1_SaverGrid.h"
#include "A6_2_WriteStep.h"

void ReadReserv(); // Œ√ƒ¿ œ–Œ÷≈—— œŒ…ƒ®“ - «¿…Ã”—‹ –≈«≈–¬Õ€Ã œ≈–≈—Œ’–¿Õ≈Õ»≈Ã!!!!!!!

void A6_0_WriteData0() {
//cout << "M6_0_WriteData0()" << endl;
	int x, y;
//! write(*,*) 'S6_WRITE_DATA0.for'
	cout << "" << endl;
//	cout << "vvod from files: ";  //<< endl;
//	cin >> nf;

	if(nf == 1) { //if (nf .eq. 1) then
		ReadReserv();
		FTime.open("Time.dat", ios::app);
		Lir1.open("Lir1.dat", ios::app);
		LirP.open("LirP.dat", ios::app);
		LirRo.open("LirRo.dat", ios::app);
		LirT1.open("LirT1.dat", ios::app);
		LirT2.open("LirT2.dat", ios::app);
		LirU1.open("LirU1.dat", ios::app);
		LirU2.open("LirU2.dat", ios::app);
		LirV1.open("LirV1.dat", ios::app);
		LirV2.open("LirV2.dat", ios::app);
	}
	
	else {
		FTime.open("Time.dat");
		Lir1.open("Lir1.dat");
		LirP.open("LirP.dat");
		LirRo.open("LirRo.dat");
		LirT1.open("LirT1.dat");
		LirT2.open("LirT2.dat");
		LirU1.open("LirU1.dat");
		LirU2.open("LirU2.dat");
		LirV1.open("LirV1.dat");
		LirV2.open("LirV2.dat");

		nf = 0;
		cout << "Otkroem vse faily v standartnom statuse (Peresapis)." << endl;
		system("pause"); //«‡‰ÂÊÍ‡ ˝Í‡Ì‡.
		A6_1_SaverGrid(); //¬ ÙÓÚ‡ÌÂ ÒÂÚÍ‡ Á‡ÔËÒ˚‚‡Î‡Ò¸ ‚ ˝ÚÓÈ ÔÓ‰ÔÓ„‡ÏÏÂ. —ÂÈ˜‡Ò ÓÚ‰ÂÎ¸ÌÓ.
		A6_2_WriteStep();
		
		//ofstream FTime("Time.dat"); FTime.precision(8);
		/*
		//!        open( 9, file = 'xp.dat' )
		open(24, file = '_Time.dat')

		open(2, file = 'LIr1.dat')
		open(3, file = 'LIrp.dat')

		open(4, file = 'LIrv1.dat')
		open(5, file = 'LIru1.dat')
		open(6, file = 'LIrro.dat')

		open(7, file = 'LIrv2.dat')
		open(8, file = 'LIru2.dat')

		open(10, file = 'LIrt1.dat')
		open(11, file = 'LIrt2.dat')


		CALL S6_1_WRITE_STEP(Time, tpaz, om, c1v, c2v, gam1, gam2)
		M6_2S_WriteStep();
	*/
	} //endif:(nf == 1).
	
	for(y=1; y<=Nyp; y++) { //do k = 1,nyp
		for(x=1; x<=Nxp; x++) { //do i = 1,nxp
			KRO1[0][y][x] = RO1[y][x];
			KRO2[0][y][x] = RO2[y][x];
			KU1[0][y][x]  = U1[y][x];
			KU2[0][y][x]  = U2[y][x];
			KV1[0][y][x]  = V1[y][x];
			KV2[0][y][x]  = V2[y][x];
			KE1[0][y][x]  = E1[y][x];
			KE2[0][y][x]  = E2[y][x];
	} //enddo
} //enddo

//2	FORMAT (500E18.6)
} //END SUBROUTINE S6_0_WRITE_DATA0


void ReadReserv() {
//	int i, k;
	int x, y;
	double a;
	cout << "OTKRYTIE' VSEH FAILY V STATUSE 'OLD'." << endl;
	system("pause"); //«‡‰ÂÊÍ‡ ˝Í‡Ì‡.
	//Rr21n.open("Rr21n.dat", ios::in);
	ifstream Rr21n("Rr21n.dat");
	Rr21n >> Time;
	TimePaz = Time * tpaz;
	Rr21n >> Nyp;
	Rr21n >> Nxp;
	for(x=1; x<=Nxp; x++) {for(y=1; y<=Nyp; y++) Rr21n >> RO1[y][x];}
	for(x=1; x<=Nxp; x++) {for(y=1; y<=Nyp; y++) Rr21n >> RO2[y][x];}
	for(x=1; x<=Nxp; x++) {for(y=1; y<=Nyp; y++) Rr21n >> U1[y][x];}
	for(x=1; x<=Nxp; x++) {for(y=1; y<=Nyp; y++) Rr21n >> U2[y][x];}
	for(x=1; x<=Nxp; x++) {for(y=1; y<=Nyp; y++) Rr21n >> V1[y][x];}
	for(x=1; x<=Nxp; x++) {for(y=1; y<=Nyp; y++) Rr21n >> V2[y][x];}
	for(x=1; x<=Nxp; x++) {for(y=1; y<=Nyp; y++) Rr21n >> E1[y][x];}
	for(x=1; x<=Nxp; x++) {for(y=1; y<=Nyp; y++) Rr21n >> E2[y][x];}
	//	Rr21n << "\n\n";

	Rr21n.close();
	//ifstream vvod("rr21n.dat");
	
	/*open(20, file = 'rr21n.dat', form = 'unformatted', status = 'old')

	//!	read(20)Time
	//!	read(20)nyp,nxp
		read(20)Time,nyp,nxp,((ro1(K,i),i=1,nxp),K=1,nyp),
	 *((ro2(K,i),i=1,nxp),K=1,nyp),
	 *((u1(K,i),i=1,nxp),K=1,nyp),
	 *((u2(K,i),i=1,nxp),K=1,nyp),
	 *((v1(K,i),i=1,nxp),K=1,nyp),
	 *((v2(K,i),i=1,nxp),K=1,nyp),
	 *((e1(K,i),i=1,nxp),K=1,nyp),
	 *((e2(K,i),i=1,nxp),K=1,nyp)
		close(20)

		TmConsol  = Time + DtConsol
		TmWriteStep = Time + DtWriteStep
		TmRestart = Time + DtRestart

		open(24, file = '_Time.dat',ACCESS='APPEND', STATUS='OLD')
		open(2, file = 'LIr1.dat',ACCESS='APPEND', STATUS='OLD')
		open(3, file = 'LIrp.dat',ACCESS='APPEND', STATUS='OLD')
		open(4, file = 'LIrv1.dat',ACCESS='APPEND', STATUS='OLD')
		open(5, file = 'LIru1.dat',ACCESS='APPEND', STATUS='OLD')
		open(6, file = 'LIrro.dat',ACCESS='APPEND', STATUS='OLD')

		open(7, file= 'LIrv2.dat',ACCESS='APPEND', STATUS='OLD')
		open(8, file= 'LIru2.dat',ACCESS='APPEND', STATUS='OLD')
		open(10,file='LIrt1.dat',ACCESS='APPEND', STATUS='OLD')
		open(11,file='LIrt2.dat',ACCESS='APPEND', STATUS='OLD')
	*/
}