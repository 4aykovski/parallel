#pragma once

void WriteHeadOneFile(ofstream&  File) {
	File << HeadFiles;
	File << Time * tpaz << endl;
}

void A6_2_WriteStep() {

//	FTime << Time * tpaz << ' ' << Time*tpaz+1 << endl;

	WriteHeadOneFile(Lir1);
	WriteHeadOneFile(LirP);
	WriteHeadOneFile(LirRo);
	WriteHeadOneFile(LirT1);
	WriteHeadOneFile(LirT2);
	WriteHeadOneFile(LirU1);
	WriteHeadOneFile(LirU2);
	WriteHeadOneFile(LirV1);
	WriteHeadOneFile(LirV2);

int x, y;
//int i, k;
//	FTime << scientific << HY[i] * 100 << endl;
	

	for(x=1; x<=Nxp; x++) { //DO I = 1, NXp
		for(y=3; y<=Nyp; y++) { //DO K = 1, NYp
LirRo << RO1[y][x] + RO2[y][x] << ' ';
Lir1 << RO1[y][x] / (RO1[y][x] + RO2[y][x] / om) << ' ';
LirP << RO1[y][x]*E1[y][x]*(gam1-1) + RO2[y][x]*E2[y][x]*(gam2-1) << ' ';
LirV1 << V1[y][x] * dL * 1000 / tpaz << ' ';
LirU1 << U1[y][x] * dL * 1000 / tpaz << ' ';

LirV2 << V2[y][x] * dL * 1000 / tpaz << ' ';
LirU2 << U2[y][x] * dL * 1000 / tpaz << ' ';
LirT1 << E1[y][x] / c1v << ' ';
LirT2 << E2[y][x] / c2v << ' ';
		}
	}

	Lir1 << endl;
	LirP << endl;
	LirRo << endl;
	LirT1 << endl;
	LirT2 << endl;
	LirU1 << endl;
	LirU2 << endl;
	LirV1 << endl;
	LirV2 << endl;
//	open(2, file = 'LIr1.dat')
//	open(3, file = 'LIrp.dat')
//	open(4, file = 'LIrv1.dat')
//	open(5, file = 'LIru1.dat')
//	open(6, file = 'LIrro.dat')
//	open(7, file = 'LIrv2.dat')
//	open(8, file = 'LIru2.dat')
//	open(10, file = 'LIrt1.dat')
//	open(11, file = 'LIrt2.dat')
}
