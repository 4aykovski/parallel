#pragma once

//MODULE M7_4_WRITE_RESTART

//SUBROUTINE S7_4_WRITE_RESTART(nxp,nyp)
void A7_39_WritRestart() {

//	int i, k;
	int x, y;
	
//!		write(*,*) 'start zapisi'
	//open(20,file='rr21n.dat',form='unformatted',status='REPLACE')
	Rr21n.open("Rr21n.dat");

//	cout << "tpaz = " << tpaz << endl;
	Rr21n << Time << "\n";
	Rr21n << Nyp << " " << Nxp << "\n";

for(x=1;x<=Nxp;x++) {for(y=1;y<=Nyp;y++) Rr21n << RO1[y][x] << ' ';} 
	Rr21n << "\n\n";
for(x=1;x<=Nxp;x++) {for(y=1;y<=Nyp;y++) Rr21n << RO2[y][x] << ' ';}
	Rr21n << "\n\n";
for(x=1;x<=Nxp;x++) {for(y=1;y<=Nyp;y++) Rr21n << U1[y][x] << ' ';}
	Rr21n << "\n\n";
for(x=1;x<=Nxp;x++) {for(y=1;y<=Nyp;y++) Rr21n << U2[y][x] << ' ';}
	Rr21n << "\n\n";
for(x=1;x<=Nxp;x++) {for(y=1;y<=Nyp;y++) Rr21n << V1[y][x] << ' ';}
	Rr21n << "\n\n";
for(x=1;x<=Nxp;x++) {for(y=1;y<=Nyp;y++) Rr21n << V2[y][x] << ' ';}
	Rr21n << "\n\n";
for(x=1;x<=Nxp;x++) {for(y=1;y<=Nyp;y++) Rr21n << E1[y][x] << ' ';}
	Rr21n << "\n\n";
for(x=1;x<=Nxp;x++) {for(y=1;y<=Nyp;y++) Rr21n << E2[y][x] << ' ';}
	
	Rr21n.close();

	TmRestart = Time + DtRestart;
//!		write(*,*) 'Final zapisi  t(mkc)=',Time*tpaz

} //END SUBROUTINE S7_4_WRITE_RESTART



  
  /*	write(20) Time, nyp, nxp, ((ro1(K, i), i = 1, nxp), K = 1, nyp),
	 *((ro2(K,i),i=1,nxp),K=1,nyp),
	 *((u1(K,i),i=1,nxp),K=1,nyp),
	 *((u2(K,i),i=1,nxp),K=1,nyp),
	 *((v1(K,i),i=1,nxp),K=1,nyp),
	 *((v2(K,i),i=1,nxp),K=1,nyp),
	 *((e1(K,i),i=1,nxp),K=1,nyp),
	 *((e2(K,i),i=1,nxp),K=1,nyp)
		close(20) */

