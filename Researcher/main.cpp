//xi defecte xi * 0.5 v=0.01

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

//Geometry of computation window and time parameters, in units of Dist*lex

#define N 81 //Size in X direction																		
#define M 81 //Size in Y direction

int n=N-1; //N-1
int m=M-1; //M-1
int l=(N-1)/2; //n/2
int L=(M-1)/2; //m/2
double Dist=0.5; //Distance between divisions in exchange lengths
double invDist=2.0;//(1/0.5);				// inversa de Dist, si cambio Dist, cambiar això!!!!!!!!!!!!!!!
int Time=40000; //Limit time of RK4
int Time2=100000000; //Demomento no usas este
double h=0.001; // Time step

//Visual parameters
double Angulo = 0*3.14159265/10; //Angulo del vector D
double Angulo2 = 0*3.14159265/10; //Angulo del corrente
double Angulo3 = 4*3.14159265/10; //Angulo SPAWN skyrmion
double RRR = 15; //Radio SPAWN
int vort = -1; // 1 for skyrmions -1 for antiskyrmions
double R=2, RR=2; //radios of skyrmions and deffects
double Ancho=0; //Ancho de la linea DEFECTO
double Position = 35; // Inicio de linea DEFECTO

//Physical parameters

double xi=1.68;//1.68;//1.290238; //DM88 parameter (xi) 0-1.68/max2
double xiInc=0.0;
double kappa=2.0;   //Anisotropy parameter (kappa)  0-2.5 (es pot baixar un poc)
double alpha=0.3; //Gilbert damping constant (alpha) 0.1-0.3 (té relació amb l'angle del skyrmion hall effect) (I té relació amb lo rápid que arribem a l'estabilitat perquè el DM es fa mes fort)
double alpha2=0;// 1/(1+alpha^2)
double hex = -1000.0; //external applied field   1 es mig Tesla aprox
int xxx= 40; //position of the dipole
int yyy= 40;
int zzz= 6;
double hexx[N][M];

//Current parameters
double sigx=0.30*(cos(Angulo2-3.14159265/2) + sin (Angulo2-3.14159265/2));
double sigy=0.30*(cos(Angulo2-3.14159265/2) - sin (Angulo2-3.14159265/2));
double sigz=0.0;

double dj=0.00;



//Definition of energys: E=Total E1=Exchange E2=DM E3=Anisotropy E4=Zeeman

double E=0;
double E1=0;
double E2=0;
double E3=0;
double E4=0;
double E0=0;
double E01=0;
double E02=0;
double E03=0;
double E04=0;
double Magnx, Magny, Magnz;
double mx[N][M], my[N][M], mz[N][M];
double mxx[N][M], myy[N][M], mzz[N][M];
double mx0[N][M], my0[N][M], mz0[N][M];
double hx[N][M], hy[N][M], hz[N][M];
double phi[N][M];
double xixi[N][M];
double mdotsig, mdoth;





//FUNCIO VALOR ABSOLUT
double Abs(double KKK){

	if (KKK>=0){
		return KKK;
	}

	else{
		return -KKK;
	}

}





// FUNCIO ANGLE QUE NO S'USA EN CARTESIANAS SUPOSO
double Angle(){

	int i, j;
	double x, y;


	for(i=0; i<=n; i++){
		for(j=0; j<=m; j++){

			x=i-l;
			y=j-L;

			if(x==0){
				phi[i][j]=2.0*atan(1.0)*((Abs(y))/(y));
			}

			else{
				phi[i][j]=atan(y/x);
			}
	
			if((x<0) & (y!=0)){
				phi[i][j]=(4.0*atan(1.0)-Abs(phi[i][j]))*(-phi[i][j]/Abs(phi[i][j]));
			}
			
			else if (x<0){
				phi[i][j]=-4.0*atan(1.0);
			}

		}

	}

	phi[l][L]=0;
	
	return 0;
}






//FUNCIO PER PRINTEAR LA MAGNETITZACIÓ en un instant en forma de matriu
void PrintM(double Matrix[N][M]){

	int i , j;

	for(j=M; j>=0; j--){
		for(i=0; i<N; i++){
	
			printf("%.5lf ", Matrix[i][j] );
			
		}
		
	printf("\n");
	
	}
	
	printf("\n");

}









//Derivatives

//First X

double dx(double v[N][M], int i, int j){						//retorna derivada respecte x de una coordenada i, j

	if(i==n){
		return (v[0][j]-v[n-1][j])*0.5*invDist;
	}

	if(i==0){
		return (v[1][j]-v[n][j])*0.5*invDist;
	}
	
	else{
		return (v[i+1][j]-v[i-1][j])*0.5*invDist;
	}
	
}


//First Y

double dy(double v[N][M], int i, int j){						//retorna derivada respecte y de una coordenada i, j

	if(j==m){
		return (v[i][0]-v[i][j-1])*0.5*invDist;
	}

	if(j==0){
		return (v[i][j+1]-v[i][m])*0.5*invDist;
	}
	
	else{
		return (v[i][j+1]-v[i][j-1])*0.5*invDist;
	}
	
}


//Second X


double dx2(double v[N][M], int i, int j){								//retorna derivada SEGONA respecte x de una coordenada i, j

	if(i==n){
		return (v[0][j]-2.0*v[i][j]+v[i-1][j])*invDist*invDist;
	}

	if(i==0){
		return (v[i+1][j]-2.0*v[i][j]+v[n][j])*invDist*invDist;
	}
	
	else{
		return (v[i+1][j]-2.0*v[i][j]+v[i-1][j])*invDist*invDist;
	}
	
}


//Second Y

double dy2(double v[N][M], int i, int j){								//retorna derivada SEGONA respecte y de una coordenada i, j

	if(j==m){
		return (v[i][0]-2.0*v[i][j]+v[i][j-1])*invDist*invDist;
	}

	if(j==0){
		return (v[i][j+1]-2.0*v[i][j]+v[i][m])*invDist*invDist;
	}
	
	else{
		return (v[i][j+1]-2.0*v[i][j]+v[i][j-1])*invDist*invDist;
	}
}






//Functions that compute the different terms of the effective field

double Exchange(double v[N][M], int i, int j){
	return (dx2(v, i , j)+dy2(v, i, j));
}

double Moriyax (int i, int j){
	return -xixi[i][j]*(dx(mz, i, j)*cos(Angulo)-dy(mz, i, j)*sin(Angulo));
}

double Moriyay (int i, int j){
	return vort*xixi[i][j]*(-dx(mz, i, j)*sin(Angulo)-dy(mz, i, j)*cos(Angulo));
}

double Moriyaz (int i, int j){
	return xixi[i][j]*(dx(mx, i, j)*cos(Angulo)-dy(mx, i, j)*sin(Angulo)+vort*dy(my, i, j)*cos(Angulo)+vort*dx(my, i, j)*sin(Angulo));
}

double Anisotropy (int i, int j){
	return kappa*mz[i][j];
}

double Zeeman (int i, int j){
	return hexx[i][j];
}







//LLG Function

double LLGi(double mi, double mj, double mk, double hi, double hj, double hk, double sigi, double sigj, double sigk){
	return (-(mj*hk-mk*hj)-alpha*(mi*mdoth-hi)-dj*(mi*mdotsig-sigi)+dj*alpha*(mj*sigk-mk*sigj))/(1+alpha*alpha);
//return -(Pi+alpha*(vj*Pk-vk*Pj)+(1+alpha*beta)*ui)/(1+alpha*alpha);
}






//Effective field for one coordinate i, j
void heff(int i, int j){

	hx[i][j]=Exchange(mx, i , j)+Moriyax(i, j);
	hy[i][j]=Exchange(my, i , j)+Moriyay(i, j);
	hz[i][j]=Exchange(mz, i , j)+Moriyaz(i, j)+Anisotropy(i, j)+Zeeman(i, j);

}






//Normalize Function
void Normalize( double vx, double vy, double vz){

	double Mod;

	Mod=sqrt(vx*vx+vy*vy+vz*vz);

	vx=vx/Mod;
	vy=vy/Mod;
	vz=vz/Mod;

}






//MAIN CICLE DEL PROGRAMA
void TimeEvolution(){
	
	double k1x[N][M], k1y[N][M], k1z[N][M];
	double k2x[N][M], k2y[N][M], k2z[N][M];
	double k3x[N][M], k3y[N][M], k3z[N][M];
	double k4x[N][M], k4y[N][M], k4z[N][M];
	double mxx, myy, mzz;
	double Px, Py, Pz, ux, uy, uz, Mod; // components of the vectorial product m x h
	char buffer [100];
	FILE * sFile;
	
	int i, tout, j, t, k=0;
	
	tout=Time/100;
	
	
	//CREACIÓ DE CARPETAS I IMPRESSIÓ DEL ARCHIUS
	
	//sprintf(buffer,"mkdir %d", Time);
	
	
	sprintf(buffer,"mkdir %d_%dx%d_kappa%.2f_alpha%.2f_hex%.2f_Vort%d_Angulo%.2f_sigs%.2f_%.2f_%.2f_dj%.2f_RR%.2f_Deffect%.2f",  Time, N, M, kappa, alpha, hex, vort, Angulo3, sigx, sigy, sigz, dj, RR,  xiInc);
	system(buffer);	
	
	for(t=0; t<=Time; t++){
				
		if((t%tout)==0){
			//sprintf(buffer,"%d/Sky%d.txt", Time, k);
			sprintf(buffer,"%d_%dx%d_kappa%.2f_alpha%.2f_hex%.2f_Vort%d_Angulo%.2f_sigs%.2f_%.2f_%.2f_dj%.2f_RR%.2f_Deffect%.2f/Sky%d.txt",  Time,  N, M, kappa, alpha, hex, vort, Angulo3, sigx, sigy, sigz, dj, RR,  xiInc, k);
			k=k+1;
			sFile = fopen (buffer , "w" );
			printf("t= %d\n", k);
	
	
			for (j=0; j<=m; j++){
				for (i=0; i<=n; i++){
				
					fprintf(sFile,"%d %d %lf %lf %lf \n", i, j, mz[i][j], mx[i][j], my[i][j]);
					
				}
				 
			}
			
			fclose(sFile);
			
		}	
			
		
		for(i=0; i<=n; i++){
			for(j=0; j<=m; j++){
		
				//Compute the effective field
				heff(i, j);
				
				mxx=mx[i][j];
				myy=my[i][j];
				mzz=mz[i][j];
				
				mdoth=mxx*hx[i][j]+myy*hy[i][j]+mzz*hz[i][j];
				mdotsig=mxx*sigx+myy*sigy+mzz*sigz;
				//First constant
				
				
				k1x[i][j]=LLGi(mxx, myy, mzz, hx[i][j], hy[i][j], hz[i][j], sigx, sigy, sigz);
				k1y[i][j]=LLGi(myy, mzz, mxx, hy[i][j], hz[i][j], hx[i][j], sigy, sigz, sigx);
				k1z[i][j]=LLGi(mzz, mxx, myy, hz[i][j], hx[i][j], hy[i][j], sigz, sigx, sigy);
				
				mxx=mx[i][j]+0.5*h*k1x[i][j];
				myy=my[i][j]+0.5*h*k1y[i][j];
				mzz=mz[i][j]+0.5*h*k1z[i][j];
				
				
				//Second Constant
				mdoth=mxx*hx[i][j]+myy*hy[i][j]+mzz*hz[i][j];
				mdotsig=mxx*sigx+myy*sigy+mzz*sigz;
				
				k2x[i][j]=LLGi(mxx, myy, mzz, hx[i][j], hy[i][j], hz[i][j], sigx, sigy, sigz);
				k2y[i][j]=LLGi(myy, mzz, mxx, hy[i][j], hz[i][j], hx[i][j], sigy, sigz, sigx);
				k2z[i][j]=LLGi(mzz, mxx, myy, hz[i][j], hx[i][j], hy[i][j], sigz, sigx, sigy);
				
				mxx=mx[i][j]+0.5*h*k2x[i][j];
				myy=my[i][j]+0.5*h*k2y[i][j];
				mzz=mz[i][j]+0.5*h*k2z[i][j];
				//Third Constant
				
				mdoth=mxx*hx[i][j]+myy*hy[i][j]+mzz*hz[i][j];
				mdotsig=mxx*sigx+myy*sigy+mzz*sigz;
				
				k3x[i][j]=LLGi(mxx, myy, mzz, hx[i][j], hy[i][j], hz[i][j], sigx, sigy, sigz);
				k3y[i][j]=LLGi(myy, mzz, mxx, hy[i][j], hz[i][j], hx[i][j], sigy, sigz, sigx);
				k3z[i][j]=LLGi(mzz, mxx, myy, hz[i][j], hx[i][j], hy[i][j], sigz, sigx, sigy);
				
				mxx=mx[i][j]+h*k3x[i][j];
				myy=my[i][j]+h*k3y[i][j];
				mzz=mz[i][j]+h*k3z[i][j];
				
				//Fourth Constant
				
				mdoth=mxx*hx[i][j]+myy*hy[i][j]+mzz*hz[i][j];
				mdotsig=mxx*sigx+myy*sigy+mzz*sigz;
				
				k4x[i][j]=LLGi(mxx, myy, mzz, hx[i][j], hy[i][j], hz[i][j], sigx, sigy, sigz);
				k4y[i][j]=LLGi(myy, mzz, mxx, hy[i][j], hz[i][j], hx[i][j], sigy, sigz, sigx);
				k4z[i][j]=LLGi(mzz, mxx, myy, hz[i][j], hx[i][j], hy[i][j], sigz, sigx, sigy);
		
			}
		}
		
		//PrintM(k1x);
		//PrintM(k1y);
		//PrintM(k1z);
		
		//printf("%d\n", t);
		
		//Update m and normalize
		
		for(i=0; i<=n; i++){
			for(j=0; j<=m; j++){
		
		
				mx[i][j]+=(h/6.0)*(k1x[i][j]+2.0*k2x[i][j]+2.0*k3x[i][j]+k4x[i][j]);
				my[i][j]+=(h/6.0)*(k1y[i][j]+2.0*k2y[i][j]+2.0*k3y[i][j]+k4y[i][j]);
				mz[i][j]+=(h/6.0)*(k1z[i][j]+2.0*k2z[i][j]+2.0*k3z[i][j]+k4z[i][j]);
				
				Mod=mx[i][j]*mx[i][j]+my[i][j]*my[i][j]+mz[i][j]*mz[i][j];
				mx[i][j]=mx[i][j]/sqrt(Mod);
				my[i][j]=my[i][j]/sqrt(Mod);
				mz[i][j]=mz[i][j]/sqrt(Mod);
		
			}
		}
			
	}

}

void Initialize(){

	int i, j;
	double x, y, xx, yy;
	double Mod;
	Angle();
	//Condicions inicials
	

	
	
	for (i=0; i<=n; i++){
		for (j=0; j<=m; j++){			
			
			xixi[i][j]=xi;
				
			if (i>=Position) {
				if (i<=Position+Ancho){
					xixi[i][j]=xi+xiInc;
				}
			}
			
			if (i>=Position+5) {
				if (i<=Position+Ancho+5){
					xixi[i][j]=xi;
				}
			}
					
			x=Dist*(i-l);
			y=Dist*(j-L);
			
			if (x*x+y*y <= R*R) {
	
				mzz[i][j]=cos((x*x+y*y)*(3.14159265/2)/(R*R));
				mxx[i][j]=-sin((x*x+y*y)*(3.14159265/2)/(R*R))*sin(phi[i][j]-3.14159265/2 + Angulo);
				myy[i][j]=sin((x*x+y*y)*(3.14159265/2)/(R*R))*cos(phi[i][j]-vort*3.14159265/2 + Angulo);
				
				hx[i][j]=0;
				hy[i][j]=0;
				hz[i][j]=0;
			
			}
			else {
				mzz[i][j]=-1;
				mxx[i][j]=0;
				myy[i][j]=0;
				
				hx[i][j]=0;
				hy[i][j]=0;
				hz[i][j]=0;
			}

			//EXTERNAL FIELD
			hexx[i][j] = hex*zzz / ( ((i-xxx)*(i-xxx)+(j-yyy)*(j-yyy)+zzz*zzz)*((i-xxx)*(i-xxx)+(j-yyy)*(j-yyy)+zzz*zzz) );
					
		}		
	}
	
	int	xxxx=(int)round(RRR*cos(Angulo3));
	int yyyy=(int)round(RRR*sin(Angulo3));
	
	for (i=0; i<n; i++){
		for (j=0; j<=m; j++){	
			mz[i][j]=mzz[i][j];
			mx[i][j]=mxx[i][j];
			my[i][j]=myy[i][j];
		}
	}
	
	for (i=20; i<=n-20; i++){
		for (j=20; j<=m-20; j++){	
		
			mz[i][j]=mzz[i-xxxx][j-yyyy];
			mx[i][j]=mxx[i-xxxx][j-yyyy];
			my[i][j]=myy[i-xxxx][j-yyyy];
			
		}
	}
	
	
}






//EXECTUABLE DEL PROGRAMA

int main(){
	
	Initialize();
	TimeEvolution();
	return 0;

}

