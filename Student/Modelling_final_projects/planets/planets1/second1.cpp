#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>


double f(double rj, double ri, double rj_bad1,double ri_bad1,double rj_bad2, double ri_bad2, double mj){
    const double pi = 4.0 * atan(1.0); //Pi
    double resultf; 
    resultf = 4*pi*pi*((rj-ri)*mj/(0.000001+pow((rj-ri)*(rj-ri)+(rj_bad1-ri_bad1)*(rj_bad1-ri_bad1)+(rj_bad2-ri_bad2)*(rj_bad2-ri_bad2),(3/2)))); 			//falta fotre els +1/2ki tamb� al denominador fent el sumatori en comptes de fent-ho aix�
    return resultf;
}

double U(double rj, double ri, double rj_bad1,double ri_bad1,double rj_bad2, double ri_bad2, double mj, double mi){
	const double pi = 4.0 * atan(1.0); //Pi
    double resultf; 
    resultf =4*pi*pi*mi*mj/(0.001+pow((rj-ri)*(rj-ri)+(rj_bad1-ri_bad1)*(rj_bad1-ri_bad1)+(rj_bad2-ri_bad2)*(rj_bad2-ri_bad2),0.5)); 			//falta fotre els +1/2ki tamb� al denominador fent el sumatori en comptes de fent-ho aix�
    return resultf;
}

double fi(double k1, double x, double y){
	double f_post; 
	f_post = x-k1*pow(x, 2)-x*y;
	return f_post;
}

double g(double k2, double k3, double x, double y){
    double g_post;
    g_post = k2*x*y-k3*y;
	return g_post;
}



int main(){
	
	//Definim les variables
	int N = 4;	//number of bodies
	int T = 126;
	int D = 3; //number of dimensions
	
	int n[N];
	double m[N];
	int t;
	int d[D];
	
	double dt = 0.01;
	double dx = 0.01;
	double R[N][D];
	double VR[N][D];
	double K1[N][D];
	double L1[N][D];
	double K2[N][D];
	double L2[N][D];
	double K3[N][D];
	double L3[N][D];
	double K4[N][D];
	double L4[N][D];
	
	
	//Imposem condicions inicials:
	
	//BLUE PESAT
    R[0][0]=1.0;  
    R[0][1]=0.0;
    R[0][2]=1.0;
	m[0]=1.0;

	VR[0][0]=0.0;
    VR[0][1]=0.0;
    VR[0][2]=0.0;

	//RED
    R[1][0]=1.0;
    R[1][1]=0.0;
 	R[1][2]=-1.0;
	m[1]=1.0;
	
	VR[1][0]=0.0;
    VR[1][1]=-1.0;
    VR[1][2]=0.0;

	//YELLOW
 	R[2][0]=-1.0;
    R[2][1]=1.0;
    R[2][2]=0.0;
	m[2]=1.0;
	
	VR[2][0]=1.0;
    VR[2][1]=0.0;
    VR[2][2]=0.0;
    
	//BLACK PLANETA
    R[3][0]=0.0;
    R[3][1]=1.0;
 	R[3][2]=1.0;
	m[3]=1.0;
	
    VR[3][0]=0.0;
    VR[3][1]=0.0;
    VR[3][2]=0.0;
	
	//CCII Energies
	double Ectotali=0;
	double Ec[N];
	for (int i=0; i<N;i++){
		Ec[i]=0.5*m[i]*(pow(VR[i][0],2)+pow(VR[i][1],2)+pow(VR[i][2],2));	
	}
	for (int l=0; l<N;l++){
		Ectotali=Ectotali+Ec[l];	
	}
	printf("Energia cinetica inicial=%lf\n",Ectotali);
	
	
	double Eptotali=0;
	double Epi[N];
	for (int k=0; k<N;k++){
		double Epib=0;
		n[0]=k;
		for (int i=1;i<N;i++){
			n[i]=(n[0]+i)%N;												
			Epib=Epib-0.5*U(R[n[i]][0],R[n[0]][0],R[n[i]][1],R[n[0]][1],R[n[i]][2],R[n[0]][2],m[n[i]],m[n[0]]);
		}
		Epi[k]=Epib;
	}
	for (int l=0; l<N;l++){
		Eptotali=Eptotali+Epi[l];	
	}
	printf("Energia potencial total=%lf\n",Eptotali);
	



	//printeamos las condiciones de contorno
	for (n[0]=0;n[0]<N;n[0]++){
		char buffer0[80]={0}; 							
		sprintf(buffer0,"animacion_cuerpo%d.txt",n[0]);
	    FILE * fpointer0 = fopen (buffer0,"w"); 
	    fprintf(fpointer0,"%lf	%lf	%lf	i\n", R[n[0]][0], R[n[0]][1], R[n[0]][2]);
	}


	//PREDATOR%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//Definim array de temperatures pel predator, array de preses, array de depredadors i array de temps (Majúscules desadimensionalitzat)
	double augment = 100.0;
	//double augment = 100000000.0;
	int augment2 = augment;
	int I = (T-1)*augment; //NO CANVIAR EL 10 A MENYS QUE ES CANVIÏ LA h DEL RK-4 PREDATOR. Multiplicar sempre per h^-1.
	double TEMP[I];
	double x[I];
	double y[I];
	double temps[I];
	double X[I];
	double Y[I];
	double TEMPS[I];

	//Condicions inicials PREDATOR
	x[1] = 500;
	y[1] = 490;
	temps[0] = 0;



	
	//EMPIEZA LO CHULO	
	for (t=1;t<T;t++){		//Avanzamos el tiempo			
		for (n[0]=0;n[0]<N;n[0]++){
			for (d[0]=0;d[0]<D;d[0]++){
				
				for (int i=1; i<N;i++){
					n[i]=(n[0]+i)%N;
				}
				
				for (int i=1; i<D;i++){
					d[i]=(d[0]+i)%D;
				}
				
				double g=0.0;
				for (int i=1; i<N;i++){
					g=g+f(R[n[i]][d[0]],R[n[0]][d[0]],R[n[i]][d[1]],R[n[0]][d[1]],R[n[i]][d[2]],R[n[0]][d[2]],m[n[i]]);
				}
				
				K1[n[0]][d[0]]=dt*VR[n[0]][d[0]];
				L1[n[0]][d[0]]=dt*g;
			}
		}	
		for (n[0]=0;n[0]<N;n[0]++){			//aqui ciclamos todas las K
			for (d[0]=0;d[0]<D;d[0]++){
				
				for (int i=1; i<N;i++){
					n[i]=(n[0]+i)%N;
				}
				
				for (int i=1; i<D;i++){
					d[i]=(d[0]+i)%D;
				}
				
				double g=0.0;
				for (int i=1; i<N;i++){
					g=g+f(R[n[i]][d[0]]+K1[n[i]][d[0]]/2, R[n[0]][d[0]]+K1[n[0]][d[0]]/2 ,R[n[i]][d[1]]+K1[n[i]][d[1]]/2, R[n[0]][d[1]]+K1[n[0]][d[1]]/2, R[n[i]][d[2]]+K1[n[i]][d[2]]/2, R[n[0]][d[2]]+K1[n[0]][d[2]]/2, m[n[i]]);	
				}
				
				K2[n[0]][d[0]]=dt*(VR[n[0]][d[0]]+L1[n[0]][d[0]]/2);
				L2[n[0]][d[0]]=dt*g; //falta veure com afegim el K1/2 al modul del denominador y no nomes adalt
			}
		}
		for (n[0]=0;n[0]<N;n[0]++){			//aqui ciclamos todas las K
			for (d[0]=0;d[0]<D;d[0]++){
				
				for (int i=1; i<N;i++){
					n[i]=(n[0]+i)%N;
				}
				
				for (int i=1; i<D;i++){
					d[i]=(d[0]+i)%D;
				}
				
				double g=0.0;
				for (int i=1; i<N;i++){
					g=g+f(R[n[i]][d[0]]+K2[n[i]][d[0]]/2, R[n[0]][d[0]]+K2[n[0]][d[0]]/2, R[n[i]][d[1]]+K2[n[i]][d[1]]/2, R[n[0]][d[1]]+K2[n[0]][d[1]]/2, R[n[i]][d[2]]+K2[n[i]][d[2]]/2, R[n[0]][d[2]]+K2[n[0]][d[2]]/2, m[n[i]]);
				}
				
				K3[n[0]][d[0]]=dt*(VR[n[0]][d[0]]+L2[n[0]][d[0]]/2);
				L3[n[0]][d[0]]=dt*g;
			}
		}
		for (n[0]=0;n[0]<N;n[0]++){			//aqui ciclamos todas las K
			for (d[0]=0;d[0]<D;d[0]++){
				
				for (int i=1; i<N;i++){
					n[i]=(n[0]+i)%N;
				}
				
				for (int i=1; i<D;i++){
					d[i]=(d[0]+i)%D;
				}
				
				double g=0.0;
				for (int i=1; i<N;i++){
					g=g+f(R[n[i]][d[0]]+K3[n[i]][d[0]], R[n[0]][d[0]]+K3[n[0]][d[0]], R[n[i]][d[1]]+K3[n[i]][d[1]], R[n[0]][d[1]]+K3[n[0]][d[1]], R[n[i]][d[2]]+K3[n[i]][d[2]], R[n[0]][d[2]]+K3[n[0]][d[2]], m[n[i]]);
				}
				
				K4[n[0]][d[0]]=dt*(VR[n[0]][d[0]]+L3[n[0]][d[0]]);
				L4[n[0]][d[0]]=dt*g;
			}
		}
		for (n[0]=0;n[0]<N;n[0]++){
			for (d[0]=0;d[0]<D;d[0]++){
				R[n[0]][d[0]]=R[n[0]][d[0]]+(K1[n[0]][d[0]]+2*K2[n[0]][d[0]]+2*K3[n[0]][d[0]]+K4[n[0]][d[0]])/6;
				VR[n[0]][d[0]]=VR[n[0]][d[0]]+(L1[n[0]][d[0]]+2*L2[n[0]][d[0]]+2*L3[n[0]][d[0]]+L4[n[0]][d[0]])/6;
			}
		}
		
		
		
		
		//Aqu� se printean los datos de animaci�n para cada t	
		int salto = 1;
		double tc= t % salto;  //PUEDES CAMBIAR CADA CUANTOS T IMPRIMES!!!!
		if (tc==0){
			for (n[0]=0;n[0]<N;n[0]++){
				char buffer[80]={0}; 							
				sprintf(buffer,"animacion_cuerpo%d.txt",n[0]);
			    FILE * fpointer = fopen (buffer,"a"); 
			    fprintf(fpointer,"%lf	%lf	%lf	i\n", R[n[0]][0], R[n[0]][1], R[n[0]][2]);	
		    }
		}
		
		
		//Energia conservada
		double Ectotal=0;
		double Ec[N];
		for (int i=0; i<N;i++){
			Ec[i]=0.5*m[i]*(pow(VR[i][0],2)+pow(VR[i][1],2)+pow(VR[i][2],2));	
		}
		for (int l=0; l<N;l++){
			Ectotal=Ectotal+Ec[l];	
		}
		
		
		double Eptotal=0;
		double Ep[N];
		for (int k=0; k<N;k++){
			double Epb=0;
			n[0]=k;
			for (int i=1;i<N;i++){
				n[i]=(n[0]+i)%N;												
				Epb=Epb-0.5*U(R[n[i]][0],R[n[0]][0],R[n[i]][1],R[n[0]][1],R[n[i]][2],R[n[0]][2],m[n[i]],m[n[0]]);
				///printf("%lf	%d %d\n",Epb,i,k);
				///printf("%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf\n \n",R[n[i]][0],R[n[0]][0],R[n[i]][1],R[n[0]][1],R[n[i]][2],R[n[0]][2],m[n[i]],m[n[0]]);
			}
			Ep[k]=Epb;
		}
		for (int l=0; l<N;l++){
			Eptotal=Eptotal+Ep[l];	
		}
			
		double Etotal=(Eptotal+Ectotal);
	    FILE * fpointerD = fopen ("Energies.dat","a"); 
	    fprintf(fpointerD,"%lf	%lf	%lf	%lf\n",t*dt, Ectotal,Eptotal,Etotal);	
		
		
		
		
		
		
		//PREDATOR%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		//Calculem distància mitja planeta a estrelles
		double r1 = pow(pow(R[3][0]-R[0][0], 2)+pow(R[3][1]-R[0][1], 2)+pow(R[3][2]-R[0][2], 2), 0.5);
		double r2 = pow(pow(R[3][0]-R[1][0], 2)+pow(R[3][1]-R[1][1], 2)+pow(R[3][2]-R[1][2], 2), 0.5);
		double r3 = pow(pow(R[3][0]-R[2][0], 2)+pow(R[3][1]-R[2][1], 2)+pow(R[3][2]-R[2][2], 2), 0.5);
		double r_mean = (m[0]*r1+m[1]*r2+m[2]*r3)/3;

		//A dins de funcions de predator, la r mitja va sense adimensionalitzar, la omega ja la deixa amb dimensió de temperatura i la sigma sense dimensió tot el numerador. Desfem adim.
		double r_mitja = r_mean*384402000; //Radi orbital Lluna = 384402000m
    	double omega=1000000; //1000000
		double Tem = omega/pow(r_mitja+0.01, 0.5);
		
		for(int j=0; j<augment2;j++){
			TEMP[augment2*(t-1)+j] = Tem; //Afegim la T a l'array. Canviarà valor cada I (augments).
		}
	}

	//RK-4 PREDATOR//
	const double pi = 4.0 * atan(1.0); //Pi 
	//double sigma_1=50;
    //double sigma_2=45;
	//double alpha=0.01/(sigma_1*pow(pi, 0.5));
    //double beta=0.09/(sigma_2*pow(pi, 0.5));
    //double a=0.5;
    //double b=0.4;
    //double l=0.8;
	double sigma_1=100000;
    double sigma_2=90000;
	double alpha=0.9/(sigma_1*pow(pi, 0.5));
    double beta=0.9/(sigma_2*pow(pi, 0.5));
    double a=0.5;
    double b=0.4;
    double l=80;

	//double T10 = 65;//Posar T estable preses
	//double T20 = 65;//Posar t estable depredadors
	double T10 = 45;//Posar T estable preses
	double T20 = 30;//Posar t estable depredadors


	for(int k=1; k<I;k++){
		if (x[k]>0 and y[k]>0){

			double k1 = (alpha*sigma_1*pow(pi, 0.5)/b)*exp(pow((TEMP[k]-T10)/(sigma_1), 2));
			double k2 = l/b;
			double k3 = (beta*sigma_2*pow(pi, 0.5)/a)*exp(pow((TEMP[k]-T20)/(sigma_2), 2));
			double h = dt/augment;

			double q1 = h*fi(k1, x[k], y[k]);
			double l1 = h*g(k2, k3, x[k], y[k]);
			double q2 = h*fi(k1, x[k]+q1/2, y[k]+l1/2);
			double l2 = h*g(k2, k3, x[k]+q1/2, y[k]+l1/2);
			double q3 = h*fi(k1, x[k]+q2/2, y[k]+l2/2);
			double l3 = h*g(k2, k3, x[k]+q2/2, y[k]+l2/2);
			double q4 = h*fi(k1, x[k]+q3, y[k]+l3);
			double l4 = h*g(k2, k3, x[k]+q3, y[k]+l3);

			x[k+1] = x[k] + (q1+2*q2+2*q3+q4)/6;
			y[k+1] = y[k] + (q1+2*q2+2*q3+q4)/6;
			temps[k+1] = temps[k]+h; //Avancem una h temporal

			//Aquí desfem canvis de variable. Ens preparem per fer plot de dades en dimensions.
			X[k] = (a*x[k])/b;
			Y[k] = (a*y[k])/b;
			TEMPS[k] = temps[k]/a;
			printf("%lf %lf %lf	%lf	%d	%d\n", TEMPS[k], TEMP[k], X[k], Y[k], k, k/10+1); // POSAR // QUAN ACABEM, NO S'HA DE VEURE RES A TERMINAL
		}	
	}

	
	//PLOT A ARXIU EXTERN DE PREDATOR NO FUNCIONA :(
	FILE * fpointerDPi = fopen ("poblacions.dat","w"); 
	fprintf(fpointerDPi, "%lf %lf %lf\n", TEMPS[0], X[0], Y[0]);
	FILE * fpointerDP = fopen ("poblacions.dat","a"); 
	for (int i=1; i<I; i++){ //imprimimos el t, X, Y de cada iteraci�n en el fichero para tener los datos
		fprintf(fpointerDP, "%lf %lf %lf	%lf	%d	%lf\n", TEMPS[i], TEMP[i], X[i], Y[i], i, i/10.0+1.0); 
	}
	
	return 0;
}
