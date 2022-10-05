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

int main(){
	
	//Definim les variables
	int N = 4;	//number of bodies
	int T = 1300;
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
    R[0][0]=0.0;  
    R[0][1]=0.0;
    R[0][2]=0.0;
	m[0]=1.0;

	VR[0][0]=0.0;
    VR[0][1]=0.0;
    VR[0][2]=-15.0;

	//RED
    R[1][0]=5.0/1;
    R[1][1]=0.0;
 	R[1][2]=0.0;
	m[1]=0.01;
	
	VR[1][0]=0.0;
    VR[1][1]=5.0;
    VR[1][2]=-15.0;

	//YELLOW
 	R[2][0]=-3.536/1;
    R[2][1]=3.536/1;
    R[2][2]=0.0;
	m[2]=0.01;
	
	VR[2][0]=-3.536;
    VR[2][1]=-3.536;
    VR[2][2]=-15.0;
    
	//BLACK PLANETA
    R[3][0]=-3.536/1;
    R[3][1]=-3.536/1;
 	R[3][2]=0;
	m[3]=0.01;
	
    VR[3][0]=3.536;
    VR[3][1]=-3.536;
    VR[3][2]=-15.0;
	
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
		int salto = 10;
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
		printf("En t=%d Ec=%lf, Ep=%lf, Et=%lf\n\n",t,Ectotal,Eptotal,Etotal);
	}
	return 0;
}
