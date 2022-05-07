
#include"Kernel.h"
#include"radix.h"
#define a 5.0
#define x0 5.0

using namespace std;

//Obtención de las densidades de cada una de las partículas en el método SPH, donde será necesario la posción para evaluar el kernel, el parámetro de suavizado para definir el kernel y un arreglo donde se guardaran los valores de la densidad, el número de iperaciones es de orden N^2.
void Densidad0(int N, double m[], double R[], double h, double D[]){
  
	for(int i=0; i<N; i++){
		D[i]=0.0;
		for(int j=0; j<N; j++){
			D[i]+=m[j]*ker(h, R[i], R[j]);
		} 
	}
}
	//We define the four steps like a function, so

void Densidad0Eff(int N,double xmax, double xmin, int xf, double m[], double R[], double hf, double D[]){
	int Rd[N];
	double h_hash=2*hf;//Only for this case
	Discretx(N, h_hash, xmin, R, Rd);//With this function, we obtain the position in discret coordinates
	//define the index max for the
	int H_x=int((xmax-xmin)/h_hash);//We define the max value for discret position. This is only for the hf=h_hash
	//Define the key value
	int key[N], keyS[N], idx[N]; //Define the key array and the key sort array and the array with the permuted index
	for(int i=0; i<N;i++){
		key[i]=Rd[i]; //Here define the value of key for each particle
		idx[i]=i;
	}
	//We need 4 arrays:
	// 1.- Index original for particles idx[]
	// 2.- New Index for sort keyS: KeyS[i]: i-> is the New Index
	// 3.-
	RadixSPH(N,key,keyS,idx); //Ordered the values of the Key and give value for the permuted index, with the noral index is the original index, so for
	
	// the array idx[i] provide the number of the particle (original): idx[i] and i is the permuted index where the key is ordered.
	int nclass=H_x;//The number of classes in this case 1D is the same that the H_x
	// Define the classes
	//POdemos hacer mass eficiente partiendo de un loop sobre todas las particulas e ir revisando cada cuando cambia un escalon
	int idxmin[nclass]={}, idxmax[nclass]={}; // the number of the class is in the index of this arrays
	bool act[nclass]={false}; //This array only say us if the cell is ocuped or not
	int numclass;
	numclass=keyS[0];
	act[numclass]=true;
	idxmin[numclass]=0;
	for(int i=1; i<N; i++){
		if(keyS[i]!=numclass){
			idxmax[numclass]=i-1;
			numclass=keyS[i];
			act[numclass]=true;
			idxmin[numclass]=i;
		}
	}
	idxmax[numclass]=N-1;
	for(int i=0; i<N; i++){
		printf("%d key sort: %d \n",i, keyS[i]);
	}
	
	for(int i=0; i<N; i++){
		D[idx[i]]=0.0;
		for(int nclass2=keyS[i]-xf;nclass2<=keyS[i]+xf;nclass2++){
			printf("El indice es %d, el KeyS: %d, entonces la nclass2 es %d \n",i, keyS[i], nclass2);
			if(act[nclass2]=true){
				for(int j=idxmin[nclass2]; j<=idxmax[nclass2];j++){
					D[idx[i]]+=m[idx[j]]*ker(hf, R[idx[i]], R[idx[j]]);
				}
			}
		}
	}
	
/*
	for(int i=0; i<N; i++){
			D[idx[i]]=0.0; // i is the permuted index, for each particle we calculate the density
			//If we know the permuted index, we can search the cell by the key and search
			if(keyS[i]==0 || keyS[i]==H_x){
				if(keyS[i]==0){
					for(int j=idxmin[keyS[i]]; j<idxmax[keyS[i]+1]; j++){
						D[idx[i]]+=m[idx[i]]*ker(hf, R[idx[i]], R[idx[j]]);
					}
				}
				if(keyS[i]==H_x){
					for(int j=idxmin[keyS[i]-1]; j<idxmax[keyS[i]]; j++){
						D[idx[i]]+=m[idx[i]]*ker(hf, R[idx[i]], R[idx[j]]);
					}
				}
			}else{
				for(int j=idxmin[keyS[i]-1]; j<idxmax[keyS[i]+1]; j++){	
					D[idx[i]]+=m[idx[i]]*ker(hf, R[idx[i]], R[idx[j]]);
				}
			}
	}
*/
}
//La derivada de la densidad se obtiene a través de la masa de cada una de las partículas como factor que multiplica la derivada del kernel y una constante q, aunque se puede omitir el factor q, como segunda discretización.
void Densidad1(int N, double m[], double R[], double h[], double D[], double Dx[]){
	double q=0.0;	
	//first discretization -----> Reproduce las derivadas del pefil analitico
	//--------->001, 002, 003, 006, 008, 009

	for(int i=0; i<N; i++){
		Dx[i]=0.0;
		for(int j=0; j<N; j++){
			Dx[i]+=m[j]*dker(h[i], R[i], R[j]);
			}
		}

	//--------->001, 002, 003, 006, 008, 009
   //second discretization
   //---------> 001, 002, 004, 005---->007
/*
  for(int i=0; i<N; i++){
    Dx[i]=0.0;
    for(int j=0; j<N; j++){
      q=(D[j]-D[i])/D[j];
      Dx[i]+=m[j]*q*dker(h[i], R[i], R[j]);
      }
  }
*/
	//------------>001, 002, 004, 005----->007
/*
  //thirth discretization
    for(int i=0; i<N; i++){
		Dx[i]=0.0;
		for(int j=0; j<N; j++){
			q=(D[i]-D[j])/D[i];
			Dx[i]+=m[j]*q*dker(h[i], R[i], R[j]);
		}
	}
*/
/*	//fourth discretization
	for(int i=0; i<N; i++){
		Dx[i]=0.0;
		for(int j=0; j<N; j++){
				q=D[i]*(1.0/D[i]+1.0/D[j]):
				Dx[i]+=m[j]*q*dker(h[i], R[i], R[j]);
			}
		}
		*/
/*
		//////////for desplacement
	
		for(int i=0; i<N; i++){
			Dx[i]=0.0;
			for(int j=0; j<N; j++){
				Dx[i]+=m[j]*dker(h[i], R[i]-1.0, R[j]);
			}
		  }
	*/	
}
//Al igual que la primera derivada, se pueden obtener con los mismos factores y la segunda derivada del Kernel.
void Densidad2(int N , double m[], double R[], double h[], double D[], double Dx[], double Dxx[]){
	double q=0.0;	
	//first discretization--------> reproduce la derivida de la dendidad analitica

  //-----------------_>003, 006, 009
  for(int i=0; i<N; i++){
    Dxx[i]=0.0;
    for(int j=0; j<N; j++){
		Dxx[i]+=m[j]*ddker(h[i], R[i], R[j]);
    }
  }
	//---------------->003, 006, 009

/*   //second discretization
  for(int i=0; i<N; i++){
    Dxx[i]=0.0;
    for(int j=0; j<N; j++){
      q=(Dx[j]-Dx[i])/Dx[j];
      Dxx[i]+=m[j]*q*dker(h[i], R[i], R[j]);
      }
  }
  */
   //thirth discretization
/*
   //----------> 001, 002, 004, 005------->007, 008
  for(int i=0; i<N; i++){
    Dxx[i]=0.0;
    for(int j=0; j<N; j++){
      q=(D[j]-D[i])/D[j];
      Dxx[i]+=m[j]*q*ddker(h[i], R[i], R[j]);
      }
  }
	//----------> 001, 002, 004, 005-------->007, 008
*/
/*  //fourth discretization
    for(int i=0; i<N; i++){
		Dxx[i]=0.0;
		for(int j=0; j<N; j++){
			q=(D[i]-D[j])/D[i];
			Dxx[i]+=m[j]*q*ddker(h[i], R[i], R[j]);
		}
	}
	
/*	//fiveth discretization
	for(int i=0; i<N; i++){
		Dx[i]=0.0;
		for(int j=0; j<N; j++){
				q=D[i]*(1.0/D[i]+1.0/D[j]):
				Dx[i]+=m[j]*q*ddker(h[i], R[i], R[j]);
			}
		}
		*/
/*
		////////for desplacement
		for(int i=0; i<N; i++){
			Dxx[i]=0.0;
			for(int j=0; j<N; j++){
				Dxx[i]+=m[j]*ddker(h[i], R[i]-1.0, R[j]);
			}
		  }
	*/	
}	

void Densidad3(int N , double m[], double R[], double h[], double D[], double Dx[], double Dxx[], double Dxxx[]){
	
	for(int i=0; i<N; i++){
		Dxxx[i]=0.0;
		for(int j=0; j<N; j++){
			Dxxx[i]+=m[j]*dddker(h[i], R[i], R[j]);
			}
		}
	 
/*
			////////for desplacement
		for(int i=0; i<N; i++){
			Dxxx[i]=0.0;
			for(int j=0; j<N; j++){
				Dxxx[i]+=m[j]*dddker(h[i], R[i]-1.0, R[j]);
			}
		  }
		  */ 
}

double DA(int Ptype, double Ri){
        double z, alt;
	if(Ptype==1){
		alt=sin(PI*Ri/a);
		z=(2.0/a)*alt*alt;
	}
	if(Ptype==2){
		z=(1.0/sqrt(PI))*exp(-Ri*Ri);
	}
	if(Ptype==3){
		z=1.0;
	}
	if(Ptype==4){
//		z=sin(Ri)*sin(Ri)/(Ri+0.1);
		alt=1.0/cosh(Ri-x0);
		z=0.25*(alt*alt);
	}
	if(Ptype==5){
		double g=1.0;
		double r0=sqrt(PI*g/4.0); 
		double R=PI*Ri/r0; 
		double rho0=PI/(4.0*r0*r0*r0);
		z=rho0*sin(R)/R;
	}
	if(Ptype==6){
//		z=(0.25)*(1.0/(cosh(Ri-5.0)*cosh(Ri-5.0)))+(0.25)*(1.0/(cosh(Ri+5.0)*cosh(Ri+5.0)));
	z=(0.25)*(1.0/(cosh(Ri+5.0)*cosh(Ri+5.0)))+(0.25)*(1.0/(cosh(Ri-5.0)*cosh(Ri-5.0)));
		}
	return z;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////Densidad 2D


