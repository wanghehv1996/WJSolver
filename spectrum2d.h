#ifndef SPECTRUM_H_
#define SPECTRUM_H_

#include<math.h>  
#include<stdlib.h>

using namespace std;  
#define PI (3.14159265358979323846)
const int __res=32;
double *__omegas;
void __init_omega_arr(){
	__omegas = new double[__res*2-1];
	for(int i=1;i<=32;i++){
		__omegas[2*i-2]=1./( pow(sin(i*PI/2/(__res+1)),2) + pow(sin(i*PI/2/(__res+1)),2) );
		if(i<32)
			__omegas[2*i-1]=1./( pow(sin(i*PI/2/(__res+1)),2) + pow(sin((i+1)*PI/2/(__res+1)),2) );
	}
}

double GetOmegaWithKM(int k, int m, int xres, int yres){
	double omega = 1./(sin(k*PI/2/(xres+1))*sin(k*PI/2/(xres+1)) + sin(m*PI/2/(yres+1))*sin(m*PI/2/(yres+1)));
	return omega;
}


double *__wkx;
double *__wmy;
double *__wkwk;
int __init_w(int xres, int yres){
	if(__wkx)
		delete[] __wkx;
	if(__wmy)
		delete[] __wmy;
	if(__wkwk)
		delete[] __wkwk;
	__wkx = new double[xres*xres];
	__wmy = new double[yres*yres];
	__wkwk = new double[xres*yres];
	
	for(int k=0;k<xres;k++)
	for(int x=0;x<xres;x++)
		__wkx[k*xres+x] =  sin(1.0*(k+1)*(x+1)*PI/(xres+1));
	
	for(int m=0;m<yres;m++)
	for(int y=0;y<yres;y++)
		__wmy[m*yres+y] =  sin(1.0*(m+1)*(y+1)*PI/(yres+1));

	for(int k=0;k<xres;k++)
	for(int m=0;m<yres;m++){
		double wkwk=0;
		for(int i=0;i<xres*yres;i++){
			int x = i%xres;
			int y = i/xres;
			double wkm = __wkx[k*xres+x] * __wmy[m*yres+y];
			wkwk+=wkm*wkm;
		}
		__wkwk[m*xres+k] = wkwk;
		std::cout<<wkwk<<'\t';
	}
}
double *project(const double *error, int xres, int yres){
	double *ck = new double[xres*yres];
	if(!__wkx || !__wmy)
		__init_w(xres,yres);

	for(int k=0;k<xres;k++)
	for(int m=0;m<yres;m++){
		double err_wkm=0;
		double wk_wk=0;
		for(int y=0;y<yres;y++)
		for(int x=0;x<xres;x++){
			// int x = i%xres;
			// int y = i/xres;
			//err*w_k
			double wkm = __wkx[k*xres+x] * __wmy[m*yres+y];
			// err_wkm += error[i] * sin(1.0*(k+1)*(x+1)*PI/(xres+1)) * sin(1.0*(m+1)*(y+1)*PI/(yres+1));
			err_wkm += error[y*xres+x] * wkm;
			//|w_k|
			// wk_wk += sin(1.0*(k+1)*(x+1)*PI/(xres+1)) * sin(1.0*(k+1)*(x+1)*PI/(xres+1)) *
			// 		sin(1.0*(m+1)*(y+1)*PI/(yres+1)) * sin(1.0*(m+1)*(y+1)*PI/(yres+1));
			// wk_wk += wkm*wkm;
		}		// std::cout<<k<<' '<<m<<' '<<wk_wk<<std::endl;
		ck[m *xres + k] = err_wkm/(__wkwk[0]);
		// std::cout<<k<<' '<<m<<' '<<wk_wk<<' '<<ck[m *xres + k]<<std::endl;
	}
	return ck;
}

double GetOmega(const double *error, int xres, int yres){
	double *ck = project(error, xres, yres);
	// double *ck = new double[xres*yres];
	
		
	// for(int k=0;k<xres;k++)
	// for(int m=0;m<yres;m++){
	// 	double err_wk=0;
	// 	double wk_wk=0;
	// 	for(int i=0;i<xres*yres;i++){
	// 		int x = i%xres;
	// 		int y = i/xres;
	// 		//err*w_k
	// 		err_wk += error[i] * sin(1.0*(k+1)*(x+1)*PI/(xres+1)) * sin(1.0*(m+1)*(y+1)*PI/(yres+1));
	// 		//|w_k|
	// 		wk_wk += sin(1.0*(k+1)*(x+1)*PI/(xres+1)) * sin(1.0*(k+1)*(x+1)*PI/(xres+1)) *
	// 				sin(1.0*(m+1)*(y+1)*PI/(yres+1)) * sin(1.0*(m+1)*(y+1)*PI/(yres+1));
	// 	}

	// 	ck[m*xres+k] = err_wk/(wk_wk);
	// }


	double maxck = 0;
	double maxk = 0;
	double maxm = 0;
	for(int k=0;k<xres;k++)
	for(int m=0;m<yres;m++){
		ck[m*xres+k] = ck[m*xres+k]<0?-ck[m*xres+k]:ck[m*xres+k];
		// std::cout<<std::endl<<maxck<<' '<<ck[m*xres+k]<<std::endl;
		// if(maxck<abs(ck[m*xres+k])){
		if(maxck<(ck[m*xres+k])){
			maxck = (ck[m*xres+k]);
			maxk = k+1;
			maxm = m+1;
			
		}
	}
	// std::cout<<ck[0*xres+0]<<' '<<ck[0*xres+1]<<' '<<ck[19*xres+19]<<std::endl;
	
	double omega = 1./( sin(maxk*PI/2/(xres+1))*sin(maxk*PI/2/(xres+1)) + sin(maxm*PI/2/(yres+1))*sin(maxm*PI/2/(yres+1)));
	std::cout<<std::endl<<maxk<<' '<<maxm<<' '<<omega<<std::endl;
	return omega;
}

double GetOmegaNeigh(const double *error, int xres, int yres){
	double omega = GetOmega(error, xres, yres);
	if(!__omegas)
		__init_omega_arr();
	for(int i=0;i<2*xres-1;i++){
		if(omega>=__omegas[i]){
			std::cout<<__omegas[i]<<' '<<omega<<std::endl;
			if(i==0)
				return __omegas[i];
			if(omega-__omegas[i] < __omegas[i-1]-omega)
				return __omegas[i];
			else
				return __omegas[i-1];
		}
	}
}

double GetBestOmega(const double *error, int xres ,int yres){
	double *ckm = project(error, xres, yres);
	// std::cout<<std::endl<<maxk<<std::endl;
	double bestsigckm = 1e20;
	double bestomega = 0;
	double bestk = 0;
	double bestm = 0;
	
	for(int k = 1;k<xres+1;k++)
	for(int m = 1;m<yres+1;m++){
		double omega = GetOmegaWithKM(k,m,xres,yres);
		double sigckm = 0;
		for(int x = 0;x<xres;x++)
		for(int y = 0;y<yres;y++){
			double lambda = 1-omega*( sin(PI*(x+1)/2/(xres+1))*sin(PI*(x+1)/2/(xres+1)) + sin(PI*(y+1)/2/(yres+1))*sin(PI*(y+1)/2/(yres+1)) );
			sigckm+=abs(lambda*ckm[y*xres+x]);
		}
		// std::cout<<k<<" "<<omega<<" makes "<<sigckm<<std::endl;
		if(sigckm<bestsigckm){
			bestomega = omega;
			bestsigckm = sigckm;
			bestk = k;
			bestm = m;
		}
	}
	
	std::cout<<"accurate k="<<bestk<<" m="<<bestm<<" omega="<<bestomega<<std::endl;
	return bestomega;	
}

double GetBestROmega(const double *error, int xres ,int yres){
	double *ckm = project(error, xres, yres);
	// std::cout<<std::endl<<maxk<<std::endl;
	double bestsigckm = 1e20;
	double bestomega = 0;
	double bestk = 0;
	double bestm = 0;
	
	for(int k = 1;k<xres+1;k++)
	for(int m = 1;m<yres+1;m++){
		double omega = GetOmegaWithKM(k,m,xres,yres);
		double sigckm = 0;
		for(int x = 0;x<xres;x++)
		for(int y = 0;y<yres;y++){
			double lambda = 1-omega*( sin(PI*(x+1)/2/(xres+1))*sin(PI*(x+1)/2/(xres+1)) + sin(PI*(y+1)/2/(yres+1))*sin(PI*(y+1)/2/(yres+1)) );
			double lambdaA = ( sin(PI*(x+1)/2/(xres+1))*sin(PI*(x+1)/2/(xres+1)) + sin(PI*(y+1)/2/(yres+1))*sin(PI*(y+1)/2/(yres+1)) );
			sigckm+=abs(lambdaA * lambda*ckm[y*xres+x]);
		}
		// std::cout<<k<<" "<<omega<<" makes "<<sigckm<<std::endl;
		if(sigckm<bestsigckm){
			bestomega = omega;
			bestsigckm = sigckm;
			bestk = k;
			bestm = m;
		}
	}
	
	std::cout<<"accurate k="<<bestk<<" m="<<bestm<<" omega="<<bestomega<<std::endl;
	return bestomega;	
}

double GetBestOmegaNeigh(const double *error, int xres ,int yres){
	double *ckm = project(error, xres, yres);
	double bestsigckm = 1e20;
	double bestomega = 0;
	if(!__omegas)
		__init_omega_arr();
	for(int i=0;i<2*xres-1;i++){
		double sigckm = 0;
		for(int x = 0;x<xres;x++)
		for(int y = 0;y<yres;y++){
			double lambda = 1-__omegas[i]*( sin(PI*(x+1)/2/(xres+1))*sin(PI*(x+1)/2/(xres+1)) + sin(PI*(y+1)/2/(yres+1))*sin(PI*(y+1)/2/(yres+1)) );
			sigckm+=abs(lambda*ckm[y*xres+x]);
		}
		if(sigckm<bestsigckm){
			bestomega = __omegas[i];
			bestsigckm = sigckm;
		}
	}
	std::cout<<"accurate omega="<<bestomega<<std::endl;
	return bestomega;	
}

double GetBestROmegaNeigh(const double *error, int xres ,int yres){
	double *ckm = project(error, xres, yres);
	double bestsigckm = 1e20;
	double bestomega = 0;
	if(!__omegas)
		__init_omega_arr();
	for(int i=0;i<2*xres-1;i++){
		double sigckm = 0;
		for(int x = 0;x<xres;x++)
		for(int y = 0;y<yres;y++){
			double lambda = 1-__omegas[i]*( sin(PI*(x+1)/2/(xres+1))*sin(PI*(x+1)/2/(xres+1)) + sin(PI*(y+1)/2/(yres+1))*sin(PI*(y+1)/2/(yres+1)) );
			double lambdaA = ( sin(PI*(x+1)/2/(xres+1))*sin(PI*(x+1)/2/(xres+1)) + sin(PI*(y+1)/2/(yres+1))*sin(PI*(y+1)/2/(yres+1)) );
			sigckm+=abs(lambdaA * lambda*ckm[y*xres+x]);
		}
		if(sigckm<bestsigckm){
			bestomega = __omegas[i];
			bestsigckm = sigckm;
		}
	}
	std::cout<<"accurate omega="<<bestomega<<std::endl;
	return bestomega;	
}
double CheckWithResidual(double *error,double *residual, int xres, int yres){
	double *rck = project(residual, xres, yres);

	for(int k=0;k<xres;k++)
	for(int m=0;m<yres;m++){
		double lambdaA = 4.0 * ( sin(1.0*(k+1)*PI/2/(xres+1))*sin(1.0*(k+1)*PI/2/(xres+1)) + 
							sin(1.0*(m+1)*PI/2/(yres+1))*sin(1.0*(m+1)*PI/2/(yres+1)) );
		rck[m*xres + k]/=lambdaA;
	}
	double *errors = new double[xres * yres];

	for(int i=0;i<xres*yres;i++){
		int x = i%xres;
		int y = i/xres;
		double err_i=0;
		for(int k=0;k<xres;k++)
		for(int m=0;m<yres;m++){
			err_i+= sin(1.0*(k+1)*(x+1)*PI/(xres+1)) * sin(1.0*(m+1)*(y+1)*PI/(yres+1)) *rck[m*xres+k];
		}
		errors[i] = err_i;
	}

	for(int i=0;i<xres*yres;i++){
		std::cout<<i<<"\t"<<errors[i]<<'\t'<<error[i]<<std::endl;
	}
	return 0;
}

//2Dimension version
double SolveWithResidual(double *result,double *residual, int xres, int yres){
	double *rck = project(residual, xres, yres);

	for(int k=0;k<xres;k++)
	for(int m=0;m<yres;m++){
		double lambdaA = 4.0 * ( sin(1.0*(k+1)*PI/2/(xres+1))*sin(1.0*(k+1)*PI/2/(xres+1)) + 
							sin(1.0*(m+1)*PI/2/(yres+1))*sin(1.0*(m+1)*PI/2/(yres+1)) );
		rck[m*xres + k]/=lambdaA;
	}
	double *errors = new double[xres * yres];

	for(int i=0;i<xres*yres;i++){
		int x = i%xres;
		int y = i/xres;
		double err_i=0;
		for(int k=0;k<xres;k++)
		for(int m=0;m<yres;m++){
			err_i+= sin(1.0*(k+1)*(x+1)*PI/(xres+1)) * sin(1.0*(m+1)*(y+1)*PI/(yres+1)) *rck[m*xres+k];
		}
		errors[i] = err_i;
	}

	for(int i=0;i<xres*yres;i++){
		result[i]+=errors[i];
	}
	delete []rck;
	delete []errors;
	return 0;
}


#endif