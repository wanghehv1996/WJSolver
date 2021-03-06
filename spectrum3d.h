#ifndef SPECTRUM_H_
#define SPECTRUM_H_

#include<math.h>  
#include<stdlib.h>

using namespace std;  
#define PI (3.14159265358979323846)
const int __res=32;
double *__omegas;
void __init_omega_arr(){
	__omegas = new double[__res*3-1];
	for(int i=1;i<=32;i++){
		__omegas[3*i-3]=1./( pow(sin(i*PI/2/(__res+1)),2) + pow(sin(i*PI/2/(__res+1)),2) + pow(sin(i*PI/2/(__res+1)),2) );
		if(i<32){
			__omegas[3*i-2]=1./( 2*pow(sin(i*PI/2/(__res+1)),2) +  pow(sin((i+1)*PI/2/(__res+1)),2) );
			__omegas[3*i-1]=1./( pow(sin(i*PI/2/(__res+1)),2) +  2*pow(sin((i+1)*PI/2/(__res+1)),2) );
		}
	}
}

double GetOmegaWithKM(int k, int m, int xres, int yres){
	double omega = 1./(sin(k*PI/2/(xres+1))*sin(k*PI/2/(xres+1)) + sin(m*PI/2/(yres+1))*sin(m*PI/2/(yres+1)));
	return omega;
}

double GetOmegaWithKMN(int k, int m,int n, int xres, int yres, int zres){
	double omega = 1./(sin(k*PI/2/(xres+1))*sin(k*PI/2/(xres+1)) + sin(m*PI/2/(yres+1))*sin(m*PI/2/(yres+1))  + sin(n*PI/2/(zres+1))*sin(n*PI/2/(zres+1)) );
	return omega;
}

double *__wkx;
double *__wmy;
double *__wnz;
double *__wkwk;
double *__lambda2k;
double *__lambda2m;
double *__lambda2n;
int __init_w(int xres, int yres, int zres){
	if(__wkx)		delete[] __wkx;
	if(__wmy)		delete[] __wmy;
	if(__wnz)		delete[] __wnz;
	if(__wkwk)		delete[] __wkwk;

	if(__lambda2k)	delete[] __lambda2k;
	if(__lambda2m)	delete[] __lambda2m;
	if(__lambda2n)	delete[] __lambda2n;

	__wkx = new double[xres*xres];
	__wmy = new double[yres*yres];
	__wnz = new double[zres*zres];

	__wkwk = new double[xres*yres*zres];

	__lambda2k = new double[xres];
	__lambda2m = new double[yres];
	__lambda2n = new double[zres];
	
	for(int k=0;k<xres;k++)
	for(int x=0;x<xres;x++)
		__wkx[k*xres+x] =  sin(1.0*(k+1)*(x+1)*PI/(xres+1));
	
	for(int m=0;m<yres;m++)
	for(int y=0;y<yres;y++)
		__wmy[m*yres+y] =  sin(1.0*(m+1)*(y+1)*PI/(yres+1));

	for(int n=0;n<zres;n++)
	for(int z=0;z<zres;z++)
		__wnz[n*zres+z] =  sin(1.0*(n+1)*(z+1)*PI/(zres+1));

	// for(int k=0;k<xres;k++)
	// for(int m=0;m<yres;m++)
	// for(int n=0;n<zres;n++){
	{
		int k=0,m=0,n=0;
		double wkwk=0;
		// for(int i=0;i<xres*yres*zres;i++){
		for(int i=0;i<xres*yres*zres;i++){
			int x = i%xres;
			int y = i%(xres*yres)/xres;
			int z = i/(xres*yres);
			double wkm = __wkx[k*xres+x] * __wmy[m*yres+y] * __wnz[n*zres+z];
			wkwk+=wkm*wkm;
		}
		__wkwk[m*xres+k] = wkwk;
		// std::cout<<wkwk<<'\t';
	}

	for(int k =0;k<xres;k++)
		__lambda2k[k] = sin((k+1)*PI/2/(xres+1))*sin((k+1)*PI/2/(xres+1));
	for(int m =0;m<yres;m++)
		__lambda2m[m] = sin((m+1)*PI/2/(yres+1))*sin((m+1)*PI/2/(yres+1));
	for(int n =0;n<zres;n++)
		__lambda2n[n] = sin((n+1)*PI/2/(zres+1))*sin((n+1)*PI/2/(zres+1));
}

double *project(const double *error, int xres, int yres, int zres){
	double *ck = new double[xres*yres*zres];
	if(!__wkx || !__wmy || !__wnz)
		__init_w(xres,yres,zres);
	for(int k=0;k<xres;k++)
	for(int m=0;m<yres;m++)
	for(int n=0;n<zres;n++){
		double err_wkmn=0;
		double wk_wk=0;
		for(int z=0;z<zres;z++)
		for(int y=0;y<yres;y++)
		for(int x=0;x<xres;x++){
			double wkmn = __wkx[k*xres+x]* __wmy[m*yres+y]* __wnz[n*zres+z];
			err_wkmn += error[z*xres*yres+y*xres+x] * wkmn;
		}
		// ck[n*yres*xres + m*xres + k] = err_wkmn/(wk_wk);
		ck[n*yres*xres + m*xres + k] = err_wkmn/(__wkwk[0]);

	}
	return ck;
}


//3Dimension
double GetOmega(const double *error, int xres, int yres, int zres){
	double *ck = project(error, xres, yres, zres);
	// // check if the project is right
	// for(int i=0;i<res;i++){
	// 	double err_i = 0;
	// 	for(int k=0;k<res;k++){
	// 		err_i+=sin(1.0*(k+1)*(i+1)*PI/(res+1))*ck[k];
	// 	}
	// 	if(err_i!=error[i])
	// 		std::cout<<"diff"<<error[i]<<"!="<<err_i<< ' ' <<err_i-error[i]<<std::endl;
	// }

	double maxck = 0;
	int maxk = 0;
	int maxm = 0;
	int maxn = 0;
	for(int k=0;k<xres;k++)
	for(int m=0;m<yres;m++)
	for(int n=0;n<zres;n++){
		ck[n*yres*xres + m*xres + k] = ck[n*yres*xres + m*xres + k]<0?-ck[n*yres*xres + m*xres + k]:ck[n*yres*xres + m*xres + k];
		if(maxck<(ck[n*yres*xres + m*xres + k])){
			maxck = (ck[n*yres*xres + m*xres + k]);

			maxk = k;
			maxm = m;
			maxn = n;
		}
	}
	
	double omega = 3./2./(__lambda2k[maxk]+__lambda2m[maxm]+__lambda2n[maxn]);

	std::cout<<std::endl<<maxk+1<<' '<<maxm+1<<' '<<maxn+1<<' '<<omega<<std::endl;
	return omega;
}


double GetOmegaNeigh(const double *error, int xres, int yres, int zres){
	double omega = GetOmega(error, xres, yres, zres);
	if(!__omegas)
		__init_omega_arr();
	for(int i=0;i<3*xres-1;i++){
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

double GetBestOmega(const double *error, int xres ,int yres, int zres){
	double *ckm = project(error, xres, yres, zres);
	// std::cout<<std::endl<<maxk<<std::endl;
	double bestsigckm = 1e20;
	double bestomega = 0;
	double bestk = 0;
	double bestm = 0;
	double bestn = 0;
	
	for(int k = 1;k<xres+1;k++)
	for(int m = 1;m<yres+1;m++)
	for(int n = 1;n<zres+1;n++){
		double omega = GetOmegaWithKMN(k,m,n,xres, yres, zres);
		double sigckm = 0;
		for(int x = 0;x<xres;x++)
		for(int y = 0;y<yres;y++)
		for(int z = 0;z<zres;z++){
			// double lambda = 1-2./3.*omega*( 
			// 	sin(PI*(x+1)/2/(xres+1))*sin(PI*(x+1)/2/(xres+1)) + 
			// 	sin(PI*(y+1)/2/(yres+1))*sin(PI*(y+1)/2/(yres+1)) +
			// 	sin(PI*(z+1)/2/(zres+1))*sin(PI*(z+1)/2/(zres+1)) );
			double lambda = 1-2./3.*omega*(__lambda2k[x]+__lambda2m[y]+__lambda2n[z]);
			sigckm+=abs(lambda*ckm[z*yres*xres + y*xres + x]);
		}
		// std::cout<<k<<" "<<omega<<" makes "<<sigckm<<std::endl;
		if(sigckm<bestsigckm){
			bestomega = omega;
			bestsigckm = sigckm;
			bestk = k;
			bestm = m;
			bestn = n;
		}
	}
	
	std::cout<<"accurate k="<<bestk<<" m="<<bestm<<" omega="<<bestomega<<std::endl;
	return bestomega;	
}

// double GetBestROmega(const double *error, int xres ,int yres){
// 	double *ckm = project(error, xres, yres);
// 	// std::cout<<std::endl<<maxk<<std::endl;
// 	double bestsigckm = 1e20;
// 	double bestomega = 0;
// 	double bestk = 0;
// 	double bestm = 0;
	
// 	for(int k = 1;k<xres+1;k++)
// 	for(int m = 1;m<yres+1;m++){
// 		double omega = GetOmegaWithKM(k,m,xres,yres);
// 		double sigckm = 0;
// 		for(int x = 0;x<xres;x++)
// 		for(int y = 0;y<yres;y++){
// 			double lambda = 1-omega*( sin(PI*(x+1)/2/(xres+1))*sin(PI*(x+1)/2/(xres+1)) + sin(PI*(y+1)/2/(yres+1))*sin(PI*(y+1)/2/(yres+1)) );
// 			double lambdaA = ( sin(PI*(x+1)/2/(xres+1))*sin(PI*(x+1)/2/(xres+1)) + sin(PI*(y+1)/2/(yres+1))*sin(PI*(y+1)/2/(yres+1)) );
// 			sigckm+=abs(lambdaA * lambda*ckm[y*xres+x]);
// 		}
// 		// std::cout<<k<<" "<<omega<<" makes "<<sigckm<<std::endl;
// 		if(sigckm<bestsigckm){
// 			bestomega = omega;
// 			bestsigckm = sigckm;
// 			bestk = k;
// 			bestm = m;
// 		}
// 	}
	
// 	std::cout<<"accurate k="<<bestk<<" m="<<bestm<<" omega="<<bestomega<<std::endl;
// 	return bestomega;	
// }

double GetBestOmegaNeigh(const double *error, int xres ,int yres, int zres){
	double *ckm = project(error, xres, yres, zres);
	double bestsigckm = 1e20;
	double bestomega = 0;
	if(!__omegas)
		__init_omega_arr();
	for(int i=0;i<3*xres-1;i++){
		double sigckm = 0;
		for(int x = 0;x<xres;x++)
		for(int y = 0;y<yres;y++)
		for(int z = 0;z<zres;z++){
			double lambda = 1-2./3.*__omegas[i]*(__lambda2k[x]+__lambda2m[y]+__lambda2n[z]);
			// double lambda = 1-__omegas[i]*( 
			// 	sin(PI*(x+1)/2/(xres+1))*sin(PI*(x+1)/2/(xres+1)) + 
			// 	sin(PI*(y+1)/2/(yres+1))*sin(PI*(y+1)/2/(yres+1)) +
			// 	sin(PI*(z+1)/2/(zres+1))*sin(PI*(z+1)/2/(zres+1)));
			sigckm+=abs(lambda*ckm[z*xres*yres + y*xres + x]);
		}
		if(sigckm<bestsigckm){
			bestomega = __omegas[i];
			bestsigckm = sigckm;
		}
	}
	std::cout<<"accurate omega="<<bestomega<<std::endl;
	return bestomega;	
}

double GetBestROmegaNeigh(const double *error, int xres ,int yres, int zres){
	double *ckm = project(error, xres, yres, zres);
	double bestsigckm = 1e20;
	double bestomega = 0;
	if(!__omegas)
		__init_omega_arr();
	for(int i=0;i<3*xres-1;i++){
		double sigckm = 0;
		for(int x = 0;x<xres;x++)
		for(int y = 0;y<yres;y++)
		for(int z = 0;z<zres;z++){
			// double lambda = 1-__omegas[i]*( sin(PI*(x+1)/2/(xres+1))*sin(PI*(x+1)/2/(xres+1)) + sin(PI*(y+1)/2/(yres+1))*sin(PI*(y+1)/2/(yres+1)) );
			double lambda = 1-2./3.*__omegas[i]*(__lambda2k[x]+__lambda2m[y]+__lambda2n[z]);
			// double lambdaA = ( sin(PI*(x+1)/2/(xres+1))*sin(PI*(x+1)/2/(xres+1)) + sin(PI*(y+1)/2/(yres+1))*sin(PI*(y+1)/2/(yres+1)) );
			double lambdaA = (__lambda2k[x]+__lambda2m[y]+__lambda2n[z]);
			sigckm+=abs(lambdaA * lambda*ckm[z*xres*yres + y*xres + x]);
		}
		if(sigckm<bestsigckm){
			bestomega = __omegas[i];
			bestsigckm = sigckm;
		}
	}
	std::cout<<"accurate omega="<<bestomega<<std::endl;
	return bestomega;	
}

double *solveDirect(const double *residual, int xres, int yres, int zres){
	double *errors = new double[xres*yres*zres];
	double *lck = new double[xres*yres*zres];
	if(!__wkx || !__wmy || !__wnz)
		__init_w(xres,yres,zres);
	for(int k=0;k<xres;k++)
	for(int m=0;m<yres;m++)
	for(int n=0;n<zres;n++){
		double err_wkmn=0;
		double wk_wk=0;
		for(int z=0;z<zres;z++)
		for(int y=0;y<yres;y++)
		for(int x=0;x<xres;x++){
			double wkmn = __wkx[k*xres+x]* __wmy[m*yres+y]* __wnz[n*zres+z];
			err_wkmn += residual[z*xres*yres+y*xres+x] * wkmn;
		}
		// ck[n*yres*xres + m*xres + k] = err_wkmn/(wk_wk);
		lck[n*yres*xres + m*xres + k] = err_wkmn/(__wkwk[0]);
		lck[n*yres*xres + m*xres + k] = 
			lck[n*yres*xres + m*xres + k]/4./(__lambda2k[k]+__lambda2m[m]+__lambda2n[n]);

	}

	for(int z=0;z<zres;z++)
	for(int y=0;y<yres;y++)
	for(int x=0;x<xres;x++){
		double err_i=0;
		for(int k=0;k<xres;k++)
		for(int m=0;m<yres;m++)
		for(int n=0;n<zres;n++){
			// err_i+= sin(1.0*(k+1)*(x+1)*PI/(xres+1)) * sin(1.0*(m+1)*(y+1)*PI/(yres+1)) *rck[m*xres+k];
			err_i+= __wkx[k*xres+x]* __wmy[m*yres+y]* __wnz[n*zres+z]*lck[n*yres*xres + m*xres + k];
		}
		errors[z*yres*xres + y*xres + x] = err_i;
	}

	return errors;
}
#endif