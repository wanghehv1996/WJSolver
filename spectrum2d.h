#ifndef SPECTRUM_H_
#define SPECTRUM_H_

#include<math.h>  
#include<stdlib.h>

using namespace std;  
// #define PI (3.14159265358979323846)
#define PI (3.1415926535897932384626433832795028841971693993751)
const int __res=64;
double *__omegas;
void __init_omega_arr(){
	__omegas = new double[__res*2-1];
	for(int i=1;i<=__res;i++){
		__omegas[2*i-2]=1./( pow(sin(i*PI/2/(__res+1)),2) + pow(sin(i*PI/2/(__res+1)),2) );
		if(i<__res)
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
double *__lambda2k;
double *__lambda2m;
int __init_w(int xres, int yres){
	if(__wkx)		delete[] __wkx;
	if(__wmy)		delete[] __wmy;
	if(__wkwk)		delete[] __wkwk;

	if(__lambda2k)	delete[] __lambda2k;
	if(__lambda2m)	delete[] __lambda2m;

	__wkx = new double[xres*xres];
	__wmy = new double[yres*yres];

	__wkwk = new double[xres*yres];

	__lambda2k = new double[xres];
	__lambda2m = new double[yres];
	
	for(int k=0;k<xres;k++)
	for(int x=0;x<xres;x++)
		__wkx[k*xres+x] =  sin(1.0*(k+1)*(x+1)*PI/(xres+1));
	
	for(int m=0;m<yres;m++)
	for(int y=0;y<yres;y++)
		__wmy[m*yres+y] =  sin(1.0*(m+1)*(y+1)*PI/(yres+1));

	// for(int k=0;k<xres;k++)
	// for(int m=0;m<yres;m++){
	{
		int k=0,m=0;
		double wkwk=0;
		for(int i=0;i<xres*yres;i++){
			int x = i%xres;
			int y = i/xres;
			// wkwk+=sin(1.0*(k+1)*(x+1)*PI/(xres+1))*sin(1.0*(k+1)*(x+1)*PI/(xres+1))*
			// 	sin(1.0*(m+1)*(y+1)*PI/(yres+1))*sin(1.0*(m+1)*(y+1)*PI/(yres+1));
			double wkm = __wkx[k*xres+x] * __wmy[m*yres+y];
			wkwk+=wkm*wkm; 
		}
		__wkwk[m*xres+k] = wkwk;
	}
	for(int k =0;k<xres;k++)
		__lambda2k[k] = sin((k+1)*PI/2/(xres+1))*sin((k+1)*PI/2/(xres+1));
	for(int m =0;m<yres;m++)
		__lambda2m[m] = sin((m+1)*PI/2/(yres+1))*sin((m+1)*PI/2/(yres+1));
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
			double wkm = __wkx[k*xres+x] * __wmy[m*yres+y];
			err_wkm += error[y*xres+x] * wkm;
		}
		ck[m *xres + k] = err_wkm/(__wkwk[0]);
		
	}
	return ck;
}

double GetOmegaResidual(const double *residual, int xres, int yres){
	double *ck = project(residual, xres, yres);

	double maxck = 0;
	int maxk = 0;
	int maxm = 0;

	for(int k=0;k<xres;k++)
	for(int m=0;m<yres;m++){
		
		double lambdaA = 4.0*(__lambda2k[k]+__lambda2k[m]);
		ck[m*xres + k]/=lambdaA;
		ck[m*xres+k] = ck[m*xres+k]<0?-ck[m*xres+k]:ck[m*xres+k];

		if(maxck<(ck[m*xres+k])){
			maxck = (ck[m*xres+k]);
			maxk = k;
			maxm = m;
			
		}
	}
	double omega = 1./(__lambda2k[maxk]+__lambda2m[maxm]);
	std::cout<<std::endl<<maxk+1<<' '<<maxm+1<<' '<<omega<<std::endl;
	return omega;
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
	int maxk = 0;
	int maxm = 0;
	for(int k=0;k<xres;k++)
	for(int m=0;m<yres;m++){
		ck[m*xres+k] = ck[m*xres+k]<0?-ck[m*xres+k]:ck[m*xres+k];
		// std::cout<<std::endl<<maxck<<' '<<ck[m*xres+k]<<std::endl;
		// if(maxck<abs(ck[m*xres+k])){
		if(maxck<(ck[m*xres+k])){
			maxck = (ck[m*xres+k]);
			maxk = k;
			maxm = m;
			
		}
	}
	// double omega = 1./( sin(maxk*PI/2/(xres+1))*sin(maxk*PI/2/(xres+1)) + sin(maxm*PI/2/(yres+1))*sin(maxm*PI/2/(yres+1)));
	double omega = 1./(__lambda2k[maxk]+__lambda2m[maxm]);
	std::cout<<std::endl<<maxk+1<<' '<<maxm+1<<' '<<maxck<<' '<<omega<<std::endl;
	return omega;
}

double *__kRough;
double *__sin2Rough;
double *__omegasRough;
double __calomega(double lamb, double k){
	int res=64;
	return (1.-lamb)/(2*sin(k*PI/2/(res+1))*sin(k*PI/2/(res+1)));
}
double __calk(double lamb, double omega){
	int res=64;
	double sinv = sqrt((1.-lamb)/omega);
	if(sinv>1)
		return res;
	std::cout<<"asing"<<sinv<<' '<<asin(sinv)<<' '<<asin(sinv)*2.*(res+1)/PI<<std::endl;
	return asin(sinv)*2.*(res+1)/PI;
}
//k < k0 use omega[0]
//k in [k0,k1] use omega[1]
void __init_omega_rough_arr(){
	int res = 64;
	if(__omegasRough)
		delete[] __omegasRough;
	if(__kRough)
		delete[] __kRough;
	if(__sin2Rough)
		delete[] __sin2Rough;
	
	__omegasRough = new double[10];
	__kRough = new double[10];
	__sin2Rough = new double[10];
	double k = 1;
	double lamb = 0;
	int i = 0;
	while(k<res){
		double omega = __calomega(lamb,k);
		
		lamb = 0.1;
		double nk = __calk(-lamb, omega);

		if(nk==res)//omega is out side
		{
			__omegasRough[i]=__calomega(0,(nk+k)/2);
			__kRough[i]=nk;
			__sin2Rough[i]=2*sin((nk)*PI/2/(res+1))*sin((nk)*PI/2/(res+1));
			k=nk;
		}
		else{
			__omegasRough[i]=omega;
			__kRough[i]=nk;
			__sin2Rough[i]=2*sin((nk)*PI/2/(res+1))*sin((nk)*PI/2/(res+1));
			k=nk;
		}
		std::cout<<omega<<' '<<k<<std::endl;
		
		i++;
	}
}

double __find_omega_rough_arr(double sin2){
	int i = 0;
	for(i;i<10;i++){
		if(__sin2Rough[i]>=sin2)
			break;
	}
	std::cout<<"find __find_omega_rough_arr"<<__kRough[i]<<' '<<__sin2Rough[i]<<' '<<__omegasRough[i]<<std::endl;
	return __omegasRough[i];
}

double GetOmegaRough(const double *error, int xres, int yres){
	double *ck = project(error, xres, yres);

	double maxck = 0;
	int maxk = 0;
	int maxm = 0;
	for(int k=0;k<xres;k++)
	for(int m=0;m<yres;m++){
		ck[m*xres+k] = ck[m*xres+k]<0?-ck[m*xres+k]:ck[m*xres+k];
		// std::cout<<std::endl<<maxck<<' '<<ck[m*xres+k]<<std::endl;
		// if(maxck<abs(ck[m*xres+k])){
		if(maxck<(ck[m*xres+k])){
			maxck = (ck[m*xres+k]);
			maxk = k;
			maxm = m;
			
		}
	}
	// double omega = 1./( sin(maxk*PI/2/(xres+1))*sin(maxk*PI/2/(xres+1)) + sin(maxm*PI/2/(yres+1))*sin(maxm*PI/2/(yres+1)));
	double sin2 = __lambda2k[maxk]+__lambda2m[maxm];
	std::cout<<"sin2="<<sin2;
	if(!__kRough || !__sin2Rough || !__omegasRough)
		__init_omega_rough_arr();
	double omega = __find_omega_rough_arr(sin2);
	// double omega = 1./(__lambda2k[maxk]+__lambda2m[maxm]);
	
	std::cout<<std::endl<<maxk+1<<' '<<maxm+1<<' '<<omega<<std::endl;
	return omega;
}

//simply from large to small
//doesnt work because of truncation error
//same in 2d
//choose a serial may make things better
//but not good as choose the largest one
//hi->lo,hi,hi,...
double i_os=__res*2-2;
int clear_lo = 0;
double GetOmegaSerial(const double *error, int xres, int yres){
	if(!__omegas)
		__init_omega_arr();
	double omega = __omegas[int(i_os)];
	if(i_os==__res)
		if(!clear_lo){
			i_os--;
			clear_lo=1;
		}else
			i_os=__res*2-2;
	else
		i_os-=1;
	if(i_os<0)
		i_os=__res*2-2;
	return omega;
}
int i_ot=0;
int don=0;
double GetOmegaTypicalChoose(const double *residual, int xres, int yres){
	if(!__omegas)
		__init_omega_arr();
	double *rck = project(residual, xres, yres);
	int maxk,maxm;
	double maxrck=0;
	for(int k=0;k<xres;k++)
	for(int m=0;m<yres;m++){
		double lambdaA = 4.0 * ( __lambda2k[k]+__lambda2m[m] );
		rck[m*xres + k]/=lambdaA;
		rck[m*xres + k] = rck[m*xres + k]<0? -rck[m*xres + k]:rck[m*xres + k];
		if(maxrck<rck[m*xres + k]){
			maxk = k;
			maxm = m ;
			maxrck = rck[m*xres + k];
		}
	}
	std::cout<<maxk<<' '<<maxm<<' '<<maxrck<<std::endl;

	double omega = 0;
	if(i_ot%2==0)
		if(!don)
			omega = __omegas[i_ot/4];
		else omega = 0.6666;
	else if(i_ot%4==1)
		omega = __omegas[2*__res-2];
	else
		omega = __omegas[64-(i_ot/2+i_ot%4+2)%(32)];
	i_ot++;
	if(i_ot>=31*4){
		don = 1;
		i_ot=0;
	}
	return omega;
}

/****************
 * Choose omega by project each line to find maxck
 ****************/
double GetOmegaLineChoose(const double *residual, int xres, int yres){
	double *ck = new double[xres];
	double maxck=0;
	int maxckk = 0,maxcky = 0;

	double *cm = new double[yres];
	double maxcm=0;
	int maxcmx = 0,maxcmm = 0;
	
	for(int y=0;y<yres;y++){
		for(int k=0;k<xres;k++){
			double res_wk=0;
			double wk_wk=1;
			for(int x=0;x<xres;x++){
				res_wk += residual[y*xres+x]*__wkx[k*xres+x];
			}
			double lambdaA=4.0*__lambda2k[k];
			ck[k] = res_wk/wk_wk/lambdaA;
			ck[k] = ck[k]<0? -ck[k]:ck[k];
			if(ck[k]>maxck){
				maxck = ck[k];
				maxckk = k;
				maxcky = y;
			}
		}
	}
	for(int x=0;x<xres;x++)
		for(int m=0;m<yres;m++){
			double res_wm=0;
			double wm_wm=1;
			for(int y=0;y<yres;y++){
				res_wm += residual[y*xres+x]*__wmy[m*yres+y];
			}
			double lambdaA=4.0*__lambda2m[m];
			cm[m] = res_wm/wm_wm/lambdaA;
			cm[m] = cm[m]<0? -cm[m]:cm[m];
			if(cm[m]>maxcm){
				maxcm = cm[m];
				maxcmm = m;
				maxcmx = x;
			}
		}
	std::cout<<"maxck"<<maxck<<" ckk="<<maxckk<<" cky="<<maxcky<<std::endl;
	std::cout<<"maxcm"<<maxcm<<" cmm="<<maxcmm<<" cmx="<<maxcmx<<std::endl;
	return 0;
}

/****************
 * Choose omega by find maxck in partial section
 ****************/
double GetOmegaPartChoose(const double *residual, int xres, int yres){
	int blocksize = 32;
	double maxck = 0;
	int maxk = 0, maxm = 0;

	for(int Y=0;Y<yres;Y+=blocksize)
	for(int X=0;X<xres;X+=blocksize){
		double *partR = new double[blocksize*blocksize];
		for (int y = 0;y<blocksize;y++)
		for (int x = 0;x<blocksize;x++){
			partR[y*blocksize+x] = residual[(Y+y)*xres+X+x];
		}

		double *ck = project(partR,blocksize,blocksize);
		// double *ck;
		for(int k=0;k<blocksize;k++)
		for(int m=0;m<blocksize;m++){
			
			double lambdaA = 4.0*(__lambda2k[k]+__lambda2k[m]);
			ck[m*blocksize + k]/=lambdaA;
			ck[m*blocksize+k] = ck[m*blocksize+k]<0?-ck[m*blocksize+k]:ck[m*blocksize+k];

			if(maxck<(ck[m*blocksize+k])){
				maxck = (ck[m*blocksize+k]);
				maxk = k;
				maxm = m;
				
			}
		}
		if(ck!=NULL)
			delete []ck;
		// std::cout<<"get ck"<<std::endl;
	}
	std::cout<<maxk<<" "<<maxm<<std::endl;
	// return 1;
	return GetOmegaWithKM((2*maxk-1),2*maxm-1,xres,yres);
}

double GetOmegaNeigh(const double *error, int xres, int yres){
	double omega = GetOmega(error, xres, yres);
	if(!__omegas)
		__init_omega_arr();
	for(int i=0;i<2*xres-1;i++){
		// std::cout<<omega<<__omegas[i];
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
	int bestk = 0;
	int bestm = 0;
	
	for(int k = 1;k<xres+1;k++)
	for(int m = 1;m<yres+1;m++){
		double omega = GetOmegaWithKM(k,m,xres,yres);
		double sigckm = 0;
		for(int x = 0;x<xres;x++)
		for(int y = 0;y<yres;y++){
			// double lambda = 1-omega*( sin(PI*(x+1)/2/(xres+1))*sin(PI*(x+1)/2/(xres+1)) + sin(PI*(y+1)/2/(yres+1))*sin(PI*(y+1)/2/(yres+1)) );
			double lambda = 1-omega*(__lambda2k[x]+__lambda2m[y]);
			if(lambda*ckm[y*xres+x]>=0)
				sigckm+=(lambda*ckm[y*xres+x]);
			else
				sigckm-= (lambda*ckm[y*xres+x]);
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
	int bestk = 0;
	int bestm = 0;
	
	for(int k = 1;k<xres+1;k++)
	for(int m = 1;m<yres+1;m++){
		double omega = GetOmegaWithKM(k,m,xres,yres);
		double sigckm = 0;
		for(int x = 0;x<xres;x++)
		for(int y = 0;y<yres;y++){
			// double lambda = 1-omega*( sin(PI*(x+1)/2/(xres+1))*sin(PI*(x+1)/2/(xres+1)) + sin(PI*(y+1)/2/(yres+1))*sin(PI*(y+1)/2/(yres+1)) );
			double lambda = 1-omega*(__lambda2k[x]+__lambda2m[y]);
			// double lambdaA = ( sin(PI*(x+1)/2/(xres+1))*sin(PI*(x+1)/2/(xres+1)) + sin(PI*(y+1)/2/(yres+1))*sin(PI*(y+1)/2/(yres+1)) );
			double lambdaA = (__lambda2k[x]+__lambda2m[y]);
			if(lambda*ckm[y*xres+x]>=0)
				sigckm+=(lambda*ckm[y*xres+x]);
			else
				sigckm-= (lambda*ckm[y*xres+x]);		}
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
			// double lambda = 1-__omegas[i]*( sin(PI*(x+1)/2/(xres+1))*sin(PI*(x+1)/2/(xres+1)) + sin(PI*(y+1)/2/(yres+1))*sin(PI*(y+1)/2/(yres+1)) );
			double lambda = 1-__omegas[i]*(__lambda2k[x]+__lambda2m[y]);
			if(lambda*ckm[y*xres+x]>=0)
				sigckm+=(lambda*ckm[y*xres+x]);
			else
				sigckm-= (lambda*ckm[y*xres+x]);
		}
		if(sigckm<bestsigckm){
			bestomega = __omegas[i];
			bestsigckm = sigckm;
		}
	}
	std::cout<<"accurate omega="<<bestomega<<"sigckm"<<bestsigckm<<std::endl;
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
			// double lambda = 1-__omegas[i]*( sin(PI*(x+1)/2/(xres+1))*sin(PI*(x+1)/2/(xres+1)) + sin(PI*(y+1)/2/(yres+1))*sin(PI*(y+1)/2/(yres+1)) );
			double lambda = 1-__omegas[i]*(__lambda2k[x]+__lambda2m[y]);
			// double lambdaA = ( sin(PI*(x+1)/2/(xres+1))*sin(PI*(x+1)/2/(xres+1)) + sin(PI*(y+1)/2/(yres+1))*sin(PI*(y+1)/2/(yres+1)) );
			double lambdaA = (__lambda2k[x]+__lambda2m[y]);
			if(lambda*ckm[y*xres+x]>=0)
				sigckm+=(lambdaA * lambda*ckm[y*xres+x]);
			else
				sigckm-= (lambdaA * lambda*ckm[y*xres+x]);
		}
		if(sigckm<bestsigckm){
			bestomega = __omegas[i];
			bestsigckm = sigckm;
		}
	}
	std::cout<<"accurate omega="<<bestomega<<"sigckm"<<bestsigckm<<std::endl;
	return bestomega;	
}

double CheckWithResidual(double *error,double *residual, int xres, int yres){
	double *rck = project(residual, xres, yres);

	for(int k=0;k<xres;k++)
	for(int m=0;m<yres;m++){
		double lambdaA = 4.0 * ( __lambda2k[k]+__lambda2m[m] );
		rck[m*xres + k]/=lambdaA;
	}

	double *ck = project(error, xres, yres);
	for(int i=0;i<xres*yres;i++){
		// std::cout<<i<<"ckck\t"<<rck[i]<<'\t'<<ck[i]<<std::endl;
	}	

	double *errors = new double[xres * yres];

	for(int i=0;i<xres*yres;i++){
		int x = i%xres;
		int y = i/xres;
		double err_i=0;
		for(int k=0;k<xres;k++)
		for(int m=0;m<yres;m++){
			// err_i = err_i+sin(1.0*(k+1)*(x+1)*PI/(xres+1)) * sin(1.0*(m+1)*(y+1)*PI/(yres+1)) *ck[m*xres+k];
			err_i += __wkx[k*xres+x]*__wmy[m*yres+y]*ck[m*xres+k];
		}
		errors[i] = err_i;
		// std::cout<<err_i<<" "<<x<<" "<<y<<" -------------"<<std::endl;
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
		double lambdaA = 4.0*(__lambda2k[k]+__lambda2k[m]);
		rck[m*xres + k]/=lambdaA;
	}
	double *errors = new double[xres * yres];

	for(int i=0;i<xres*yres;i++){
		int x = i%xres;
		int y = i/xres;
		double err_i=0;
		for(int k=0;k<xres;k++)
		for(int m=0;m<yres;m++){
			err_i+= __wkx[k*xres+x] * __wmy[m*yres+y] *rck[m*xres+k];
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