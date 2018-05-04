#ifndef SPECTRUM_H_
#define SPECTRUM_H_

#include <eigen3/Eigen/IterativeLinearSolvers>
#include <eigen3/Eigen/Eigen>
#include <eigen3/Eigen/LU>
#include <eigen3/Eigen/QR>

#include<math.h>  
#include<stdlib.h>

using namespace std;  
#define PI (3.14159265358979323846)  

double *project(const double *error, int res){
	double *ck = new double[res];
	for(int k=0;k<res;k++){
		double err_wk=0;
		double wk_wk=0;
		for(int i=0;i<res;i++){
			//err*w_k
			err_wk += error[i]*sin(1.0*(k+1)*(i+1)*PI/(res+1));
			//|w_k|
			wk_wk += sin(1.0*(k+1)*(i+1)*PI/(res+1))*sin(1.0*(k+1)*(i+1)*PI/(res+1));
		}
		ck[k] = err_wk/(wk_wk);
	}
	return ck;
}

double *project(const double *error, int xres, int yres){
	double *ck = new double[xres*yres];
	for(int k=0;k<xres;k++)
	for(int m=0;m<yres;m++){
		double err_wkm=0;
		double wk_wk=0;
		for(int i=0;i<xres*yres;i++){
			int x = i%xres;
			int y = i/xres;
			//err*w_k
			err_wkm += error[i] * sin(1.0*(k+1)*(x+1)*PI/(xres+1)) * sin(1.0*(m+1)*(y+1)*PI/(yres+1));
			//|w_k|
			wk_wk += sin(1.0*(k+1)*(x+1)*PI/(xres+1)) * sin(1.0*(k+1)*(x+1)*PI/(xres+1)) *
					sin(1.0*(m+1)*(y+1)*PI/(yres+1)) * sin(1.0*(m+1)*(y+1)*PI/(yres+1));
		}
		ck[m *xres + k] = err_wkm/(wk_wk);

	}
	return ck;
}

double GetOmega(int k, int res){
	// std::cout<<k<<std::endl;
	if(k>=res+1||k<=0)
		return 1;
	double omega = 1./(2*sin(k*PI/2/(res+1))*sin(k*PI/2/(res+1)));
	return omega;
}

double GetOmega(const double *error, int res){
	double *ck = new double[res];
	for(int k=0;k<res;k++){
		double err_wk=0;
		double wk_wk=0;
		for(int i=0;i<res;i++){
			//err*w_k
			err_wk += error[i]*sin(1.0*(k+1)*(i+1)*PI/(res+1));
			//|w_k|
			wk_wk += sin(1.0*(k+1)*(i+1)*PI/(res+1))*sin(1.0*(k+1)*(i+1)*PI/(res+1));
		}
		ck[k] = err_wk/(wk_wk);
	}
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
	double maxk = 0;
	for(int k=0;k<res;k++){
		// std::cout<<ck[k]<<' ';
		if(maxck<abs(ck[k])){
			maxck = abs(ck[k]);
			maxk = k+1;
		}
	}
	// std::cout<<std::endl<<maxk<<std::endl;
	double omega = 1./(2*sin(maxk*PI/2/(res+1))*sin(maxk*PI/2/(res+1)));
	return omega;
}

//2Dimension
double GetOmega(const double *error, int xres, int yres){
	double *ck = new double[xres*yres];
	
		
	for(int k=0;k<xres;k++)
	for(int m=0;m<yres;m++){
		double err_wk=0;
		double wk_wk=0;
		for(int i=0;i<xres*yres;i++){
			int x = i%xres;
			int y = i/xres;
			//err*w_k
			err_wk += error[i] * sin(1.0*(k+1)*(x+1)*PI/(xres+1)) * sin(1.0*(m+1)*(y+1)*PI/(yres+1));
			//|w_k|
			wk_wk += sin(1.0*(k+1)*(x+1)*PI/(xres+1)) * sin(1.0*(k+1)*(x+1)*PI/(xres+1)) *
					sin(1.0*(m+1)*(y+1)*PI/(yres+1)) * sin(1.0*(m+1)*(y+1)*PI/(yres+1));
		}

		ck[m*xres+k] = err_wk/(wk_wk);
	}
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
	double maxk = 0;
	double maxm = 0;
	for(int k=0;k<xres;k++)
	for(int m=0;m<yres;m++){
		if(maxck<abs(ck[m*xres+k])){
			maxck = abs(ck[m*xres+k]);
			maxk = k+1;
			maxm = m+1;
		}
	}
	// std::cout<<ck[0*xres+0]<<' '<<ck[0*xres+1]<<' '<<ck[19*xres+19]<<std::endl;
	
	double omega = 1./( sin(maxk*PI/2/(xres+1))*sin(maxk*PI/2/(xres+1)) + sin(maxm*PI/2/(yres+1))*sin(maxm*PI/2/(yres+1)));
	std::cout<<std::endl<<maxk<<' '<<maxm<<' '<<omega<<std::endl;
	return omega;
}

double GetOmegaLow(double *error, int res){
	double *ck = project(error, res);
	double maxck = 0;
	double maxk = 0;
	for(int k=0;k<=(res)/2;k++){
		// std::cout<<ck[k]<<' ';
		if(maxck<abs(ck[k])){
			maxck = abs(ck[k]);
			maxk = k+1;
		}
	}
	// std::cout<<std::endl<<maxk<<std::endl;
	double omega = 1./(2*sin(maxk*PI/2/(res+1))*sin(maxk*PI/2/(res+1)));
	return omega;
}
double GetOmegaHigh(double *error, int res){
	double *ck = project(error, res);
	double maxck = 0;
	double maxk = 0;
	for(int k=(res)/2;k<(res);k++){
		// std::cout<<ck[k]<<' ';
		if(maxck<abs(ck[k])){
			maxck = abs(ck[k]);
			maxk = k+1;
		}
	}
	// std::cout<<std::endl<<maxk<<std::endl;
	double omega = 1./(2*sin(maxk*PI/2/(res+1))*sin(maxk*PI/2/(res+1)));
	return omega;
}

double k1 = 32/2;
double k2 = 32/2;

double GetOmegaSmoothK(double *error, int res){
	double *ck = project(error, res);
	double maxck = 0;
	double maxk = 0;

	for(int k=0;k<(res);k++){
		
		if(maxck<abs(ck[k])){
			maxck = abs(ck[k]);
			maxk = k+1;
		}
	}
	

	double choosek = maxk*0.9+k2*0.1;
		std::cout<<maxk<<' '<<k2<<' '<<choosek<<std::endl;
	double omega = 1./(2*sin(choosek*PI/2/(res+1))*sin(choosek*PI/2/(res+1)));

	k1=k2;
	k2=maxk;
	return omega;
}

double w1 = 1;
double w2 = 1;
double GetOmegaSmoothW(double *error, int res){
	double *ck = project(error, res);
	double maxck = 0;
	double maxk = 0;

	for(int k=0;k<(res);k++){
		// std::cout<<ck[k]<<' ';
		if(maxck<abs(ck[k])){
			maxck = abs(ck[k]);
			maxk = k+1;
		}
	}
	
	double omega = 1./(2*sin(maxk*PI/2/(res+1))*sin(maxk*PI/2/(res+1)));
	double chooseomega = omega*0.9+w2*0.1;
	std::cout<<omega<<' '<<w2<<' '<<chooseomega<<std::endl;
	w1=w2;
	w2=omega;
	return chooseomega;
}


double GetOmegaSmooth(double *error, int res){
	double *ck = project(error, res);
	double maxck = 0;
	double maxk = 0;

	double maxck2 = 0;
	double maxk2 = 0;
	for(int k=0;k<(res);k++){
		// std::cout<<ck[k]<<' ';
		if(maxck<abs(ck[k])){
			maxck = abs(ck[k]);
			maxk = k+1;
		}
	}
	for(int k=0;k<(res);k++){
		// std::cout<<ck[k]<<' ';
		if(maxck2<abs(ck[k]) && k!=maxk-1){
			maxck2 = abs(ck[k]);
			maxk2 = k+1;
		}
	}
	double maxkavg = (maxk+maxk2)/2;
	double omega = 1./(2*sin(maxkavg*PI/2/(res+1))*sin(maxkavg*PI/2/(res+1)));
	// w1=w2;
	// w2=omega;
	return omega;
}

double GetOmegaWithRand(double *error, int res){
	double *ck = project(error, res);
	double maxck = 0;
	double maxk = 0;

	for(int k=0;k<(res);k++){
		// std::cout<<ck[k]<<' ';
		if(maxck<abs(ck[k])){
			maxck = abs(ck[k]);
			maxk = k+1;
		}
	}
	// std::cout<<maxk;
	// maxk=maxk+(rand()%9-4);
	// if(maxk<1)
	// 	maxk=1;
	// if(maxk>=res)
	// 	maxk=res;
	// std::cout<<' '<<maxk<<std::endl;

	double omega = 1./(2*sin(maxk*PI/2/(res+1))*sin(maxk*PI/2/(res+1)));
	// std::cout<<omega;
	omega*=((1.0*(rand()%20))+90.0)/100.0;
	// std::cout<<' '<<omega<<' '<<((1.0*(rand()%10))+95.0)/100.0<<std::endl;
	return omega;
}

double Check(double *error, double *residual, int res){
	double *rck = project(residual, res);
	for(int k=0;k<res;k++){
		// std::cout<<residual[k]<<' '<<error[k]<<std::endl;
	}
	double *ck = project(error, res);
	for(int k=0;k<res;k++){
		double lambdaA=4.0*sin(1.0*(k+1)*PI/2/(res+1))*sin(1.0*(k+1)*PI/2/(res+1));
		
		std::cout<<k<<"\t"<<ck[k]*lambdaA/rck[k]<<"rck="<<rck[k]<<";ldack="<<ck[k]*lambdaA<<";ck="<<ck[k]<<std::endl;
	}
	return 0;
}

double CheckWithResidual(double *error,double *residual, int res){
	double *rck = project(residual, res);
	for(int k=0;k<res;k++){
		double lambdaA=4.0*sin(1.0*(k+1)*PI/2/(res+1))*sin(1.0*(k+1)*PI/2/(res+1));
		rck[k]/=lambdaA;
	}
	double *errors = new double[res];

	for(int i=0;i<res;i++){
		double err_i=0;
		for(int k=0;k<res;k++){
			err_i+= sin(1.0*(k+1)*(i+1)*PI/(res+1))*rck[k];
		}
		errors[i] = err_i;
	}

	for(int i=0;i<res;i++){
		std::cout<<i<<"\t"<<errors[i]<<'\t'<<error[i]<<std::endl;
	}
	return 0;
}

//2Dimension version
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

double SolveWithResidual(double *result,double *residual, int res){
	double *rck = project(residual, res);
	for(int k=0;k<res;k++){
		double lambdaA=4.0*sin(1.0*(k+1)*PI/2/(res+1))*sin(1.0*(k+1)*PI/2/(res+1));
		rck[k]/=lambdaA;
	}
	double *errors = new double[res];

	for(int i=0;i<res;i++){
		double err_i=0;
		for(int k=0;k<res;k++){
			err_i+= sin(1.0*(k+1)*(i+1)*PI/(res+1))*rck[k];
		}
		errors[i] = err_i;
	}

	for(int i=0;i<res;i++){
		result[i]+=errors[i];
	}
	delete []rck;
	delete []errors;
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