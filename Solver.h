#ifndef SOLVER_H_
#define SOLVER_H_
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <string.h>
#include <iomanip>

#define SQRTPI 1.77245385091
#define PI (3.14159265358979323846)
#define PRECISION 1e-6

#define PB_BOUND_Dirichlet	1
#define PB_BOUND_Neumann	2



class Solver{
public:
	Solver(){_xres=0; _yres=0; _btype=PB_BOUND_Dirichlet; _mask=NULL; _writeflag=false;};
	Solver(int xres, int yres, int btype):_xres(xres), _yres(yres), _btype(btype), _mask(NULL),_writeflag(false){};
	
	void SetWriteFlag(bool flag){_writeflag = flag;};
	void SetMask(int *mask){_mask = mask;};
	void SetPrecision(double precision){_precision = precision;};
	void SetMaxIter(int maxiter){_maxiter = maxiter;};
	void WriteLog(std::string filename, std::string log){
		std::ofstream ofs(filename.c_str(),std::ios::app);
		if(ofs)
			ofs<<log;
		ofs.close();
	};
	virtual double* Solve(double* b){};

protected:
	//Eigen::SparseMatrix<double> _A;
	int *_mask;
	int _xres;
	int _yres;
	int _btype;
	int _writeflag;
	double _precision = PRECISION;
	int _maxiter = 100;
	int GetIndex(int x, int y){
		if(x<0||x>_xres-1)	return -1;
		if(y<0||y>_yres-1)	return -1;
		return y*_xres+x;
	}

	double GetFromArray(double *arr, int x, int y){
		int ind = GetIndex(x,y);

		return (ind>=0 && arr!=NULL)?arr[ind]:0;
	}

	double GetNeighbCoef(int x, int y){
		if(x<0||x>_xres-1)	return 0;
		if(y<0||y>_yres-1)	return 0;
		if(_mask)
			if(_mask[GetIndex(x,y)]==PB_BOUND_Dirichlet||_mask[GetIndex(x,y)]==PB_BOUND_Neumann)return 0;
		return -1;
	}

	double GetCenterCoef(int x, int y){
		int center = -4;
		//suppose we use dirichlet boundary outside
		int xminusIndex = GetIndex(x-1, y);
		if(xminusIndex>=0 && _mask)	center = (_mask[xminusIndex]==PB_BOUND_Neumann)?center+1:center;
		else				center = (_btype==PB_BOUND_Neumann)?center+1:center;
		int xplusIndex = GetIndex(x+1, y);
		if(xplusIndex<=_xres-1 && _mask)	center = (_mask[xplusIndex]==PB_BOUND_Neumann)?center+1:center;
		else				center = (_btype==PB_BOUND_Neumann)?center+1:center;
		int yminusIndex = GetIndex(x, y-1);
		if(yminusIndex>=0 && _mask)	center = (_mask[yminusIndex]==PB_BOUND_Neumann)?center+1:center;
		else				center = (_btype==PB_BOUND_Neumann)?center+1:center;
		int yplusIndex = GetIndex(x, y+1);
		if(yplusIndex<=_yres-1 && _mask)	center = (_mask[yplusIndex]==PB_BOUND_Neumann)?center+1:center;
		else				center = (_btype==PB_BOUND_Neumann)?center+1:center;

		return -center;
		return center;

	}

public:
	int writeIter(const double *field, int itr, std::string prefix){
		char buffer[256];
		sprintf(buffer,"%04i", itr);
		std::string number = std::string(buffer);
		std::string filename = prefix+number+".dat";

		std::ofstream fout(filename.c_str(), std::ios::binary);
		if(!fout){
			std::cout<<"WriteIter failed: not such file \'"<<filename<<"\'"<<std::endl;
			return 0;
		}
		for(int i=0;i<_xres*_yres;i++)
			fout.write((char*)&(field[i]), sizeof(double));
		fout.close();
		return 0;
	}
};

class JacobiSolver : public Solver{
public:
	JacobiSolver(int xres, int yres, int btype){_xres=xres; _yres=yres; _btype=btype; _mask=NULL;};
	virtual double* Solve(double* b){


		double diff = 100;
		int iter = 0;
		int size = _xres*_yres;
		double *result = new double[size];
		double *resulttmp = new double[size];
		for(int i=0;i<size;i++)
			result[i]=0,resulttmp[i]=0;
		while(iter<=_maxiter && diff>_precision){
			diff=0;
			for(int i = 0;i<size;i++){
				int x = i%_xres;
				int y = i/_xres;
			
				double all = 
					GetNeighbCoef(x-1,y)*GetFromArray(result,x-1,y)+
					GetNeighbCoef(x+1,y)*GetFromArray(result,x+1,y)+
					GetNeighbCoef(x,y-1)*GetFromArray(result,x,y-1)+
					GetNeighbCoef(x,y+1)*GetFromArray(result,x,y+1);
				resulttmp[i] = (b[i]- all)/GetCenterCoef(x,y);
			
				diff+= (b[i]-(GetCenterCoef(x,y)*result[i]+all))*(b[i]-(GetCenterCoef(x,y)*result[i]+all));
			}
			//write out each iteration result
			if(_writeflag)
				writeIter(resulttmp, iter, "./jacobi_");

			std::cout<<iter<<" iteration: diff="<<diff<<std::endl;
			double *swap = result;
			result = resulttmp;
			resulttmp = swap;
			
			iter++;

		}
		delete[] resulttmp;
		return result;
	};
};

class WeightedJacobiSolver : public Solver{
public:
	WeightedJacobiSolver(int xres, int yres, int btype){_xres=xres; _yres=yres; _btype=btype; _mask=NULL;};
	WeightedJacobiSolver(int xres, int yres, int btype, double w){_w=w; _xres=xres; _yres=yres; _btype=btype; _mask=NULL;};
	void SetW(double w){_w=w;}
	void SetWriteError(bool b){_writeerror=b;}
	virtual double* Solve(double* b){
		std::ofstream ofs;
		if(_writeerror){
			ofs.open("./jacobi_err.csv",std::ios::app);
			if(!ofs){
				std::cout<<"cannot open jacobi_err"<<std::endl;
				return NULL;
			}
			ofs<<'\n'<<_w<<',';
		}

		double diff = 100;
		int iter = 0;
		int size = _xres*_yres;
		double *result = new double[size];
		double *resulttmp = new double[size];
		for(int i=0;i<size;i++)
			result[i]=0,resulttmp[i]=0;

		while(iter<=_maxiter && diff>_precision){
			diff=0;
			for(int i = 0;i<size;i++){
				int x = i%_xres;
				int y = i/_xres;
			
				double all = 
					GetNeighbCoef(x-1,y)*GetFromArray(result,x-1,y)+
					GetNeighbCoef(x+1,y)*GetFromArray(result,x+1,y)+
					GetNeighbCoef(x,y-1)*GetFromArray(result,x,y-1)+
					GetNeighbCoef(x,y+1)*GetFromArray(result,x,y+1)+
					GetCenterCoef(x,y)*GetFromArray(result,x,y);
				//x'=x-w*(Ax)*D_inv+w*b*D_inv
				resulttmp[i] = GetFromArray(result,x,y) - _w*all/GetCenterCoef(x,y) + _w*b[i]/GetCenterCoef(x,y);
			
				diff+= (b[i]-all)*(b[i]-all);
			}
			//write out each iteration result
			if(_writeflag)
				writeIter(resulttmp, iter, "./wjacobi_");

			std::cout<<iter<<" iteration: diff="<<diff<<std::endl;
			if(_writeerror){
				ofs<<diff<<',';
			}
			double *swap = result;
			result = resulttmp;
			resulttmp = swap;
			
			iter++;

		}
		// ofs.close();
		delete[] resulttmp;
		return result;
	};
private:
	double _w;
	bool _writeerror = false;
};

class GaussSeidelSolver : public Solver{
public:
	GaussSeidelSolver(int xres, int yres, int btype){_xres=xres; _yres=yres; _btype=btype; _mask=NULL;};
	virtual double* Solve(double* b){


		double diff = 100;
		int iter = 0;
		int size = _xres*_yres;
		double *result = new double[size];
		for(int i=0;i<size;i++)
			result[i]=0;
		while(iter<=_maxiter && diff>_precision){
			diff=0;
			for(int i = 0;i<size;i++){
				double origin = result[i];
				int x = i%_xres;
				int y = i/_xres;

				double all = 
					GetNeighbCoef(x-1,y)*GetFromArray(result,x-1,y)+
					GetNeighbCoef(x+1,y)*GetFromArray(result,x+1,y)+
					GetNeighbCoef(x,y-1)*GetFromArray(result,x,y-1)+
					GetNeighbCoef(x,y+1)*GetFromArray(result,x,y+1);
				diff+= (b[i]-(GetCenterCoef(x,y)*result[i]+all))*(b[i]-(GetCenterCoef(x,y)*result[i]+all));
				result[i] = (b[i]- all)/GetCenterCoef(x,y);
			}
			if(_writeflag)
				writeIter(result, iter, "./GaussSeidel_");
			std::cout<<iter<<" iteration: diff="<<diff<<std::endl;
			
			iter++;

		}
		return result;
	};
};

class SORSolver : public Solver{
public:
	SORSolver(int xres, int yres, int btype){_w=1; _xres=xres; _yres=yres; _btype=btype; _mask=NULL;};
	SORSolver(int xres, int yres, int btype, double w){_w=w; _xres=xres; _yres=yres; _btype=btype; _mask=NULL;};
	void SetW(double w){_w=w;}
	void SetWriteError(bool b){_writeerror=b;}
	virtual double* Solve(double* b){
		std::ofstream ofs;
		if(_writeerror){
			ofs.open("./sor_err.csv",std::ios::app);
			if(!ofs){
				std::cout<<"cannot open sor_err"<<std::endl;
				return NULL;
			}
			ofs<<'\n'<<_w<<',';
		}
		double diff = 100;
		int iter = 0;
		int size = _xres*_yres;
		double *result = new double[size];
		for(int i=0;i<size;i++)
			result[i]=0;
		while(iter<=_maxiter && diff>_precision){
			diff=0;
			for(int i = 0;i<size;i++){
				double origin = result[i];
				int x = i%_xres;
				int y = i/_xres;

				double all = 
					GetNeighbCoef(x-1,y)*GetFromArray(result,x-1,y)+
					GetNeighbCoef(x+1,y)*GetFromArray(result,x+1,y)+
					GetNeighbCoef(x,y-1)*GetFromArray(result,x,y-1)+
					GetNeighbCoef(x,y+1)*GetFromArray(result,x,y+1);

				diff+= (b[i]-(GetCenterCoef(x,y)*result[i]+all))*(b[i]-(GetCenterCoef(x,y)*result[i]+all));

				result[i] = (1-_w)*origin + _w*(b[i]- all)/GetCenterCoef(x,y);
			}
			if(_writeflag)
				writeIter(result, iter, "./SOR_");
			if(_writeerror){
					ofs<<diff<<',';
				}
			std::cout<<iter<<" iteration: diff="<<diff<<std::endl;
			
			iter++;

		}
		return result;
	};
private:
	double _w;
	double _writeerror=false;
};

class CGSolver : public Solver{
public:
	CGSolver(int xres, int yres, int btype){_xres=xres; _yres=yres; _btype=btype; _mask=NULL;};
	//CGSolver(int xres, int yres, int btype, double w){_w=w; _xres=xres; _yres=yres; _btype=btype; _mask=NULL;};
	//void SetW(double w){_w=w;}
	virtual double* Solve(double* b){
		double diff = 100;
		int iter = 0;
		int size = _xres*_yres;
		double *result = new double[size];

		double *xk = new double[size];
		double *rk = new double[size];
		double *vk = new double[size];
		double *Avk = new double[size];
		for(int i=0;i<size;i++)
			result[i]=0,xk[i]=0,rk[i]=0,vk[i]=0,Avk[i]=0;
		double alphak, betak;
		for(int i=0;i<size;i++)
			rk[i]=b[i],vk[i]=b[i];
		double rkm_rkm;
		double rk_rk;
		for(int i = 0;i<size;i++)
			rk_rk += rk[i]*rk[i];


		while(iter<=_maxiter && diff>_precision){
			// diff=0;
			//get A*vk
			for(int i = 0;i<size;i++){
				int x = i%_xres;
				int y = i/_xres;
				Avk[i] = 
					GetNeighbCoef(x-1,y)*GetFromArray(vk,x-1,y)+
					GetNeighbCoef(x+1,y)*GetFromArray(vk,x+1,y)+
					GetNeighbCoef(x,y-1)*GetFromArray(vk,x,y-1)+
					GetNeighbCoef(x,y+1)*GetFromArray(vk,x,y+1)+
					GetCenterCoef(x,y)*GetFromArray(vk,x,y);
			}
			//get alphak
			
			rkm_rkm = rk_rk;
			double vk_Avk = 0;
			for(int i = 0;i<size;i++){
				vk_Avk+=vk[i]*Avk[i];
			}
			alphak = rkm_rkm/vk_Avk;
			//get new xk
			for(int i = 0;i<size;i++){
				xk[i] = xk[i]+alphak*vk[i];
				rk[i] = rk[i]-alphak*Avk[i];
			}

			//get new betak
			rk_rk = 0;
			for(int i = 0;i<size;i++){
				rk_rk += rk[i]*rk[i];
			}
			betak = rk_rk/rkm_rkm;

			diff = rk_rk;
			//get vk
			for(int i = 0;i<size;i++){
				vk[i] = rk[i]+betak*vk[i];
			}
			if(_writeflag)
			{
				char buffer[256];
				sprintf(buffer,"%04i", iter);
				std::string number = std::string(buffer);
				std::string filename = "./CG/CG_"+number+".dat";

				std::ofstream fout(filename.c_str(), std::ios::binary);
				for(int i=0;i<_xres*_yres;i++)
					fout.write((char*)&(xk[i]), sizeof(double));
				fout.close();
			}

			std::cout<<iter<<" iteration: diff="<<diff<<std::endl;
			iter++;
		}
		return xk;
	};
};

class ICPCGSolver : public Solver{
public:
	ICPCGSolver(int xres, int yres, int btype){_xres=xres; _yres=yres; _btype=btype; _mask=NULL;};
	//CGSolver(int xres, int yres, int btype, double w){_w=w; _xres=xres; _yres=yres; _btype=btype; _mask=NULL;};
	virtual double* Solve(double* b){
		incompleteCholesky();
		double diff = 100;
		int iter = 0;
		int size = _xres*_yres;
		double *result = new double[size];

		double *xk = new double[size];
		double *rk = new double[size];
		double *vk = new double[size];
		double *Avk = new double[size];
		for(int i=0;i<size;i++)
			result[i]=0,xk[i]=0,rk[i]=0,vk[i]=0,Avk[i]=0;
		double alphak, betak;
		for(int i=0;i<size;i++)
			rk[i]=b[i];
		double rkm_rkm = 0;
		double rk_rk = 0;
		for(int i = 0;i<size;i++)
			rk_rk += rk[i]*rk[i];

		double rkm_zkm = 0;
		double rk_zk = 0;
		double *zk = new double[size];
		for(int i=0;i<size;i++)
			zk[i]=0;
		applyPrecond(zk, rk);
		for(int i = 0;i<size;i++){
			vk[i]=zk[i];
		}
		for(int i = 0;i<size;i++)
			rk_zk += rk[i]*zk[i];

		while(iter<=_maxiter && diff>_precision){
			// diff=0;
			//get A*vk
			for(int i = 0;i<size;i++){
				int x = i%_xres;
				int y = i/_xres;
				Avk[i] = 
					GetNeighbCoef(x-1,y)*GetFromArray(vk,x-1,y)+
					GetNeighbCoef(x+1,y)*GetFromArray(vk,x+1,y)+
					GetNeighbCoef(x,y-1)*GetFromArray(vk,x,y-1)+
					GetNeighbCoef(x,y+1)*GetFromArray(vk,x,y+1)+
					GetCenterCoef(x,y)*GetFromArray(vk,x,y);
			}
			//get alphak
			
			rkm_rkm = rk_rk;
			//-------------------------
			rkm_zkm = rk_zk;
			//-------------------------
			double vk_Avk = 0;
			for(int i = 0;i<size;i++){
				vk_Avk+=vk[i]*Avk[i];
			}
			//alphak = rkm_rkm/vk_Avk;
			//-------------------------
			alphak = rkm_zkm/vk_Avk;
			//-------------------------
			//get new xk
			for(int i = 0;i<size;i++){
				xk[i] = xk[i]+alphak*vk[i];
				rk[i] = rk[i]-alphak*Avk[i];
			}
			//get new betak
			rk_rk = 0;
			for(int i = 0;i<size;i++){
				rk_rk += rk[i]*rk[i];
			}
			//-------------------------
			applyPrecond(zk, rk);
			rk_zk=0;
			for(int i = 0;i<size;i++)
				rk_zk += rk[i]*zk[i];
			//-------------------------
			//betak = rk_rk/rkm_rkm;
			//-------------------------
			betak = rk_zk/rkm_zkm;
			//-------------------------

			diff = rk_rk;
			// diff = rk_zk;

			//get vk
			for(int i = 0;i<size;i++){
				// vk[i] = rk[i]+betak*vk[i];
			//-------------------------
				vk[i] = zk[i]+betak*vk[i];
			//-------------------------
			}
			if(_writeflag)
				writeIter(xk, iter, "./ICPCG_");
			std::cout<<iter<<" iteration: diff="<<diff<<std::endl;
			
			iter++;

		}
		return xk;
	};
private:
	double *Lii, *Lxp1, *Lyp1;


	double* incompleteCholesky(){
		int size = _xres*_yres;
		Lii = new double[size];
		Lxp1 = new double[size];
		Lyp1 = new double[size];
		
		double *L = new double[size*size];
		for(int i = 0;i<size;i++){
			int x = i%_xres;
			int y = i/_xres;
			Lii[i] = GetCenterCoef(x,y);
			Lxp1[i] = GetNeighbCoef(x+1,y);
			Lyp1[i] = GetNeighbCoef(x,y+1);
		}
		for(int i = 0;i<size;i++){
			Lii[i] = sqrt(Lii[i]);
			Lxp1[i]*= 1.0/Lii[i];
			Lyp1[i]*= 1.0/Lii[i];

			if(i<size-1)
				Lii[i]-= Lxp1[i]*Lxp1[i];
			if(i<size-_xres)
				Lii[i]-= Lyp1[i]*Lyp1[i];
		}
		for(int i = 0;i<size;i++){
			Lii[i] = 1/Lii[i];
		}

	}

	double *applyPrecond(double *dst, double *var){
		int size = _xres*_yres;
		for(int i = 0;i<size;i++){
			int x = i%_xres;
			int y = i/_xres;
			dst[i] = var[i];
			if(i>0)
				dst[i] -= dst[i-1] * Lxp1[i-1];
			if(i>_xres)
				dst[i] -= dst[i-_xres] * Lyp1[i - _xres];
			dst[i]*=Lii[i];

		}

		for(int i = size-1;i>=0;i--){
			int x = i%_xres;
			int y = i/_xres;
			if(i<size-1)
				dst[i] -= dst[i+1] * Lxp1[i];
			if(i<size-_xres)
				dst[i] -= dst[i+_xres] * Lyp1[i];
			dst[i]*=Lii[i];
		}
		return dst;
	}

	double GetM(int x, int y){
		//return m[y*_xres+x];
		
		if(x==y){
			double res=Lii[x]*Lii[x];
			if(x>0)
				res+=Lxp1[x-1]*Lxp1[x-1];;
			if(x>=_xres)
				res+=Lyp1[x-_xres]*Lyp1[x-_xres];
			//std::cout<<res<<'\t';
			return res;
		}
		else if(x+1==y){
			double res = Lii[x]*Lxp1[x];
			//std::cout<<res<<'\t';
			return res;
		}
		else if(x==y+1){
			double res = Lii[y]*Lxp1[y];
			//std::cout<<res<<'\t';
			return res;
		}
		else if(x+_xres==y){
			double res = Lii[x]*Lyp1[x];
			//std::cout<<res<<'\t';
			return res;
		}else if(y+_xres==x){
			double res = Lii[y]*Lyp1[y];
			//std::cout<<res<<'\t';
			return res;
		}
		else if(x+_xres-1==y&&x>0){
			double res = Lxp1[x-1]*Lyp1[x-1];
			//std::cout<<res<<'\t';
			return res;
		}else if(y+_xres-1==x&&y>0){
			double res = Lxp1[y-1]*Lyp1[y-1];
			//std::cout<<res<<'\t';
			return res;
		}else return 0;//std::cout<<0<<'\t';
	}
};








double* incompleteCholesky(int res)
{	
	int num = res*res;
	double *A = new double[num*num];
	for(int i=0;i<num;i++){
		A[i*num+i]=4;
		if(i%res>0)			A[i*num+i-1]=-1;
		if(i%res < res-1)	A[i*num+i+1]=-1;
		if(i/res>0)			A[i*num+i-res]=-1;
		if(i/res < res-1)	A[i*num+i+res]=-1;
	}

	// std::cout<<"A="<<std::setprecision(4)<<std::endl;
	// for(int y=0;y<num;y++){
	// 	for(int x=0;x<num;x++)
	// 		std::cout<<A[y*num+x]<<'\t';
	// 	std::cout<<std::endl;
	// }


	double* L= new double[num*num];
	

	double temp;
	for (int i = 0; i < num; i++)
	{
		temp = 0;
		for (int j = 0; j <= i-1; j++)
		{
			temp = temp + L[i*num + j]*L[i*num + j];
		}
		L[i*num + i] = sqrt(A[i*num + i] - temp);
		
		for (int j = i+1; j < num; j++)
		{
			if (A[i*num + j]!=0 && L[i*num + i]!=0)
			{
				temp = 0;
				for (int k = 0; k <= i-1; k++)
				{
					temp = temp + L[i*num + k]*L[j*num + k];
				}
				L[j*num + i] = (A[i*num + j] - temp)/L[i*num + i];
			}
		}
	}

	std::cout<<"L="<<std::endl;
	for(int y=0;y<num;y++){
		for(int x=0;x<num;x++)
			std::cout<<L[y*num+x]<<'\t';
		std::cout<<std::endl;
	}

	double* Lr= new double[num*num];
	for (int i = 0; i < num; i++)
	{

		Lr[i*num+i] = 1.0/L[i*num+i];
		for(int j = 0; j < i; j++){
			if(L[i*num+j]<0.00001)
				continue;
			for(int k = j; k < i; k++)
				Lr[i*num+j]+= L[i*num+k]*Lr[k*num+j];
			Lr[i*num+j]*=(-Lr[i*num+i]);
		}
	}	
	// std::cout<<"LR="<<std::endl;
	// for(int y=0;y<num;y++){
	// 	for(int x=0;x<num;x++)
	// 		if(Lr[y*num+x]>0.001)
	// 		std::cout<<1<<'\t';
	// 		else
	// 			std::cout<<Lr[y*num+x]<<'\t';
	// 	std::cout<<std::endl;
	// }

	double* LrL= new double[num*num];
	for (int i = 0; i < num; i++)
	{
		for(int j = 0; j < num; j++){
			LrL[i*num+j]=0;
			for(int k = 0; k < num; k++)
				LrL[i*num+j]+= Lr[i*num+k]*L[k*num+j];
		}
	}	
	// std::cout<<"LRL="<<std::endl;
	// for(int y=0;y<num;y++){
	// 	for(int x=0;x<num;x++)
	// 		std::cout<<LrL[y*num+x]<<'\t';
	// 	std::cout<<std::endl;
	// }



	double* M = new double[num*num];
	for (int i = 0; i < num; i++)
	{
		for (int j = 0; j < num; j++)
		{
			M[i*num + j]  = 0;
			for (int k = 0; k < num; k++)
			{
				M[i*num + j] = M[i*num + j] + Lr[i*num + k]*Lr[j*num + k];
			}
		}
	}

	std::cout<<"M="<<std::endl;
	for(int y=0;y<num;y++){
		for(int x=0;x<num;x++)
			std::cout<<M[y*num+x]<<'\t';
		std::cout<<std::endl;
	}

	double* MA = new double[num*num];
	for(int y=0;y<num;y++){
		for(int x=0;x<num;x++){
			double mai=0;
			for(int i=0;i<num;i++){
				//std::cout<<M[y*num+x]*A[]<<'\t';
				mai+=M[y*num+i]*A[i*num+x];
			}
			std::cout<<mai<<'\t';
		}
		std::cout<<std::endl;
	}	
	return M;
}
#endif