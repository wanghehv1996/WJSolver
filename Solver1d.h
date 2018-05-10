#ifndef SOLVER_1D_H_
#define SOLVER_1D_H_
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include <iomanip>

#include "spectrum.h"
#include "function.h"

#define SQRTPI 1.77245385091
#define PI (3.14159265358979323846)
#define PRECISION 1e-6

#define PB_BOUND_Dirichlet	1
#define PB_BOUND_Neumann	2

class Solver1d{
public:
	Solver1d(){_xres=0; _btype=PB_BOUND_Dirichlet; _mask=NULL; _writeflag=false;};
	Solver1d(int xres, int btype):_xres(xres), _btype(btype), _mask(NULL),_writeflag(false){};
	
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
	int _btype;
	int _writeflag;
	double _precision = PRECISION;
	int _maxiter = 100;
	int GetIndex(int x){
		if(x<0||x>_xres-1)	return -1;
		return x;
	}

	double GetFromArray(const double *arr, int x){
		if(x>=0&&x<_xres)
			return arr[x];
		else
			return 0;
		int ind = GetIndex(x);
		return (x>=0)?arr[x]:0;
	}

	double GetNeighbCoef(int x){
		if(x<0||x>_xres-1)	return 0;
		if(_mask)
			if(_mask[GetIndex(x)]==PB_BOUND_Dirichlet||_mask[GetIndex(x)]==PB_BOUND_Neumann)return 0;
		return -1;
	}

	double GetCenterCoef(int x){
		int center = -2;
		//suppose we use dirichlet boundary outside
		int xminusIndex = GetIndex(x-1);
		if(xminusIndex>=0 && _mask)	center = (_mask[xminusIndex]==PB_BOUND_Neumann)?center+1:center;
		else				center = (_btype==PB_BOUND_Neumann)?center+1:center;
		int xplusIndex = GetIndex(x+1);
		if(xplusIndex<=_xres-1 && _mask)	center = (_mask[xplusIndex]==PB_BOUND_Neumann)?center+1:center;
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
		for(int i=0;i<_xres;i++)
			fout.write((char*)&(field[i]), sizeof(double));
		fout.close();
		return 0;
	}
};


class WeightedJacobiSolver1d : public Solver1d{
public:
	WeightedJacobiSolver1d(int xres, int btype){_xres=xres; _btype=btype; _mask=NULL;};
	WeightedJacobiSolver1d(int xres, int btype, double w){_w=w; _xres=xres; _btype=btype; _mask=NULL;};
	void SetW(double w){_w=w;}
	int solven = 0;

	int writeTrainData(const double *error, const double *residual, const double *ck, std::string filename){
		std::ofstream fout(filename.c_str(), std::ios::app | std::ios::binary);
		if(!fout){
			std::cout<<"WriteIter failed: not such file \'"<<filename<<"\'"<<std::endl;
			return 0;
		}
		for(int i=0;i<_xres;i++)
			fout.write((char*)&(error[i]), sizeof(double));
		for(int i=0;i<_xres;i++)
			fout.write((char*)&(residual[i]), sizeof(double));
		for(int i=0;i<_xres;i++)
			fout.write((char*)&(ck[i]), sizeof(double));
		fout.write((char*)&(_w), sizeof(double));
		
		fout.close();

		return 0;
	}

	virtual double* Solve(double* b){

		int iter = 0;
		int size = _xres;
		double *result = new double[size];
		double *resulttmp = new double[size];

		double sq_residual = 100;
		double sq_error = 0;
		double *residual = new double[size];
		double *error = new double[size];
		//1. Set Initial Guess
		for(int i=0;i<size;i++){
			// result[i]=function1d_1wave(i),resulttmp[i]=0;
			// result[i]=function1d_allzero(i),resulttmp[i]=0;
			// result[i]=function1d_2wave(i),resulttmp[i]=0;
			result[i]=function1d_random(i),resulttmp[i]=0;
			// result[i]=function1d_low(i),resulttmp[i]=0;
			// result[i]=function1d_wave3(i),resulttmp[i]=0;
		}

		while(iter<=_maxiter && sq_residual>_precision){
			sq_residual=0;
			sq_error =0;
			for(int i=0; i<size; i++){
				error[i] = function1d_allzero(i)-result[i];
				// error[i] = function1d_1wave(i)-result[i];
				sq_error += error[i]*error[i];
			}
			// if(iter==0)
			_w = GetOmega(error,size);
			// _w = 0.666666;
			// _w = GetOmega(iter+1,size);
			// _w = (iter%2!=0)?GetOmegaLow(result,size):GetOmegaHigh(result,size);
			// _w = GetOmegaSmoothW(result,size);
			// if(iter==0)
				// initRand();
			// _w = GetOmegaWithRand(result,size);
			// if(_w>1)
			// 	_w=1;

			//2. Do Iteration
			for(int i = 0;i<size;i++){
				int x = i;
				double all = 
					GetNeighbCoef(x-1)*GetFromArray(result,x-1)+
					GetNeighbCoef(x+1)*GetFromArray(result,x+1)+
					GetCenterCoef(x)*GetFromArray(result,x);
				//x'=x-w*(Ax)*D_inv+w*b*D_inv
				resulttmp[i] = GetFromArray(result,x) - _w*all/GetCenterCoef(x) + _w*b[i]/GetCenterCoef(x);				
			
				residual[i] = b[i]-all;
				sq_residual += residual[i]*residual[i];
				// std::cout<<i<<" , "<<iter<<" , "<<residual[i]<<" , "<<error[i]<<" , "<<resulttmp[i]<<std::endl;
			}
			
			// Check(error, residual, size);
			// SolveWithResidual(result, residual, size);
			//3. Record the result
			char buffer[256];
			sprintf(buffer,"%04i", solven);
			std::string number = std::string(buffer);
			// writeTrainData(error, residual, project(result,size), "./traindata/wjacobi_data"+number+".dat");
			writeTrainData(error, residual, project(result,size), "./traindata/wjacobi_databest.dat");

			// Cout Brief Information
			std::cout<<"| "<<iter<<" | "<<_w<<" | "<<sq_residual<<" | "<<sq_error<<" | "<<std::endl;
			// std::cout<<_w<<' '<<err<<' '<<sq_residual<<std::endl;

			double *swap = result;
			result = resulttmp;
			resulttmp = swap;
			
			iter++;
		}
		if(iter>_maxiter)
			std::cout<<"----------no converge"<<sq_residual<<std::endl;
		else
			std::cout<<iter<<":converge in "<<sq_residual<<std::endl;
		delete[] resulttmp;
		return result;
	};
private:
	double _w;
};

#endif