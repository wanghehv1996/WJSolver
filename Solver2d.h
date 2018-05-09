#ifndef SOLVER_2D_H_
#define SOLVER_2D_H_

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include <iomanip>

#include "spectrum2d.h"
#include "function.h"

#define SQRTPI 1.77245385091
#define PI (3.14159265358979323846)
#define PRECISION 1e-6

#define PB_BOUND_Dirichlet	1
#define PB_BOUND_Neumann	2



class Solver2d{
public:
	Solver2d(){_xres=0; _yres=0; _btype=PB_BOUND_Dirichlet; _mask=NULL; _writeflag=false;};
	Solver2d(int xres, int yres, int btype):_xres(xres), _yres(yres), _btype(btype), _mask(NULL),_writeflag(false){};
	
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


class WeightedJacobiSolver2d : public Solver2d{
public:
	WeightedJacobiSolver2d(int xres, int yres, int btype){_xres=xres; _yres=yres; _btype=btype; _mask=NULL;};
	WeightedJacobiSolver2d(int xres, int yres, int btype, double w){_w=w; _xres=xres; _yres=yres; _btype=btype; _mask=NULL;};
	void SetW(double w){_w=w;}

	int writeTrainData(const double *error, const double *residual, const double *ck, std::string filename){
		std::ofstream fout(filename.c_str(), std::ios::app | std::ios::binary);
		if(!fout){
			std::cout<<"WriteIter failed: not such file \'"<<filename<<"\'"<<std::endl;
			return 0;
		}
		for(int i=0; i < _xres * _yres; i++)
			fout.write((char*)&(error[i]), sizeof(double));
		for(int i=0; i < _xres * _yres; i++)
			fout.write((char*)&(residual[i]), sizeof(double));
		for(int i=0; i < _xres * _yres; i++)
			fout.write((char*)&(ck[i]), sizeof(double));
		fout.write((char*)&(_w), sizeof(double));
		
		fout.close();

		return 0;
	}

	int writeTrainDataRCk(const double *residual, const double *ck, std::string filename){
		std::ofstream fout(filename.c_str(), std::ios::app | std::ios::binary);
		if(!fout){
			std::cout<<"WriteIter failed: not such file \'"<<filename<<"\'"<<std::endl;
			return 0;
		}
		for(int i=0; i < _xres * _yres; i++)
			fout.write((char*)&(residual[i]), sizeof(double));
		for(int i=0; i < _xres * _yres; i++)
			fout.write((char*)&(ck[i]), sizeof(double));
		
		fout.close();

		return 0;
	}

	virtual double* Solve(double* b){

		int iter = 0;
		int size = _xres*_yres;
		double *result = new double[size];
		double *resulttmp = new double[size];

		double sq_residual = 100;
		double sq_error = 0;
		double *residual = new double[size];
		double *error = new double[size];
		
		//1. Set Initial Guess
		for(int y=0;y<_yres;y++)
		for(int x=0;x<_xres;x++){
			// result[i]=function1d_1wave(i),resulttmp[i]=0;
			// result[i]=function1d_allzero(i),resulttmp[i]=0;
			// result[i]=function1d_2wave(i),resulttmp[i]=0;
			result[ GetIndex(x,y) ] = function2d_random(x,y),resulttmp[ GetIndex(x, y) ]=0;
			// result[ GetIndex(x,y) ] = function2d_waves(x,y),resulttmp[ GetIndex(x, y) ]=0;
			if(_mask)
				if(_mask[GetIndex(x,y)]==PB_BOUND_Dirichlet||_mask[GetIndex(x,y)]==PB_BOUND_Neumann)
					result[ GetIndex(x,y) ] =0;
		}

		while(iter<=_maxiter && sq_residual>_precision){
			sq_residual=0;
			sq_error =0;
			for(int i=0; i<size; i++){
				int x = i%_xres;
				int y = i/_xres;
				error[i] = function2d_allzero(x,y)-result[i];
				// error[i] = function1d_1wave(i)-result[i];
				sq_error += error[i]*error[i];

				double all = 
					GetNeighbCoef(x-1,y)*GetFromArray(result,x-1,y)+
					GetNeighbCoef(x+1,y)*GetFromArray(result,x+1,y)+
					GetNeighbCoef(x,y-1)*GetFromArray(result,x,y-1)+
					GetNeighbCoef(x,y+1)*GetFromArray(result,x,y+1)+
					GetCenterCoef(x,y)*GetFromArray(result,x,y);
				residual[i] = b[i]-all;
			}

			// _w = GetOmegaNeigh(error,_xres, _yres);
			// _w = GetBestOmega(error,_xres, _yres);
			// _w = GetBestROmega(error,_xres, _yres);
			// _w = GetBestOmegaNeigh(error,_xres, _yres);
			// _w = GetBestROmegaNeigh(error,_xres, _yres);
			// _w = GetOmega(error,_xres, _yres);
			// _w=0.666666;
			// _w = GetOmegaResidual(residual,_xres,_yres);
			_w = GetOmegaRough(residual,_xres,_yres);


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
				if(GetCenterCoef(x,y))
					resulttmp[i] = GetFromArray(result,x,y) - _w*all/GetCenterCoef(x,y) + _w*b[i]/GetCenterCoef(x,y);
				residual[i] = b[i]-all;
				sq_residual += residual[i]*residual[i];
			}
			// CheckWithResidual(error, residual, _xres, _yres);
			// SolveWithResidual(result, residual, _xres, _yres);

			//3. Record the result
			writeTrainData(error, residual, project(error, _xres, _yres), "./traindata2/wjacobi_data_rough.dat");
			// writeTrainDataRCk(residual, project(error, _xres, _yres), "./traindata2/wjacobi_data.dat");
			
			// Cout Brief Information
			std::cout<<"| "<<iter<<" | "<<_w<<" | "<<sq_residual<<" | "<<sq_error<<" | "<<std::endl;
			std::cout<<"------------"<<std::endl;
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