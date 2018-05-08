#ifndef SOLVER_3D_H_
#define SOLVER_3D_H_

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include <iomanip>

#include "spectrum3d.h"
#include "function.h"

#define SQRTPI 1.77245385091
#define PI (3.14159265358979323846)
#define PRECISION 1e-6

#define PB_BOUND_Dirichlet	1
#define PB_BOUND_Neumann	2

class Solver3d{
public:
	Solver3d(){_xres=0; _yres=0; _zres=0; _btype=PB_BOUND_Dirichlet; _mask=NULL; _writeflag=false;};
	Solver3d(int xres, int yres, int zres, int btype):_xres(xres), _yres(yres), _zres(zres), _btype(btype), _mask(NULL),_writeflag(false){};
	
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
	int _zres;
	int _btype;
	int _writeflag;
	double _precision = PRECISION;
	int _maxiter = 100;
	int GetIndex(int x, int y, int z){
		if(x<0||x>_xres-1)	return -1;
		if(y<0||y>_yres-1)	return -1;
		if(z<0||z>_zres-1)	return -1;
		return z*_xres*_yres + y*_xres + x;
	}

	double GetFromArray(const double *arr, int x, int y, int z){
		int ind = GetIndex(x,y,z);
		return (ind>=0 && arr!=NULL)?arr[ind]:0;
	}

	double GetNeighbCoef(int x, int y, int z){
		if(x<0||x>_xres-1)	return 0;
		if(y<0||y>_yres-1)	return 0;
		if(z<0||z>_zres-1)	return 0;
		if(_mask)
			if(_mask[GetIndex(x,y,z)]==PB_BOUND_Dirichlet||_mask[GetIndex(x,y,z)]==PB_BOUND_Neumann)return 0;
		return -1;
	}

	double GetCenterCoef(int x, int y, int z){
		int center = -6;
		//suppose we use dirichlet boundary outside
		int xminusIndex = GetIndex(x-1, y, z);
		if(xminusIndex>=0 && _mask)	center = (_mask[xminusIndex]==PB_BOUND_Neumann)?center+1:center;
		else				center = (_btype==PB_BOUND_Neumann)?center+1:center;
		int xplusIndex = GetIndex(x+1, y, z);
		if(xplusIndex<=_xres-1 && _mask)	center = (_mask[xplusIndex]==PB_BOUND_Neumann)?center+1:center;
		else				center = (_btype==PB_BOUND_Neumann)?center+1:center;
		
		int yminusIndex = GetIndex(x, y-1, z);
		if(yminusIndex>=0 && _mask)	center = (_mask[yminusIndex]==PB_BOUND_Neumann)?center+1:center;
		else				center = (_btype==PB_BOUND_Neumann)?center+1:center;
		int yplusIndex = GetIndex(x, y+1, z);
		if(yplusIndex<=_yres-1 && _mask)	center = (_mask[yplusIndex]==PB_BOUND_Neumann)?center+1:center;
		else				center = (_btype==PB_BOUND_Neumann)?center+1:center;

		int zminusIndex = GetIndex(x, y, z-1);
		if(zminusIndex>=0 && _mask)	center = (_mask[zminusIndex]==PB_BOUND_Neumann)?center+1:center;
		else				center = (_btype==PB_BOUND_Neumann)?center+1:center;
		int zplusIndex = GetIndex(x, y, z+1);
		if(zplusIndex<=_yres-1 && _mask)	center = (_mask[zplusIndex]==PB_BOUND_Neumann)?center+1:center;
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
		for(int i=0;i<_xres*_yres*_zres;i++)
			fout.write((char*)&(field[i]), sizeof(double));
		fout.close();
		return 0;
	}
};


class WeightedJacobiSolver3d : public Solver3d{
public:
	WeightedJacobiSolver3d(int xres, int yres, int zres, int btype){_xres=xres; _yres=yres; _zres=zres; _btype=btype; _mask=NULL;};
	WeightedJacobiSolver3d(int xres, int yres, int zres, int btype, double w){_w=w; _xres=xres; _yres=yres;_zres=zres; _btype=btype; _mask=NULL;};
	void SetW(double w){_w=w;}

	int writeTrainData(const double *error, const double *residual, const double *ck, std::string filename){
		std::ofstream fout(filename.c_str(), std::ios::app | std::ios::binary);
		if(!fout){
			std::cout<<"WriteIter failed: not such file \'"<<filename<<"\'"<<std::endl;
			return 0;
		}
		for(int i=0; i < _xres * _yres * _zres; i++)
			fout.write((char*)&(error[i]), sizeof(double));
		for(int i=0; i < _xres * _yres * _zres; i++)
			fout.write((char*)&(residual[i]), sizeof(double));
		for(int i=0; i < _xres * _yres * _zres; i++)
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
		for(int i=0; i < _xres * _yres * _zres; i++)
			fout.write((char*)&(residual[i]), sizeof(double));
		for(int i=0; i < _xres * _yres * _zres; i++)
			fout.write((char*)&(ck[i]), sizeof(double));
		
		fout.close();

		return 0;
	}

	virtual double* Solve(double* b){

		int iter = 0;
		int size = _xres*_yres*_zres;
		double *result = new double[size];
		double *resulttmp = new double[size];

		double sq_residual = 100;
		double sq_error = 0;
		double *residual = new double[size];
		double *error = new double[size];
		
		//1. Set Initial Guess
		for(int z=0;z<_zres;z++)
		for(int y=0;y<_yres;y++)
		for(int x=0;x<_xres;x++){
			result[ GetIndex(x,y,z) ] = function3d_random(x,y,z),resulttmp[ GetIndex(x, y, z) ]=0;
			
			if(_mask)
				if(_mask[GetIndex(x,y,z)]==PB_BOUND_Dirichlet||_mask[GetIndex(x,y,z)]==PB_BOUND_Neumann)
					result[ GetIndex(x,y,z) ] =0;
		}

		while(iter<=_maxiter && sq_residual>_precision){
			sq_residual=0;
			sq_error =0;
			for(int i=0; i<size; i++){
				int x = i%_xres;
				int y = i%(_xres*_yres)/_xres;
				int z = i/(_xres*_yres);
				error[i] = function3d_allzero(x,y,z)-result[i];
				// error[i] = function1d_1wave(i)-result[i];
				sq_error += error[i]*error[i];
			}
			std::cout<<size <<'-'<<sq_error<<std::endl;

			// _w = GetOmegaNeigh(error,_xres, _yres);
			// _w = GetBestOmega(error,_xres, _yres);
			// _w = GetBestROmega(error,_xres, _yres);
			// _w = GetBestOmegaNeigh(error,_xres, _yres);
			// _w = GetBestROmegaNeigh(error,_xres, _yres);
			_w = GetOmega(error,_xres, _yres,_zres);
			// _w=0.666666;
			std::cout<<"get"<<std::endl;

			for(int i = 0;i<size;i++){
				int x = i%_xres;
				int y = i%(_xres*_yres)/_xres;
				int z = i/(_xres*_yres);
			
				double all = 
					GetNeighbCoef(x-1,y,z)*GetFromArray(result,x-1,y,z)+
					GetNeighbCoef(x+1,y,z)*GetFromArray(result,x+1,y,z)+
					GetNeighbCoef(x,y-1,z)*GetFromArray(result,x,y-1,z)+
					GetNeighbCoef(x,y+1,z)*GetFromArray(result,x,y+1,z)+
					GetNeighbCoef(x,y,z-1)*GetFromArray(result,x,y,z-1)+
					GetNeighbCoef(x,y,z+1)*GetFromArray(result,x,y,z+1)+
					GetCenterCoef(x,y,z)*GetFromArray(result,x,y,z);

				//x'=x-w*(Ax)*D_inv+w*b*D_inv
				if(GetCenterCoef(x,y,z))
					resulttmp[i] = GetFromArray(result,x,y,z) - _w*all/GetCenterCoef(x,y,z) + _w*b[i]/GetCenterCoef(x,y,z);
				residual[i] = b[i]-all;
				sq_residual += residual[i]*residual[i];
			}
			std::cout<<"iter"<<std::endl;
			// CheckWithResidual(error, residual, _xres, _yres);
			// SolveWithResidual(result, residual, _xres, _yres);

			//3. Record the result
			// writeTrainData(error, residual, project(error, _xres, _yres, _zres), "./traindata3/wjacobi_data_2div3.dat");
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