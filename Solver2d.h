#ifndef SOLVER_2D_H_
#define SOLVER_2D_H_

#include "Solver.h"
#include "function.h"
#include "spectrum2d.h"

class WeightedJacobiSolver2d : public Solver{
public:
	WeightedJacobiSolver2d(int xres, int yres, int btype){_xres=xres; _yres=yres; _btype=btype; _mask=NULL;};
	WeightedJacobiSolver2d(int xres, int yres, int btype, double w){_w=w; _xres=xres; _yres=yres; _btype=btype; _mask=NULL;};
	void SetW(double w){_w=w;}

	double GetNeighbCoefDot(int x, int y){
		if(x<0||x>_xres-1)	return 0;
		if(y<0||y>_yres-1)	return 0;
		if(_mask)
			if(_mask[GetIndex(x,y)]==PB_BOUND_Dirichlet||_mask[GetIndex(x,y)]==PB_BOUND_Neumann)return 0;

		return -1;
	}

	double GetCenterCoefDot(int x, int y){
		int center = -4;
		if(_mask)
			if(_mask[GetIndex(x,y)]==PB_BOUND_Dirichlet||_mask[GetIndex(x,y)])
				return 1;
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
			
			// result[i]=function1d_low(i),resulttmp[i]=0;
			// result[i]=function1d_wave3(i),resulttmp[i]=0;
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
			}

			// _w = GetOmegaNeigh(error,_xres, _yres);
			// _w = GetBestOmega(error,_xres, _yres);
			// _w = GetBestROmega(error,_xres, _yres);
			// _w = GetBestOmegaNeigh(error,_xres, _yres);
			// _w = GetBestROmegaNeigh(error,_xres, _yres);
			_w = GetOmega(error,_xres, _yres);
			// _w=0.666666;


			for(int i = 0;i<size;i++){
				int x = i%_xres;
				int y = i/_xres;
			
				double all = 
					GetNeighbCoefDot(x-1,y)*GetFromArray(result,x-1,y)+
					GetNeighbCoefDot(x+1,y)*GetFromArray(result,x+1,y)+
					GetNeighbCoefDot(x,y-1)*GetFromArray(result,x,y-1)+
					GetNeighbCoefDot(x,y+1)*GetFromArray(result,x,y+1)+
					GetCenterCoefDot(x,y)*GetFromArray(result,x,y);

				//x'=x-w*(Ax)*D_inv+w*b*D_inv
				if(GetCenterCoefDot(x,y))
					resulttmp[i] = GetFromArray(result,x,y) - _w*all/GetCenterCoefDot(x,y) + _w*b[i]/GetCenterCoefDot(x,y);
				residual[i] = b[i]-all;
				sq_residual += residual[i]*residual[i];
			}
			// CheckWithResidual(error, residual, _xres, _yres);
			// SolveWithResidual(result, residual, _xres, _yres);

			//3. Record the result
			// writeTrainData(error, residual, project(error, _xres, _yres), "./traindata2/wjacobi_data_64best.dat");
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