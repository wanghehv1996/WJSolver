#include "Solver2d.h"
// #include "Solver1d.h"
#include <iostream>
#include <fstream>
#include<time.h>

double initRand(){
	srand((unsigned)time(NULL));
}


int main(){

	int xres = 32;//32;//32;
	int yres = 32;//32;
	/*
	 * Try diff solver
	 */
	// // Solver* solver = new ICPCGSolver(xres, yres, PB_BOUND_Dirichlet);
	// Solver* solver = new CGSolver(xres, yres, PB_BOUND_Dirichlet);
	// Solver* solver = new JacobiSolver(xres, yres, PB_BOUND_Dirichlet);
	// Solver* solver = new GaussSeidelSolver(xres, yres, PB_BOUND_Dirichlet);
	// Solver* solver = new SORSolver(xres, yres, PB_BOUND_Dirichlet, 1.25);
	
	// solver->SetWriteFlag(true);
	// solver->SetPrecision(1e-10);
	// solver->SetMaxIter(100);
	// double* force = new double [xres*yres];
	// for(int y = 0; y < yres; y++)
	// 	for(int x = 0; x < xres; x++){
	// 		force[y*xres+x] = divfunction(x,y);
	// 		//std::cout<<force[y*xres+x];
	// 	}
	// double *result = solver->Solve(force);

	// std::ofstream ofs;
	// ofs.open("result.txt");
	// for(int y = 0; y < yres; y++){
	// 	for(int x = 0; x < xres; x++){
	// 		ofs<<result[y*xres+x]<<' ';
	// 	}
	// 	ofs<<std::endl;
	// }
	// ofs.close();
	// Py_Init();
	// Py_ShowMat(result, xres, yres);
	// Py_End();

	/*
	 *
	 *
	initRand();
	system("rm ./traindata/wjacobi_data.dat");
	WeightedJacobiSolver1d* solver1d = new WeightedJacobiSolver1d(xres, PB_BOUND_Dirichlet, 1);
	solver1d->SetPrecision(1e-10);
	solver1d->SetMaxIter(100);//(49);
	double* force = new double [xres];
	for(int x = 0; x < xres; x++){
		// force[x] = divfunction1d_random(x);
		// force[x] = divfunction1d_1wave(x);
		force[x] = divfunction1d_allzero(x);
		// force[x] = divfunction1d_3wave(x);
	}
	solver1d->Solve(force);
	 *
	 *
	 */

	initRand();
	WeightedJacobiSolver2d* solver2d = new WeightedJacobiSolver2d(xres, yres, PB_BOUND_Dirichlet, 1);
	solver2d->SetPrecision(1e-10);
	solver2d->SetMaxIter(100);//(49);

	double* force = new double [xres * yres];
	for(int y = 0; y < yres; y++)
	for(int x = 0; x < xres; x++){
		// force[y*xres + x] = divfunction2d_random(x,y);
		force[y*xres + x] = divfunction2d_allzero(x,y);
	}
	for(int i=0;i<100;i++)
		solver2d->Solve(force);

	/*
	 * Generate data
	 */
	// WeightedJacobiSolver1d* solver1d = new WeightedJacobiSolver1d(xres, PB_BOUND_Dirichlet, 1);
	// solver1d->SetPrecision(1e-10);
	// solver1d->SetMaxIter(31);//(49);
	// initRand();
	// solver1d->solven=1;
	// double* force = new double [xres];
	// for(int x = 0; x < xres; x++){
	// 	force[x] = divfunction1d_allzero(x);
	// }
	// for(int i =0;i<200;i++){
	// 	solver1d->Solve(force);
		
	// }
}

