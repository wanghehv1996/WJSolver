#ifndef _READ_H_
#define _READ_H_

#include<fstream>
#include<iostream>

double *readfile(std::string filename, int xres, int yres){
	std::ifstream fin(filename, std::ios::binary);
	
	double *mat = new double[xres*yres];

	double nNum;
	for(int i=0;i<xres*yres;i++){
		fin.read((char*)&nNum, sizeof(double));
		mat[i] = nNum;
	}
	return nNum;
}

#endif