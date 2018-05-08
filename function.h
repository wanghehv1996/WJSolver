#ifndef FUNCTION_H_
#define FUNCTION_H_
#include <math.h>
#include <stdlib.h>
//------------------------------
// 1 DIMENSION
//------------------------------ 

//the divergence function of 1 wave, whose cycle is 1.5*2=3
//because we use the 1.5 as the k, your grids resolution should be divisible with 3 (here we use 33)
double divfunction1d_1wave(int x){
	double cycle = 1.5;
	double maxres = 33;
	//z = sin(x/cycle)
	double xx = x+1;
	// Here instead of using f''(x), we use exactly f(x+1)+f(x-1)-2f(x). That will be more accurate
	double div = sin((xx-1.)*PI/cycle)+sin((xx+1.)*PI/cycle)-2*sin(xx*PI/cycle);

	// if(x==0)
	// 	div -= sin(-1.*PI/cycle);
	// else if(x==maxres-1)
	// 	div -= sin(maxres*PI/cycle);
	return -div;
}

//the divergence function of all zero
double divfunction1d_allzero(int x){
	return 0;
}

//the divergence function of 3 wave
double divfunction1d_3wave(int x){
	double cycle = 2;//k=N/2
	double cycle1 = 1.5;//k=N/1.5
	// double cycle2 = 8;//k=N/8
	double amp = 1;
	double amp1 = 2;
	// double amp2 = 0.5;

	double maxres = 36;
	//z = sin(x/cycle)
	double xx = x;
	// Here instead of using f''(x), we use exactly f(x+1)+f(x-1)-2f(x). That will be more accurate
	double div = amp*( sin((xx-1.)*PI/cycle) + sin((xx+1.)*PI/cycle) - 2*sin(xx*PI/cycle) );
	div += amp1*( sin((xx-1.)*PI/cycle1) + sin((xx+1.)*PI/cycle1) - 2*sin(xx*PI/cycle1) );
	// div += amp2*( sin((xx-1.)*PI/cycle2) + sin((xx+1.)*PI/cycle2) - 2*sin(xx*PI/cycle2) );

	if(x==0)
		div -= amp * sin(-1.*PI/cycle) + amp1 * sin(-1.*PI/cycle1)/* + amp2 * sin(-1.*PI/cycle2)*/;
	else if(x==maxres-1)
		div -= amp * sin(maxres*PI/cycle) + amp1 * sin(maxres*PI/cycle1)/* + amp2 * sin(maxres*PI/cycle2)*/;
	return -div;
}

//the divergence function of 1 wave, whose cycle is 1.5*2=3
double function1d_1wave(int x){
	double cycle1 = 1;
	double cycle2 = 2;
	double cycle3 = 1.5;
	//z = sin(x/cycle)
	double xx = x+1;
	// return sin(xx*PI/cycle1)+sin(xx*PI/cycle2)+sin(xx*PI/cycle3);
	return sin(xx*PI/cycle3);
}

//the divergence function of all zero
double function1d_allzero(int x){
	return 0;
}

//the divergence function of 3 wave
double function1d_2wave(int x){
	double cycle = 2;
	double cycle1 = 1.5;
	double amp = 1;
	double amp1 = 2;
	double xx = x+1;
	return amp * sin(xx*PI/cycle) + amp1 * sin(xx*PI/cycle1);
}

double function1d_random(int x){
	double random =  1.0*rand()/double(RAND_MAX);
	return random*2.0-1.0;
	// double xx = x+1;
	// return (xx-10)*(xx-10)+100.0/xx+log(xx)*10+10;
}

double function1d_complex(int x){
	double xx = x+1;
	return (xx-10)*(xx-10)+100.0/xx+log(xx)*10+10;
}

double function1d_low(int x){
	// 1
	double cycle = 1.8;
	double xx = x;
	return sin((xx+2)*PI/cycle);

	// 2
	// double xres = 32;
	// double k1 = 4, k2 = 28;
	// double phase1 = 2, phase2 = 8;
	// double amp1 = 1, amp2 = 1;
	// double xx = x;
	// return amp1 * sin((xx+phase1)*k1*PI/xres) + amp2 * sin((xx+phase2)*k2*PI/xres);

	// 3
	// double xres = 32;
	// double k1 = 6, k2 = 12;
	// double phase1 = 4, phase2 = 2;
	// double amp1 = 1, amp2 = 1;
	// double xx = x;
	// return amp1 * sin((xx+phase1)*k1*PI/xres) + amp2 * sin((xx+phase2)*k2*PI/xres);

	// 4
	// double xres = 32;
	// double k1 = 20, k2 = 30;
	// double phase1 = 16, phase2 = 9;
	// double amp1 = 1, amp2 = 1;
	// double xx = x;
	// return amp1 * sin((xx+phase1)*k1*PI/xres) + amp2 * sin((xx+phase2)*k2*PI/xres);
}

double function1d_wave3(int x){
	double xx = x+1;
	return sin((xx)*PI/3)+sin((xx)*PI/1.5);
}

double function1d_minmax(int x){
	double xx = x+1;
	return sin((xx)*PI/33)+2*sin((32*xx)*PI/33);
}

//------------------------------
// 2 DIMENSION
//------------------------------ 
double divfunction(int x,int y){
	double cycle = 8;
	double maxres = 32;
	//z = sin(x/4)*sin(y/4)
	double xx = x;
	double yy = y;
	double div = -2.*(PI/cycle)*(PI/cycle)*sin(xx*PI/cycle)*sin(yy*PI/cycle);
	if(x==0)
		div -= sin(-1.*PI/cycle)*sin(yy*PI/cycle);
	else if(x==maxres-1)
		div -= sin(maxres*PI/cycle)*sin(yy*PI/cycle);

	if(y==0)
		div -= sin(-1.*PI/cycle)*sin(xx*PI/cycle);
	else if(y==maxres-1)
		div -= sin(maxres*PI/cycle)*sin(xx*PI/cycle);
	return -div;
}


double complexdivfunction(int x, int y){
	double maxres = 32;
	double cycle = 4;
	double xx = x;
	double yy = y;
	double div = -2.*(PI/cycle)*(PI/cycle)*sin(xx*PI/cycle)*sin(yy*PI/cycle);
	if(x==0)
		div -= sin(-1.*PI/cycle)*sin(yy*PI/cycle);
	else if(x==maxres-1)
		div -= sin(maxres*PI/cycle)*sin(yy*PI/cycle);

	if(y==0)
		div -= sin(-1.*PI/cycle)*sin(xx*PI/cycle);
	else if(y==maxres-1)
		div -= sin(maxres*PI/cycle)*sin(xx*PI/cycle);


	return -div;
}

double divfunction2d_allzero(int x, int y){
	return 0;
}

double function2d_allzero(int x, int y){
	return 0;
}

double function2d_waves(int x, int y){
	double result = 0;
	result += 2*sin((1.+x)*2*PI/33)*sin((1.+y)*2*PI/33);
	result += 2*sin((1.+x)*5*PI/33)*sin((1.+y)*6*PI/33);
	result += 5*sin((1.+x)*20*PI/33)*sin((1.+y)*20*PI/33);
	return result;
}

double function2d_random(int x, int y){
	double random =  1.0*rand()/double(RAND_MAX);
	return random*2.0-1.0;
}


//------------------------------
// 3 DIMENSION
//------------------------------ 
double divfunction3d_allzero(int x, int y, int z){
	return 0;
}

double function3d_allzero(int x, int y, int z){
	return 0;
}

double function3d_random(int x, int y, int z){
	double random =  1.0*rand()/double(RAND_MAX);
	return random*2.0-1.0;
}
#endif