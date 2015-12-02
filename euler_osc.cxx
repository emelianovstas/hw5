#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>

using namespace std;

void Eulerforward (const int N, double* x, const double dt);
void Eulerbackward (const int N, double* y, const double dt);
void OUT (const int N, double* x,double* y, const double dt);

int main() {
	const double dt= M_PI/10;
	const int N = 20*M_PI/dt;
	double* x= new double [2*N];
	double* y= new double [2*N];

x[0]=1;
x[N]=0;
y[0]=1;
y[N]=0;

	Eulerforward(N,x,dt);
	Eulerbackward(N,y,dt);
	OUT(N,x,y,dt);

delete[] x;
delete[] y;

	return 0;
}

void Eulerforward (const int N, double* x, const double dt)
{
	for (int i=0; i<(N-1); i++)
		{
			x[i+1]=x[i]+dt*x[N+i];
			x[i+1+N]=x[i+N]-dt*x[i];
		}
}

void Eulerbackward (const int N, double* y, const double dt)
{

	for (int i=0; i<(N-1); i++)
		{
			y[i+1]=(1/(dt*dt+1))*y[i]+dt*y[N+i];
			y[i+1+N]=(1/(dt*dt+1))*y[i+N]-dt*y[i];
		}
}

void OUT (const int N, double* x, double* y, const double dt)
{
ofstream output("output.txt");
	for (int i=0; i<(N-1); i++)
		{
			output<<i*dt<<"\t"<<x[i]<<"\t"<<y[i]<<endl;
		}
	output.close();
}
