#include <iostream>
#include <cmath>
#include <cassert>

// better function object (with some inline maybe?), pass function by reference or use template
double fun(double x) {
	return sqrt(x)*log(x);
}

double integrate(double a, double b, unsigned n)
{	
	assert(n>0);
	double S = 0.0;			// the sum
	double dx = (b-a)/n;	// interval size
	double x = 0.0;			// x* in the formula
	double x1 = a;			// x_i
	double x2 = a+dx;		// x_(i+1)
	for (int i = 0; i < n; ++i)
	{
		x = (x1+x2) / 2;
		S += fun(x)*dx;
		x1 = x2;
		x2 += dx;
	}
	return S;
} 

int main()
{
	typedef double bound_t;

	bound_t a = 1;
	bound_t b = 4;

	std::cout << integrate(a,b,1000) << std::endl;
}