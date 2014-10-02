#include <iostream>
#include <cmath>
#include <cassert>
#include <functional>


double func(double x) {
	return std::sqrt(x)*std::log(x);
}

// maybe use function object or template instead
inline double integrate(std::function<double(double)> fun, double a, double b, unsigned n)
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
	unsigned bins = 100;

	std::cout << integrate(func, a,b,bins) << std::endl;
}