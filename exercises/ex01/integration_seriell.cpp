#include <iostream>
#include <cmath>
#include <cassert>
#include <functional>

inline double integrate(std::function<double(double)> fun, double a, double b, unsigned n)
{	
	assert(n>0);
	double sum = 0.0;			// the sum
	double dx = (b-a)/n;	// interval size
	double x = 0.0;			// x* in the formula
	double x1 = a;			// x_i
	double x2 = a+dx;		// x_(i+1)
	for (int i = 0; i < n; ++i)
	{
		x = (x1+x2) / 2;
		sum += fun(x)*dx;
		x1 = x2;
		x2 += dx;
	}
	return sum;
} 

int main()
{
	typedef double bound_t;
	bound_t a = 1;
	bound_t b = 4;
	unsigned bins = 1000;

	// lambda function
	auto f = [] (double x) { return std::sqrt(x)*std::log(x); };

	std::cout << integrate(f,a,b,bins) << std::endl;
}