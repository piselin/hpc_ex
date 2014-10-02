#include <cmath>
#include <cassert>
#include <vector>
#include <iostream>
#include <thread>
#include <mutex>

void integrate(
	std::function<double(double)> fun, 
	std::pair<double, std::mutex>& result,
 	double a, double b, unsigned n)
{	
	assert(n>0);
	double sum = 0.0;		// the sum
	double dx = (b-a)/n;	// interval double
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
	// accuire the lock, it will be released onced the scope is left
	std::lock_guard<std::mutex> l (result.second);
 	result.first += sum;
}

int main()
{
	typedef double bound_t;
	bound_t a = 1;
	bound_t b = 4;
	unsigned const bins = 50;

	// setup of the threads
	std::size_t const nthreads = std::max(1u, std::thread::hardware_concurrency());
 	std::vector<std::thread> threads(nthreads);

 	// protect the sum with a mutex
	std::pair<double, std::mutex> result;
	result.first=0.;

	// the function which we integrate
	auto f = [] (double x) { return std::sqrt(x)*std::log(x); };

	// figure out how to distribute the bins
	int bpt = bins / nthreads; // bins per thread
	int extrabins = bins%nthreads;
	
	double dx = (double)(b-a)/bins;
	
	bound_t lowerb = a;
	bound_t upperb = a;

	for (unsigned i = 0; i < nthreads-1; ++i)
	{
		upperb += bpt*dx;
		// std::cout << "loop " << i << std::endl;
		// std::cout << "integrate from " << lowerb << " to " << upperb << std::endl;
		threads[i] = std::thread(integrate,f,std::ref(result),lowerb,upperb,bpt);
		lowerb = upperb;
	}

	upperb += (bpt+extrabins)*dx;
	threads[nthreads-1] = std::thread(integrate,f,std::ref(result),lowerb,upperb,bpt);
  
  	for (std::thread& t : threads)
    	t.join();

    std::cout << result.first << std::endl;
	return 0;
}