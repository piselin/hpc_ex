#include <cmath>
#include <cassert>
#include <vector>
#include <iostream>
#include <thread>
#include <mutex>
#include <iomanip>
#include "timer.hpp"

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
	// acquire the lock, it will be released onced the scope is left
	std::lock_guard<std::mutex> l (result.second);
 	result.first += sum;
}

int main()
{
	std::setprecision(18);
	typedef double bound_t;
	bound_t a = 1;
	bound_t b = 4;
	unsigned const bins = 1000;

	timer myt;

	// the function which we integrate
	auto f = [] (double x) { return std::sqrt(x)*std::log(x); };

	// figure out how many threads we have
	std::size_t const maxthreads = std::max(1u, std::thread::hardware_concurrency());

	// do experiment for each number of possible threads
	for (int nthreads = 1; nthreads <= maxthreads; ++nthreads)
	{
		// setup of the threads
	 	std::vector<std::thread> threads(nthreads);
	 	// protect the sum with a mutex
		std::pair<double, std::mutex> result;
		
		// figure out how to distribute the bins
		int bpt = bins / nthreads; // bins per thread
		int extrabins = bins%nthreads;
		double dx = (double)(b-a)/bins;

		// ------- exmeriment --------
		myt.start();
		double elapsed = 0.;
		unsigned nruns = 100; // 
		for (int k = 0; k < nruns; ++k)
		{
			// make sure that the result does not get added up
			result.first=0.;
			// reset the interval bounds
			bound_t lowerb = a;
			bound_t upperb = a;
	
			for (unsigned i = 0; i < nthreads-1; ++i)
			{
				upperb += bpt*dx;
				threads[i] = std::thread(integrate,f,std::ref(result),lowerb,upperb,bpt);
				lowerb = upperb;
			}

			upperb += (bpt+extrabins)*dx;
			threads[nthreads-1] = std::thread(integrate,f,std::ref(result),lowerb,upperb,bpt);
		  
		  	for (std::thread& t : threads)
		    	t.join();

		    myt.stop();
		    // ------ end of experiment ------
		    elapsed += myt.get_timing();
		}
		elapsed = elapsed / nruns;
		
	    std::cout << result.first << " " << nthreads << " " << elapsed << std::endl;
	}

	return 0;
}