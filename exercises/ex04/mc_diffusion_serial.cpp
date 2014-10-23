#include <iostream>
#include <random>
#include <functional> // needed for std::bind
#include <cmath>

typedef double value_t; // computations are done on the real numbers

class Diffusion
{
public:
	Diffusion(const value_t x0, const value_t y0, const value_t d, const int N)
	: x0_(x0), y0_(y0), d_(d), x_(x0), y_(y0), N_(N)
	{}

	value_t mc_scheme()
	{
		for (unsigned i = 0; i < N_; ++i)
		{
			rnd_walk();
			rho += g(x_,y_);
		}
		rho = (double) rho/N_;
		return rho;
	}

	
private:
	const value_t x0_, y0_; // starting point
	value_t x_,y_; // current position
	value_t xc, yc;
	value_t rho = 0.0; // the solution
	const value_t d_; // step size
	const value_t seed_ = 42;
	const value_t omega_l = 0.0; // lower bound of domain
	const value_t omega_u = 1.0; // upper bound of domain
	int cnt = 0;
	const int N_; // number of sample points

	std::function<double()> real_rand = std::bind(
		std::uniform_real_distribution<double>(0,1), std::mt19937(seed_));

	// function on the boundary
	value_t g(value_t x, value_t y) 
	{
		return x;
	}

	// one random step
	void step() {
		value_t a = real_rand(); // get new random number
		x_ += d_*cos(2*M_PI*a);
		y_ += d_*sin(2*M_PI*a);
		cnt++;
	}
	// one full random walk starting from (x0,y0) until the domain is left
	void rnd_walk() {
		// make sure everything is ready for a new independent random walk
		bool inside_of_domain = true;
		// reset starting point
		x_ = x0_; 
		y_ = y0_; 
		// reset counter
		cnt = 0;
	
		// start the walk
		while(inside_of_domain)
		{
			// check if we are inside of domain omega
			if (x_ > omega_l && x_ < omega_u &&
				y_ > omega_l && y_ < omega_u) 
			{
				step();
			} else
			{
				//std::cout << "left domain at (" << x_ << "," <<y_<< ") after " << cnt << " steps" << std::endl;
				inside_of_domain = false;
			}		
		}
		// find the point (xc, yc) where we left the domain
		xc = x_;
		yc = y_;
	}
};

int main()
{
	//**************************************************************************
	// std::mt19937 mt_rand(42);
	// std::uniform_real_distribution<double> dis(0,1);
	// dis(mt_rand) // to generate rnd number

	// this is more c++11 style! to generate rnd number just call real_rand()
	//auto real_rand = std::bind(std::uniform_real_distribution<double>(0,1), std::mt19937(42));
	//***************************************************************************
	

	const value_t x0 = 0.3;
	const value_t y0 = 0.4;
	const value_t d = 0.01;
	const value_t s = 0.3;
	value_t r = 0.0;
	int N = 1;

	for (int i = 1; i < 19; ++i)
	{
		N *= 2;
		Diffusion system(x0, y0, d, N);
		r = system.mc_scheme();
		std::cout << N << " " << r << " " << std::abs(s-r) <<std::endl;
	}
}