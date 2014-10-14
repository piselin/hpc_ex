#include <iostream>
#include <vector>
#include <fstream>
#include "timer.hpp"

typedef double value_type;
typedef std::size_t size_type;

class Diffusion2D {
public:
	Diffusion2D(const value_type D, 
				const value_type L, 
				const size_type N,
				const value_type dt)
	: D_(D), L_(L), N_(N), dt_(dt)
	{
		dx_ = L_ / (N_ - 1);
		alpha_ = dt_ *D_ / (dx_*dx_);

		// initialize vectors, we have homogeneous dirichlet bc so it's usefull to initialize everything with 0
		rho_.resize(N_*N_, 0.);
		rho_tmp.resize(N_*N_, 0.);

		initialize_density();
	}

	void advance() 
	{
		for (size_type i = 1; i < N_-1; ++i) // ignore the boundary, there we have dirchlet b.c.
		{
			for (size_type j = 1; j < N_-1; ++j)
			{
				// implementation of equation (7) from my report	
				rho_tmp[i*N_+j] = alpha_*(
					rho_[i*N_+1+j]
					+rho_[i*N_-1+j]
					+rho_[(i+1)*N_+j]
					+rho_[(i-1)*N_+j]
					-4*rho_[i*N_+j])
				+rho_[i*N_+j];
			}
		}
		// one timestep is done, now we can exchange the old and the new solution
		using std::swap;
		swap(rho_tmp,rho_);
	}

	
	// write the content of the vector into a file
	void write_density(std::string const& filename) const 
	{
		std::ofstream out_file(filename, std::ios::out);
		for(size_type i = 0; i < N_; ++i)
		{
			for (size_type j = 0; j < N_; ++j)
			{
				out_file << rho_[i*N_+j] << " ";
			}
			out_file << std::endl;
		}
		out_file.close();
	}

private:
	value_type D_; // diffusion constant
	value_type L_;
	value_type dx_; // distance between two grid points
	value_type dt_; // timestep
	value_type alpha_; // constant factor, see formula (6)
	size_type N_; // number of gridpoints per dimension (omega)

	std::vector<value_type> rho_; // current solution
	std::vector<value_type> rho_tmp; // temporary vector to store stuff 

	// the initial condition (t = 0) is given in equation (3). 
	// the initial densitity is 1 on the square N/4 * N/4 and 0 everywhere else
	// FIXME: attention, the rest of the vector outside of the inner square remains 0. This is not very portable to other problems.
	void initialize_density() 
	{
		for(size_type i = N_/4; i < N_*3/4; i++)
		{
			for(size_type j = N_/4; j < N_*3/4; j++)
			{
				rho_[i*N_+j] = 1;
			}
		}
	}
};

int main() 
{
	const value_type D = 1.0;
	const value_type L = 2.0;
	const size_type N = 1024;
	const value_type dt = 0.00000001;

	Diffusion2D system(D, L, N, dt);
	system.write_density("density.0.dat");

	const value_type tmax = 10000 * dt;
    value_type time = 0;
    
    timer t;
    
    t.start();
    while (time < tmax) {
        system.advance();
        time += dt;
    }
    t.stop();
    
    std::cout << "Timing : " << N << " " << t.get_timing() << std::endl;
    
    system.write_density("density_serial.dat");


	return 0;
}