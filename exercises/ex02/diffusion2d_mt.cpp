#include <iostream>
#include <vector>
#include <fstream>
#include <thread>
#include <mutex>
#include "timer.hpp"

using namespace std;

const unsigned N = 128; // the number of gridpoints per dimension (omega)
const double d = 1; // diffusion constant
const double dt = 0.001; // timestep
const double dx = (double)2/N; // distance between two grid points
const unsigned max_time = 1; // the amount of time we want to simulate

// this is needed for plotting
const void print_matrix(const std::vector<double>& u, std::ofstream& out)
{
	for(int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			out << u[i+N*j] << " ";
		}
		out << endl;
	}
}

void solve_eq(std::vector<double> &u, const unsigned nthreads)
{
	const double nsteps = max_time /  dt; // figure out how many time steps we have
	std::vector<double> temp = u;
	const double alpha = dt/(dx*dx)*d; // the coefficient is always the same

	std::vector<std::thread> threads(nthreads);


	int nslices = (N-2)/nthreads; // Fixme (pi) what if this division yields remainder?
	std::cout << "slices = " << (double) (N-2)/nthreads << std::endl;

	for (int k = 0; k < nsteps; ++k) // each loop is a time step
	{

		/*
		multi threading has to happen inside here, as in the explicit euler scheme
		each timestep is enterly dependent on the one before
		*/
		for (unsigned thr = 0; thr < nthreads; thr++)
		{
			threads[thr] = std::thread([&, thr]() {

				for (int i = thr*nslices+1; i < thr*nslices + nslices - 1; ++i)
				{
					//cout << thr << ": i = " << i << endl; 
					for (int j = 1; j < N-1; ++j)
					{
						// implementation of equation (7) from my report	
						temp[i+N*j] = alpha*(u[i+1+N*j]+u[i-1+N*j]+u[i+N*(j+1)]+u[i+N*(j-1)]-4*u[i+N*j])+u[i+N*j];
					}
				}
			});
		}
	
		for (std::thread& thr : threads)
			thr.join();
		u = temp; // one timestep is done, now we can owerwrite the old solution
	}
}

// the initial condition (t = 0) is given in equation (3). 
// the initial densitity is 1 on the square N/4 * N/4 and 0 everywhere else
void set_initial_condition(std::vector<double> &u)
{
	// fixme (pi) is there a better way with itterators?
	for(unsigned i = N/4; i < N*3/4; i++)
	{
		for(unsigned j = N/4; j < N*3/4; j++)
		{
			u[i+N*j] = 1;
		}
	}
}

int main(int argc, char *argv[])
{

	unsigned int nthreads = 1;
	if (argc > 1) 
		nthreads = atoi(argv[1]);

	//std::cout << "#Threads: " << nthreads << endl;

	std::vector<double> u(N*N, 0); //we have homogeneous dirichlet boundary conditions
	ofstream out("initial_mt.out");
	set_initial_condition(u);
	print_matrix(u,out);
	solve_eq(u, nthreads);
	ofstream out2("solution_mt.out");
	print_matrix(u, out2);
	return 0;
}