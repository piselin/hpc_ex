#include <iostream>
#include <vector>
#include <fstream>

using namespace std;

const unsigned N = 256; // the number of gridpoints per dimension (omega)
const double d = 1; // diffusion constant
const double dt = 0.000001; // timestep
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

void solve_eq(std::vector<double> &u)
{
	const double nsteps = max_time /  dt; // figure out how many time steps we have
	std::vector<double> temp = u;
	const double alpha = dt/(dx*dx)*d; // the coefficient is always the same

	for (int k = 0; k < nsteps; ++k) // each loop is a time step
	{
		for (int i = 1; i < N-1; ++i) // ignore the boundary, there we have dirchlet b.c.
		{
			for (int j = 1; j < N-1; ++j)
			{
				// implementation of equation (7) from my report	
				temp[i+N*j] = alpha*(u[i+1+N*j]+u[i-1+N*j]+u[i+N*(j+1)]+u[i+N*(j-1)]-4*u[i+N*j])+u[i+N*j];
			}
		}
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

int main()
{
	// set up the matrix
	std::vector<double> u(N*N, 0); //we have homogeneous dirichlet boundary conditions
	ofstream out("initial.out");
	set_initial_condition(u);
	print_matrix(u,out);

	//print_matrix(u);
	//cout << u[72] << endl;
	solve_eq(u);
	ofstream out2("solution.out");
	print_matrix(u, out2);


// for (auto x:m)
	// 	std::cout << x << std::endl;

	return 0;
}