#include <iostream>
#include <math.h>
#include <iomanip>
#include <vector>
#include <stdlib.h>
#include <time.h>
#include <algorithm>
using namespace std;


// calculate the mean value of a vector
double VecMean(vector <double> x) {
	
	int n = x.size();
	double sum = 0;

	for (int i = 0; i < n; i++)
		sum += x[i];

	double mean = sum / n;
	
	return mean;
}

// N(0, 1) Density Function
double f(double x) {

	double pi = 4.0 * atan(1.0);
	
	return exp(-x * x * 0.5) / sqrt(2 * pi);
}

// Boole's Law
bool Bool(double Startpoint, double Endpoint, int n) {
	vector <double> X(n + 1, 0.0);
	vector <double> Y(n + 1, 0.0);
	double delta_x = (Endpoint - Startpoint) / double(n);
	double sum = 0;

	for (int i = 0; i <= n; i++) {
		X[i] = Startpoint + i * delta_x;
		Y[i] = f(X[i]);
	}

	for (int i = 0; i < (n - 1) / 4; i++) {
		int idx = 4 * i;
		sum += (1 / 45) * (14 * Y[idx] + 64 * Y[idx + 1] + 24 * Y[idx + 2] + 64 * Y[idx + 3] + 14 * Y[idx + 4] * delta_x);
	}

	return sum;
}

// N(0, 1) cdf by Boole's Rule
double N(double x) { return Bool(-10.0, x, 240); }

// Black Scholes Closed Form Solution
double BS(double S, double K, double v, double T, double r, double q, char PutCall) {
	
	double d1 = (log(S / K) + (r - q + 0.5 * v * v) * T) / v / sqrt(T);
	double d2 = d1 - v * sqrt(T);
	double Call = S * exp(-q * T) * N(d1) - K * exp(-r * T) * N(d2);

	if (PutCall == 'C')
		return Call;
	else
		return Call - S * exp(-q * T) + K * exp(-r * T);
}

void Run(double S, double K, double v, double T, double r, double q, char PutCall, int Nsims) {

	vector <double> ST(Nsims, 0.0);      // Initialise Terminal Prices
	vector <double> ST_K(Nsims, 0.0);    // Initialise Call Payoff
	vector <double> K_ST(Nsims, 0.0);    // Initialise Put Payoff

	double u1, u2, Z;
	double pi = 3.141592653589793;

	cout << "Running Simulation... " << endl;

	for (int i = 0; i < Nsims; i++) {
		// Independant Uniform RVs
		u1 = ((double)rand() / ((double)(RAND_MAX)+(double)(1)));
		u2 = ((double)rand() / ((double)(RAND_MAX)+(double)(1)));

		// Floor u1 to avoid error with log function
		u1 = max(u1, 1.0e-10);

		// Box-Muller Transformation
		Z = sqrt(-2.0 * log(u1)) * sin(2 * pi * u2);
		ST[i] = S * exp((r - q - 0.5 * v * v) * T + v * sqrt(T) * Z);
		ST_K[i] = max(ST[i] - K, 0.0);
		K_ST[i] = max(K - ST[i], 0.0);
	}

	// Discounted Simulated Price Averaged
	double BSCallSim = exp(-r * T) * VecMean(ST_K);
	double BSPutSim = exp(-r * T) * VecMean(K_ST);

	// Closed Form Solution
	double BSCall = BS(S, K, v, T, r, q, 'C');
	double BSPut = BS(S, K, v, T, r, q, 'P');

	// Error
	double CallError = BSCall - BSCallSim;
	double PutError = BSPut - BSPutSim;

	//cout << setprecision(4) << fixed;
	cout << "Using " << Nsims << " Simulations..." << endl;
	cout << " " << endl;
	cout << "Method CallPrice PutPrice " << endl;
	cout << "Simulation: " << BSCallSim << " " << BSPutSim << endl;
	cout << "Closed Form: " << BSCall << " " << BSPut << endl;
	cout << "Error: " << CallError << PutError << endl;
	cout << "-----------------------------------" << endl;
	cout << " " << endl;
	system("PAUSE");

}


int main()
{
	srand(time(0));            // Set random seed for PSRNG
	double S = 100.0;          // Spot Price
	double K = 100.0;          // Strike Price
	double T = 1;              // Maturity (Years)
	double r = 0.05;           // Interest Rate
	double q = 0;              // Dividend Yeild
	double v = 0.2;            // Volatility (%)
	int Nsims = 1e7;           // Number of Simulation Pathways

	Run(S, K, v, T, r, q, 'C', Nsims);

	return 0;
}
