#include <cplex.h>
#include <ctype.h>
#include <fstream>
#include <stdlib.h>
#include <string.h>
#include <iostream>

using namespace std;

//	Define global constants 
#define NAMELEN 50

class GenLP {
public:
	GenLP();
	GenLP(CPXENVptr the_env);
	virtual ~GenLP();

	double MIPScalarproduct(double* array1, double* array2, int beg, int end);
	/**
     * This function computes the scalar product of two vectors given by the elements of
     * two arrays with same beginning and ending index specified
     * @param array1 first array
     * @param array2 second array
     * @param beg beginning index for both arrays
     * @param end ending index for both arrays
     * @return returns the scalar product of the two arrays
     */

	int LPSetupLP(char* fname);
	/**
	 * Sets up the LP program: reads an MPS dual file
	 * @param filename name of the dual MPS file
	 * Changes the problem, adds columns e and f
	 * @return 0 if memory allocation is success, otherwise return a nonzero integer
	 */

	int LPSetreqdsoln(const double* reqdx);
	/**
	 * This method sets the required solution x_d obtaned from MIP
	 * @param reqdx the required solution from MIP
	 */

	int LPAddcutlp(const double* solnx);
	/**
	 * This method adds a cut d^k x^d >= d^k x^k or d^k x^d >= d^k x^k_sigma
	 * @param solnx the incumbent solution x^k/x^k_sigma from CutMIP/RegMIP
	 */

	int LPsolvelp();
	/**
	 * This method solves the Cost generator LP
	 * All relevant cuts need to be added by calling LPAddcutlp(double* solnx)
     */

	int LPSetcurrsolns();
	/**
	 * This method sets the solutions to y^k, e^k, and f^k after solving Cost generator LP
	 * Then d^k = c^k - e^k + f^k is set
	 */

	double* LPGetobjcoeffd();
	/**
	 * This method returns d^k  the current objective coefficient solution to the Cost generator LP
	 */

	void LPPrintobjd();
	/**
	 * This method printsthe incumbent solution d
	 */

private:
	CPXENVptr env;	 // CPLEX environment
	CPXLPptr  lp;    // Pointer to CPLEX LP obj
	int curr_ncols;   // Current LP number of columns (variables)
	int curr_nrows;	 // Current LP number of rows (constraints)
	int ini_nrows;
	int ini_ncols;
	int status, solnstat;
	double* objcoeffd; //Objective coefficient d^k to be passed to MIP
	double* incum_solnyef; //Solutions to y,e, and f
	double* reqd_soln;
	double* inicoeffc;
};
