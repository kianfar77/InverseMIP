

#include <cplex.h>
#include <ctype.h>
#include <fstream>
#include <stdlib.h>
#include <string.h>
#include <iostream>

using namespace std;

//	Define global constants 
#define NAMELEN 50
#define INFBOUND 1e+20

class ProbMIP {

public:
	ProbMIP();
	ProbMIP(CPXENVptr env);
	virtual ~ProbMIP();

	double MIPScalarproduct(const double* array1, const double* array2, int beg, int end);
	/**
	 * This function computes the scalar product of two vectors given by the elements of
	 * two arrays with same beginning and ending index specified
	 * @param array1 first array
	 * @param array2 second array
	 * @param beg beginning index for both arrays
	 * @param end ending index for both arrays
	 * @return returns the scalar product of the two arrays
	 */

	double MIPL1Norm(const double* array1, const double* array2, int len);
	/**
	 * This function computes the L1 norm ||array1 - array2|| given the lenght of the arrays
	 * @param array1 first array
	 * @param array2 second array
	 * @param len lenght of the arrays
	 * @return returns the L1 norm of the difference of two arrays
	 */

	int MIPSetupMIP(char* filename);
	/**
	 * Sets up the MIP program: reads an MPS file
	 * @param filename  name of the MPS file
	 * Set the problem up as a maximization problem
	 * Calls MIPMallocRegMIP() to allocate memory of RegMIP arrays
	 * @return 0 if memory allocation is success, otherwise return a nonzero integer
	 */

	int MIPWriteProblem(char* newname);
	/**
	 * Writes the MIP prblem object as an MPS file 
	 * @param filename  name of the MPS file written
	 */

	int MIPWriteDualRelaxed();
	/**
	 * Write the dual relaxed problem as filename_dual.mps
	 * The initial coefficient c must be set before calling this method
	 * @return 0 if write is success, otherwise return a nonzero integer
	 */

	int MIPMallocRegMIP();
	/**
	 * Dynamically allocate memory for required member variables
	 * mipsigma_nrows and mipsigma_ncols must be set after the regularized MIP problem is set
	 */

	int MIPSolveMIP();
	/**
	 * This method calls CPLEX to solve the MIP problem and sets
	 * the optimal objective value and current optimal solution
	 * @return 0 on success, otherwise returns a nonzero value
	 */
	int MIPSetupCutMIP(char* newfname);
	/**
	 * This method sets up the MIP Cut problem 
	 */

	int MIPSetupRegMIP();
	/**
	 * This method sets up the Regularized MIP problem
	 * Converts the MIP problem read via MIPSetupMIP()
	 * MIP problem must be read using MIPSetupMIP(char* filename)
	 * Initial objective coefficient must be set before solving
	 * the regularized MIP by calling MIPSetupObjcoeffsigma()
	 * Required solution x_d should be obtained
	 * Adds the necessary additional columns and rows to get the regularized problem
	 */

	int MIPSetSigma(const double* solnx);
	/**
	 * Set/Update the value of sigma for the regularized problem
	 * Must be called before MIP Reg is solved
	 * @param solnx The incumbent solution x^k to cut mip 
	 */

	int MIPSetupReqdSolnandIniObj(int seed);
	/**
	 * This method gets the required solution x^d and initial cost vector c
	 * @param seed srand() value to be used
	 * Needs to be called before dual is written using MIPWriteDualRelaxed()
	 */

	int MIPGetIniObj(const int* arraylb, const int* arrayub, int seed);
	/**
	 * This method gets the initial cost vector c, having set the value of x^d
	 * @param arraylb the lower boud on coefficient values
	 * @param arraylb the upper boud on coefficient values
	 * @param seed srand() value to be used
	 * Needs to be called before dual is written using MIPWriteDualRelaxed()
	 */

	int MIPSetCurrobjsigma(const double* objd);
	/**
	 * This method setsup the current objective coefficients for the regularized MIP
	 * when an objective array is passed as an argument
	 * @param objd The current objective coefficient obtained from the Cost Generator LP 
	 * mipsigma_ncols must be set before calling this method
	 * Value of sigma must be set/updated if necessary using MIPSetSigma()
	 * obj coeff = (d^k - sigma) e - (d^k + sigma) f
	 * Values are obtained from the cost generator LP
	 */

	int MIPSetCurrobjsigma();
	/**
	 * This method setsup the current objective coefficients for the regularized MIP
	 * when no objective array is passed as an argument
	 * The incumbent objective is used as to get the objective
	 * mipsigma_ncols must be set before calling this method
	 * Value of sigma must be set/updated if necessary using MIPSetSigma()
	 * obj coeff = (d^k - sigma) e - (d^k + sigma) f
	 * Values are obtained from the cost generator LP
	 */

	int MIPSetCurrobj(const double* objd);
	/**
	 * Gets the objective coefficient d^k for the upcoming iteration
	 * The GenLP should provide the current objective coefficient except for the initial iteration
	 * For Regularized MIP MIPSetupObjcoeffsigma() must be called
	 * after this method to set the objective coefficients for regularized MIP
	 */

	double MIPRegsetIncumsolnxandObjval();
	/**
	 * Sets the incumbant solution after solving Regularized MIP
	 * x^k = x^d + e^k - f^k; objval = objval_regMIP + d^k x^d
	 * Return diff = d^k x^k - sigma ||x^k -x ^d|| - d^k x^d = objval_regMIP
	 */

	int MIPCutsetIncumsolnxandObjval();
	/**
	 * Sets the incumbant solution after solving Cut MIP
	 */

	int MIPGetSoln();
	/**
	 * Sets and prints the temp solution 
	 */

	int MIPSetCutSoln();
	/**
	 * This method sets the current solution x^k to the incumbant solution array
	 */

	int MIPSetRegSoln();
	/**
	 * This method sets the curent solutions y^k, e^k, and f^k
	 * Then the incumbent array is set ti x^k_sigma = x^d + e^k - f^k
	 */

	double* MIPGetIncumsolnx();
	/**
	 * Returns the incumbent solution set after solving Regularized MIP
	 */

	static double* MIPGetreqdsoln();
	/**
	 * Returns the required solution x^d
	*/

	static double* MIPGetiniobjcoeff();
	/**
	 * Returns the initial objective coefficient c
	*/

	double MIPChecktermcondn();
	/**
	 * Compare the obj value of MIPCut/MipReg with c^k x^d
	 * Return the difference
	 */

	static int MIPGetncols();
	/**
	 * Returns mip_ncols
	 */

private:
	static char* filename;
	static char* dualname;
	static double* reqd_soln;
	static double* ini_objcoeff;
	static double* orig_objcoeff; //Cost vector in the file read
	static int mip_ncols;   // Current MIP number of columns (variables)
	static int mip_nrows;	 // Current MIP number of rows (constraints)

	CPXENVptr env;	 // CPLEX environment
	CPXLPptr  mip;    // Pointer to CPLEX MIP obj
	int status;
	int solnstat;
	double sigma;
	double objval, objshift;
	double* curr_objcoeff;
	double* incum_solnx; //Solution x^k to be passed to cost generator LP
	double* opt_soln;
	int mipsigma_nrows;
	int mipsigma_ncols;
	int mipsigma_ncolswithx;
	double* curr_objcoeffsigma;
	double* incum_solnxsigma;
};