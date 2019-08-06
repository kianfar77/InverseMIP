/* --------------------------------------------------------------------------
* Implementation of Cutting Plane Algorithm and Regularized Cutting Plane algorithm 
* to solve Inverse Mixed Integer Programs
*
* Author: Vishnu Vijayaraghavan
* Date: May 10, 2018
* Revised: Jun 25, 2019
* Version 12.9.1
* 
* --------------------------------------------------------------------------
*
*  Main.cpp - Reading in and optimizing a problem 
*  The MIP problem to be used is specified as command line argument during execution 
*
*  #define SOLVE_MIPCUT to solve Cutting Plane Algorithm
*  #define SOLVE_MIPREG to solve Regularized Cutting Plae Algorithm 
*
*  The seed to be used to generate random instances of the Inverse MIP problem
*  shoule be set as the value of the variable seed
* --------------------------------------------------------------------------*/

#include <cplex.h>
#include <ctype.h>
#include<ctime>
#include <fstream>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include "GenLP.h"
#include "ProbMIP.h"

using namespace std;

#define NONZERO_LB 1e-6	    // Accuracy level for defining zero
#define NAMELEN 50

//#define SOLVE_MIPCUT
#define SOLVE_MIPREG

int
main(int argc, char* argv[]) {
	CPXENVptr env = NULL;
	int status = 0, seed = 96;
	int mip_numcols;
	double objshift;
	double* solnx, * solnd, * reqdx, * inicoeff;
	char* filename, * newfname, * dualfname;

	filename = new char[NAMELEN];
	newfname = new char[NAMELEN];
	dualfname = new char[NAMELEN];
	
	strcpy_s(filename, NAMELEN, argv[1]);
	strcpy_s(newfname, NAMELEN, filename); strcat_s(newfname, NAMELEN, "_new");
	strcpy_s(dualfname, NAMELEN, filename); strcat_s(dualfname, NAMELEN, "_dual");

	/* Initialize the CPLEX environment */
	/*This routine returns a pointer to a CPLEX environment.*/
	env = CPXopenCPLEX(&status);
	if (env == NULL) {
		char errmsg[1024];
		std::cerr << endl << "Algorithm: Could not open CPLEX environment." << endl;
		CPXgeterrorstring(env, status, errmsg);
		std::cerr << errmsg << endl;
		std::cerr << "Terminating ..." << endl;
		exit(status);
	}
	// Turn on/off output to the screen
	status = CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_OFF);
	if (status != 0) {
		std::cerr << endl << "CPLEX Error: " << status << endl;
		std::cerr << "Algorithm: Failure to turn on screen indicator!" << endl;
		std::cerr << "Terminating ..." << endl;
		exit(status);
	}
	// Set PRESOLVE either ON or OFF
	status = CPXsetintparam(env, CPX_PARAM_PREIND, CPX_OFF);
	if (status) {
		std::cerr << endl << "CPLEX Error: " << status << endl;
		std::cerr << "Failure to turn on/off CPLEX presolve" << endl;
		std::cerr << "Terminating ..." << endl;
		exit(status);
	}

  // Initialize the LP/MIP objects (always after setting the CPLEX environment!)
	ProbMIP CutMIP(env);
	ProbMIP RegMIP(env);
	GenLP CostLP(env);

	status = RegMIP.MIPSetupMIP(filename);
	if (status != 0) {
		std::cerr << endl << "CPLEX Error: " << status << endl;
		std::cerr << "Algorithm: Failure to setup MIP CUT!" << endl;
		std::cerr << "Terminating ..." << endl;
		exit(status);
	}

	status = RegMIP.MIPSetupReqdSolnandIniObj(seed);
	if (status != 0) {
		std::cerr << endl << "CPLEX Error: " << status << endl;
		std::cerr << "Algorithm: Failure to get required solution and initial objective coefficients for MIP CUT!" << endl;
		std::cerr << "Terminating ..." << endl;
		exit(status);
	}

	status = RegMIP.MIPWriteProblem(newfname);
	if (status != 0) {
		std::cerr << endl << "CPLEX Error: " << status << endl;
		std::cerr << "Algorithm: Failure to write the MIP problem!" << endl;
		std::cerr << "Terminating ..." << endl;
		exit(status);
	}

	status = RegMIP.MIPWriteDualRelaxed();
	if (status != 0) {
		std::cerr << endl << "CPLEX Error: " << status << endl;
		std::cerr << "Algorithm: Failure to write the dual LP relaxation problem!" << endl;
		std::cerr << "Terminating ..." << endl;
		exit(status);
	}

	status = RegMIP.MIPSetupRegMIP();
	if (status != 0) {
		std::cerr << endl << "CPLEX Error: " << status << endl;
		std::cerr << "Algorithm: Failure to setup the Regukarized MIP problem!" << endl;
		std::cerr << "Terminating ..." << endl;
		exit(status);
	}

	//status = RegMIP.MIPWritetempProblem("fileorigreg");

	//Setting up the Reg MIP problem
	status = RegMIP.MIPMallocRegMIP();
	if (status != 0) {
		std::cerr << endl << "CPLEX Error: " << status << endl;
		std::cerr << "Algorithm: Failure to allocate memory for MIP CUT arrays!" << endl;
		std::cerr << "Terminating ..." << endl;
		exit(status);
	}

	//Setting up the  CUT MIP problem
	status = CutMIP.MIPSetupCutMIP(newfname);
	if (status != 0) {
		std::cerr << endl << "CPLEX Error: " << status << endl;
		std::cerr << "Algorithm: Failure to set up the MIP CUT problem!" << endl;
		std::cerr << "Terminating ..." << endl;
		exit(status);
	}

	//status = CutMIP.MIPWritetempProblem("tempcut");

	//Setting up the LP Cost Generator problem
	status = CostLP.LPSetupLP(dualfname);
	if (status != 0) {
		std::cerr << endl << "CPLEX Error: " << status << endl;
		std::cerr << "Algorithm: Failure to set up the MIP CUT problem!" << endl;
		std::cerr << "Terminating ..." << endl;
		exit(status);
	}

	mip_numcols = ProbMIP::MIPGetncols();

	solnx = new double[mip_numcols];
	solnd = new double[mip_numcols];
	inicoeff = new double[mip_numcols];

	reqdx = ProbMIP::MIPGetreqdsoln();
	
	std::cout << "reqdx " << endl;
	for (int j = 0; j < mip_numcols; j++)
		cout << reqdx[j] << " ";
	std::cout << endl;

	status = CostLP.LPSetreqdsoln(reqdx);
	if (status != 0) {
		std::cerr << endl << "CPLEX Error: " << status << endl;
		std::cerr << "Algorithm: Failure to set required solution for LP!" << endl;
		std::cerr << "Terminating ..." << endl;
		exit(status);
	}

	inicoeff = ProbMIP::MIPGetiniobjcoeff();

	clock_t start, end;
	int iter = 0;
	double difference;
	 
	start = clock();

 #ifdef SOLVE_MIPCUT
	//Start algorithm
	cout << "**************MIP CUT algorithm starts*****************" << endl;
	do {

		iter++;
		status = CutMIP.MIPSolveMIP();
		if (status != 0) {
			std::cerr << endl << "CPLEX Error: " << status << endl;
			std::cerr << "Algorithm: Failure to solve MIP CUT problem!" << endl;
			std::cerr << "Iteration: " << iter << endl;
			std::cerr << "Terminating ..." << endl;
			exit(status);
		}

		status = CutMIP.MIPCutsetIncumsolnxandObjval();
		if (status != 0) {
			std::cerr << endl << "CPLEX Error: " << status << endl;
			std::cerr << "Algorithm: Failure to set incumbent solution for MIP CUT problem!" << endl;
			std::cerr << "Terminating ..." << endl;
			exit(status);
		}

		difference = CutMIP.MIPChecktermcondn();
		cout << "difference - MIPCUT: " << difference << endl;

		solnx = CutMIP.MIPGetIncumsolnx();

		status = CostLP.LPAddcutlp(solnx);
		if (status != 0) {
			std::cerr << endl << "CPLEX Error: " << status << endl;
			std::cerr << "Algorithm: Failure to add cut to Cost Generator LP!" << endl;
			std::cerr << "Terminating ..." << endl;
			exit(status);
		}

		status = CostLP.LPsolvelp();
		if (status != 0) {
			std::cerr << endl << "CPLEX Error: " << status << endl;
			std::cerr << "Algorithm: Failure to solve Cost Generator LP problem!" << endl;
			std::cerr << "Iteration: " << iter << endl;
			std::cerr << "Terminating ..." << endl;
			exit(status);
		}

		status = CostLP.LPSetcurrsolns();
		if (status != 0) {
			std::cerr << endl << "CPLEX Error: " << status << endl;
			std::cerr << "Algorithm: Failure to set solutions to Cost Generator LP!" << endl;
			std::cerr << "Terminating ..." << endl;
			exit(status);
		}

		solnd = CostLP.LPGetobjcoeffd();

		status = CutMIP.MIPSetCurrobj(solnd);
		if (status != 0) {
			std::cerr << endl << "CPLEX Error: " << status << endl;
			std::cerr << "Algorithm: Failure to set current objective for MIP Cut!" << endl;
			std::cerr << "Terminating ..." << endl;
			exit(status);
		}
	cout << "************************************" << endl;
	} while (difference > NONZERO_LB);
#endif

#ifdef SOLVE_MIPREG

	cout << "**************MIP REG algorithm starts*****************" << endl;

	iter++;
	std::cout << "Iteration : " << iter << endl;

	status = CutMIP.MIPSolveMIP();
	if (status != 0) {
		std::cerr << endl << "CPLEX Error: " << status << endl;
		std::cerr << "Algorithm: Failure to solve MIP CUT problem!" << endl;
		std::cerr << "Iteration: " << iter << endl;
		std::cerr << "Terminating ..." << endl;
		exit(status);
	}

	status = CutMIP.MIPCutsetIncumsolnxandObjval();
	if (status != 0) {
		std::cerr << endl << "CPLEX Error: " << status << endl;
		std::cerr << "Algorithm: Failure to set incumbent solution for MIP CUT problem!" << endl;
		std::cerr << "Terminating ..." << endl;
		exit(status);
	}

	difference = CutMIP.MIPChecktermcondn();
	cout << "difference - Initial MIPCUT: " << difference << endl;

	solnx = CutMIP.MIPGetIncumsolnx();

	status = CostLP.LPAddcutlp(solnx);
	if (status != 0) {
		std::cerr << endl << "CPLEX Error: " << status << endl;
		std::cerr << "Algorithm: Failure to add cut to Cost Generator LP!" << endl;
		std::cerr << "Terminating ..." << endl;
		exit(status);
	}

	status = RegMIP.MIPSetSigma(solnx);
	if (status != 0) {
		std::cerr << endl << "CPLEX Error: " << status << endl;
		std::cerr << "Algorithm: Failure to set sigma for Regularized MIP!" << endl;
		std::cerr << "Terminating ..." << endl;
		exit(status);
	}

	status = RegMIP.MIPSetCurrobjsigma();
	if (status != 0) {
		std::cerr << endl << "CPLEX Error: " << status << endl;
		std::cerr << "Algorithm: Failure to change the objective coefficients for Regularized MIP!" << endl;
		std::cerr << "Terminating ..." << endl;
		exit(status);
	}


	status = RegMIP.MIPSolveMIP();
	if (status != 0) {
		std::cerr << endl << "CPLEX Error: " << status << endl;
		std::cerr << "Algorithm: Failure to change the objective coefficients for Regularized MIP!" << endl;
		std::cerr << "Terminating ..." << endl;
		exit(status);
	}

	difference = RegMIP.MIPRegsetIncumsolnxandObjval();
	cout << "difference - MIPREG: " << difference << endl;

	if (difference < NONZERO_LB) {
		status = RegMIP.MIPSetSigma(NULL);
		if (status != 0) {
			std::cerr << endl << "CPLEX Error: " << status << endl;
			std::cerr << "Algorithm: Failure to set sigma for Regularized MIP!" << endl;
			std::cerr << "Terminating ..." << endl;
			exit(status);
		}
	}
	else {
		solnx = RegMIP.MIPGetIncumsolnx();
		status = CostLP.LPAddcutlp(solnx);
		if (status != 0) {
			std::cerr << endl << "CPLEX Error: " << status << endl;
			std::cerr << "Algorithm: Failure to add cut to Cost Generator LP!" << endl;
			std::cerr << "Terminating ..." << endl;
			exit(status);
		}
	}

	std::cout << "*************************************" << endl;

	do {
		iter++;

		std::cout << "Iteration: " << iter << endl;

		status = CostLP.LPsolvelp();
		if (status != 0) {
			std::cerr << endl << "CPLEX Error: " << status << endl;
			std::cerr << "Algorithm: Failure to solve Cost Generator LP problem!" << endl;
			std::cerr << "Iteration: " << iter << endl;
			std::cerr << "Terminating ..." << endl;
			exit(status);
		}

		status = CostLP.LPSetcurrsolns();
		if (status != 0) {
			std::cerr << endl << "CPLEX Error: " << status << endl;
			std::cerr << "Algorithm: Failure to set solutions to Cost Generator LP!" << endl;
			std::cerr << "Terminating ..." << endl;
			exit(status);
		}

		solnd = CostLP.LPGetobjcoeffd();

		status = RegMIP.MIPSetCurrobjsigma(solnd);
		if (status != 0) {
			std::cerr << endl << "CPLEX Error: " << status << endl;
			std::cerr << "Algorithm: Failure to uodate objective coefficients for Reg MIP!" << endl;
			std::cerr << "Terminating ..." << endl;
			exit(status);
		}

		status = RegMIP.MIPSolveMIP();
		if (status != 0) {
			std::cerr << endl << "CPLEX Error: " << status << endl;
			std::cerr << "Algorithm: Failure to change the objective coefficients for Regularized MIP!" << endl;
			std::cerr << "Terminating ..." << endl;
			exit(status);
		}

		difference = RegMIP.MIPRegsetIncumsolnxandObjval();
		cout << "difference - MIPREG: " << difference << endl;

		solnx = RegMIP.MIPGetIncumsolnx();
		
		if (difference > NONZERO_LB) {
			status = CostLP.LPAddcutlp(solnx);
			if (status != 0) {
				std::cerr << endl << "CPLEX Error: " << status << endl;
				std::cerr << "Algorithm: Failure to add cut to Cost Generator LP!" << endl;
				std::cerr << "Terminating ..." << endl;
				exit(status);
			}
		}
		
		if (difference < NONZERO_LB) {
			status = CutMIP.MIPSetCurrobj(solnd);
			if (status != 0) {
				std::cerr << endl << "CPLEX Error: " << status << endl;
				std::cerr << "Algorithm: Failure to set current objective for MIP Cut!" << endl;
				std::cerr << "Terminating ..." << endl;
				exit(status);
			}

			status = CutMIP.MIPSolveMIP();
			if (status != 0) {
				std::cerr << endl << "CPLEX Error: " << status << endl;
				std::cerr << "Algorithm: Failure to solve MIP CUT problem!" << endl;
				std::cerr << "Iteration: " << iter << endl;
				std::cerr << "Terminating ..." << endl;
				exit(status);
			}

			status = CutMIP.MIPCutsetIncumsolnxandObjval();
			if (status != 0) {
				std::cerr << endl << "CPLEX Error: " << status << endl;
				std::cerr << "Algorithm: Failure to set incumbent solution for MIP CUT problem!" << endl;
				std::cerr << "Terminating ..." << endl;
				exit(status);
			}

			difference = CutMIP.MIPChecktermcondn();
			cout << "difference - MIPCUT: " << difference << endl;

			solnx = CutMIP.MIPGetIncumsolnx();
			
			if (difference > NONZERO_LB) {

				status = RegMIP.MIPSetSigma(solnx);
				if (status != 0) {
					std::cerr << endl << "CPLEX Error: " << status << endl;
					std::cerr << "Algorithm: Failure to update sigma for Regularized MIP!" << endl;
					std::cerr << "Terminating ..." << endl;
					exit(status);
				}

				status = CostLP.LPAddcutlp(solnx);
			    if (status != 0) {
				    std::cerr << endl << "CPLEX Error: " << status << endl;
				    std::cerr << "Algorithm: Failure to add cut to Cost Generator LP!" << endl;
				    std::cerr << "Terminating ..." << endl;
				    exit(status);
			    }
			}
		}

	std::cout << "************************************" << endl;
	} while(difference > NONZERO_LB);
#endif

	//difference > NONZERO_LB && difference < - NONZERO_LB
	std::cout << "Number of iterations: " << iter << endl;

	end = clock();
	std::cout << "CPU time for execution: "
		<< (double)(end - start) / CLOCKS_PER_SEC << " seconds." << endl;


	std::cout << "Optimal d is: " << endl;
	for (int j = 0; j < mip_numcols; j++)
		std::cout << solnd[j] << " ";
	cout << endl;

	objshift = 0;

	for (int j = 0; j < mip_numcols; j++) {
		objshift += fabs(solnd[j] - inicoeff[j]);
	}
	
	cout << "Delta C value is: " << objshift << endl;

	delete[] filename;
	delete[] newfname;
	delete[] dualfname;
	delete[] solnx;
	delete[] solnd;
	delete[] inicoeff;

	return(status);

} /*END main*/










