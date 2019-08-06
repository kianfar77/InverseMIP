
#include "ProbMIP.h"

//Static members

char* ProbMIP::filename = NULL;
char* ProbMIP::dualname = NULL;
double* ProbMIP::reqd_soln = NULL;
double* ProbMIP::ini_objcoeff = NULL;
double* ProbMIP::orig_objcoeff = NULL;
int ProbMIP::mip_ncols = 0;
int ProbMIP::mip_nrows = 0;

//Construction/Destruction

ProbMIP::ProbMIP() {
	status = 0;
	solnstat = 0;
	objval = 0;
	objshift = 0;
	mipsigma_nrows = 0;
	mipsigma_ncols = 0;
	mipsigma_ncolswithx = 0;
	sigma = CPX_INFBOUND;
	env = NULL;
	mip = NULL;
	incum_solnx = NULL;
	incum_solnxsigma = NULL;
	opt_soln = NULL;
	curr_objcoeff = NULL;
	curr_objcoeffsigma = NULL;
}

ProbMIP::ProbMIP(CPXENVptr the_env) {
	status = 0;
	solnstat = 0;
	objval = 0;
	objshift = 0;
	mipsigma_nrows = 0;
	mipsigma_ncols = 0;
	mipsigma_ncolswithx = 0;
	status = 0;
	sigma = CPX_INFBOUND;
	env = the_env;
	mip = NULL;
	incum_solnx = NULL;
	incum_solnxsigma = NULL;
	opt_soln = NULL;
	curr_objcoeff = NULL;
	curr_objcoeffsigma = NULL;
}

ProbMIP::~ProbMIP() {

}

//Methods

double ProbMIP::MIPScalarproduct(const double* array1, const double* array2, int beg, int end) {
	int j;
	double prod = 0;
	for (j = beg; j <= end; j++) {
		prod += array1[j] * array2[j];
	}
	return prod;
}

double ProbMIP::MIPL1Norm(const double* array1, const double* array2, int len) {
	int j;
	double lsum = 0;
	for (j = 0; j < len; j++) {
		lsum += fabs(array1[j] - array2[j]);
	}
	return lsum;
}

int ProbMIP::MIPSetupMIP(char* fname) {

	int objsen;

	filename = new char[NAMELEN];
	dualname = new char[NAMELEN];
	if (filename == NULL || dualname == NULL) {
		cerr << "Memory allocation failure for filename arrays." << endl;
		return(1);
	}
	strcpy_s(filename, NAMELEN, fname);
	strcpy_s(dualname, NAMELEN, fname);strcat_s(dualname, NAMELEN, "_dual");

	/* Create the problem object in the CPLEX environment, using the filename as the problem name */
	mip = CPXcreateprob(env, &status, filename);
	if (mip == NULL) {
		cerr << "Failed to create MIP.\n" << endl;
		return(status);
	}

	/* Now read the file, and copy the data into the created mip object*/
	status = CPXreadcopyprob(env, mip, filename, "MPS");
	if (status) {
		cerr << "Failed to read and copy the problem data." << endl;
		return(1);
	}

	//Get the number of rows and columns
	mip_nrows = CPXgetnumrows(env, mip);
	mip_ncols = CPXgetnumcols(env, mip);

	//Allocating memory for arrays that will be populated later
	opt_soln = new double[mip_ncols];
	reqd_soln = new double[mip_ncols];
	incum_solnx = new double[mip_ncols];
	curr_objcoeff = new double[mip_ncols];
	orig_objcoeff = new double[mip_ncols];
	ini_objcoeff = new double[mip_ncols];
	if (opt_soln == NULL || reqd_soln == NULL || incum_solnx == NULL || orig_objcoeff == NULL) {
		cerr << "Memory allocation failure for objective coefficient array." << endl;
		return(1);
	}

	objsen = CPXgetobjsen(env, mip);
	if (status) {
		cerr << "Failed to get the objective sense." << endl;
		return(status);
	}

	status = CPXgetobj(env, mip, orig_objcoeff, 0, mip_ncols - 1);
	if (status) {
		cerr << "Failed to get the objective coefficients." << endl;
		return(status);
	}

	//Change minimization problem to maximization
	if (objsen == CPX_MIN) {
		CPXchgobjsen(env, mip, CPX_MAX);
		if (status) {
			cerr << "Failed to change the objective sense." << endl;
			return(status);
		}

		for (int j = 0; j < mip_ncols; j++) {
			orig_objcoeff[j] = -orig_objcoeff[j];
		}
	}

	return(0);
}

int ProbMIP::MIPSetupCutMIP(char* newfname) {
	
	/* Create the problem object in the CPLEX environment, using the filename as the problem name */
	mip = CPXcreateprob(env, &status, filename);
	if (mip == NULL) {
		cerr << "Failed to create MIP." << endl;
		return(status);
	}

	/* Now read the file, and copy the data into the created mip object*/
	status = CPXreadcopyprob(env, mip, newfname, "MPS");
	if (status) {
		cerr << "Failed to read and copy the problem data." << endl;
		return(1);
	}

	if (mip_ncols == 0) {
		cerr << "Number of columns wrongly set." << endl;
		return(1);
	}

	//Allocating memory for arrays that will be populated later
	opt_soln = new double[mip_ncols];
	incum_solnx = new double[mip_ncols];
	curr_objcoeff = new double[mip_ncols];
	if (opt_soln == NULL || incum_solnx == NULL || curr_objcoeff == NULL) {
		cerr << "Memory allocation failure for objective coefficient array." << endl;
		return(1);
	}

	for (int j = 0; j < mip_ncols; j++) {
		curr_objcoeff[j] = ini_objcoeff[j];
	}

	return(0);
}

int ProbMIP::MIPMallocRegMIP() {

	mipsigma_nrows = CPXgetnumcols(env, mip);
	mipsigma_ncolswithx = CPXgetnumcols(env, mip);
	mipsigma_ncols = CPXgetnumcols(env, mip) - mip_ncols;


	curr_objcoeffsigma = new double[mipsigma_ncolswithx];
	incum_solnxsigma = new double[mipsigma_ncolswithx];
	if (curr_objcoeffsigma == NULL || incum_solnxsigma == NULL) {
		cerr << "Failure to allocate memory to MIP Sigma problem data, bailing out ..." << endl;
		return(1);
	}
	return(0);
}

int ProbMIP::MIPWriteProblem(char* newname) {
	
	status = CPXwriteprob(env, mip, newname, "MPS");
	if (status) {
		cerr << "Failure to CPLEX call CPXwriteprob (.), error code: " << status << endl;
		return(status);
	}

	return(0);
}

int ProbMIP::MIPWriteDualRelaxed() {

	char* ctype;
	ctype = new char[mip_ncols];

	//Get the type of variables (continuous, integer, binary)
	status = CPXgetctype(env, mip, ctype, 0, mip_ncols - 1);
	if (status) {
		cerr << "Failure to CPLEX call CPXgetctype (.), error code: " << status << endl;
		return(status);
	}

	/*Change the MIP proble to its LP relaxed version*/
	status = CPXchgprobtype(env, mip, CPXPROB_LP);
	if (status) {
		cerr << "Failure to CPLEX call CPXchgprobtype (.), error code: " << status << endl;
		return(status);
	}

	status = CPXdualwrite(env, mip, dualname, &objshift);
	if (status) {
		cerr << "CPXdualwrite(.) failed" << endl;
		return(status);
	}

	if (objshift != 0) {
		cerr << "Dual problem needs modification. Objective shift: " << objshift << endl;
		return(1);
	}

	/*Change the problem back to MIP*/
	status = CPXchgprobtype(env, mip, CPXPROB_MILP);
	if (status) {
		cerr << "Failure to CPLEX call CPXchgprobtype (.), error code: " << status << endl;
		return(status);
	}

	status = CPXcopyctype(env, mip, ctype);
	if (status) {
		cerr << "Failure to CPLEX call CPXcopyctype(.), error code: " << status << endl;
		return(status);
	}

	delete[] ctype;

	return(0);
}

int ProbMIP::MIPGetSoln() {
	double* solntemp;
	int len;

	len = CPXgetnumcols(env, mip);
	solntemp = new double[len];

	status = CPXgetx(env, mip, solntemp, 0, len - 1);
	 
	delete[] solntemp;

	return(0);
}

int ProbMIP::MIPSolveMIP() {

	status = CPXmipopt(env, mip);
	if (status) {
		cerr << "Failed to optimize MIP." << endl;
		return(status);
	}

	solnstat = CPXgetstat(env, mip);
	if (status) {
		cerr << "Failure to CPLEX call CPXgetstat (.), error code: " << status << endl;
		return(status);
	}
	return(0);
}

int ProbMIP::MIPSetupRegMIP() {

	int nrows_B, ncols_B, nzcntA, nzcntB, rcmatspaceA, surplusA;
	double* cmatvalA, * cmatvalB, * rmatvalA;
	int* cmatbegA, * cmatindA, * cmatcntA, * rmatbegA, * rmatindA, * rmatcntA, * cmatbegB, * cmatindB, * cmatcntB;
	double* lbA, * ubA, * lbB, * ubB, * rangeAB;
	double* rhsA, * rhsB;
	double* objB;
	char* senseAB, * ctypeA;

	rcmatspaceA = mip_nrows * mip_ncols;

	cmatbegA = new int[mip_ncols];
	cmatindA = new int[rcmatspaceA];
	cmatcntA = new int[mip_ncols];
	cmatvalA = new double[rcmatspaceA];

	if (cmatbegA == NULL || cmatcntA == NULL || cmatindA == NULL || cmatvalA == NULL) {
		cerr << "Memory allocation failure for matrix A columns." << endl;
		return 1;
	}

	rmatbegA = new int[mip_nrows];
	rmatindA = new int[rcmatspaceA];
	rmatcntA = new int[mip_nrows];
	rmatvalA = new double[rcmatspaceA];
	if (rmatbegA == NULL || rmatcntA == NULL || rmatindA == NULL || rmatvalA == NULL) {
		cerr << "Memory allocation failure for matrix A rows." << endl;
		return 1;
	}

	//Get the constraint matrix A 
	status = CPXgetrows(env, mip, &nzcntA, rmatbegA, rmatindA, rmatvalA, rcmatspaceA, &surplusA, 0, mip_nrows - 1);
	if (status) {
		cerr << "Failure to CPLEX call CPXgetcols (.), error code: " << status << endl;
		return(status);
	}

	status = CPXgetcols(env, mip, &nzcntA, cmatbegA, cmatindA, cmatvalA, rcmatspaceA, &surplusA, 0, mip_ncols - 1);
	if (status) {
		cerr << "Failure to CPLEX call CPXgetcols (.), error code: " << status << endl;
		return(status);
	}

	lbA = new double[mip_ncols];
	ubA = new double[mip_nrows];
	rhsA = new double[mip_nrows];
	senseAB = new char[mip_nrows];
	rangeAB = new double[mip_nrows];
	ctypeA = new char[mip_ncols];
	if (lbA == NULL || ubA == NULL || rhsA == NULL || senseAB == NULL || rangeAB == NULL) {
		cerr << "Memory allocation failure for ub, lb, sense, rhs, and range arrays in A" << endl;
		return(1);
	}

	//Get the sense of constraints in A
	status = CPXgetsense(env, mip, senseAB, 0, mip_nrows - 1);
	if (status) {
		cerr << "Failure to CPLEX call CPXgetsense (.), error code: " << status << endl;
		return(status);
	}

	//Get the type of variables (continuous, integer, binary)
	status = CPXgetctype(env, mip, ctypeA, 0, mip_ncols - 1);
	if (status) {
		cerr << "Failure to CPLEX call CPXgetctype (.), error code: " << status << endl;
		return(status);
	}

	//Get the rhs of A
	status = CPXgetrhs(env, mip, rhsA, 0, mip_nrows - 1);
	if (status) {
		cerr << "Failure to CPLEX call CPXgetrhs(.), error code: " << status << endl;
		return(status);
	}

	//Get the range value of ranged constraints in A
	status = CPXgetrngval(env, mip, rangeAB, 0, mip_nrows - 1);
	if (status) {
		cerr << "Failure to CPLEX call CPXgetrhs(.), error code: " << status << endl;
		return(status);
	}

	//Get the lower bound on variables x
	status = CPXgetlb(env, mip, lbA, 0, mip_ncols - 1);
	if (status) {
		cerr << "Failure to CPLEX call CPXgetlb(.), error code: " << status << endl;
		return(status);
	}

	//Get the upper bound on variales x
	status = CPXgetub(env, mip, ubA, 0, mip_ncols - 1);
	if (status) {
		cerr << "Failure to CPLEX call CPXgetlb(.), error code: " << status << endl;
		return(status);
	}

	for (int j = 0; j < mip_ncols - 1; j++)
		cmatcntA[j] = cmatbegA[j + 1] - cmatbegA[j];
	cmatcntA[mip_ncols - 1] = nzcntA - cmatbegA[mip_ncols - 1];

	for (int i = 0; i < mip_nrows - 1; i++)
		rmatcntA[i] = rmatbegA[i + 1] - rmatbegA[i];
	rmatcntA[mip_nrows - 1] = nzcntA - rmatbegA[mip_nrows - 1];

	/*******************  Get the rmat vectors for constraint matrice B = A|-A *******************************/

	nrows_B = mip_nrows;
	ncols_B = 2 * mip_ncols;
	nzcntB = 2 * nzcntA;

	//cout << "numcolsB " << ncols_B << endl;
	//cout << "numrowsB " << nrows_B << endl;

	cmatbegB = new int[ncols_B];
	cmatindB = new int[nzcntB];
	cmatcntB = new int[ncols_B];
	cmatvalB = new double[nzcntB];

	if (cmatbegB == NULL || cmatcntB == NULL || cmatindB == NULL || cmatvalB == NULL) {
		cerr << "Memory allocation failure for matrix B." << endl;
		return 1;
	}

	for (int j = 0; j < mip_ncols; j++) {
		cmatbegB[j] = cmatbegA[j];
		cmatcntB[j] = cmatcntA[j];
	}

	for (int j = 0; j < mip_ncols; j++) {
		cmatcntB[j + mip_ncols] = cmatcntA[j];
		cmatbegB[j + mip_ncols] = cmatbegB[j+mip_ncols-1] + cmatcntB[j + mip_ncols - 1];
	}

	for (int i = 0; i < nzcntA; i++) {
		cmatindB[i] = cmatindA[i];
		cmatvalB[i] = cmatvalA[i];
	}

	for (int i = 0; i < nzcntA; i++) {
		cmatindB[i + nzcntA] = cmatindA[i];
		cmatvalB[i + nzcntA] = -cmatvalA[i];
	}

	/*************** Get the bounds, temporary objective coefficients, cytpes for the variables e and f (x = x_d + e - f) and rhs for the regularized MIP problem ********************/

	objB = new double[ncols_B];
	lbB = new double[ncols_B];
	ubB = new double[ncols_B];
	rhsB = new double[nrows_B];
	if (objB == NULL || lbB == NULL || ubB == NULL) {
		cerr << "Failed to allocate memory to ctype array" << endl;
		return(1);
	}

	int l = 0;

	for (int i = 0; i < nrows_B; i++) {
		rhsB[i] = rhsA[i];
		for (int j = 0; j < rmatcntA[i]; j++) {
			rhsB[i] -= rmatvalA[l] * reqd_soln[rmatindA[l]];
			l++;
		}
	}

	for (int j = 0; j < ncols_B / 2; j++)
		objB[j] = curr_objcoeff[j];
	for (int j = ncols_B / 2; j < ncols_B; j++)
		objB[j] = -curr_objcoeff[j - mip_ncols];

	for (int j = 0; j < ncols_B; j++) {
		lbB[j] = 0;
		ubB[j] = CPX_INFBOUND;
	}

	//Copy the prolem to the mip object
	status = CPXcopylp(env, mip, ncols_B, nrows_B, CPX_MAX, objB, rhsB, senseAB, cmatbegB, cmatcntB, cmatindB, cmatvalB, lbB, ubB, rangeAB);
	if (status) {
		cerr << "Failure to CPLEX call CPXcopylpwnames(.), error code: " << status << endl;
		return(status);
	}

	/************************* Adding constraints based on lb <= x <= ub **********************/

	int* rmatbegtemp1, * rmatindtemp, * rmatindtemp1;
	double* rmatvaltemp, * rhstemp, * rmatvaltemp1;
	char* sensetemp;

	rmatbegtemp1 = new int[1];
	rmatbegtemp1[0] = 0;
	rmatvaltemp = new double[2];
	rmatvaltemp[0] = 1;
	rmatvaltemp[1] = -1;
	rmatindtemp = new int[2];
	sensetemp = new char[1];
	rhstemp = new double[1];

    //Adding lb constraints
	for (int j = 0; j < mip_ncols; j++) {
		if (lbA[j] > -INFBOUND) {
			rmatindtemp[0] = j;
			rmatindtemp[1] = j + mip_ncols;
			sensetemp[0] = 'G';
			rhstemp[0] = lbA[j] - reqd_soln[j];
			status = CPXaddrows(env, mip, 0, 1, 2, rhstemp, sensetemp, rmatbegtemp1, rmatindtemp, rmatvaltemp, NULL, NULL);
			if (status) {
				cerr << " Failure to CPLEX call CPXaddrows(.): Failed to add e - f >= lb - xd  constraints." << endl;
				cerr << " CPLEX error code = " << status << endl;
			}
		}
	}
	
	//Adding ub constraints
	for (int j = 0; j < mip_ncols; j++) {
		if (ubA[j] < INFBOUND) {
			rmatindtemp[0] = j;
			rmatindtemp[1] = j + mip_ncols;
			sensetemp[0] = 'L';
			rhstemp[0] = ubA[j] - reqd_soln[j];
			status = CPXaddrows(env, mip, 0, 1, 2, rhstemp, sensetemp, rmatbegtemp1, rmatindtemp, rmatvaltemp, NULL, NULL);
			if (status) {
				cerr << " Failure to CPLEX call CPXaddrows(.): Failed to add e - f >= lb - xd  constraints." << endl;
				cerr << " CPLEX error code = " << status << endl;
			}

		}
	}

	//*******************Adding variables x = x^d + e - f*********************************

	status = CPXnewcols(env, mip, mip_ncols, NULL, lbA, ubA, ctypeA, NULL);
	if (status) {
		cerr << " Failure to CPLEX call CPXnewcols(.): Failed to add x variables." << endl;
		cerr << " CPLEX error code = " << status << endl;
	}

	rmatindtemp1 = new int[3];
	rmatvaltemp1 = new double[3];

	rmatvaltemp1[0] = -1;
	rmatvaltemp1[1] =  1;
	rmatvaltemp1[2] =  1;
	sensetemp[0] = 'E';

	for (int j = 0; j < mip_ncols; j++) {
		rmatindtemp1[0] = j;
		rmatindtemp1[1] = j + mip_ncols;
		rmatindtemp1[2] = j + 2*mip_ncols;
		rhstemp[0] = reqd_soln[j];
		status = CPXaddrows(env, mip, 0, 1, 3, rhstemp, sensetemp, rmatbegtemp1, rmatindtemp1, rmatvaltemp1, NULL, NULL);
		if (status) {
			cerr << " Failure to CPLEX call CPXaddrows(.): Failed to add x - x^d = e - f constraints." << endl;
			cerr << " CPLEX error code = " << status << endl;
		}
	}

	status = MIPMallocRegMIP();
	if (status) {
		cerr << "Failure to call MIPMallocRegMIP(), error code: " << status << endl;
		return(status);
	}

	delete[]  cmatbegA;
	delete[]  cmatindA;
	delete[]  cmatcntA;
	delete[]  cmatvalA;

	delete[]  rmatbegA;
	delete[]  rmatindA;
	delete[]  rmatcntA;
	delete[]  rmatvalA;

	delete[]  lbA;
	delete[]  ubA;
	delete[]  rhsA;
	delete[]  senseAB;
	delete[]  rangeAB;
	delete[]  ctypeA;

	delete[]  cmatbegB;
	delete[]  cmatindB;
	delete[]  cmatcntB;
	delete[]  cmatvalB;

	delete[]  objB;
	delete[]  lbB;
	delete[]  ubB;
	delete[]  rhsB;

	delete[]  rmatbegtemp1;
	delete[]  rmatvaltemp;
	delete[]  rmatindtemp;
	delete[]  sensetemp;
	delete[]  rhstemp;

	delete[]  rmatindtemp1;
	delete[]  rmatvaltemp1;

	return(0);
}

int ProbMIP::MIPSetSigma(const double* solnx) {

	if (solnx == NULL)
		sigma = sigma / 2;
	else {
		for (int j = 0; j < mip_ncols; j++)
			incum_solnx[j] = solnx[j];

		double numr1, numr2, denr;
		int tempsigma;

		numr1 = MIPScalarproduct(curr_objcoeff, incum_solnx, 0, mip_ncols - 1);
		numr2 = MIPScalarproduct(curr_objcoeff, reqd_soln, 0, mip_ncols - 1);

		denr = MIPL1Norm(incum_solnx, reqd_soln, mip_ncols);

		if (denr == 0) {
			cerr << " Failure to update/initialize sigma - zero denomenator" << endl;
			return(1);
		}

		tempsigma = (numr1 - numr2) / denr;

		if (tempsigma < 0) {
			cerr << " Failure to update/initialize sigma - nonzero value" << endl;
			return(1);
		}

		if (tempsigma < sigma / 2)
			sigma = tempsigma;
		else
			sigma = sigma / 2;
	}

	return(0);
}

int ProbMIP::MIPSetCutSoln() {
	
	status = CPXgetx(env, mip, incum_solnx, 0, mip_ncols - 1);
	if (status) {
		cerr << " Failure to CPLEX call CPXgetx(.)." << endl;
		cerr << " CPLEX error code = " << status << endl;
		return(status);
	}
	return(0);
}

int ProbMIP::MIPSetRegSoln() {

	status = CPXgetx(env, mip, incum_solnxsigma, 0, mipsigma_ncolswithx - 1);
	if (status) {
		cerr << " Failure to CPLEX call CPXgetx(.)." << endl;
		cerr << " CPLEX error code = " << status << endl;
		return(status);
	}

	if (mipsigma_ncolswithx != 3 * mip_ncols) {
		cerr << "The solutions arrays are of wrong size" << endl;
		return(1);
	}

	for (int j = 0; j < mip_ncols; j++)
		incum_solnx[j] = reqd_soln[j] + incum_solnxsigma[j] - incum_solnxsigma[j + mip_ncols];
	return(0);
}

int ProbMIP::MIPSetupReqdSolnandIniObj(int seed) {

	int status, rseed, cnt = 0;

	int* randlb, * randub, * indices;
	double* randobj;

	randlb = new int[mip_ncols];
	randub = new int[mip_ncols];
	indices = new int[mip_ncols];
	randobj = new double[mip_ncols];

	for (int j = 0; j < mip_ncols; j++) {
		randlb[j] = (int)orig_objcoeff[j] - 2*fabs(orig_objcoeff[j]);
		randub[j] = (int)orig_objcoeff[j] + 2*fabs(orig_objcoeff[j]);
		indices[j] = j;
	}

	do {
		rseed = seed + cnt;
		srand(rseed);
		for (int j = 0; j < mip_ncols; j++) {
			randobj[j] = randlb[j] + (rand() % (randub[j] - randlb[j] + 1));
		}

		status = CPXchgobj(env, mip, mip_ncols, indices, randobj);
		if (status) {
			cerr << " Failure to CPLEX call CPXchgobj(.)." << endl;
			cerr << " CPLEX error code = " << status << endl;
			return(status);
		}
		status = MIPSolveMIP();
		if (status) {
			cerr << " Failure to CPLEX call MIPSolveMIP(.)." << endl;
			cerr << " CPLEX error code = " << status << endl;
			return(status);
		}

		if (solnstat == CPXMIP_OPTIMAL) {
			cnt = 10;
			status = CPXgetx(env, mip, reqd_soln, 0, mip_ncols - 1);
			if (status) {
				cerr << " Failure to CPLEX call CPXgetx(.)." << endl;
				cerr << " CPLEX error code = " << status << endl;
				return(status);
			}
		}
		cnt++;
	} while (cnt < 10);

	if (cnt == 11) {
		status = MIPGetIniObj(randlb, randub, rseed + 255);
		if (status) {
			cerr << " Failure to CPLEX call MIPGetIniObj(.)." << endl;
			cerr << " CPLEX error code = " << status << endl;
			return(status);
		}
	}
	else {
		cerr << "All random objective coefficients failed to give a bounded MIP problem." << endl;
		return(1);
	}

	delete[] randlb;
	delete[] randub;
	delete[] indices;
	delete[] randobj;

	return(0);
}

int ProbMIP::MIPGetIniObj(const int* arraylb, const int* arrayub, int seed) {

	int* indices, rseed, cnt = 0;
	double* randobj, * tempsoln;
	double tempobjval, prod;

	indices = new int[mip_ncols];
	randobj = new double[mip_ncols];
	tempsoln = new double[mip_ncols];

	for (int j = 0; j < mip_ncols; j++) {
		indices[j] = j;
	}

	do {
		rseed = seed + cnt;
		srand(rseed);
		for (int j = 0; j < mip_ncols; j++) {
			randobj[j] = arraylb[j] + (rand() % (arrayub[j] - arraylb[j] + 1));
			//cout << "rand: " << randobj[j] << endl;
		}

		status = CPXchgobj(env, mip, mip_ncols, indices, randobj);
		if (status) {
			cerr << " Failure to CPLEX call CPXchgobj(.)." << endl;
			cerr << " CPLEX error code = " << status << endl;
			return(status);
		}

		status = MIPSolveMIP();
		if (status) {
			cerr << " Failure to CPLEX call MIPSolveMIP(.)." << endl;
			cerr << " CPLEX error code = " << status << endl;
			return(status);
		}

		if (solnstat == CPXMIP_OPTIMAL) {
			status = CPXgetx(env, mip, tempsoln, 0, mip_ncols - 1);
			if (status) {
				cerr << " Failure to CPLEX call CPXgetx(.)." << endl;
				cerr << " CPLEX error code = " << status << endl;
				return(status);
			}

			prod = MIPScalarproduct(reqd_soln, randobj, 0, mip_ncols - 1);
			status = CPXgetobjval(env, mip, &tempobjval);
			if (status) {
				cerr << " Failure to CPLEX call CPXgetobjval(.)." << endl;
				cerr << " CPLEX error code = " << status << endl;
				return(status);
			}
			if (prod != tempobjval)
				cnt = 10;
		}
		cnt++;

	} while (cnt < 10);

	if (cnt == 11) {
		for (int j = 0; j < mip_ncols; j++) {
			ini_objcoeff[j] = randobj[j];
			curr_objcoeff[j] = randobj[j];
		}
	}
	else {
		cout << "All random objective coefficients failed to give a bounded MIP problem." << endl;
		return(1);
	}
	delete[] indices;
	delete[] randobj;
	delete[] tempsoln;

	return(0);
}

int ProbMIP::MIPSetCurrobj(const double* objd) {

	int* indices, objsen;
	indices = new int[mip_ncols];

	for (int j = 0; j < mip_ncols; j++) {
		curr_objcoeff[j] = objd[j];
		indices[j] = j;
	}

	status = CPXchgobj(env, mip, mip_ncols, indices, curr_objcoeff);
	if (status) {
		cerr << " Failure to CPLEX call CPXchgobj(.)." << endl;
		cerr << " CPLEX error code = " << status << endl;
		return(status);
	}

	objsen = CPXgetobjsen(env, mip);

	//Change minimization problem to maximization
	if (objsen == CPX_MIN) {
		CPXchgobjsen(env, mip, CPX_MAX);
		if (status) {
			cerr << "Failed to change the objective sense." << endl;
			return(status);
		}
	}

	delete[] indices;

	return(0);
}

double ProbMIP::MIPRegsetIncumsolnxandObjval() {

	double objvaltemp;

	status = CPXgetx(env, mip, incum_solnxsigma, 0, mipsigma_ncolswithx - 1);
	if (status) {
		cerr << "Failed to call CPXgetx(.)." << endl;
		return(1);
	}

	for (int j = 0; j < mip_ncols; j++) {
		incum_solnx[j] = reqd_soln[j] + incum_solnxsigma[j] - incum_solnxsigma[j + mip_ncols];
	}

	status = CPXgetobjval(env, mip, &objvaltemp);
	if (status) {
		cerr << "Failed to call CPXgetobjval(.)." << endl;
		return(1);
	}

	objval = MIPScalarproduct(curr_objcoeff, reqd_soln, 0, mip_ncols - 1) + objvaltemp;

	return(objvaltemp);
}

int ProbMIP::MIPCutsetIncumsolnxandObjval() {
	
	status = CPXgetx(env, mip, incum_solnx, 0, mip_ncols - 1);
	if (status) {
		cerr << "Failed to call CPXgetx(.)." << endl;
		return(1);
	}

	status = CPXgetobjval(env, mip, &objval);
	if (status) {
		cerr << "Failed to call CPXgetobjval(.)." << endl;
		return(1);
	}
	
	return(0);
}

double* ProbMIP::MIPGetIncumsolnx() {
	return incum_solnx;
}

double* ProbMIP::MIPGetiniobjcoeff() {
	return ini_objcoeff;
}

double* ProbMIP::MIPGetreqdsoln() {
	return reqd_soln;
}

int ProbMIP::MIPSetCurrobjsigma() {
	
	int* indices, objsen;
	indices = new int[mipsigma_ncolswithx];

	for (int j = 0; j < mipsigma_ncolswithx; j++)
		indices[j] = j;

	for (int j = 0; j < mip_ncols; j++)
		curr_objcoeffsigma[j] = curr_objcoeff[j] - sigma;

	for (int j = 0; j < mip_ncols; j++)
		curr_objcoeffsigma[j + mip_ncols] = -curr_objcoeff[j] - sigma;

	for (int j = 0; j < mip_ncols; j++)
		curr_objcoeffsigma[j + 2*mip_ncols] = 0;

	status = CPXchgobj(env, mip, mipsigma_ncolswithx, indices, curr_objcoeffsigma);
	if (status) {
		cerr << " Failure to CPLEX call CPXchgobj(.)." << endl;
		cerr << " CPLEX error code = " << status << endl;
		return(status);
	}

	objsen = CPXgetobjsen(env, mip);

	//Change minimization problem to maximization
	if (objsen == CPX_MIN) {
		CPXchgobjsen(env, mip, CPX_MAX);
		if (status) {
			cerr << "Failed to change the objective sense." << endl;
			return(status);
		}
	}

	delete[] indices;

	return(0);
}

int ProbMIP::MIPSetCurrobjsigma(const double* objd) {

	int* indices, objsen;;
	indices = new int[mipsigma_ncolswithx];

	for (int j = 0; j < mip_ncols; j++) {
		curr_objcoeff[j] = objd[j];
	}

	for (int j = 0; j < mipsigma_ncolswithx; j++)
		indices[j] = j;

	for (int j = 0; j < mip_ncols; j++)
		curr_objcoeffsigma[j] = curr_objcoeff[j] - sigma;

	for (int j = 0; j < mip_ncols; j++)
		curr_objcoeffsigma[j + mip_ncols] = -curr_objcoeff[j] - sigma;

	for (int j = 0; j < mip_ncols; j++)
		curr_objcoeffsigma[j + 2*mip_ncols] = 0;

	status = CPXchgobj(env, mip, mipsigma_ncolswithx, indices, curr_objcoeffsigma);
	if (status) {
		cerr << " Failure to CPLEX call CPXchgobj(.)." << endl;
		cerr << " CPLEX error code = " << status << endl;
		return(status);
	}

	objsen = CPXgetobjsen(env, mip);

	//Change minimization problem to maximization
	if (objsen == CPX_MIN) {
		CPXchgobjsen(env, mip, CPX_MAX);
		if (status) {
			cerr << "Failed to change the objective sense." << endl;
			return(status);
		}
	}

	delete[] indices;

	return(0);
}

int ProbMIP::MIPGetncols() {
	
	return(mip_ncols);
}

double ProbMIP::MIPChecktermcondn() {
	
	double diff, difftemp1, difftemp2;

	difftemp1 = MIPScalarproduct(reqd_soln, curr_objcoeff, 0, mip_ncols - 1);
	difftemp2 = MIPScalarproduct(incum_solnx, curr_objcoeff, 0, mip_ncols - 1);

	diff = difftemp2 - difftemp1;

	return(diff);	
}
