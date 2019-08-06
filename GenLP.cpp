
#include "GenLP.h"

// Construction/Destruction

GenLP::GenLP() {
	status = 0;
	solnstat = 0; 
	curr_ncols = 0;
	curr_nrows = 0;
	ini_nrows = 0;
	ini_ncols = 0;

	env = NULL;
	lp = NULL;
	objcoeffd = NULL;
	incum_solnyef = NULL;
	inicoeffc = NULL;
	reqd_soln = NULL;
}

GenLP::GenLP(CPXENVptr the_env) {
	status = 0;
	solnstat = 0;
	curr_ncols = 0;
	curr_nrows = 0;
	ini_nrows = 0;
	ini_ncols = 0;

	env = the_env;
	lp = NULL;
	objcoeffd = NULL;
	incum_solnyef = NULL;
	inicoeffc = NULL;
	reqd_soln = NULL; 
}

GenLP::~GenLP() {
}

//Methods

double GenLP::MIPScalarproduct(double* array1, double* array2, int beg, int end) {
	int j;
	double prod = 0;
	for (j = beg; j <= end; j++) {
		prod += array1[j] * array2[j];
	}
	return prod;
}

int GenLP::LPSetupLP(char* fname) {
	lp = CPXcreateprob(env, &status, fname);
	if (lp == NULL) {
		cerr << "Failed to create MIP." << endl;
	}

	/* Now read the file, and copy the data into the created lp object*/
	status = CPXreadcopyprob(env, lp, fname, "MPS");
	if (status) {
		cerr << "Failed to read and copy the problem data." << endl;
		return(status);
	}

	//Get the number of rows and columns
	ini_nrows = CPXgetnumrows(env, lp);
	ini_ncols = CPXgetnumcols(env, lp);

	//Setting the parameters to add new columns e and f
	int ccnt = 2 * ini_nrows;
	double* obj, * cmatval, * lb, * ub;
	int* cmatbeg, * cmatind, * indices;

	objcoeffd = new double[ini_nrows];
	inicoeffc = new double[ini_nrows];
	reqd_soln = new double[ini_nrows];

	status = CPXgetrhs(env, lp, inicoeffc, 0, ini_nrows - 1);
	if (status) {
		cerr << "Failed to read and copy the problem data." << endl;
		return(status);
	}

	cmatval = new double[ccnt];
	cmatind = new int[ccnt];
	cmatbeg = new int[ccnt];
	lb = new double[ccnt];
	ub = new double[ccnt];

	if (cmatval == NULL || cmatind == NULL || cmatbeg == NULL) {
		cerr << "Failed to allocate memory to new column arrays" << endl;
		return(1);
	}

	for (int j = 0; j < ccnt; j++) {
		cmatbeg[j] = j;
		lb[j] = 0;
		ub[j] = CPX_INFBOUND;
	}

	for (int j = 0; j < ini_nrows; j++) {
		cmatval[j] = 1;
		cmatind[j] = j;
	}

	for (int j = ini_nrows; j < ccnt; j++) {
		cmatval[j] = -1;
		cmatind[j] = j - ini_nrows;
	}

	status = CPXaddcols(env, lp, ccnt, ccnt, NULL, cmatbeg, cmatind, cmatval, lb, ub, NULL);
	if (status) {
		cerr << "Failure to CPLEX call CPXaddcols(.), error code: " << status << endl;
		return(status);
	}

	//Change the objective sense to minimization problem 
	CPXchgobjsen(env, lp, CPX_MIN);
	if (status) {
		cerr << "Failure to CPLEX call CPXchgobjsen(.), error code: " << status << endl;
		return(status);
	}

	//Change the objective coefficients
	curr_nrows = CPXgetnumrows(env, lp);
	curr_ncols = CPXgetnumcols(env, lp);

	if (curr_ncols != ccnt + ini_ncols) {
		cerr << "Failure while creating the cost generator LP, error code: " << status << endl;
		status = 1;
		return(status);
	}

	incum_solnyef = new double[curr_ncols];

	obj = new double[curr_ncols];
	indices = new int[curr_ncols];
	for (int j = 0; j < ini_ncols; j++) {
		obj[j] = 0;
		indices[j] = j;
	}

	for (int j = ini_ncols; j < curr_ncols; j++) {
		obj[j] = 1;
		indices[j] = j;
	}

	status = CPXchgobj(env, lp, curr_ncols, indices, obj);
	if (status) {
		cerr << "Failure to CPLEX call CPXchgobj(.), error code: " << status << endl;
		return(status);
	}

	delete[] cmatval;
	delete[] cmatind;
	delete[] cmatbeg;
	delete[] lb;
	delete[] ub;
	delete[] obj;
	delete[] indices;

	return(0);
}

int GenLP::LPSetreqdsoln(const double* reqdx) {

	for (int j = 0; j < ini_nrows; j++) {
		reqd_soln[j] = reqdx[j];
	}
	return(0);
}

int GenLP::LPAddcutlp(const double* solnx) {

	double* temparray, * rmatval;
	char* sense;
	double* rhs;
	int nzcnt = 0;
	int arrlen = 2 * ini_nrows;
	int* rmatbeg, * rmatind;

	temparray = new double[ini_nrows];
	sense = new char[1];
	rmatbeg = new int[1];
	rhs = new double[1];
	rmatval = new double[arrlen];
	rmatind = new int[arrlen];

	for (int j = 0; j < ini_nrows; j++) {
		temparray[j] = solnx[j] - reqd_soln[j];
	}

	for (int j = 0; j < ini_nrows; j++) {
		if (temparray[j] != 0) {
			rmatind[nzcnt] = j + ini_ncols;
			rmatval[nzcnt++] = temparray[j];
		}
	}

	for (int j = 0; j < nzcnt; j++) {
		rmatind[j + nzcnt] = rmatind[j] + ini_nrows;
		rmatval[j + nzcnt] = -rmatval[j];
	}

	rhs[0] = MIPScalarproduct(inicoeffc, temparray, 0, ini_nrows - 1);
	sense[0] = 'G';
	rmatbeg[0] = 0;
	nzcnt = 2 * nzcnt;

	status = CPXaddrows(env, lp, 0, 1, nzcnt, rhs, sense, rmatbeg, rmatind, rmatval, NULL, NULL);
	if (status) {
		cerr << "Failure to call CPXaddrows(.)." << endl;
		return(status);
	}

	curr_nrows = CPXgetnumrows(env, lp);

	delete[] temparray;
	delete[] sense;
	delete[] rmatbeg;
	delete[] rhs;
	delete[] rmatval;
	delete[] rmatind;

	return(0);
}

int GenLP::LPsolvelp() {

	status = CPXlpopt(env, lp);
	if (status) {
		cerr << "Failed to optimize LP." << endl;
		return(status);
	}

	solnstat = CPXgetstat(env, lp);
	if (status) {
		cerr << "Failure to CPLEX call CPXgetstat (.), error code: " << status << endl;
		return(status);
	}
	return(0);
}

int GenLP::LPSetcurrsolns() {

	status = CPXgetx(env, lp, incum_solnyef, 0, curr_ncols - 1);
	if (status) {
		cerr << "Failed to call CPXgetx(.)." << endl;
		return(1);
	}

	if (curr_ncols != ini_ncols + 2*ini_nrows) {
		cerr << "Solution array size mismatch." << status << endl;
		return(1);
	}

	for (int j = 0; j < ini_nrows; j++) {
		objcoeffd[j] = inicoeffc[j] - incum_solnyef[j+ini_ncols] + incum_solnyef[j+ini_ncols+ini_nrows];
	}

	return(0);
}

double* GenLP::LPGetobjcoeffd() {

	return objcoeffd;
}

void GenLP::LPPrintobjd() {

	std::cout << "The objective coeffcient d^k is: " << endl;
	for (int j = 0; j < ini_nrows; j++) {
		std::cout << j << ": " << objcoeffd[j] << " ";
	}
	cout << endl;
}
