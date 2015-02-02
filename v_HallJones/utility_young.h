#define MAX(a,b) (((a)>(b))?(a):(b))
#define MIN(a,b) (((a)<(b))?(a):(b))

inline
double GS(double y, double tau0, double tau1, double tau2)
{
	double bracket1 = pow(y, -tau1) + tau2;
	double bracket2 = y - pow(bracket1, -1 / tau1);
	return tau0 * bracket2;
}

inline
double GSp(double y, double tau0, double tau1, double tau2)
{
	double bracket1 = 1 + tau2 * pow(y, tau1);
	double bracket2 = 1 - pow(bracket1, -1 / tau1 - 1);
	return tau0 * bracket2;
}

/** define user function here:
 *
 * Please following the two steps to write user function:
 * 1. In first step, take x[] and data[] input, output futureState[][] which will be used in interpolation
 * 2. In second step, take x[], data[], vFuture[] and gFuture[] as input, output v, grad[]
 *  */
USER_FUNC_HEAD
{
	/**
	 * preprocessor, do not edit
	 */
	USER_FUNC_PRE;

	/**
	 * step 1:
	 * input: x[], data[]
	 * output: stateFuture[][] of length [nvec, stateDim]
	 */
	// private input ebp
	int POPCOUNT = 0;
	// input stack goes from top to bottom
#define POP(var) double var = data[POPCOUNT++]
#define POPARRAY(var, length) double* var = &data[POPCOUNT]; POPCOUNT+=length

	// global input ebp
	int GPOPCOUNT = 0;
	// global input stack goes from top to bottom
#define GPOP(var) double var = g_shared_data[GPOPCOUNT++]
#define GPOPARRAY(var, length) double* var = &g_shared_data[GPOPCOUNT]; GPOPCOUNT+=length

	// user input ebp
	int UPOPCOUNT = 0;
	// user stack goes from bottom to top
#define UPOP(var) double var = data[MAXDATA-(++UPOPCOUNT)]

	// get global data
	// utility
	GPOP(beta); // discount factor
	GPOP(lambda); // comsumption leisure weight
	GPOP(zeta); // consumption weight
	GPOP(rho); // elasticity of health
	GPOP(nu); // risk aversion
	GPOP(d); // flow utility of alive
	// survival
	GPOP(psi1);
	GPOP(psi2);
	// health production function
	GPOP(A_h); // production level
	GPOP(theta); // production curvature
	// tax
	GPOP(tau_ss);
	GPOP(tau0_gs);
	GPOP(tau1_gs);
	GPOP(tau2_gs);

	// individual data
	POP(budget);
	POP(hl);
	POP(hl_after_derpreciation);
	POP(wage);
	POP(kinc);
	POP(gamma); // coinsurance rate
	POP(chi); // deductible

	// controls
	double kp = x[0];
	double m = x[1];
	double n = x[2];

	// labor income
	double ninc = wage * n *(1 - tau_ss);

	// tax
	double tax = GS(ninc + kinc, tau0_gs, tau1_gs, tau2_gs);
	double taxp = GSp(ninc + kinc, tau0_gs, tau1_gs, tau2_gs);
	double dinc_dn = wage*(1 - tau_ss)*(1 - taxp);

	//
	double c = budget - (kp + m - (1 - gamma)*MAX(m - chi, 0)) + (ninc - tax);
	double dc_dkp = -1;
	double dc_dm = -((m > chi) ? gamma : 1);
	double dc_dn = dinc_dn;
	// c = MAX(c, exp(-1e20));

	// 
	double l = 1 - n;
	double dl_dn = -1;
	// l = MAX(l, exp(-1e20));

	// h
	double incr_h = A_h*pow(m, theta);
	// double incr_h = hl*A_h*pow(m, theta);
	double h = hl_after_derpreciation + incr_h;
	double dh_dm = theta*incr_h / m;

	double inf = 1e20;

    // uncomment to shutdown labor leisure choice
    l = 1;
    dl_dn = 0;

	// u
    /*
	double cl_comp = (c > 0 && l > 0) ? pow(c, zeta)*pow(l, 1 - zeta) : -inf;
	double dcl_comp_dc = (c > 0) ? zeta*cl_comp / c : inf;
	double dcl_comp_dl = (l > 0) ? (1 - zeta)*cl_comp / l : inf;
	double u_base = lambda*pow(cl_comp, rho) + (1 - lambda)*pow(h, rho);
	double u_base_power = pow(u_base, (1 - nu) / rho);
	double u_base_power_minus_1 = u_base_power / u_base;
	double u = u_base_power / (1 - nu);
	double du_dh = (1 - lambda)*pow(h, rho - 1)*u_base_power_minus_1;
	double du_dcl_comp = u_base_power_minus_1*lambda*pow(cl_comp, rho - 1);
	double du_dc = du_dcl_comp * dcl_comp_dc;
	double du_dl = du_dcl_comp * dcl_comp_dl;
    */

	double cl_comp = (c > 0 && l > 0) ? pow(c, zeta)*pow(l, 1 - zeta) : -inf;
	double dcl_comp_dc = (c > 0) ? zeta*cl_comp / c : inf;
	double dcl_comp_dl = (l > 0) ? (1 - zeta)*cl_comp / l : -inf;
	double u = pow(cl_comp, 1 - nu) / (1 - nu) + lambda*pow(h, 1 - rho) / (1 - rho);
	double du_dh = lambda*pow(h, -rho);
	double du_dcl_comp = pow(cl_comp, -nu);
	double du_dc = du_dcl_comp * dcl_comp_dc;
    double du_dl = du_dcl_comp * dcl_comp_dl;
	

	// psi
    /*
    double psi = 1 - 1 / (pow(h, psi1) + psi2);
    double dpsi_dh = 1 / pow((pow(h, psi1) + psi2), 2)*psi1*pow(h, psi1 - 1);
    */
    /*
    double psi = psi2 - pow(h, -psi1);
    double dpsi_dh = psi1*pow(h, -psi1 - 1);
    */
    /*
	// psi
	double psi = 1 - 1 / h;
	double dpsi_dh = 1 / h / h;
    */
    double pow_part = pow(h, psi2);
    double exp_part = exp(-psi1*pow_part);
    double psi = 1 - exp_part;
    double dpsi_dh = psi1*psi2*(pow_part / h)*exp_part;
	/*
	// psi
	double psi = h;
	double dpsi_dh = 1;
	*/


	// future State
	stateFuture[0][0] = kp;
	stateFuture[0][1] = h;

	/**
	 * interpolate future state, do not edit
	 */
	USER_FUNC_INTERP;

	/**
	 * step 2:
	 * input: x[], data[], vFuture[], gFuture[][]
	 * output:
	 *    v: scalar value eavaluated at x
	 *    grad[] of size controlDim: gradient of v evaluated at x
	 *    cons[] of size nonlin: evaluation of constraint
	 *    consgrad[][] of size [nonlin, controlDim]: jacobian of cons
	 */
	if (f) {
		v = d + u + beta*psi*vFuture[0];
	}

	if (grad) {
		double dcontv_dh = beta*(dpsi_dh*vFuture[0] + psi*gFuture[0][1]);
		grad[0] = du_dc*dc_dkp + beta*psi*gFuture[0][0];
		grad[1] = du_dc*dc_dm + (du_dh + dcontv_dh)*dh_dm;
		grad[2] = du_dc*dc_dn + du_dl*dl_dn;
	}

	if (cons) {
		cons[0] = c;
		cons[1] = l;
	}

	if (consgradRaw) {
		consgrad[0][0] = dc_dkp;
		consgrad[0][1] = dc_dm;
		consgrad[0][2] = dc_dn;

		consgrad[1][0] = 0;
		consgrad[1][1] = 0;
		consgrad[1][2] = dl_dn;
	}
	USER_FUNC_RETURN;
}

#undef MIN
#undef MAX
