// This is a personal academic project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: http://www.viva64.com
#include <boost/config.hpp>
#ifdef BOOST_MSVC
#pragma warning(disable : 4996)
#endif

#include "include/CTR.hpp"
#include <execution>
#include <iostream>
#include <fstream> // headerfor writing data files

// overloaded class constructor
CTR::CTR(const std::array<std::shared_ptr<Tube>, 3UL> &Tb, blaze::StaticVector<double, 6UL> &q, double Tol, mathOp::rootFindingMethod method) : m_accuracy(Tol), m_method(method), m_Tubes(Tb), m_beta(blaze::subvector<0UL, 3UL>(q)), m_q(q)
{
	this->m_eps = 5e-4;
	this->m_theta_0 = {0, q[4UL] - q[3UL], q[5UL] - q[4UL]};
	this->m_e3 = {0, 0, 1};
	this->m_segment = std::make_unique<Segment>(this->m_Tubes, this->m_beta);
	this->m_r_0 = 0.0;
	this->m_h_0 = {1.0, 0.0, 0.0, 0.0};

	this->m_stateEquations = std::make_unique<ODESystem>();
	this->m_stateObserver = std::make_unique<Observer>(this->m_y, this->m_s);
}

// copy constructor
CTR::CTR(const CTR &rhs) : m_accuracy(rhs.m_accuracy), m_eps(rhs.m_eps), m_method(rhs.m_method), m_Tubes(rhs.m_Tubes),
						   m_beta(rhs.m_beta), m_q(rhs.m_q), m_theta_0(rhs.m_theta_0), m_e3(rhs.m_e3),
						   m_r_0(rhs.m_r_0), m_h_0(rhs.m_h_0), m_y(rhs.m_y), m_s(rhs.m_s)
{
	this->m_segment = std::make_unique<Segment>(this->m_Tubes, this->m_beta);
	this->m_stateEquations = std::make_unique<ODESystem>();
	this->m_stateObserver = std::make_unique<Observer>(this->m_y, this->m_s);
}

// move constructor
CTR::CTR(CTR &&rhs) noexcept
{
	// handling self assignment
	if (this != &rhs)
	{
		this->m_accuracy = rhs.m_accuracy;
		this->m_eps = rhs.m_eps;
		this->m_method = rhs.m_method;
		this->m_Tubes = std::move(rhs.m_Tubes);
		this->m_beta = std::move(rhs.m_beta);
		this->m_q = std::move(rhs.m_q);
		this->m_theta_0 = std::move(rhs.m_theta_0);
		this->m_e3 = std::move(rhs.m_e3);
		this->m_segment = std::move(rhs.m_segment);
		this->m_r_0 = std::move(rhs.m_r_0);
		this->m_h_0 = std::move(rhs.m_h_0);
		this->m_y = std::move(rhs.m_y);
		this->m_s = std::move(rhs.m_s);
		this->m_stateEquations = std::move(rhs.m_stateEquations);
		this->m_stateObserver = std::move(rhs.m_stateObserver);
	}
}

// CTR destructor
CTR::~CTR()
{
	// nothing to be done as smart pointers are being used
}

// copy assignment operator
CTR &CTR::operator=(const CTR &rhs)
{
	// handling self assignment
	if (this != &rhs)
	{
		this->m_accuracy = rhs.m_accuracy;
		this->m_eps = rhs.m_eps;
		this->m_method = rhs.m_method;
		this->m_Tubes = rhs.m_Tubes;
		this->m_beta = rhs.m_beta;
		this->m_q = rhs.m_q;
		this->m_theta_0 = rhs.m_theta_0;
		this->m_e3 = rhs.m_e3;
		this->m_segment = std::make_unique<Segment>(this->m_Tubes, this->m_beta);
		this->m_r_0 = rhs.m_r_0;
		this->m_h_0 = rhs.m_h_0;
		this->m_y = rhs.m_y;
		this->m_s = rhs.m_s;
		this->m_stateEquations = std::make_unique<ODESystem>();
		this->m_stateObserver = std::make_unique<Observer>(this->m_y, this->m_s);
	}

	return *this;
}

// move assignment operator
CTR &CTR::operator=(CTR &&rhs) noexcept
{
	// handling self assignment
	if (this != &rhs)
	{
		this->m_accuracy = rhs.m_accuracy;
		this->m_eps = rhs.m_eps;
		this->m_method = rhs.m_method;
		this->m_Tubes = std::move(rhs.m_Tubes);
		this->m_beta = std::move(rhs.m_beta);
		this->m_q = std::move(rhs.m_q);
		this->m_theta_0 = std::move(rhs.m_theta_0);
		this->m_e3 = std::move(rhs.m_e3);
		this->m_segment = std::move(rhs.m_segment);
		this->m_r_0 = std::move(rhs.m_r_0);
		this->m_h_0 = std::move(rhs.m_h_0);
		this->m_y = std::move(rhs.m_y);
		this->m_s = std::move(rhs.m_s);
		this->m_stateEquations = std::move(rhs.m_stateEquations);
		this->m_stateObserver = std::move(rhs.m_stateObserver);
	}

	return *this;
}

// function that resets the initial parameters for the ODESolver
void CTR::reset(const blaze::StaticVector<double, 5UL> &initGuess)
{
	vec3d uz_0 = {initGuess[2UL], initGuess[3UL], initGuess[4UL]};
	// alpha1_0 =  alpha_1 - beta_1 * uz_1(0)
	double alpha1_0 = this->m_q[3UL] - this->m_beta[0UL] * uz_0[0UL];

	// clearing the observer's containers
	if (!this->m_s.empty())
	{
		this->m_y.clear();
		this->m_y.reserve(300UL);
		this->m_s.clear();
		this->m_s.reserve(300UL);
	}

	// theta_i(0) = alpha_1 - alpha_i - (beta_i * uz_i(0) - beta_1 * uz_1(0))
	this->m_theta_0 = {0.0, this->m_q[4UL] - this->m_beta[1UL] * uz_0[1UL] - alpha1_0, this->m_q[5UL] - this->m_beta[2UL] * uz_0[2UL] - alpha1_0};

	// origin of the local frame at s = 0
	this->m_r_0 = 0.0;
	// transforming proximal orientation to quaternion representation
	mathOp::euler2Quaternion(0.0, alpha1_0, 0.0, this->m_h_0);
}

// function that solves (integrates) the CTR ode (state) equations
blaze::StaticVector<double, 5UL> CTR::ODESolver(const blaze::StaticVector<double, 5UL> &initGuess)
{
	// initGuess = [u1_x(0) u1_y(0) u1_z(0) u2_z(0) u3_z(0)]  --> vector of initial guesses for solving the BVP
	this->reset(initGuess); // resets CTR parameters and variables for a new iteration of the ode-solver

	// retrieving the bending & torsional stiffness and precurvatures in all segments of the CTR in the current configuration
	auto tuple = this->m_segment->returnParameters();

	blaze::HybridMatrix<double, 3UL, 18UL, blaze::columnMajor> EI(std::get<0UL>(tuple)), GJ(std::get<1UL>(tuple)), U_x(std::get<2UL>(tuple)), U_y(std::get<3UL>(tuple));
	std::vector<double> S(std::get<4UL>(tuple));

	// ##################################################### NUMERICAL METHODS FOR ODE INTEGRATION #####################################################

	// ********************************  8-th ORDER ADAPTIVE ADAMS-BASHFORTH-MOULTON STEPPER ********************************
	boost::numeric::odeint::adaptive_adams_bashforth_moulton<8UL, state_type, double, state_type, double,
															 boost::numeric::odeint::vector_space_algebra, boost::numeric::odeint::default_operations, boost::numeric::odeint::initially_resizer>
		abm8_stepper;

	// ********************************  4-th ORDER CLASSIC RUNGE-KUTTA STEPPER ********************************
	// typedef boost::numeric::odeint::runge_kutta4_classic<state_type, double, state_type, double, boost::numeric::odeint::vector_space_algebra,
	// 	boost::numeric::odeint::default_operations, boost::numeric::odeint::initially_resizer> rk4_stepper;

	// ********************************  5-th ORDER CASH-KARP STEPPER ********************************
	// typedef boost::numeric::odeint::runge_kutta_cash_karp54<state_type, double, state_type, double, boost::numeric::odeint::vector_space_algebra,
	// 	boost::numeric::odeint::default_operations, boost::numeric::odeint::initially_resizer> rkk54_stepper;

	// ********************************  5-th ORDER DORMAND-PRINCE RUNGE-KUTTA ********************************
	// typedef boost::numeric::odeint::runge_kutta_dopri5<state_type, double, state_type, double, boost::numeric::odeint::vector_space_algebra> rkd5_stepper;

	// ********************************  BULIRSCH-STOER STEPPER ********************************
	// typedef boost::numeric::odeint::bulirsch_stoer<state_type, double, state_type, double, boost::numeric::odeint::vector_space_algebra,
	// 	boost::numeric::odeint::default_operations, boost::numeric::odeint::initially_resizer> blstr_stepper;

	// ******************************** RUUNGE-KUTTA-FEHLBERG (RKF78) STEPPER ********************************
	// typedef boost::numeric::odeint::runge_kutta_fehlberg78<state_type, double, state_type, double, boost::numeric::odeint::vector_space_algebra,
	// 	boost::numeric::odeint::default_operations, boost::numeric::odeint::initially_resizer> rk78_stepper;

	// #################################################################################################################################################

	// start and end points, in terms of arc-length s, of each CTR segment and initial step-size for integration (ds)
	double s_start, s_end, ds;

	// instantiating variables needed for enforcing the continuity of internal moments across CTR segments
	vec3d u1, u2, u3, u1_next, u1_ast, u2_ast, u3_ast, u1_ast_next, u2_ast_next, u3_ast_next;
	blaze::DiagonalMatrix<blaze::StaticMatrix<double, 3UL, 3UL, blaze::rowMajor>> K1, K2, K3, K1_next, K2_next, K3_next;
	blaze::StaticMatrix<double, 3UL, 3UL> R_theta_2, R_theta_3;
	double d_theta2, d_theta3;

	// initial guess for the BC (twist curvatures)
	blaze::StaticVector<double, 2UL> u1_xy_BC = { initGuess[0UL], initGuess[1UL] };
	blaze::StaticVector<double, 3UL> theta_BC(this->m_theta_0), uz_BC = {initGuess[2UL], initGuess[3UL], initGuess[4UL]};

	// instantiating the vector of initial conditions for solving the state equations (20 x 1)
	state_type y_0;	

	// iterating through the tube segments comprising the CTR
	size_t len_seg = S.size() - 1;
	for (size_t seg = 0; seg < len_seg; ++seg)
	{
		// 1st element of y computes the curvature of the first (innermost) tube along the x direction
		// 2nd element of y computes the curvature of the first (innermost) tube along the y direction
		// next 3 elements of y are the torsional curvatures for the three tubes, e.g., y = [u1_z  u2_z  u3_z]
		// next 3 elements of y are twist angles, theta_i = [theta_1  theta_2  theta_3]
		// last 7 elements are r(position) and h(quaternion-orientations) of the local frame, respectively at each arc-length s

		// initializing the initial conditions vector (15 x 1)
		y_0 = {u1_xy_BC[0UL], u1_xy_BC[1UL],
			   uz_BC[0UL], uz_BC[1UL], uz_BC[2UL],
			   theta_BC[0UL], theta_BC[1UL], theta_BC[2UL],
			   m_r_0[0UL], m_r_0[1UL], m_r_0[2UL],
			   m_h_0[0UL], m_h_0[1UL], m_h_0[2UL], m_h_0[3UL]};

		// specifying the interval of integration (in terms of tube segment arc-lengths)
		s_start = S[seg];
		s_end = S[seg + 1];
		ds = (s_end - s_start) / 20; // 20 points per segment

		// passing the tube parameters in the segment to the state equation method
		this->m_stateEquations->setEquationParameters(blaze::column(U_x, seg), blaze::column(U_y, seg), blaze::column(EI, seg), blaze::column(GJ, seg));

		// ##################################################### NUMERICAL INTEGRATION #####################################################
		// Employs the selected stepper (Numerical method) and integrates the system of ODEs along the segment considered
		boost::numeric::odeint::integrate_adaptive(abm8_stepper, *this->m_stateEquations, y_0, s_start, s_end, ds, *this->m_stateObserver);

		// rate of twist of tube 2 at the end of the i - th segment (y0 is changed to the approximative solution f(x,t') at the end of integration)
		d_theta2 = y_0[3UL] - y_0[2UL];
		// rate of twist of tube 3 at the end of the i - th segment
		d_theta3 = y_0[4UL] - y_0[2UL];

		// setting the boundary(initial conditions) for the next segment
		uz_BC = blaze::subvector<2UL, 3UL>(y_0);
		theta_BC = blaze::subvector<5UL, 3UL>(y_0);
		m_r_0 = blaze::subvector<8UL, 3UL>(y_0);
		m_h_0 = blaze::subvector<11UL, 4UL>(y_0);
		u1 = blaze::subvector<0UL, 3UL>(y_0);

		if (seg < len_seg - 1)
		{
			// stiffness matrices in the current segment
			blaze::diagonal(K1) = { EI(0UL, seg), EI(0UL, seg), GJ(0UL, seg) };
			blaze::diagonal(K2) = { EI(1UL, seg), EI(1UL, seg), GJ(1UL, seg) };
			blaze::diagonal(K3) = { EI(2UL, seg), EI(2UL, seg), GJ(2UL, seg) };
			// precurvature vectors in the current segment
			u1_ast = {U_x(0, seg), U_y(0, seg), 0.0};
			u2_ast = {U_x(1, seg), U_y(1, seg), 0.0};
			u3_ast = {U_x(2, seg), U_y(2, seg), 0.0};
			// stiffness matrices in the next segment
			blaze::diagonal(K1_next) = { EI(0UL, seg + 1), EI(0UL, seg + 1), GJ(0UL, seg + 1) };
			blaze::diagonal(K2_next) = { EI(1UL, seg + 1), EI(1UL, seg + 1), GJ(1UL, seg + 1) };
			blaze::diagonal(K3_next) = { EI(2UL, seg + 1), EI(2UL, seg + 1), GJ(2UL, seg + 1) };
			// precurvature vectors in the next segment
			u1_ast_next = {U_x(0, seg + 1), U_y(0, seg + 1), 0.0};
			u2_ast_next = {U_x(1, seg + 1), U_y(1, seg + 1), 0.0};
			u3_ast_next = {U_x(2, seg + 1), U_y(2, seg + 1), 0.0};

			// orientation of tube 2's local frame at the end of the i-th segment
			R_theta_2 = mathOp::rotz(theta_BC[1UL]);
			// orientation of tube 3's local frame at the end of the i-th segment
			R_theta_3 = mathOp::rotz(theta_BC[2UL]);

			// curvatures of tubes 2 and 3 at the end of the current segment
			u2 = mathOp::transposePreMultiply(R_theta_2, u1) + d_theta2 * m_e3;
			u3 = mathOp::transposePreMultiply(R_theta_3, u1) + d_theta3 * m_e3;

			// estimating curvature of next segment through continuity of internal bending moments
			// * * * # --- # ====>> For reference, check: Computing Jacobians and Compliance Matrices(Rucker, 2011) pg. 5 <<==== # --- # * * *
			u1_next = blaze::inv(K1_next + K2_next + K3_next) * (K1 * (u1 - u1_ast) + R_theta_2 * K2 * (u2 - u2_ast) + R_theta_3 * K3 * (u3 - u3_ast) + K1_next * u1_ast_next + R_theta_2 * K2_next * (u2_ast_next - d_theta2 * m_e3) + R_theta_3 * K3_next * (u3_ast_next - d_theta3 * m_e3));

			// effectively setting the initial condition for u1_dot along the x, y directions
			u1_xy_BC = blaze::subvector<0UL, 2UL>(u1_next);
		}
		else
		{
			// precurvature of the innermost tube at the distal end
			u1_ast = {U_x(0, seg), U_y(0, seg), 0.0};
		}
	}

	//
	//	****************  #####  -------------- ___ DISTAL BOUNDARY CONDITIONS ___ --------------  #####  ****************
	//			1) internal moment at the tip of tube 1 must equal the external moment applied
	//			   Namely, R1K1(u1(L1) - u1 * (L1)) = Wm = > u1(l1) - u1 * (l1)-K1\R1'Wm = 0
	//
	//			2) at the distal ends of the remaining tubes, the axial component of the internal moments must equal zero
	//			   Namely, ui, z(Li) - ui, z* (Li) = 0

	// Residue vector due to infringment of the distal boundary conditions
	blaze::StaticVector<double, 5UL> Residue;

	// Residue = [ u1(0) - U_x(0,end), u2(1) - U_y(0,end), u1(2), 0, 0 ]
	Residue = {u1[0UL] - u1_ast[0UL],
			   u1[1UL] - u1_ast[1UL],
			   u1[2UL],
			   0.0,
			   0.0};

	// lambda function that finds the u_z curvatures at the distal ends of tubes 2 and 3
	auto computeResidue = [&](double distalEnd, size_t index)
	{
		// must use some tolerance when comparing floating points
		auto itt = std::find_if(this->m_s.begin(), this->m_s.end(), [&](double x)
								{ return (std::abs(x - distalEnd) <= 1e-7) ? true : false; }); // finds where tube ends (with a 0.0001mm tolerance)

		auto id = std::distance(this->m_s.begin(), itt);
		// ui_z at the distal end of the i-th tube
		Residue[2 + index] = this->m_y[id][2 + index];
	};

	// Computing the Residues associated to the twist curvatures (rate of twist) at the distal ends of tubes 2 and 3
	vec3d distEnd(this->m_segment->getDistEnd()); // arc-lengths at which the distal ends of the tubes are currently

	computeResidue(distEnd[1UL], 1UL);
	computeResidue(distEnd[2UL], 2UL);

	return Residue;
}

// function that computes the finite-differences Jacobian for solving the BVP
blaze::StaticMatrix<double, 5UL, 5UL> CTR::jac_BVP(const blaze::StaticVector<double, 5UL> &initGuess, const blaze::StaticVector<double, 5UL> &residue)
{
	blaze::StaticMatrix<double, 5UL, 5UL, blaze::columnMajor> jac_bvp;
	//  perturbed initial conditions for computing the partial derivative
	blaze::StaticVector<double, 5UL> initGuessPerturbed(initGuess), residuePerturbed;

	for (size_t iter = 0UL; iter < 5UL; ++iter)
	{
		initGuessPerturbed[iter] += this->m_eps;
		// perturbed residue
		residuePerturbed = this->ODESolver(initGuessPerturbed);
		// building the finite-differences Residue jacobian
		blaze::column(jac_bvp, iter) = (residuePerturbed - residue) / this->m_eps;
		// restoring the original value of the array
		initGuessPerturbed[iter] = initGuess[iter];
	}

	return jac_bvp;
}

// function that computes the finite-differences Jacobian wrt actuation inputs
blaze::StaticMatrix<double, 3UL, 6UL> CTR::jacobian(const blaze::StaticVector<double, 5UL> &initGuess, const blaze::StaticVector<double, 3UL> &tipPos)
{
	blaze::StaticMatrix<double, 3UL, 6UL, blaze::columnMajor> jac;
	blaze::StaticVector<double, 6UL> q_Original(this->m_q), q_Perturbed(this->m_q);

	for (size_t iter = 0UL; iter <= 5UL; ++iter)
	{
		q_Perturbed[iter] += this->m_eps;
		this->setConfiguration(q_Perturbed);

		// recalculates the CTR transition points and segments for beta actuation
		if (iter <= 3UL)
			m_segment->recalculateSegments(this->m_Tubes, this->m_beta);

		this->ODESolver(initGuess);

		// computing the tip position of the perturbed CTR
		blaze::column(jac, iter) = (this->getTipPos() - tipPos) / this->m_eps;
		// restoring the original CTR joint values
		q_Perturbed[iter] = q_Original[iter];
	}
	// sets the joint values to their original, unperturbed configuration values
	this->setConfiguration(q_Original);
	// this->ODESolver(initGuess);

	return jac;
}

// function that implements Powell's Dog Leg Method (Nonlinear root-finding method for solving the BVP)
bool CTR::PowellDogLeg(blaze::StaticVector<double, 5UL> &initGuess)
{
	bool found;
	size_t k = 0UL;
	const size_t k_max = 300UL;
	double alpha, beta, delta, eps1, eps2, rho, c;
	blaze::StaticVector<double, 5UL> g, f, f_new, x_new, h_sd, h_gn, h_dl;
	blaze::StaticMatrix<double, 5UL, 5UL> J;

	// initializing parameters
	delta = 1;
	eps1 = eps2 = 1e-22;

	f = this->ODESolver(initGuess);
	J = this->jac_BVP(initGuess, f);
	g = blaze::trans(J) * f;

	// checking if the initial guess satisfies the BVP without the need of any further refinement
	found = ((blaze::linfNorm(f) <= this->m_accuracy) || (blaze::linfNorm(g) <= eps1)) ? true : false;

	while (!found && (k < k_max))
	{
		k++;

		alpha = blaze::sqrNorm(g) / blaze::sqrNorm(J * g);
		h_sd = -alpha * g;			 // steepest descend (this is a direction, not a step!)
		h_gn = -mathOp::pInv(J) * f; // Gauss-Newton step (Least Square solution)

		// two candidates for the step to take from this point, a = alpha*h_sd & b = h_gn

		// computing the dog leg direction
		if (blaze::norm(h_gn) <= delta)
			h_dl = h_gn;
		else
		{
			if (blaze::norm(h_sd) >= delta)
				h_dl = delta * blaze::normalize(h_sd);
			else
			{
				c = blaze::trans(h_sd) * (h_gn - h_sd);

				if (c <= 0.0)
				{
					beta = (-c + sqrt(c * c + blaze::sqrNorm(h_gn - h_sd) * (delta * delta - blaze::sqrNorm(h_sd)))) / blaze::sqrNorm(h_gn - h_sd);
				}
				else
				{
					beta = (delta * delta - blaze::sqrNorm(h_sd)) / (c + sqrt(c * c + blaze::sqrNorm(h_gn - h_sd) * (delta * delta - blaze::sqrNorm(h_sd))));
				}

				h_dl = h_sd + beta * (h_gn - h_sd); // Dog Leg step
			}
		}

		if (blaze::norm(h_dl) <= eps2 * (blaze::norm(initGuess) + eps2))
			found = true;
		else
		{
			x_new = initGuess + h_dl;
			f_new = this->ODESolver(x_new);
			rho = (blaze::sqrNorm(f) - blaze::sqrNorm(f_new)) / (0.5 * blaze::trans(h_dl) * ((delta * h_dl) - g));

			if (rho > 0.0)
			{
				initGuess = std::move(x_new);
				f = std::move(f_new);
				J = this->jac_BVP(initGuess, f);
				g = blaze::trans(J) * f;

				if ((blaze::linfNorm(f) <= this->m_accuracy) || (blaze::linfNorm(g) <= eps1))
					found = true;
			}

			if (rho > 0.75)
				delta = std::max(delta, 3 * blaze::norm(h_dl));
			else
			{
				if (rho < 0.25)
					delta *= 0.5;
			}

			if (delta < eps2 * (blaze::norm(initGuess) + eps2))
				found = true;
		}
	}

	return found;
}

// function that implements the Levenberg-Marquardt Method (Nonlinear root-finding method for solving the BVP)
bool CTR::Levenberg_Marquardt(blaze::StaticVector<double, 5UL> &initGuess)
{
	size_t k = 0UL;
	const size_t k_max = 300UL;
	blaze::StaticVector<double, 5UL> h, g, f, f_new;
	blaze::StaticMatrix<double, 5UL, 5UL> J, A;
	blaze::IdentityMatrix<double> I(5UL);
	double rho, nu = 2, mu, tau = 1e-3, e1 = 1e-18, e2 = 1e-25;
	bool found;

	// computing the residue and residue Jacobian associated to initGuess
	f = this->ODESolver(initGuess);
	J = this->jac_BVP(initGuess, f);
	A = blaze::trans(J) * J;
	g = blaze::trans(J) * f;
	found = (blaze::linfNorm(g) <= e1) ? true : false;
	mu = tau * blaze::max(blaze::diagonal(A));

	// starting the iterative minimization loop
	while ((!found) && (k < k_max))
	{
		k++;
		blaze::solve(blaze::declsym(A + (mu * I)), h, -g);

		f_new = this->ODESolver(initGuess + h);
		rho = (blaze::sqrNorm(f) - blaze::sqrNorm(f_new)) / (0.5 * blaze::trans(h) * ((mu * h) - g));

		if (rho > 0.0)
		{
			// accept the decrease in the function
			initGuess += h;
			// computing the residue Jacobian at the new initial guess
			J = this->jac_BVP(initGuess, f_new);
			A = blaze::trans(J) * J;
			f = std::move(f_new);
			g = blaze::trans(J) * f;
			found = (blaze::linfNorm(g) <= e1) ? true : false;
			mu = mu * std::max(0.33333333, 1 - blaze::pow(2 * rho - 1, 3));
			nu = 2;
		}
		else
		{
			mu = mu * nu;
			nu = 2 * nu;
		}

		// checking if the tolerance has been satisfied
		if (blaze::linfNorm(f) <= this->m_accuracy)
			found = true;
	}

	return found;
}

// function that implements the Broyden (Nonlinear root-finding method for solving the BVP)
bool CTR::Broyden(blaze::StaticVector<double, 5UL> &initGuess)
{
	// found: returns true (false) when the root-finding method converges (does not converge) within k_max iterations
	bool found;

	// initial Hessian matrix --> computed via finite differences
	blaze::StaticMatrix<double, 5UL, 5UL> JacInv, JacInvNew;

	// setting up and starting my handmadeBFGS method
	blaze::StaticVector<double, 5UL> F, Fold, X, Xold, deltaX, deltaF; // staticVectors are automatically initialized to 0

	// Residue yielded by the initial guess for the CTR BVP
	F = this->ODESolver(initGuess); // F(x_k)	: residue
	X = std::move(initGuess);		// x_k		: initial guess
	JacInvNew = JacInv = mathOp::pInv(this->jac_BVP(X, F));

	// checking if the initial guess already satisfies the BVP
	found = (blaze::linfNorm(F) <= this->m_accuracy) ? true : false;

	size_t k = 0UL;
	const size_t k_max = 300UL;
	while (!found && (k < k_max))
	{
		k++;

		deltaX = X - Xold; // dX := x_k - x_k-1
		deltaF = F - Fold; // dF := F(x_k) - F(x_k-1)

		JacInv = std::move(JacInvNew);
		if ((blaze::norm(deltaX) > 0.0) && (blaze::norm(deltaF) > 0.0))
			JacInvNew = JacInv + ((deltaX - JacInv * deltaF) / (blaze::trans(deltaX) * JacInv * deltaF)) * blaze::trans(deltaX) * JacInv;
		else
			JacInvNew = JacInv;

		Xold = std::move(X);
		Fold = std::move(F);

		// update the initial guess
		X = Xold - JacInv * F;
		F = this->ODESolver(X);

		while (blaze::isnan(F))
		{
			if (blaze::isnan(X))
				X = 0.0;
			else
				X /= blaze::max(blaze::abs(X));

			F = this->ODESolver(X);
			JacInv = JacInvNew = mathOp::pInv(this->jac_BVP(X, F));
			Xold = std::move(X);
			X = Xold - JacInv * F;
		}

		if (k % 10 == 0.0)
		{
			JacInv = JacInvNew = mathOp::pInv(this->jac_BVP(X, F));
			X = Xold - JacInv * F;
		}

		if (blaze::linfNorm(F) <= this->m_accuracy)
			found = true;
	}

	initGuess = std::move(X);
	return found;
}

// function that implements Broyden's Nonlinear root-finding method for solving the BVP
bool CTR::Broyden_II(blaze::StaticVector<double, 5UL> &initGuess)
{
	bool found;

	// initial Hessian matrix --> computed via finite differences
	blaze::StaticMatrix<double, 5UL, 5UL> Jac, JacNew;

	// setting up and starting my handmadeBFGS method
	blaze::StaticVector<double, 5UL> F, Fold, X, Xold, deltaX, deltaF; // staticVectors are automatically initialized to 0

	// Residue yielded by the initial guess for the CTR BVP
	F = this->ODESolver(initGuess); // F(x_k)	: residue
	X = std::move(initGuess);		// x_k		: initial guess
	JacNew = this->jac_BVP(initGuess, F);

	// checking if the initial guess already satisfies the BVP
	found = (blaze::linfNorm(F) <= this->m_accuracy) ? true : false;

	size_t k = 0UL;
	const size_t k_max = 300UL;
	while (!found && (k < k_max))
	{
		k++;

		deltaX = X - Xold; // dX := x_k - x_k-1
		deltaF = F - Fold; // dF := F(x_k) - F(x_k-1)

		Jac = std::move(JacNew);
		if (blaze::sqrNorm(deltaX) > 0.0)
			JacNew = Jac + blaze::sqrNorm(deltaX) * (deltaF - (Jac * deltaX)) * blaze::trans(deltaX);
		else
			JacNew = Jac;

		Xold = std::move(X);
		Fold = std::move(F);

		// update the initial guess
		X = Xold - mathOp::pInv(Jac) * F;
		F = this->ODESolver(X);

		while (blaze::isnan(F))
		{
			if (blaze::isnan(X))
				X = 0.0;
			else
				X /= blaze::max(blaze::abs(X));

			F = this->ODESolver(X);
			JacNew = this->jac_BVP(X, F);
			Xold = std::move(X);
			X = Xold - mathOp::pInv(Jac) * F;
		}

		if (k % 10 == 0.0)
		{
			JacNew = this->jac_BVP(X, F);
			X = Xold - mathOp::pInv(Jac) * F;
		}

		if (blaze::linfNorm(F) <= this->m_accuracy)
			found = true;
	}

	initGuess = std::move(X);
	return found;
}

// function that implements the Newton-Raphson method (Nonlinear root-finding method for solving the BVP)
bool CTR::Newton_Raphson(blaze::StaticVector<double, 5UL> &initGuess)
{
	bool found;
	// setting up and starting my handmade Newton-Raphson method
	blaze::StaticVector<double, 5UL> Residue, Residue_new, d_Residue, dGuess; // staticVectors are automatically initialized to 0

	// Residue of the unperturbed initial guess for the CTR
	Residue = this->ODESolver(initGuess);

	found = (blaze::linfNorm(Residue) <= this->m_accuracy) ? true : false;

	//  Weighing matrices for adjusting the initial guess iteratively (Implementing a PD regulator)
	blaze::DiagonalMatrix<blaze::StaticMatrix<double, 5UL, 5UL, blaze::rowMajor>> Kp, Kd;
	blaze::StaticMatrix<double, 5UL, 5UL, blaze::columnMajor> jac_bvp;
	blaze::diagonal(Kp) = 0.3;	// 0.45 | 0.6  | 0.3
	blaze::diagonal(Kd) = 2e-3; // 3e-3 | 5e-3 | 2e-3

	size_t k = 0UL;
	const size_t k_max = 300UL;

	// starting iterations for adjusting the initial guess "u_guess ~ initGuess"
	while (!found && (k < k_max))
	{
		k++;
		jac_bvp = this->jac_BVP(initGuess, Residue);
		// error equation(globally asymptotically stable)
		dGuess = mathOp::pInv(jac_bvp) * (Kp * Residue + Kd * d_Residue);
		// updating the initial guess(weighted negative gradient of the cost function)
		initGuess -= dGuess;
		// computing the new cost associated to the newly readjusted initial guess
		Residue_new = this->ODESolver(initGuess);

		// Checking if the Jacobian has large elements beyond machine precision
		while (blaze::isnan(Residue_new))
		{
			if (blaze::isnan(initGuess))
				initGuess = 0.0;
			else
				initGuess /= blaze::max(blaze::abs(initGuess));

			Residue_new = this->ODESolver(initGuess);
			d_Residue = Residue_new - Residue;
			Residue = std::move(Residue_new);
			continue;
		}

		// cost variation due to initial guess refinement
		d_Residue = Residue_new - Residue;
		// updating the cost
		Residue = std::move(Residue_new);

		if (blaze::linfNorm(Residue) <= this->m_accuracy)
			found = true;
	}

	return found;
}

// function that implements the Modified, globally convergent Newton-Raphson method (Nonlinear root-finding method for solving the BVP)
bool CTR::Modified_Newton_Raphson(blaze::StaticVector<double, 5UL> &initGuess)
{
	/*
		Algorithm extracted from page 309 of Introduction to Numerical Analysis 3rd edition by Josef Stoer & Roland Bulirsch
	*/

	bool found;
	// computes the residue associated to the initial guess
	blaze::StaticVector<double, 5UL> f(this->ODESolver(initGuess)), d;
	blaze::StaticVector<double, 5UL, blaze::rowVector> Dh;
	blaze::StaticMatrix<double, 5UL, 5UL> D, D_inv;
	double h, h_0, lambda, gamma, improvementFactor, d_norm, Dh_norm;
	size_t j = 0UL, k = 0UL;
	const size_t k_max = 500UL;
	std::vector<double> h_k; // vector to store all h_k's
	h_k.reserve(k_max);
	std::vector<double>::iterator result; // iterator for determining the min element in the vector h_k

	found = (blaze::linfNorm(f) <= this->m_accuracy) ? true : false;

	while (!found && (k < k_max))
	{
		k++;
		// computing the residue Jacobian
		D = this->jac_BVP(initGuess, f);
		D_inv = mathOp::pInv(D);

		// search direction (directional derivative)
		d = D_inv * f;
		gamma = 1 / (blaze::norm(D_inv) * blaze::norm(D)); // gamma := 1/cond(Df)
		h_0 = blaze::sqrNorm(f);						   // h := f'f
		// Dh := D(f'f) = 2f'Df
		Dh = 2 * blaze::trans(f) * D;
		d_norm = blaze::norm(d);
		Dh_norm = blaze::norm(Dh);

		while (true)
		{
			f = this->ODESolver(initGuess - blaze::pow(0.5, j) * d);
			// std::cout << "Modified_Newton_Raphson -- j = : " << j << " | residue = " << blaze::trans(f);
			while (blaze::isnan(f))
			{
				j++;
				f = this->ODESolver(initGuess - blaze::pow(0.5, j) * d);
			}
			h = blaze::sqrNorm(f);
			improvementFactor = blaze::pow(0.5, j) * 0.25 * gamma * d_norm * Dh_norm;
			// storig the value of h_k to determine step size posteriorly
			h_k.push_back(h);

			if (h <= (h_0 - improvementFactor))
				break;
			else
				j++;
		}

		// determining the value of the step-size lambda
		result = std::min_element(h_k.begin(), h_k.end());
		// retrieving the minimum h_k
		lambda = blaze::pow(0.5, std::distance(h_k.begin(), result));
		initGuess -= lambda * d;
		h_k.clear();

		// resets the exponent variable j
		j = 0UL;

		// compute the residue associated to the newly refined initGuess
		f = this->ODESolver(initGuess);

		// checking the terminating condition
		if (blaze::linfNorm(f) <= this->m_accuracy)
		{
			found = true;
		}
		else
		{
			if (blaze::isnan(f))
			{
				if (blaze::isnan(initGuess))
					initGuess = 0.0;
				else
					initGuess /= blaze::max(blaze::abs(initGuess));

				// recompute the residue of the readjusted initial guess
				f = this->ODESolver(initGuess);
			}
		}
	}

	if (!found)
	{
		initGuess[0UL] = initGuess[1UL] = initGuess[4UL] = 0.0;
		found = this->PowellDogLeg(initGuess);

		if (!found)
		{
			initGuess[0UL] = initGuess[1UL] = initGuess[4UL] = 0.0;
			found = this->Newton_Raphson(initGuess);
		}

		if (!found)
		{
			initGuess[0UL] = initGuess[1UL] = initGuess[4UL] = 0.0;
			found = this->Levenberg_Marquardt(initGuess);
		}
	}

	return found;
}

// function that implements the CTR actuation for any inputs joint values q
bool CTR::actuate_CTR(blaze::StaticVector<double, 5UL> &initGuess, const blaze::StaticVector<double, 6UL> &q_input)
{
	// boolean flag for indicating convergence (1: zero found | 0: zero not found)
	bool found = false;

	// updating the CTR joints for desired input values
	this->setConfiguration(q_input);

	// recalculates the CTR transition points and segments
	m_segment->recalculateSegments(this->m_Tubes, this->m_beta);

	// initial guess for proximal boundary condition--[u_x(0) u_y(0) u1_z(0) u2_z(0) u3_z(0)]
	switch (this->m_method)
	{
	case mathOp::rootFindingMethod::NEWTON_RAPHSON:
		found = this->Newton_Raphson(initGuess);
		break;
	case mathOp::rootFindingMethod::LEVENBERG_MARQUARDT:
		found = this->Levenberg_Marquardt(initGuess);
		break;
	case mathOp::rootFindingMethod::POWELL_DOG_LEG:
		found = this->PowellDogLeg(initGuess);
		break;
	case mathOp::rootFindingMethod::MODIFIED_NEWTON_RAPHSON:
		found = this->Modified_Newton_Raphson(initGuess);
		break;
	case mathOp::rootFindingMethod::BROYDEN:
		found = this->Broyden(initGuess);
		break;
	case mathOp::rootFindingMethod::BROYDEN_II:
		found = this->Broyden_II(initGuess);
		break;
	}

	return found;
}

// function that implements the position control ==> returns [u_0, Jac, q_min, timeout]
std::tuple<blaze::StaticMatrix<double, 3UL, 6UL>, blaze::StaticVector<double, 6UL>, bool> CTR::posCTRL(blaze::StaticVector<double, 5UL> &initGuess, const vec3d &target, const double posTol)
{
	double t = 0.5;											 		// time base
	double minError = 1e3;										 		// minimum distance to target
	bool status;												 		// status = TRUE (FALSE) indicates convergence (lack thereof)
	blaze::StaticMatrix<double, 3UL, 6UL, blaze::columnMajor> J; 		// Jacobian matrix
	blaze::StaticMatrix<double, 6UL, 3UL, blaze::columnMajor> J_inv;	// Jacobian pseudoinverse
	blaze::IdentityMatrix<double, blaze::rowMajor> I(6UL); 				// 6 x 6 Identity matrix
	std::ofstream file_IK;

	file_IK.open("IK_data.dat", std::ios::out | std::ios::trunc);

	if (file_IK.is_open())
		file_IK << "Format: joint values q = (beta_1[meters], beta_2[meters], beta_3[meters], alpha_1[rad], alpha_2[rad], alpha_3[rad]) | tip position p = (x[meters], y[meters], z[meters])\n";
	else
		std::cout << "The IK data file failed to open!" << std::endl;

	// proportional gain for control
	blaze::DiagonalMatrix<blaze::StaticMatrix<double, 3UL, 3UL>> Kp, Kd, Ki;
	blaze::diagonal(Kp) = 1.0;
	blaze::diagonal(Ki) = 0.05;
	blaze::diagonal(Kd) = 0.5;

	// Capturing the CTR's current joint configuration
	blaze::StaticVector<double, 6UL> dqdt, q_min(this->m_q), q(this->m_q);
	// Capturing the proximal BC for the minimum distance
	blaze::StaticVector<double, 5UL> initGuessMin(initGuess);
	// Calculate the CTR Jacobian in the present configuration and retrieves convergence status
	status = this->actuate_CTR(initGuess, q);

	if (!status)
		return std::make_tuple(J, q, status);

	// Capturing the position of the end-effector
	vec3d x_CTR, tipError, d_tipError, last_tipError, int_tipError;

	x_CTR = this->getTipPos();
	// Current position error
	int_tipError = last_tipError = tipError = target - x_CTR;

	// Euclidean distance to target
	double dist2Tgt = blaze::norm(tipError);

	if (dist2Tgt < minError)
	{
		minError = dist2Tgt;
		q_min = q;

		if (dist2Tgt <= posTol)
			return std::make_tuple(J, q_min, status);
	}

	// function to implement actuators sigularity avoidance
	blaze::StaticVector<double, 6UL> f;
	// clearance between linear actuators
	double Clr = 5e-3, deltaBar = 0.0;
	// lengths of straight sections of the CTR tubes
	vec3d ls, L;
	L = {this->m_Tubes[0UL]->getTubeLength(), this->m_Tubes[1UL]->getTubeLength(), this->m_Tubes[2UL]->getTubeLength()};
	ls = {this->m_Tubes[0UL]->getStraightLen(), this->m_Tubes[1UL]->getStraightLen(), this->m_Tubes[2UL]->getStraightLen()};
	// lower and upper bounds on prismatic joint limits
	vec3d betaMax, betaMin;

	size_t N_itr = 0UL;				// iterations counter
	const size_t maxIter = 3500UL;  // maximum admissible number of iterations in the position control loop

	// parameters for local optimization (joint limits avoidance)
	double ke = 15;

	// pointers to the revolute joints of the CTR
	auto revJoints = blaze::subvector<3UL, 3UL>(q);

	auto f1 = blaze::subvector<0UL, 3UL>(f);

	// position control loop
	while ((dist2Tgt > posTol) && (N_itr < maxIter))
	{
		// incrementing the number of iterations
		N_itr++;

		// compute the Jacobian in the present configuration
		J = this->jacobian(initGuess, x_CTR);
		// Pseudo-inverse of Jacobian for resolving CTR joint motion rates
		J_inv = mathOp::pInv(J);

		// Nullspace control (collision and actuation limits)
		betaMin = {std::max(-ls[0UL] + deltaBar, L[1UL] + this->m_beta[1UL] - L[0UL] + deltaBar),
				   std::max({-ls[1UL] + deltaBar, L[2UL] + this->m_beta[2UL] - L[1UL] + deltaBar, this->m_beta[0UL] + Clr}),
				   std::max(-ls[2UL] + deltaBar, this->m_beta[1UL] + Clr)};

		betaMax = {m_beta[1UL] - Clr,
				   std::min(L[0UL] + this->m_beta[0UL] - deltaBar - L[1UL], this->m_beta[2UL] - Clr),
				   std::min(L[1UL] + this->m_beta[1UL] - deltaBar - L[2UL], -deltaBar)};

		// penalty function for local optimization (actuator collision avoidance)
		// Had to add an infinitesimal (1e-10) to the denominador (betaMax - betaMin + 1e-12) to avoid blow-up when the actuator interval degenerates to a point
		// f1 = blaze::pow(blaze::abs((betaMax + betaMin - 2 * this->m_beta) / (betaMax - betaMin + 1e-10)), ke) * blaze::sign(this->m_beta - (betaMax + betaMin) * 0.5);

		// Inverse kinematics
		f1 = 1e-4 * blaze::pow(blaze::abs((betaMax + betaMin - 2 * this->m_beta) / (betaMax - betaMin + 1e-10)), ke) * blaze::sign(this->m_beta - (betaMax + betaMin) * 0.5);
		// resolved rates -- Nullspacec local optimization (joint limit avoidance)
		dqdt = (J_inv * (Kp * tipError + Kd * d_tipError + Ki * int_tipError) + (I - blaze::trans(J_inv * J)) * (-f)) * t;
		// dqdt = J_inv * Kp * tipError * t;

		auto q_dot = [&] { // rescaling joint variables for limit avoidance
			for (size_t i = 0UL; i < 3UL; ++i)
			{
				if (this->m_beta[i] + dqdt[i] > betaMax[i])
					dqdt[i] = (betaMax[i] - this->m_beta[i]) * 0.5;

				if (this->m_beta[i] + dqdt[i] < betaMin[i])
					dqdt[i] = (betaMin[i] - this->m_beta[i]) * 0.5;
			}
			return dqdt;
		};

		dqdt = q_dot();

		// updating the CTR joints->q: [beta, theta]
		q += dqdt;

		// wrapping the actuation angles to the [-Pi,Pi) interval
		revJoints = blaze::map(revJoints, [](double theta)
							   {
				static constexpr double TWO_PI = 2 * M_PI;
				double wrappedAngle;

				// wrappedAngle = fmod(theta, TWO_PI);
				wrappedAngle = remainder(theta, TWO_PI);

				return wrappedAngle; });

		// actuate the CTR to new configuration and retrieve execution timeout status
		status = this->actuate_CTR(initGuess, q);

		// interrupts the loop execution if actuation fails
		if (!status)
		{
			initGuess[0UL] = initGuess[1UL] = initGuess[4UL] = 0.0;
			status = this->actuate_CTR(initGuess, q);

			if (!status)
			{
				/*std::cerr << "############################################# posCTR() #############################################" << std::endl;
				std::cerr << "		posCTRL() ==> Nonlinear root-finder did not converge when solving the BVP problem!" << std::endl;
				std::cerr << "############################################# posCTR() #############################################" << std::endl << std::endl;*/
				this->actuate_CTR(initGuessMin, q_min);
				return std::make_tuple(J, q_min, status);
			}
		}

		// CTR tip position after control adjustment
		x_CTR = this->getTipPos();

		// writes the joint values and tip position to the data file
		file_IK << "q = (" << q[0UL] << ", " << q[1UL] << ", " << q[2UL] << ", " << q[3UL] << ", " << q[4UL] << ", " << q[5UL] << ") \t | \t p = (" << x_CTR[0UL] << ", " << x_CTR[1UL] << ", " << x_CTR[2UL] << ")\n";

		tipError = target - x_CTR;
		d_tipError = tipError - last_tipError;
		last_tipError = tipError;
		int_tipError += tipError;
		dist2Tgt = blaze::norm(tipError);

		if (dist2Tgt < minError)
		{
			minError = dist2Tgt;
			q_min = q;
			initGuessMin = initGuess;
		}

		// stops the control loop when the position update becomes significantly small
		if (blaze::norm(dqdt) <= 1e-6)
		{
			std::cout << "Exited out of position control loop due small incremental threshold!" << std::endl;
			return std::make_tuple(J, q_min, status);
		}

		// std::cout << "Tip error after adustment: " << dist2Tgt << " Tip Position: " << blaze::trans(x_CTR);
	}

	// Actuating the CTR to the configuration which yields the minimum position error
	initGuess = std::move(initGuessMin);
	this->actuate_CTR(initGuess, q_min);

	// std::cout << "CTR Inverse Kinematics ended in " << N_itr << " iterations. Minimum error: " << minError << std::endl;

	return std::make_tuple(J, q_min, status);
}

// function that returns the Vector of tubes comprising the CTR
std::array<std::shared_ptr<Tube>, 3UL> CTR::getTubes()
{
	return this->m_Tubes;
}

// function that returns the current linear joint values of the CTR
vec3d CTR::getBeta()
{
	return this->m_beta;
}

// function that returns the current joint values of the CTR
blaze::StaticVector<double, 6UL> CTR::getConfiguration()
{
	return this->m_q;
}

// function that returns the position of the CTR tip
vec3d CTR::getTipPos()
{
	vec3d pos;
	if (!this->m_y.empty())
		pos = blaze::subvector<8UL, 3UL>(this->m_y.back());

	return pos;
}

// function that returns the arc-lenghts at each tube's distal end
vec3d CTR::getDistalEnds()
{
	return this->m_segment->getDistalEnds();
}

// function that returns the individual tube shapes
std::tuple<blaze::HybridMatrix<double, 3UL, 200UL, blaze::columnMajor>, blaze::HybridMatrix<double, 3UL, 200UL, blaze::columnMajor>, blaze::HybridMatrix<double, 3UL, 200UL, blaze::columnMajor>> CTR::getTubeShapes()
{
	blaze::HybridMatrix<double, 3UL, 200UL, blaze::columnMajor> Tb_1, Tb_2, Tb_3;
	// arc-lengths at the distal ends of each tube
	vec3d distal_idx, distalEnds(this->m_segment->getDistalEnds());

	// lambda retrieves the shape of the ith tube, transforms to patient coord frame and rescales to mm
	auto tubeShape = [&](size_t tube_index, blaze::HybridMatrix<double, 3UL, 200UL, blaze::columnMajor> &M)
	{
		blaze::StaticVector<double, 4UL> v;
		double element;
		size_t numOfPoints;

		// find the index in the arc-length vector at which each tube ends
		element = distalEnds[tube_index];
		std::vector<double>::iterator it = std::find_if(this->m_s.begin(), this->m_s.end(), [&element](double x)
														{ return (std::abs(x - element) <= 1e-7) ? true : false; }); // finds where tube ends (0.0001mm tolerance)

		numOfPoints = std::distance(this->m_s.begin(), it);

		// resizes the columns of the hybrid matrix M to the exact number of points
		M.resize(3UL, numOfPoints);

		for (size_t col = 0UL; col < numOfPoints; ++col)
		{
			v = {this->m_y[col][8UL], this->m_y[col][9UL], this->m_y[col][10UL], 1.0};
			blaze::column(M, col) = blaze::subvector<0UL, 3UL>(v);
		}
	};

	tubeShape(0UL, Tb_1);
	tubeShape(1UL, Tb_2);
	tubeShape(2UL, Tb_3);

	// returns the tuple containing the shape of the tubes
	return std::make_tuple(Tb_1, Tb_2, Tb_3);
}

// function that returns a vector with the CTR shape
std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> CTR::getShape()
{
	std::vector<double> r_x, r_y, r_z;
	r_x.reserve(300UL);
	r_y.reserve(300UL);
	r_y.reserve(300UL);

	if (this->m_y.size() > 0UL)
	{
		for (size_t i = 0UL; i < this->m_y.size(); ++i)
		{
			r_x.push_back(this->m_y[i][8UL]);
			r_y.push_back(this->m_y[i][9UL]);
			r_z.push_back(this->m_y[i][10UL]);
		}
	}

	/*printf("s = \n");
	for (double d : m_s)
		printf("%lf, ", d);
	printf("\n\n");*/

	return std::make_tuple(r_x, r_y, r_z);
}

// setter method for setting the actuation joint values (without actuating the CTR) <--> used for computing the Jacobian
void CTR::setConfiguration(const blaze::StaticVector<double, 6UL> &q)
{
	this->m_q = q;
	this->m_beta = blaze::subvector<0UL, 3UL>(this->m_q);
}

// function that sets which method to use for solving the BVP
void CTR::setBVPMethod(mathOp::rootFindingMethod mthd)
{
	this->m_method = mthd;
}