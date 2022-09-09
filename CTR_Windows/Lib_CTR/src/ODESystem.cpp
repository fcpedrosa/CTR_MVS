// This is a personal academic project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: http://www.viva64.com
#include "include/ODESystem.hpp"

// move constructor
ODESystem::ODESystem(ODESystem &&rhs) noexcept
{
	// handling self assignment
	if (this != &rhs)
	{
		this->m_u_ast_x = std::move(rhs.m_u_ast_x);
		this->m_u_ast_y = std::move(rhs.m_u_ast_y);
		this->m_EI = std::move(rhs.m_EI);
		this->m_GJ = std::move(rhs.m_GJ);
		this->m_e3 = std::move(rhs.m_e3);
	}
}

// ODESystem destructor
ODESystem::~ODESystem()
{
	// nothing to be done
}

// Copy assignment operator
ODESystem &ODESystem::operator=(const ODESystem &rhs)
{
	// handling self assignment
	if (this != &rhs)
	{
		this->m_u_ast_x = rhs.m_u_ast_x;
		this->m_u_ast_y = rhs.m_u_ast_y;
		this->m_EI = rhs.m_EI;
		this->m_GJ = rhs.m_GJ;
		this->m_e3 = rhs.m_e3;
	}
	return *this;
}

// move assignment operator
ODESystem &ODESystem::operator=(ODESystem &&rhs) noexcept
{
	// handling self assignment
	if (this != &rhs)
	{
		this->m_u_ast_x = std::move(rhs.m_u_ast_x);
		this->m_u_ast_y = std::move(rhs.m_u_ast_y);
		this->m_EI = std::move(rhs.m_EI);
		this->m_GJ = std::move(rhs.m_GJ);
		this->m_e3 = std::move(rhs.m_e3);
	}
	return *this;
}

// functor that implements the system of ODEs governing a three-tube CTR
void ODESystem::operator()(const state_type &y, state_type &dyds, const double s)
{
	// 1st element of y computes the curvature of the first (innermost) tube along the x direction
	// 2nd element of y computes the curvature of the first (innermost) tube along the y direction
	// next 3 elements of y are the torsional curvatures for the three tubes, e.g., y = [u1_z  u2_z  u3_z]
	// next 3 elements of y are twist angles, theta_i = [theta_1  theta_2  theta_3]
	// last 7 elements are r(position) and h(quaternion-orientations) of the local frame, respectively at each arc-length s

	double theta_2 = y[6UL], dtheta_2 = y[3UL] - y[2UL];
	double theta_3 = y[7UL], dtheta_3 = y[4UL] - y[2UL];

	// implementing curvature equation u_i = transpose(R_z(theta_i))*u_1 + \dot{theta_i}*e3
	vec3d u1, u2, u3;
	u1 = blaze::subvector<0UL, 3UL>(y);
	u2 = mathOp::transposePreMultiply(mathOp::rotz(theta_2), u1) + dtheta_2 * m_e3;
	u3 = mathOp::transposePreMultiply(mathOp::rotz(theta_3), u1) + dtheta_3 * m_e3;

	// estimating the twist curvatures (uz_i) and twist angles (theta_i)
	auto computeTwists = [&](size_t idx, const blaze::StaticVector<double, 3UL>& u) {
		if (m_GJ[idx] != 0.0)
		{
			// uz_i = ( (E_i * I_i) / (G_i * J_i) ) * (ux_i * uy_ast - uy_i * ux_ast)
			dyds[2 + idx] = (m_EI[idx] / m_GJ[idx]) * (u[0UL] * m_u_ast_y[idx] - u[1UL] * m_u_ast_x[idx]);
			// dtheta_i = uz_i - uz_1
			dyds[5 + idx] = u[2UL] - u1[2UL];
		}
		else
			dyds[2 + idx] = dyds[5 + idx] = 0.0;
	};

	computeTwists(0UL, u1);
	computeTwists(1UL, u2);
	computeTwists(2UL, u3);

	// estimating curvature of the first tube along the x and y directions
	blaze::DiagonalMatrix<blaze::StaticMatrix<double, 3UL, 3UL, blaze::rowMajor>> K_inv, Ki;
	K_inv(0UL, 0UL) = K_inv(1UL, 1UL) = 1 / blaze::sum(m_EI);
	K_inv(2UL, 2UL) = 1 / blaze::sum(m_GJ);

	// du_i = R(theta_i) * (K_i * d_Theta_i * d_R(theta_i)' * u_1 + u_i^ * K_i * (u_i - u_i_ast)
	vec3d du_1, du_2, du_3, ui_ast, Du;

	// partial component due to tube 1
	Ki(0UL, 0UL) = Ki(1UL, 1UL) = m_EI[0UL];
	Ki(2UL, 2UL) = m_GJ[0UL];
	ui_ast = {m_u_ast_x[0UL], m_u_ast_y[0UL], 0.0};
	du_1 = mathOp::rotz(y[5UL]) * mathOp::hatPreMultiply(u1, Ki) * (u1 - ui_ast);

	// partial component due to tube 2
	Ki(0UL, 0UL) = Ki(1UL, 1UL) = m_EI[1UL];
	Ki(2UL, 2UL) = m_GJ[1UL];
	ui_ast = {m_u_ast_x[1UL], m_u_ast_y[1UL], 0.0};
	du_2 = mathOp::rotz(y[6UL]) * (Ki * dyds[6UL] * mathOp::rotz_dot_transpose(y[6UL]) * u1 + mathOp::hatPreMultiply(u2, Ki) * (u2 - ui_ast));

	// partial component due to tube 3
	Ki(0UL, 0UL) = Ki(1UL, 1UL) = m_EI[2UL];
	Ki(2UL, 2UL) = m_GJ[2UL];
	ui_ast = {m_u_ast_x[2UL], m_u_ast_y[2UL], 0.0};
	du_3 = mathOp::rotz(y[7UL]) * (Ki * dyds[7UL] * mathOp::rotz_dot_transpose(y[7UL]) * u1 + mathOp::hatPreMultiply(u3, Ki) * (u3 - ui_ast));

	// Equilibrium bending curvature along x, y
	Du = - K_inv * (du_1 + du_2 + du_3);

	// curvature of tube 1 along the x and y directions
	blaze::subvector<0UL, 2UL>(dyds) = blaze::subvector<0UL, 2UL>(Du);

	// spatial derivative of the quaternion representation h_dot
	blaze::subvector<11UL, 4UL>(dyds) = mathOp::quaternionDiff(u1, blaze::subvector<11UL, 4UL>(y));

	// R (orientation) of the local frame at arc-length s
	blaze::StaticMatrix<double, 3UL, 3UL> R;
	mathOp::getSO3(blaze::subvector<11UL, 4UL>(y), R);

	// calculating r_dot
	blaze::subvector<8UL, 3UL>(dyds) = blaze::column<2UL>(R); // r_dot = R * e3
}

// setter method for updating the parameters forward kinematics computation
void ODESystem::setEquationParameters(const vec3d &u_ast_x, const vec3d &u_ast_y, const vec3d &EI, const vec3d &GJ)
{
	this->m_u_ast_x = std::move(u_ast_x);
	this->m_u_ast_y = std::move(u_ast_y);
	this->m_EI = std::move(EI);
	this->m_GJ = std::move(GJ);
}