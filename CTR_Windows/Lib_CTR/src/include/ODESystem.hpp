#pragma once

#include <blaze/Math.h>
#include "mathOperations.hpp"

typedef blaze::StaticVector<double, 3UL> vec3d;
typedef blaze::StaticVector<double, 15UL> state_type;

class ODESystem
{
private:
	vec3d m_u_ast_x, m_u_ast_y, m_EI, m_GJ, m_e3;
	// matrices for the derivative propagation approach

public:
	// default constructor
	ODESystem() : m_u_ast_x(0.0), m_u_ast_y(0.0), m_EI(0.0), m_GJ(0.0)
	{
		m_e3 = {0, 0, 1};
	}

	// overloaded constructor
	ODESystem(const vec3d &u_ast_x, const vec3d &u_ast_y, const vec3d &EI, const vec3d &GJ) : m_u_ast_x(u_ast_x), m_u_ast_y(u_ast_y), m_EI(EI), m_GJ(GJ)
	{
		m_e3 = {0, 0, 1};
	}

	// copy constructor
	ODESystem(const ODESystem &rhs) : m_u_ast_x(rhs.m_u_ast_x), m_u_ast_y(rhs.m_u_ast_y), m_EI(rhs.m_EI), m_GJ(rhs.m_GJ), m_e3(rhs.m_e3) {}

	// move constructor
	ODESystem(ODESystem &&rhs) noexcept;

	// ODESystem destructor
	~ODESystem();

	// Copy assignment operator
	ODESystem &operator=(const ODESystem &rhs);

	// move assignment operator
	ODESystem &operator=(ODESystem &&rhs) noexcept;

	// functor that implements the system of ODEs governing a three-tube CTR
	void operator()(const state_type &y, state_type &dyds, const double s);

	// setter method for updating the parameters for forward kinematics computation
	void setEquationParameters(const vec3d &u_ast_x, const vec3d &u_ast_y, const vec3d &EI, const vec3d &GJ);
};