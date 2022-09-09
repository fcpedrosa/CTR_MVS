#pragma once

#include <vector>
#include <iostream>
#include <blaze/Math.h>

typedef blaze::StaticVector<double, 3UL> vec3d;

class Tube
{
private:
	static const double pi_64;
	static const double pi_32;
	double m_OD;
	double m_ID;
	double m_E;
	double m_I;
	double m_G;
	double m_J;
	blaze::DiagonalMatrix<blaze::StaticMatrix<double, 3UL, 3UL, blaze::rowMajor>> m_K;
	double m_ls;
	double m_lc;
	vec3d m_u_ast;

public:
	// default constructor
	Tube();

	// Tube constructor
	Tube(double OD, double ID, double E, double G, double ls, double lc, const vec3d &u_ast);

	// Tube desctructor
	~Tube();

	// copy constructor
	Tube(const Tube &rhs);

	// move constructor
	Tube(Tube &&rhs) noexcept;

	// Copy assignment operator
	Tube &operator=(const Tube &rhs);

	// move assignment operator
	Tube &operator=(Tube &&rhs) noexcept;

	// get method for retrieving the tube parameters
	std::tuple<double, double, double, double, double, vec3d> getTubeParameters();

	// get method for retrieving the tube overall length
	double getTubeLength();

	// get method for retrieving the tube precurvature vector
	vec3d get_u_ast();

	// get method for retrieving the "scalar" curvature along x or y directions
	double get_u_ast(const size_t id);

	// set method for updating the tube precurvature vector
	void set_u_ast(const vec3d &u_ast);

	// set method for updadting the tube's precurvature along x or y directions
	void set_u_ast(const size_t id, const double u);

	// get method for retrieving the length of the straight section
	double getStraightLen();

	// get method for retrieving the length of the tube curved section
	double getCurvLen();

	// set method for updating the length of the straight section
	void setStraightLen(double ls);

	// set method for updating the legnth of the curved section
	void setCurvLen(double lc);

	// get method for retrieving the stiffness matrix
	blaze::DiagonalMatrix<blaze::StaticMatrix<double, 3UL, 3UL, blaze::rowMajor>> getK_Matrix();

	// get method for retrieving the ith entry of the main diagonal
	double getK(int i);
};