// This is a personal academic project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: http://www.viva64.com

#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif

#include <iostream>
#include "include/CTR.hpp"
#include <limits>
#include <chrono>
#include <vector>
#include <array>


int main()
{
	constexpr double inf = std::numeric_limits<double>::infinity();

	//  # # # # # # # # ---- Properties of Nitinol Tubes ---- # # # # # # # #
	double E = 58.0 * 1e9; // Young's modulus GPa
	double G = 21.5 * 1e9; // Shear modulus GPa

	// Precurvature radii for the tubes
	double R1 = 0.04; // (4cm curvature radius)
	double R2 = 0.1;  // (10 cm curvature radius)
	double R3 = inf;  // (infinite curvature radius)

	// -- ** -- Precurvature vectors (for curved portions of the tubes) -- ** -- [u_x* u_y* 0]
	vec3d u1, u2, u3;
	u1 = { 1 / R1, 0.0, 0.0 };
	u2 = { 1 / R2, 0.0, 0.0 };
	u3 = { 1 / R3, 0.0, 0.0 };

	vec3d ls, lc, OD, ID;
	// --** --Lengths of the tubes' straight sections (meters) -- ** --
	ls = { 190e-3, 120e-3, 100e-3 }; // 190, 120, 100

	// --** --Lengths of the tubes' curved sections (meters) -- ** --
	lc = { 60e-3, 80e-3, 0.0 }; // 60, 80, 0;

	// --** --Outer and Inner diameters of the tubes (meters)--** --
	OD = { 0.92e-3, 1.10e-3, 1.40e-3 };
	ID = { 0.80e-3, 0.97e-3, 1.20e-3 };

	// # # # # # ---- Instantiating the three Tube objects ---- # # # # #
	std::shared_ptr<Tube> T1 = std::make_shared<Tube>(OD[0UL], ID[0UL], E, G, ls[0UL], lc[0UL], u1); // innermost tube
	std::shared_ptr<Tube> T2 = std::make_shared<Tube>(OD[1UL], ID[1UL], E, G, ls[1UL], lc[1UL], u2); // intermediate tube
	std::shared_ptr<Tube> T3 = std::make_shared<Tube>(OD[2UL], ID[2UL], E, G, ls[2UL], lc[2UL], u3); // outermost tube

	// instantiating an array of smart pointers to CTR component tubes
	std::array<std::shared_ptr<Tube>, 3UL> Tb = { T1, T2, T3 };

	// initial joint actuation values "home position" - q = [Beta Alpha]
	vec3d Beta_0, Alpha_0;
	Beta_0 = { -50e-3, -20e-3, -10e-3 }; // -115, -100, -80
	Alpha_0 = { mathOp::deg2Rad(0.0), mathOp::deg2Rad(0.0), mathOp::deg2Rad(0.0) };

	blaze::StaticVector<double, 6UL> q_0;
	blaze::subvector<0UL, 3UL>(q_0) = Beta_0;
	blaze::subvector<3UL, 3UL>(q_0) = Alpha_0;

	// Determining the accuracy of BVP solutions
	double Tol = 1e-6;

	// tolerance for position control
	double pos_tol = 5e-4;

	// # # # # # ---- Instantiating the CTR object ---- # # # # #
	CTR robot(Tb, q_0, Tol, mathOp::rootFindingMethod::MODIFIED_NEWTON_RAPHSON);

	// initial guess for the BVP
	blaze::StaticVector<double, 5UL> initGuess;

	/*
	============================================================================================================
							YOU DON'T NEED TO CHANGE ANYTHING UP TO THIS POINT
	============================================================================================================
	*/

	// ************************** Actuating the CTR and solving the corresponding BVP **************************
	// Actuating the CTR and solving the corresponding BVP
	// q_0 = { -0.125, -0.1, -0.08, mathOp::deg2Rad(120), 0.0, 0.0 };	

	// Running a FK example => pass the initial guess vector and the actuation inputs (joint values [q_0])
	robot.actuate_CTR(initGuess, q_0);	
	std::cout << "CTR_robot Tip position is: \n" << robot.getTipPos() << std::endl;

	// Running an IK example: => pass the initial guess vector, a target position [meters], and a tolerance for residual position error [meters]
	blaze::StaticVector<double, 3UL> target = { -0.04, 0.023, 0.18 }; // in meters
	robot.posCTRL(initGuess, target, pos_tol);

	double pos_error = blaze::norm(target - robot.getTipPos());
	std::cout << "\n\nTarget position: " << blaze::trans(target);
	std::cout << "Robot's Tip position after IK: " << blaze::trans(robot.getTipPos());
	std::cout << "Residual position error: " << pos_error << std::endl;


	// passing a new target
	target = { 0.024, -0.056, 0.148 }; // in meters
	robot.posCTRL(initGuess, target, pos_tol);

	pos_error = blaze::norm(target - robot.getTipPos());
	std::cout << "\n\nTarget position: " << blaze::trans(target);
	std::cout << "Robot's Tip position after IK: " << blaze::trans(robot.getTipPos());
	std::cout << "Residual position error: " << pos_error << std::endl;


	// passing a new target
	target = { 0.05, 0.034, 0.195 }; // in meters
	robot.posCTRL(initGuess, target, pos_tol);

	pos_error = blaze::norm(target - robot.getTipPos());
	std::cout << "\n\nTarget position: " << blaze::trans(target);
	std::cout << "Robot's Tip position after IK: " << blaze::trans(robot.getTipPos());
	std::cout << "Residual position error: " << pos_error << std::endl;

	return 0;
}