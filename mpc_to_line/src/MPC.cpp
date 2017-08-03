#include "MPC.h"
#include <math.h>
#include <cppad/cppad.hpp>
#include "FG_eval.h"



//
// MPC class definition
//

MPC::MPC() {}
MPC::~MPC() {}

vector<double> MPC::Solve(Eigen::VectorXd x0, Eigen::VectorXd coeffs) 
{
    typedef CPPAD_TESTVECTOR(double) Dvector;

    const double x =      x0[0];
    const double y =      x0[1];
    const double psi =    x0[2];
    const double v =      x0[3];
    const double ctErr =  x0[4];
    const double psiErr = x0[5];

    // number of independent variables
    // numTimeSteps -> numTimeSteps - 1 actuations
    const size_t numIndepVars =     numTimeSteps * dimState             // numTimeSteps * 6 state variables
                                +   (numTimeSteps - 1) * dimActuator;   // (numTimeSteps - 1) * 2 control variables

    // Initial value of the independent variables.
    // Should be 0 except for the initial values.
    Dvector indepVars(numIndepVars);
    for (int i(0); i < numIndepVars; ++i)
    {
        indepVars[i] = 0.0;
    }

    // Set the initial variable values
    indepVars[startIdxX] =      x;
    indepVars[startIdxY] =      y;
    indepVars[startIdxPsi] =    psi;
    indepVars[startIdxV] =      v;
    indepVars[startIdxCTErr] =  ctErr;
    indepVars[startIdxPsiErr] = psiErr;

    cout << "numTimeSteps: " << numTimeSteps << ", numIndepVars: " << numIndepVars 
        << "\nindepVars: " << indepVars << "\n";

    // Lower and upper limits for state vector x0
    Dvector lowerBoundsOfVars(numIndepVars);
    Dvector upperBoundsOfVars(numIndepVars);

    // Set all non-actuators upper and lowerlimits
    // to the max negative and positive values.
    // x, y, psi, v, ctErr, psiErr
    for (size_t i(0); i < startIdxDelta; ++i) 
    {
        lowerBoundsOfVars[i] = -1.0e19;
        upperBoundsOfVars[i] = 1.0e19;
    }

    // The upper and lower limits of delta (steering angle) are set 
    // to -25 and 25 degrees (values in radians).
    // NOTE: Feel free to change this to something else.
    for (size_t i(startIdxDelta); i < startIdxA; ++i)
    {
        lowerBoundsOfVars[i] = -0.436332;
        upperBoundsOfVars[i] = 0.436332;
    }

    // Acceleration/deceleration upper and lower limits.
    // NOTE: Feel free to change this to something else.
    for (size_t i(startIdxA); i < numIndepVars; ++i)
    {
        lowerBoundsOfVars[i] = -1.0;
        upperBoundsOfVars[i] = 1.0;
    }

    // Number of constraints
    const size_t numConstraints(numTimeSteps * dimState);

    // Lower and upper limits for constraints
    // All of these should be 0 except the initial
    // state indices.
    Dvector lowerBoundsOfConstraints(numConstraints);
    Dvector upperBoundsOfConstraints(numConstraints);
    for (size_t i(0); i < numConstraints; ++i)
    {
        lowerBoundsOfConstraints[i] = 0;
        upperBoundsOfConstraints[i] = 0;
    }

    lowerBoundsOfConstraints[startIdxX] =       x;
    lowerBoundsOfConstraints[startIdxY] =       y;
    lowerBoundsOfConstraints[startIdxPsi] =     psi;
    lowerBoundsOfConstraints[startIdxV] =       v;
    lowerBoundsOfConstraints[startIdxCTErr] =   ctErr;
    lowerBoundsOfConstraints[startIdxPsiErr] =  psiErr;

    upperBoundsOfConstraints[startIdxX] =       x;
    upperBoundsOfConstraints[startIdxY] =       y;
    upperBoundsOfConstraints[startIdxPsi] =     psi;
    upperBoundsOfConstraints[startIdxV] =       v;
    upperBoundsOfConstraints[startIdxCTErr] =   ctErr;
    upperBoundsOfConstraints[startIdxPsiErr] =  psiErr;

    cout << "numConstraints: " << numConstraints << "\n";
    cout << "lowerBoundsOfConstraints: " << lowerBoundsOfConstraints << "\n";
    cout << "upperBoundsOfConstraints: " << upperBoundsOfConstraints << "\n";

    // Object that computes objective and constraints
    FG_eval fg_eval(coeffs);

    // options
    std::string options;
    options += "Integer print_level  0\n";
    options += "Sparse  true        forward\n";
    options += "Sparse  true        reverse\n";

    // place to return solution
    CppAD::ipopt::solve_result<Dvector> solution;

    // solve the problem
    CppAD::ipopt::solve<Dvector, FG_eval>(
        options, 
        indepVars, lowerBoundsOfVars, upperBoundsOfVars, 
        lowerBoundsOfConstraints, upperBoundsOfConstraints, 
        fg_eval, solution);

    // Check some of the solution values
    //bool ok = true;
    //ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

    auto cost = solution.obj_value;
    std::cout << "Cost " << cost << std::endl;

    cout << "solution size: " << solution.x.size() << ", solution: " << solution.x << "\n";
    
    return { solution.x[startIdxX + 1], solution.x[startIdxY + 1],
            solution.x[startIdxPsi + 1], solution.x[startIdxV + 1],
            solution.x[startIdxCTErr + 1], solution.x[startIdxPsiErr + 1],
            solution.x[startIdxDelta], solution.x[startIdxA] };
}
