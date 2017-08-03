#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "matplotlibcpp.h"
#include "MPC.h"


namespace plt = matplotlibcpp;
using namespace std;


//
// Helper functions to fit and evaluate polynomials.
//

// Evaluate a polynomial.
double polyeval(Eigen::VectorXd coeffs, double x)
{
    double result = 0.0;
    for (int i = 0; i < coeffs.size(); i++)
    {
        result += coeffs[i] * pow(x, i);
    }
    return result;
}

// Fit a polynomial, adapted from
// https://github.com/JuliaMath/Polynomials.jl/blob/master/src/Polynomials.jl#L676-L716
Eigen::VectorXd polyfit(
    Eigen::VectorXd xvals,
    Eigen::VectorXd yvals,
    int order)
{
    assert(xvals.size() == yvals.size());
    assert(order >= 1 && order <= xvals.size() - 1);

    Eigen::MatrixXd A(xvals.size(), order + 1);

    for (int i = 0; i < xvals.size(); i++)
    {
        A(i, 0) = 1.0;
    }

    for (int j = 0; j < xvals.size(); j++)
    {
        for (int i = 0; i < order; i++)
        {
            A(j, i + 1) = A(j, i) * xvals(j);
        }
    }

    auto Q = A.householderQr();
    auto result = Q.solve(yvals);

    return result;
}

int main()
{
    MPC mpc;

    Eigen::VectorXd ptsx(2);
    Eigen::VectorXd ptsy(2);
    ptsx << -100, 100;
    ptsy << -1, -1;

    // TODO: fit a polynomial to the above x and y coordinates
    // The polynomial is fitted to a straight line so a polynomial with
    // order 1 is sufficient.
    const int orderOfPolynomial(1);
    const Eigen::VectorXd fittedPolyCoeffs = polyfit(ptsx, ptsy, orderOfPolynomial);
    cout << "fitted polynomial coeffs: [" << fittedPolyCoeffs << "]\n";

    // NOTE: free feel to play around with these
    const double x = -1;
    const double y = 10;
    const double psi = 0;
    const double v = 10;

    // TODO: calculate the cross track error
    // The cross track error is calculated by evaluating at polynomial
    // at x, f(x) and subtracting y.
    const double ctErr = polyeval(fittedPolyCoeffs, x) - y;

    // TODO: calculate the orientation error
    // Due to the sign starting at 0, the orientation error is -f'(x).
    // derivative of coeffs[0] + coeffs[1] * x -> coeffs[1]
    const double psiErr = psi - atan(fittedPolyCoeffs[1]);

    Eigen::VectorXd state(6);
    state << x, y, psi, v, ctErr, psiErr;

    std::vector<double> valuesX = { state[0] };
    std::vector<double> valuesY = { state[1] };
    std::vector<double> valuesPsi = { state[2] };
    std::vector<double> valuesV = { state[3] };
    std::vector<double> valuesCTErr = { state[4] };
    std::vector<double> valuesPsiErr = { state[5] };
    std::vector<double> valuesDelta = {};
    std::vector<double> valuesA = {};

    //const int numIterations = 50;
    const int numIterations(30);
    for (size_t i(0); i < numIterations; ++i)
    {
        std::cout << "Iteration " << i << std::endl;

        const vector<double> vars = mpc.Solve(state, fittedPolyCoeffs);

        // store calculated values
        valuesX.push_back(vars[0]);
        valuesY.push_back(vars[1]);
        valuesPsi.push_back(vars[2]);
        valuesV.push_back(vars[3]);
        valuesCTErr.push_back(vars[4]);
        valuesPsiErr.push_back(vars[5]);

        valuesDelta.push_back(vars[6]);
        valuesA.push_back(vars[7]);

        // update the state
        state << vars[0], vars[1], vars[2], vars[3], vars[4], vars[5];

        std::cout << "x = " << vars[0] << std::endl;
        std::cout << "y = " << vars[1] << std::endl;
        std::cout << "psi = " << vars[2] << std::endl;
        std::cout << "v = " << vars[3] << std::endl;
        std::cout << "cte = " << vars[4] << std::endl;
        std::cout << "epsi = " << vars[5] << std::endl;
        std::cout << "delta = " << vars[6] << std::endl;
        std::cout << "a = " << vars[7] << std::endl;
        std::cout << std::endl;
    }

    // Plot values
    // NOTE: feel free to play around with this.
    // It's useful for debugging!
    plt::subplot(3, 1, 1);
    plt::title("CTE");
    plt::plot(valuesCTErr);
    
    plt::subplot(3, 1, 2);
    plt::title("Delta (Radians)");
    plt::plot(valuesDelta);
    
    plt::subplot(3, 1, 3);
    plt::title("Velocity");
    plt::plot(valuesV);

    plt::show();
}
