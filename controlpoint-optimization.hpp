#ifndef CONTROLPOINTOPTIMIZATION_HPP
#define CONTROLPOINTOPTIMIZATION_HPP
#include <gsl/gsl_multimin.h>
#include "splines.hpp"
#include "pose.h"
#include <alglib/optimization.h>
#include "bayesopt/bayesopt.h"
#include "bayesopt/parameters.h"
#include "utils/displaygp.hpp"

using namespace alglib;
namespace Optimization {
struct OptParams {
    Pose start, end;
    double vls, vrs, vle, vre;
    // n = number of control points
    int n;
    OptParams(): start(), end(), vls(0), vrs(0), vle(0), vre(0), n(0){}
    OptParams(Pose start, Pose end, double vls, double vrs, double vle, double vre, int n): start(start),
        end(end), vls(vls), vrs(vrs), vle(vle), vre(vre), n(n) {}
};

// optimization function
double f_cubicnCP(const gsl_vector* x, void * params);
double f_cubicnCP(unsigned int n, const double *x, double *gradient, void *func_data);
Trajectory *cubicSplinenCPOptimization(Pose start, Pose end, double vls, double vrs, double vle, double vre, int n, std::string fileid);

// Single step optimization
Trajectory* SingleStepOptimization(Pose start, Pose end, double vls, double vrs, double vle, double vre, int n, std::string fileid);

//OptParams *gparams;
class f_cubicnCP_eval: public bayesopt::ContinuousModel{
public:
    f_cubicnCP_eval(bayesopt::Parameters param):
       ContinuousModel(2, param)
     {qDebug() << "Class is constructed ";}
     double evaluateSample( const vectord& query )
     {
         if (query.size() != 2)
          {
            std::cout << "WARNING: This only works for 2D inputs." << std::endl
                      << "WARNING: Using only first two components." << std::endl;
          }
         double x[100];
             for (size_t i = 0; i < query.size(); ++i)
               {
                 x[i] = query(i);
               }
         return f_cubicnCP(query.size(), x, NULL, NULL);
     }
     bool checkReachability( const vectord& query )
     {return true;}
}; //end of class f_cubicnCP

}
#endif // CONTROLPOINTOPTIMIZATION_HPP

