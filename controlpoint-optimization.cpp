#include "controlpoint-optimization.hpp"
#include "splines.hpp"
#include <vector>
using namespace std;
#include "drawable.h"
#include "collision-checking.h"
#include <QFile>
#include <QTextStream>
#include <QString>
#include <alglib/optimization.h>
#include "bayesopt/bayesopt.h"
#include "bayesopt/parameters.h"
#include "dataset.hpp"
#include "utils/displaygp.hpp"
#include "utils/testfunctions.hpp"
#include "bopt_state.hpp"
#include "prob_distribution.hpp"

#include <iostream>
#include <fstream>
#include <sstream>
#include <chrono>
#include <sys/time.h>

using namespace alglib;

extern RenderArea *gRenderArea;
//double Optimization::f_cubic2CP(const gsl_vector *x, void *params_)
//{
//    OptParams *params = static_cast<OptParams*>(params_);
//    Pose cp1(gsl_vector_get(x,0)*fieldXConvert, gsl_vector_get(x,1)*fieldXConvert, 0);
//    Pose cp2(gsl_vector_get(x,2)*fieldXConvert, gsl_vector_get(x,3)*fieldXConvert, 0);
//    vector<Pose> midPoints;
//    midPoints.push_back(cp1);
//    midPoints.push_back(cp2);
//    CubicSpline *p = new CubicSpline(params->start, params->end, midPoints);
//    SplineTrajectory *st = new SplineTrajectory(p, params->vls, params->vrs, params->vle, params->vre);
//    double time = st->totalTime();
//    // use maxk also as a cost.
//    // try using ONLY maxk as cost function?
//    //double maxk = p->maxk();
////    delete st;
//    //vector<pair<double,float> > mp = p->lmaxk();
//    //p->lmaxk();
//    return time;
//    //return maxk;
//}


//Trajectory *Optimization::cubicSpline2CPOptimization(Pose start, Pose end, double vls, double vrs, double vle, double vre)
//{
//    OptParams params(start, end, vls, vrs, vle, vre);
//    const gsl_multimin_fminimizer_type *T =
//    gsl_multimin_fminimizer_nmsimplex2;
//    gsl_multimin_fminimizer *s = NULL;
//    gsl_vector *ss, *x;
//    gsl_multimin_function minex_func;

//    size_t iter = 0;
//    int status;
//    double size;

//    int n = 4;  // number of optimization params. here they are 2 CPs (x,y) in cm
//    /* Starting point */
//    // set some good starting values for the 2 CPs
////    Pose cp1((start.x()*2+end.x())*1/3., (start.y()*2+end.y())*1/3., 0);
////    Pose cp2((start.x()+2*end.x())*1/3., (start.y()+2*end.y())*1/3., 0);
//    Pose cp1(1000,1000,0);
//    Pose cp2(-1000,-1000,0);
//    x = gsl_vector_alloc (n);
//    gsl_vector_set(x, 0, cp1.x()/fieldXConvert);
//    gsl_vector_set(x, 1, cp1.y()/fieldXConvert);
//    gsl_vector_set(x, 2, cp2.x()/fieldXConvert);
//    gsl_vector_set(x, 3, cp2.y()/fieldXConvert);

//    /* Set initial step sizes to 100 */
//    ss = gsl_vector_alloc (n);
//    gsl_vector_set_all (ss, 100.0);

//    /* Initialize method and iterate */
//    minex_func.n = n;
//    minex_func.f = f_cubic2CP;
//    minex_func.params = &params;

//    s = gsl_multimin_fminimizer_alloc (T, n);
//    gsl_multimin_fminimizer_set (s, &minex_func, x, ss);

//    do
//    {
//      iter++;
//      status = gsl_multimin_fminimizer_iterate(s);

//      if (status)
//        break;

//      size = gsl_multimin_fminimizer_size (s);
//      status = gsl_multimin_test_size (size, 1e-2);

//      if (status == GSL_SUCCESS)
//        {
//          printf ("converged to minimum at\n");
//        }

//      printf ("%5d %10.3e %10.3e f() = %7.3f size = %.3f\n",
//              iter,
//              gsl_vector_get (s->x, 0),
//              gsl_vector_get (s->x, 1),
//              s->fval, size);
//    }
//    while (status == GSL_CONTINUE && iter < 100);

//    // make the trajectory now
//    SplineTrajectory *st;
//    {
//        Pose cp1(gsl_vector_get(s->x,0)*fieldXConvert, gsl_vector_get(s->x,1)*fieldXConvert, 0);
//        Pose cp2(gsl_vector_get(s->x,2)*fieldXConvert, gsl_vector_get(s->x,3)*fieldXConvert, 0);
//        static PointDrawable *pt1 = NULL, *pt2 = NULL;
//        if (pt1)
//            delete pt1;
//        if (pt2)
//            delete pt2;
//        pt1 = new PointDrawable(QPointF(cp1.x(), cp1.y()), gRenderArea);
//        pt2 = new PointDrawable(QPointF(cp2.x(), cp2.y()), gRenderArea);
//        vector<Pose> midPoints;
//        midPoints.push_back(cp1);
//        midPoints.push_back(cp2);
//        CubicSpline *p = new CubicSpline(start, end, midPoints);
//        st = new SplineTrajectory(p, vls, vrs, vle, vre);
////        double maxk_u, maxk;
////        maxk = p->maxk(&maxk_u);
////        qDebug() << "maxk = " << maxk << ", maxk_u = " << maxk_u;
//    }
//    gsl_vector_free(x);
//    gsl_vector_free(ss);
//    gsl_multimin_fminimizer_free (s);

//    return st;
//}


// make n cp optimizer

double Optimization::f_cubicnCP(const gsl_vector *x, void *params_)
{
    OptParams *params = static_cast<OptParams*>(params_);
    int n = params->n;
    std::vector<Pose> cps;
    for (int i = 0; i < n; i++) {
        cps.push_back(Pose(gsl_vector_get(x, 2*i)*fieldXConvert, gsl_vector_get(x, 2*i+1)*fieldXConvert, 0));
    }
    CubicSpline *p = new CubicSpline(params->start, params->end, cps);

    // check if collides with wall. if so, make the score = 1.5 x time score
    using CollisionChecking::LineSegment;
    vector<LineSegment> ls;
    ls.push_back(LineSegment(-HALF_FIELD_MAXX/fieldXConvert, -HALF_FIELD_MAXY/fieldXConvert, HALF_FIELD_MAXX/fieldXConvert, -HALF_FIELD_MAXY/fieldXConvert));
    ls.push_back(LineSegment(-HALF_FIELD_MAXX/fieldXConvert, HALF_FIELD_MAXY/fieldXConvert, HALF_FIELD_MAXX/fieldXConvert, HALF_FIELD_MAXY/fieldXConvert));
    ls.push_back(LineSegment(-HALF_FIELD_MAXX/fieldXConvert, -HALF_FIELD_MAXY/fieldXConvert, -HALF_FIELD_MAXX/fieldXConvert, +HALF_FIELD_MAXY/fieldXConvert));
    ls.push_back(LineSegment(HALF_FIELD_MAXX/fieldXConvert, -HALF_FIELD_MAXY/fieldXConvert, HALF_FIELD_MAXX/fieldXConvert, +HALF_FIELD_MAXY/fieldXConvert));
    bool collides_flag = false;
    for (int i = 0; i < ls.size(); i++) {
        vector<Pose> collisions = CollisionChecking::cubicSpline_LineSegmentIntersection(*p, ls[i]);
        if (collisions.size()) {
            collides_flag = true;
            break;
        }
    }
    SplineTrajectory *st = new SplineTrajectory(p, params->vls, params->vrs, params->vle, params->vre);
    double time = st->totalTime();
    // use maxk also as a cost.
    // try using ONLY maxk as cost function?
    //double maxk = p->maxk();
//    delete st;
    //vector<pair<double,float> > mp = p->lmaxk();
    //p->lmaxk();
    if (collides_flag)
        time *= 3;
    return time;
    //return maxk;
}

Optimization::OptParams *gparams = NULL;
double Optimization::f_cubicnCP(unsigned int n1, const double *x, double *gradient, void *func_data){
    OptParams *params = static_cast<OptParams*>(gparams);
    int n = params->n;
    std::vector<Pose> cps;
    for (int i = 0; i < n; i++) {
        cps.push_back(Pose(x[2*i]*fieldXConvert, x[2*i+1]*fieldXConvert, 0));
//        qDebug() << x[2*i]*fieldXConvert << " " <<  x[2*i+1]*fieldXConvert << endl;
    }
    CubicSpline *p = new CubicSpline(params->start, params->end, cps);

    // check if collides with wall. if so, make the score = 1.5 x time score
    using CollisionChecking::LineSegment;
    vector<LineSegment> ls;
    ls.push_back(LineSegment(-HALF_FIELD_MAXX/fieldXConvert, -HALF_FIELD_MAXY/fieldXConvert, HALF_FIELD_MAXX/fieldXConvert, -HALF_FIELD_MAXY/fieldXConvert));
    ls.push_back(LineSegment(-HALF_FIELD_MAXX/fieldXConvert, HALF_FIELD_MAXY/fieldXConvert, HALF_FIELD_MAXX/fieldXConvert, HALF_FIELD_MAXY/fieldXConvert));
    ls.push_back(LineSegment(-HALF_FIELD_MAXX/fieldXConvert, -HALF_FIELD_MAXY/fieldXConvert, -HALF_FIELD_MAXX/fieldXConvert, +HALF_FIELD_MAXY/fieldXConvert));
    ls.push_back(LineSegment(HALF_FIELD_MAXX/fieldXConvert, -HALF_FIELD_MAXY/fieldXConvert, HALF_FIELD_MAXX/fieldXConvert, +HALF_FIELD_MAXY/fieldXConvert));
    bool collides_flag = false;
    for (int i = 0; i < ls.size(); i++) {
        vector<Pose> collisions = CollisionChecking::cubicSpline_LineSegmentIntersection(*p, ls[i]);
        if (collisions.size()) {
            collides_flag = true;
            break;
        }
    }
    SplineTrajectory *st = new SplineTrajectory(p, params->vls, params->vrs, params->vle, params->vre);
    double time = st->totalTime();
//    qDebug() << "time:" << time;
//    if (collides_flag)
//        time *= 3;
    return time;
}

Trajectory* Optimization::SingleStepOptimization(Pose start, Pose end, double vls, double vrs, double vle, double vre, int n, std::string fileid) {

    OptParams params(start, end, vls, vrs, vle, vre, n);
    gparams = &params;

    int max_runs = 150;
    QString filename = QString::fromStdString(fileid);

    bayesopt::Parameters bparams = initialize_parameters_to_default();
    bparams.n_iterations = 150;
    bparams.verbose_level = 4;
    bparams.l_type = L_MCMC;
    bparams.log_filename = filename.toStdString().c_str();

    f_cubicnCP_eval *opt(new f_cubicnCP_eval(bparams));
    vectord lowerBounds(2);
    vectord upperBounds(2);
    lowerBounds(0) = -2000/fieldXConvert; lowerBounds(1) = -2000/fieldXConvert;
    upperBounds(0) = 2000/fieldXConvert; upperBounds(1) = 2000/fieldXConvert;
    opt->setBoundingBox(lowerBounds,upperBounds);

    std::string points_logger("/home/kv/Desktop/std_new_run_8.log");
    std::fstream f, f2;
    opt->initializeOptimization();
    for (size_t run = 0; run < max_runs; ++run) {
        opt->stepOptimization();
//        const double res = opt->getData()->getValueAtMinimum();
//        f.open(points_logger, std::ofstream::app);
//        f << res << "\n";
//        f.close();
    }

    qDebug() << "Start: " << start.x() << " " << start.y() << "\n";
    qDebug() << "End: " << end.x() << " " << end.y() << "\n";
//    qDebug() << "Result: " << final_result[0]*fieldXConvert << " " << final_result[1] * fieldXConvert << "\n";

    bayesopt::BOptState state;
    bayesopt::Parameters uparams;
    opt->saveOptimization(state);
    state.saveToFile("/home/kv/Desktop/state3.log");
    qDebug() << "Saved the state !\n";

    qDebug() << "Loading the saved state...\n";
    state.loadFromFile("/home/kv/Desktop/state3.log", uparams);
    qDebug() << "Done!\n";

    // Running from the saved parameters
//    opt = new f_cubicnCP_eval(uparams);
//    opt->setBoundingBox(lowerBounds,upperBounds);

//    qDebug() << "Running from the saved paramters\n";
//    std::chrono::steady_clock::time_point t1, t2;
//    opt->initializeOptimization();
    for (size_t run = 0; run < max_runs; ++run) {
//        qDebug() << "Step Optimization #" << run << "\n";
//        t1 = std::chrono::steady_clock::now();
//        opt->stepOptimization();
//        t2 = std::chrono::steady_clock::now();
//        const double res = opt->getData()->getValueAtMinimum();
    }

//        f.open(points_logger, std::ofstream::app);
//        f << res << " ";
//        f << std::chrono::duration_cast<std::chrono::nanoseconds>(t2-t1).count() << "\n";
//    vectord final_result = opt->getFinalResult();
//    qDebug() << "Result: " << final_result[0]*fieldXConvert << " " << final_result[1] * fieldXConvert << "\n";


    // Scale the actual final result
    vectord final_result = opt->getFinalResult();

    // Generate the 2d grid
    std::vector<double> x, y;
    if (start.x() < end.x()) {
       x = linspace(start.x()/fieldXConvert, end.x()/fieldXConvert, 50);
    } else {
        x = linspace(end.x()/fieldXConvert, start.x()/fieldXConvert, 50);
    }
    if (start.y() < end.y()) {
        y = linspace(start.y()/fieldXConvert, end.y()/fieldXConvert, 50);
    } else {
        y = linspace(end.y()/fieldXConvert, start.y()/fieldXConvert, 50);
    }

    vectord query(2);
    for (size_t i = 0; i < 50; ++i) {
        for (size_t j = 0; j < 50; ++j) {
            query[0] = x[i];
            query[1] = y[j];

            bayesopt::ProbabilityDistribution* pd = opt->getPrediction(query);
            double criterion = -opt->evaluateCriteria(query);
            double sample_evlt = opt->evaluateSample(query);
            double _mean = pd->getMean();
            double stddev = pd->getStd();

            qDebug() << ">>" << query[0] << " "<< query[1] << "\n";
            qDebug() << "mean is " << _mean << " " << "E: " << sample_evlt << endl;
            f2.open(points_logger, std::fstream::app);
            f2 << query[0] << " " << query[1] << " " << _mean << " "<< stddev << " " << criterion << " " << sample_evlt << "\n";
            f2.close();
        }
    }

    query[0] = final_result[0] * fieldXConvert;
    query[1] = final_result[1] * fieldXConvert;
    bayesopt::ProbabilityDistribution* pd = opt->getPrediction(query);
    double criterion = -opt->evaluateCriteria(query);
    double sample_evlt = opt->evaluateSample(query);
    f2.open(points_logger, std::fstream::app);
    double _mean = pd->getMean();
    double stddev = pd->getStd();
    f2 << "-------------------------------------------\n";
    f2 << query[0] << " " << query[1] << " " << _mean << " "<< stddev << " " << criterion << " " << sample_evlt << "\n";
    f2.close();

    // Generate the trajectory
    SplineTrajectory *st;
    {
        vector<Pose> cpsf;
        for (int i = 0; i < n; i++) {
            cpsf.push_back(Pose(final_result[2*i]*fieldXConvert, final_result[2*i+1]*fieldXConvert, 0));
        }
        static vector<PointDrawable*> pts;
        for (int i = 0; i < pts.size(); i++) {
            if (pts[i])
                delete pts[i];
        }
        pts.clear();
        for (int i = 0; i < n; i++) {
            PointDrawable *pt = new PointDrawable(QPointF(cpsf[i].x(), cpsf[i].y()), gRenderArea);
            pts.push_back(pt);
        }
        CubicSpline *p = new CubicSpline(start, end, cpsf);
        st = new SplineTrajectory(p, vls, vrs, vle, vre);
    }

    return st;
}

Trajectory *Optimization::cubicSplinenCPOptimization(Pose start, Pose end, double vls, double vrs, double vle, double vre, int n, std::string fileid)
{
    assert(n >= 0 && n <= 5);
    OptParams params(start, end, vls, vrs, vle, vre, n);
    gparams = &params;
    qDebug() << "after params";
    /* Starting point */
    // set some values for the control points
    double cps[2*n] = {0};
    for (int i = 0; i < n; i++) {
        double x = rand()/(double)RAND_MAX*1000.*((rand()%2)*2-1);
        double y = rand()/(double)RAND_MAX*1000.*((rand()%2)*2-1);
        cps[2*i] = x/fieldXConvert; cps[2*i+1] = y/fieldXConvert;
    }

    QString filename = QString::fromStdString(fileid);
//    QFile file(filename);
//    file.open(QIODevice::WriteOnly| QIODevice::Text);
//    QTextStream stream(&file);

//    double x[4] = {};
//    const double bndl[2*n] = {-1000,-1000,-1000,-1000};
//    const double bndu[2*n] = {1000,1000,1000,1000};
    const double bndl[2*n] = {-100,-100};
    const double bndu[2*n] = {100,100};
    double minf=0;

    bopt_params bparams = initialize_parameters_to_default();
    bparams.n_iterations = 40;
    bparams.verbose_level = 4;
    set_log_file(&bparams, filename.toStdString().c_str());
    set_learning(&bparams,"L_MCMC");
    int k= bayes_optimization(2*n,f_cubicnCP,NULL,bndl, bndu, cps, &minf, bparams);

//    bayesopt::Parameters bparams = initialize_parameters_to_default();
////    qDebug() << "Initialised params";
//    bparams.n_iterations = 150;
//    bparams.verbose_level = 4;
//    bparams.l_type = L_MCMC;
//    bparams.log_filename = filename.toStdString().c_str();

//    f_cubicnCP_eval optimizer(bparams);
//    //Define bounds and prepare result.
//    vectord result(2*n);
//    vectord lowerBounds(2);
//    vectord upperBounds(2);
//    lowerBounds(0) = -1000;lowerBounds(1) = -1000;
//    upperBounds(0) = 1000;upperBounds(1) = 1000;
//    //Set the bounds. This is optional. Default is [0,1]
//    //Only required because we are doing continuous optimization
//    optimizer.setBoundingBox(lowerBounds,upperBounds);
//    //Collect the result in bestPoint
//    optimizer.optimize(result);

//    stream << minf;
    // make the trajectory now
    SplineTrajectory *st;
    {
        vector<Pose> cpsf;
        for (int i = 0; i < n; i++) {
            cpsf.push_back(Pose(cps[2*i]*fieldXConvert, cps[2*i+1]*fieldXConvert, 0));
        }
        static vector<PointDrawable*> pts;
        for (int i = 0; i < pts.size(); i++) {
            if (pts[i])
                delete pts[i];
        }
        pts.clear();
        for (int i = 0; i < n; i++) {
            PointDrawable *pt = new PointDrawable(QPointF(cpsf[i].x(), cpsf[i].y()), gRenderArea);
            pts.push_back(pt);
        }
        CubicSpline *p = new CubicSpline(start, end, cpsf);
        st = new SplineTrajectory(p, vls, vrs, vle, vre);
    }

    return st;
}
