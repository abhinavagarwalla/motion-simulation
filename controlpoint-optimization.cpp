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
#include "bayesopt/bayesopt/bayesopt.h"
#include "bayesopt/bayesopt/parameters.h"
#include "bayesopt/utils/displaygp.hpp"

#include <iostream>
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
    if (collides_flag)
        time *= 3;
    return time;
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

    QString filename = "/home/abhinav/Desktop/pathplanner_extras/log/" + QString::fromStdString(fileid) + ".log";
//    QFile file(filename);
//    file.open(QIODevice::WriteOnly| QIODevice::Text);
//    QTextStream stream(&file);

//    double x[4] = {};
//    const double bndl[2*n] = {-1000,-1000,-1000,-1000};
//    const double bndu[2*n] = {1000,1000,1000,1000};
//    const double bndl[2*n] = {-1000,-1000};
//    const double bndu[2*n] = {1000,1000};
    double minf=0;

//    bopt_params bparams = initialize_parameters_to_default();
//    bparams.n_iterations = 150;
//    bparams.verbose_level = 4;
//    set_log_file(&bparams, filename.toStdString().c_str());
//    set_learning(&bparams,"L_MCMC");
//    int k= bayes_optimization(2*n,f_cubicnCP,NULL,bndl, bndu, cps, &minf, bparams);

    bayesopt::Parameters bparams = initialize_parameters_to_default();
//    qDebug() << "Initialised params";
    bparams.n_iterations = 150;
    bparams.verbose_level = 4;
    bparams.l_type = L_MCMC;
    bparams.log_filename = filename.toStdString().c_str();

    f_cubicnCP_eval optimizer(bparams);
    //Define bounds and prepare result.
    vectord result(2);
    vectord lowerBounds(2);
    vectord upperBounds(2);
    lowerBounds(0) = -1000;lowerBounds(1) = -1000;
    upperBounds(0) = 1000;upperBounds(1) = 1000;
    //Set the bounds. This is optional. Default is [0,1]
    //Only required because we are doing continuous optimization
    optimizer.setBoundingBox(lowerBounds,upperBounds);
    //Collect the result in bestPoint
    optimizer.optimize(result);

//    stream << minf;
    // make the trajectory now
    SplineTrajectory *st;
    {
        vector<Pose> cpsf;
        for (int i = 0; i < n; i++) {
            cpsf.push_back(Pose(cps[2*i], cps[2*i+1], 0));
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
