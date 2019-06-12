
#include "ceres/ceres.h"
#include "ceres/rotation.h"

#include "mex.h"
using ceres::AutoDiffCostFunction;
using ceres::CostFunction;
using ceres::Problem;
using ceres::Solver;
using ceres::Solve;



template <typename T>
        void MultMatrixMy(T * m1, T * m2, T * res){
    res[0] = m1[0] * m2[0] + m1[3] * m2[1] + m1[6] * m2[2];
    res[1] = m1[1] * m2[0] + m1[4] * m2[1] + m1[7] * m2[2];
    res[2] = m1[2] * m2[0] + m1[5] * m2[1] + m1[8] * m2[2];
    res[3] = m1[0] * m2[3] + m1[3] * m2[4] + m1[6] * m2[5];
    res[4] = m1[1] * m2[3] + m1[4] * m2[4] + m1[7] * m2[5];
    res[5] = m1[2] * m2[3] + m1[5] * m2[4] + m1[8] * m2[5];
    res[6] = m1[0] * m2[6] + m1[3] * m2[7] + m1[6] * m2[8];
    res[7] = m1[1] * m2[6] + m1[4] * m2[7] + m1[7] * m2[8];
    res[8] = m1[2] * m2[6] + m1[5] * m2[7] + m1[8] * m2[8];
}
struct myReprojectionError {
    myReprojectionError(double observed_x, double observed_y, double * Porig, double * Xorig,double blockNum)
    : _observed_x(observed_x), _observed_y(observed_y), _Porig(Porig), _Xorig(Xorig), _blockNum(blockNum){}
    
    template <typename T>
            bool operator()(const T* const P, const T* const X,
            T* residuals) const {
        // camera[0,1,2] are the angle-axis rotation.
        
        T Pnow[12];
        T Xnow[4];
        T projection[3];
        for (int i = 0; i < 12; i++){
            Pnow[i] = P[i] + _Porig[i];
        }
        for (int i = 0; i < 3; i++){
            Xnow[i] = X[i] + _Xorig[i];
        }
        Xnow[3] = T(1.0);
        projection[0] = Pnow[0] * Xnow[0] + Pnow[3] * Xnow[1] + Pnow[6] * Xnow[2] + Pnow[9] * Xnow[3];
        projection[1] = Pnow[1] * Xnow[0] + Pnow[4] * Xnow[1] + Pnow[7] * Xnow[2] + Pnow[10] * Xnow[3];
        projection[2] = Pnow[2] * Xnow[0] + Pnow[5] * Xnow[1] + Pnow[8] * Xnow[2] + Pnow[11] * Xnow[3];
        
        projection[0] = projection[0] / projection[2];
        projection[1] = projection[1] / projection[2];
        
        residuals[0] = (projection[0] - _observed_x) / _blockNum;
        residuals[1] = (projection[1] - _observed_y) / _blockNum;
        
        
        return true;
    }
    
    static ceres::CostFunction* CreateMyRep(const double observed_x,
            const double observed_y, double * Porig, double * Xorig, double  blocksNum) {
        return (new ceres::AutoDiffCostFunction<myReprojectionError, 2, 12,3>(
                new myReprojectionError(observed_x, observed_y, Porig, Xorig, blocksNum)));
    }
    
    double _observed_x;
    double _observed_y;
    
    double _blockNum;
    double * _Porig;
    double * _Xorig;
    
    
};






struct ExponentialResidual {
    ExponentialResidual(double x, double y)
    : x_(x), y_(y) {}
    
    template <typename T> bool operator()(const T* const m,
            const T* const c,
            T* residual) const {
        residual[0] = y_ - exp(m[0] * x_ + c[0]);
        return true;
    }
    
private:
    const double x_;
    const double y_;
};


/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    double multiplier;              /* input scalar */
    double *inMatrix;               /* 1xN input matrix */
    size_t ncols;                   /* size of matrix */
    double *outMatrix;              /* output matrix */
    
    
    
    /* get the value of the scalar input  */
    multiplier = mxGetScalar(prhs[0]);
    
    /* create a pointer to the real data in the input matrix  */
    double * Xs; double *xs;double *Ps;double *camPointmap;
    Xs = mxGetPr(prhs[1]);
    camPointmap=mxGetPr(prhs[2]);
    Ps=mxGetPr(prhs[3]);
    xs=mxGetPr(prhs[4]);
    int  xsrows= mxGetM(prhs[4]);
    
    int nrows = mxGetM(prhs[1]);
    ncols = mxGetN(prhs[1]);
    
    
    
    /* create the output matrix */
    plhs[0] = mxCreateDoubleMatrix((mwSize)nrows,(mwSize)ncols,mxREAL);
    double *  Xsu= mxGetPr(plhs[0]);
    for(int ii=0;ii<nrows*ncols;ii++){
        Xsu[ii]=0;
    }
    ncols = mxGetN(prhs[3]);
    nrows = mxGetM(prhs[3]);
    plhs[1] = mxCreateDoubleMatrix((mwSize)nrows,(mwSize)ncols,mxREAL);
    
    
    double *  Psu= mxGetPr(plhs[1]);
    for(int ii=0;ii<nrows*ncols;ii++){
        Psu[ii]=0;
    }
    ceres::Problem problem;
    for (int i = 0; i < xsrows; i++){
        
        
        int camIndex = int(camPointmap[i]);
        int point3DIndex = int(camPointmap[i+xsrows]);
        
        ceres::CostFunction* cost_function =
                myReprojectionError::CreateMyRep(xs[i], xs[i+xsrows], Ps+12*(camIndex - 1),Xs+3*(point3DIndex - 1),1);
        
        ceres::LossFunction * loss_function = new ceres::HuberLoss(0.1);
        problem.AddResidualBlock(cost_function,
                loss_function,
                Psu+12*(camIndex - 1),Xsu+3*(point3DIndex - 1));
    }
    
    
    ceres::Solver::Options options;
    options.function_tolerance = 0.0001;
    options.max_num_iterations = 100;
    options.num_threads = 24;
    
    options.linear_solver_type = ceres::DENSE_SCHUR;
    options.minimizer_progress_to_stdout = true;
    
    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);
    mexPrintf(summary.FullReport().c_str());
    
    
    
    
    
}
