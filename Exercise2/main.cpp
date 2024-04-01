#include <iostream>
#include <cmath>

#define dim 2
using namespace std;


struct QR{
    double Qt[dim][dim], R[dim][dim];
};

struct LU{
    double L[dim][dim], U[dim][dim];
};

double norm(double v[dim]){
    return sqrt(v[0]*v[0] + v[1]*v[1]);
}
double scalar_product(double v1[dim], double v2[dim]){
    double scalar = 0;
    for (int j = 0; j<2; j++){
        scalar += v1[j]*v2[j];
    }
    return scalar;
}

double* projection_a_on_u(double a[dim],double u[dim]){


    double scalar = scalar_product(a,u)/norm(u);
    static double projection[dim] = {scalar*u[0], scalar*u[1]};
    return projection;
}

double rel_err_aRb(double a[dim], double b[dim]){
    double diff[dim];
    for (int i = 0; i<2; i++){
        diff[i] = a[i] - b[i];
    }
    return norm(diff)/norm(b);

}


QR QR_alg(double A[dim][dim]){
    QR qr;
    double a1[dim] = {A[0][0], A[1][0]};
    double a2[dim] = {A[0][1], A[1][1]};
    double n1 = norm(a1);
    double e1[dim] = {a1[0]/n1, a1[1]/n1};


    double* proj_a2 = projection_a_on_u(a2, e1);
    double u2[dim];
    for (int j = 0; j<dim; j++){
        u2[j] = a2[j] - proj_a2[j];
    }

    double n2 = norm(u2);

    double e2[dim] = {u2[0]/n2, u2[1]/n2};

    qr.R[0][0] = scalar_product(e1,a1);
    qr.R[1][0] = 0;
    qr.R[0][1] = scalar_product(e1,a2);
    qr.R[1][1] = scalar_product(e2,a2);


    for (int i = 0; i<2; i++){
       qr.Qt[0][i] = e1[i];
    }
    for (int j = 0; j<dim; j++){
        qr.Qt[1][j] = e2[j];
    }

    return qr;
}

double* resolve_QR(QR qr, double b[dim]){
    static double y[dim];
    static double x[dim];
    for (int i = 0; i<2; i++){
        double e[dim];
        for (int j = 0; j<2; j++){
            e[j] = qr.Qt[i][j];
        }
        y[i] = scalar_product(e, b);
    }

    x[1] = y[1]/(qr.R[1][1]);
    x[0] = (y[0] - x[1]*qr.R[0][1])/qr.R[0][0];

    return x;
}

LU LU_alg(double A[dim][dim]){
    LU LU;
    LU.L[0][0] = 1.;
    LU.L[0][1] = 0;
    LU.L[1][1] = 1;
    LU.U[1][0] = 0;

    LU.U[0][0] = A[0][0];
    LU.U[0][1] = A[0][1];
    LU.L[1][0] = A[1][0]/(LU.U[0][0]);
    LU.U[1][1] = A[1][1] - (LU.L[1][0])*(LU.U[0][1]);
    return LU;
}

double* resolve_LU(double L[dim][dim], double U[dim][dim], double b[dim]){
    static double y[dim];
    static double x[dim];

    y[0] = b[0];
    y[1] = b[1] - (L[1][0])*y[0];

    x[1] = y[1]/(U[1][1]);
    x[0] = (b[0] - (U[0][1])*x[1])/(U[0][0]);
    return x;
}


int main()
{
    double x[dim] = {-1,-1};

    double A [dim][dim] = {{5.547001962252291e-01, -3.770900990025203e-02}, {8.320502943378437e-01, -9.992887623566787e-01}};
    double b[dim] ={ -5.169911863249772e-01, 1.672384680188350e-01};

    QR QR1 = QR_alg(A);
    double* solution_QR = resolve_QR(QR1, b);
    double err = rel_err_aRb(solution_QR,x);
    cout <<"Solution of Matrix 1, QR algorithm:"<< endl
         << '[' << solution_QR[0] << ',' << solution_QR[1] << ']'<< endl
         << "relative error: "<< err<< endl;

    LU LU = LU_alg(A);
    double* solution_LU = resolve_LU(LU.L, LU.U, b);
    double err_LU = rel_err_aRb(solution_LU,x);
    cout <<"Solution of Matrix 1, LU algorithm:"<< endl
         << '[' << solution_LU[0] << ',' << solution_LU[1] << ']'<< endl
         << "relative error: "<< err_LU<< endl;


    double A1 [dim][dim] = {{5.547001962252291e-01, -5.540607316466765e-01}, {8.320502943378437e-01, -8.324762492991313e-01}};
    double b1[dim] =  { -6.394645785530173e-04, 4.259549612877223e-04};


    QR QR2 = QR_alg(A1);
    double* solution_QR2 = resolve_QR(QR2, b1);
    double err1 = rel_err_aRb(solution_QR2,x);
    cout <<"Solution of Matrix 1, QR algorithm:"<< endl
         << '[' << solution_QR2[0] << ',' << solution_QR2[1] << ']'<< endl
         << "relative error: "<< err1<< endl;

    LU = LU_alg(A1);
    solution_LU = resolve_LU(LU.L, LU.U, b1);
    err_LU = rel_err_aRb(solution_LU,x);
    cout <<"Solution of Matrix 2, LU algorithm:"<< endl
         << '[' << solution_LU[0] << ',' << solution_LU[1] << ']'<< endl
         << "relative error: "<< err_LU<< endl;

    double A2 [dim][dim] = {{5.547001962252291e-01, -5.547001955851905e-01}, {8.320502943378437e-01, -8.320502947645361e-01}};
    double b2[dim] ={ -6.400391328043042e-10, 4.266924591433963e-10};


    QR1 = QR_alg(A2);
    solution_QR = resolve_QR(QR1, b2);
    err = rel_err_aRb(solution_QR,x);
    cout <<"Solution of Matrix 3, QR algorithm:"<< endl
         << '[' << solution_QR[0] << ',' << solution_QR[1] << ']'<< endl
         << "relative error: "<< err<<endl;

    LU = LU_alg(A2);
    solution_LU = resolve_LU(LU.L, LU.U, b2);
    err_LU = rel_err_aRb(solution_LU,x);
    cout <<"Solution of Matrix 3, LU algorithm:"<< endl
         << '[' << solution_LU[0] << ',' << solution_LU[1] << ']'<< endl
         << "relative error: "<< err_LU<< endl;






  return 0;
}


/*

*/
