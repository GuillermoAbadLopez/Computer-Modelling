#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

////////////////////////////////////////////////////////////////////////////////
// INITIAL VARIABLES
////////////////////////////////////////////////////////////////////////////////

#define Nx 400
#define Ny 300
#define M 500
double Lx = 20;
double Ly = 15;
double dt = 0.003;
int iterations = 100;
double a1 = 24;
double a2 = 25;

double xs[Nx][Ny];
double ys[Nx][Ny];
double Vs[Nx][Ny];
double complex bs[Nx][Ny];
double complex Psi_old[Nx][Ny];
double complex Psi_new[Nx][Ny];
double complex R[Nx][Ny][2];

////////////////////////////////////////////////////////////////////////////////
// INITIAL FUNCTIONS
////////////////////////////////////////////////////////////////////////////////

void write_positions(double dx, double dy, char name[64]){
    FILE *arxiu;
    arxiu = fopen(name, "w");
    for (int j = 0; j < Ny; j++){
        for (int i = 0; i < Nx; i++){
            fprintf(arxiu, "%lf\t", i*dx);
        }
    }

    fprintf(arxiu, "\n");

    for (int j = 0; j < Ny; j++){
        for (int i = 0; i < Nx; i++){
            fprintf(arxiu, "%lf\t", j*dy);
        }
    }

    fclose(arxiu);
    return;
}

void set_V(double Vs[Nx][Ny], double (*V)(double, double), double dx, double dy){
    for (int i = 0; i < Nx; i++){
        for (int j = 0; j < Ny; j++){
            Vs[i][j] = V(i*dx, j*dy);
        }
    }
    return;
}

void write_V(double Vs[Nx][Ny], char name[64]){
    FILE *arxiu;
    arxiu = fopen(name, "w");
    for (int j = 0; j < Ny; j++){
        for (int i = 0; i < Nx; i++){
            fprintf(arxiu, "%lf\t", Vs[i][j]);
        }
    }

    fclose(arxiu);
    return;
}

void set_psi(double complex po[Nx][Ny], double complex (*Psi)(double, double), double dx, double dy){
    for (int i = 0; i < Nx; i++){
        for (int j = 0; j < Ny; j++){
            po[i][j] = Psi(i*dx, j*dy);
        }
    }
    return;
}

void write_psi(double complex po[Nx][Ny], char name[64]){
    FILE *arxiu;
    arxiu = fopen(name, "w");
    for (int j = 0; j < Ny; j++){
        for (int i = 0; i < Nx; i++){
            fprintf(arxiu, "%lf\t%lf\t", creal(po[i][j]), cimag(po[i][j]));
        }
    }

    fclose(arxiu);
    return;
}

////////////////////////////////////////////////////////////////////////////////
// CRANK-NICHOLSON METHOD WITH REFLECTIVE BOUNDARY CONDITIONS
////////////////////////////////////////////////////////////////////////////////

void CN_solver(double complex R[Nx][Ny][2], double complex Psi_old[Nx][Ny], double complex Psi_new[Nx][Ny], double Vs[Nx][Ny], double dx, double dy, double m, char name[64]){

    double ax = dt/(4*m*dx*dx);
    double ay = dt/(4*m*dy*dy);
    double complex bs[Nx][Ny];
    for (int i = 0; i < Nx; i++){
        for (int j = 0; j < Ny; j++){
            bs[i][j] = 1 + (2*ax + 2*ay + dt*Vs[i][j]/2.0)*I;
        }
    }


    // Initial conditions
    for (int i = 0; i < Nx; i++){
        for (int j = 0; j < Ny; j++){
            R[i][j][0] = Psi_old[i][j];
        }
    }


    // Computation
    FILE *arxiu;
    arxiu = fopen(name, "w");

    for (int n = 0; n < M; n++){
        for (int k = 0; k < iterations; k++){
            for (int i = 1; i < Nx - 1; i++){
                for (int j = 1; j < Ny - 1; j++){
                    Psi_new[i][j] = conj(bs[i][j])/bs[i][j]*R[i][j][0];
                    Psi_new[i][j] = Psi_new[i][j] + ax*I/(bs[i][j])*(Psi_old[i+1][j] + Psi_old[i-1][j] + R[i+1][j][0] + R[i-1][j][0]);
                    Psi_new[i][j] = Psi_new[i][j] + ay*I/(bs[i][j])*(Psi_old[i][j+1] + Psi_old[i][j-1] + R[i][j+1][0] + R[i][j-1][0]);

                    Psi_old[i][j] = Psi_new[i][j];
                    R[i][j][1] = Psi_new[i][j];
                }
            }
        }

        if (n % 4 == 0){

            for (int j = 0; j < Ny; j++){
                for (int i = 0; i < Nx; i++){
                    fprintf(arxiu, "%lf\t%lf\t", creal(R[i][j][0]), cimag(R[i][j][0]));
                }
            }
            fprintf(arxiu, "\n");
        }

        for (int i = 0; i < Nx; i++){
            for (int j = 0; j < Ny; j++){
                R[i][j][0] = R[i][j][1];
            }
        }
    }

    fclose(arxiu);

    return;
}

////////////////////////////////////////////////////////////////////////////////
// CRANK-NICHOLSON METHOD WITH ABSORBING BOUNDARY CONDITIONS
////////////////////////////////////////////////////////////////////////////////

void CN_solver_ABC(double complex R[Nx][Ny][2], double complex Psi_old[Nx][Ny], double complex Psi_new[Nx][Ny], double Vs[Nx][Ny], double dx, double dy, double m, char name[64]){

    // Parameters for the interior
    double ax = dt/(4*m*dx*dx);
    double ay = dt/(4*m*dy*dy);
    double complex bs[Nx][Ny];
    for (int i = 0; i < Nx; i++){
        for (int j = 0; j < Ny; j++){
            bs[i][j] = 1 + (2*ax + 2*ay + dt*Vs[i][j]/2.0)*I;
        }
    }


    // Parameters for the boundary
    double g1 = (sqrt(2.0*m*a2) - sqrt(2.0*m*a1))/(a2 - a1);
    double g2 = (a2*sqrt(2.0*m*a1) - a1*sqrt(2*m*a2))/(a2 - a1);

    double complex el_1s_x[Nx];
    double complex el_2s_x[Nx];
    double complex er_1s_x[Nx];
    double complex er_2s_x[Nx];

    double complex el_1s_y[Ny];
    double complex el_2s_y[Ny];
    double complex er_1s_y[Ny];
    double complex er_2s_y[Ny];

    for (int i = 0; i < Nx; i++){
        // Conditions at y = 0
        el_1s_x[i] = 1.0 - 2.0*dt/(g1*dx) + I*dt*(g2/g1 - Vs[i][0]);
        el_2s_x[i] = 1.0 + 2.0*dt/(g1/dx) + I*dt*(g2/g1 - Vs[i][0]);
        // Conditions at y = Ly
        er_1s_x[i] = 1.0 + 2.0*dt/(g1*dx) + I*dt*(g2/g1 - Vs[i][Ny-1]);
        er_2s_x[i] = 1.0 - 2.0*dt/(g1*dx) + I*dt*(g2/g1 - Vs[i][Ny-1]);
    }

    for (int j = 0; j < Ny; j++){
        // Conditions at x = 0
        el_1s_y[j] = 1.0 - 2.0*dt/(g1*dx) + I*dt*(g2/g1 - Vs[0][j]);
        el_2s_y[j] = 1.0 + 2.0*dt/(g1/dx) + I*dt*(g2/g1 - Vs[0][j]);
        // Conditions at x = Lx
        er_1s_y[j] = 1.0 + 2.0*dt/(g1*dx) + I*dt*(g2/g1 - Vs[Nx-1][j]);
        er_2s_y[j] = 1.0 - 2.0*dt/(g1*dx) + I*dt*(g2/g1 - Vs[Nx-1][j]);
    }


    // Initial conditions
    for (int i = 0; i < Nx; i++){
        for (int j = 0; j < Ny; j++){
            R[i][j][0] = Psi_old[i][j];
        }
    }


    // Computation
    FILE *arxiu;
    arxiu = fopen(name, "w");

    for (int n = 0; n < M; n++){
        for (int k = 0; k < iterations; k++){


            // Left boundary conditions
            for (int i = 0; i < Nx; i++){
                Psi_new[i][0] = 0.0 - Psi_old[i][1] + el_1s_x[i]*R[i][0][0] + el_2s_x[i]*R[i][1][0];
                Psi_old[i][0] = Psi_new[i][0];
                R[i][0][1] = Psi_new[i][0];
            }
            for (int j = 0; j < Ny; j++){
                Psi_new[0][j] = 0.0 - Psi_old[1][j] + el_1s_y[j]*R[0][j][0] + el_2s_y[j]*R[1][j][0];
                Psi_old[0][j] = Psi_new[0][j];
                R[0][j][1] = Psi_new[0][j];
            }


            // Interior
            for (int i = 1; i < Nx - 1; i++){
                for (int j = 1; j < Ny - 1; j++){
                    Psi_new[i][j] = conj(bs[i][j])/bs[i][j]*R[i][j][0];
                    Psi_new[i][j] = Psi_new[i][j] + ax*I/(bs[i][j])*(Psi_old[i+1][j] + Psi_old[i-1][j] + R[i+1][j][0] + R[i-1][j][0]);
                    Psi_new[i][j] = Psi_new[i][j] + ay*I/(bs[i][j])*(Psi_old[i][j+1] + Psi_old[i][j-1] + R[i][j+1][0] + R[i][j-1][0]);

                    Psi_old[i][j] = Psi_new[i][j];
                    R[i][j][1] = Psi_new[i][j];
                }
            }

            // Right boundary conditions
            for (int i = 0; i < Nx; i++){
                Psi_new[i][Ny-1] = 0.0 - Psi_old[i][Ny-2] + er_1s_x[i]*R[i][Ny-2][0] + er_2s_x[i]*R[i][Ny-1][0];
                Psi_old[i][Ny-1] = Psi_new[i][Ny-1];
                R[i][Ny-1][1] = Psi_new[i][Ny-1];
            }
            for (int j = 0; j < Ny; j++){
                Psi_new[Nx-1][j] = 0.0 - Psi_old[Nx-2][j] + er_1s_y[j]*R[Nx-2][j][0] + er_2s_y[j]*R[Nx-1][j][0];
                Psi_old[Nx-1][j] = Psi_new[Nx-1][j];
                R[Nx-1][j][1] = Psi_new[Nx-1][j];
            }
        }

        // Writing the result
        if (n % 4 == 0){

            for (int j = 0; j < Ny; j++){
                for (int i = 0; i < Nx; i++){
                    fprintf(arxiu, "%lf\t%lf\t", creal(R[i][j][0]), cimag(R[i][j][0]));
                }
            }
            fprintf(arxiu, "\n");
        }

        for (int i = 0; i < Nx; i++){
            for (int j = 0; j < Ny; j++){
                R[i][j][0] = R[i][j][1];
            }
        }
    }

    fclose(arxiu);

    return;
}

////////////////////////////////////////////////////////////////////////////////
// CRANK-NICHOLSON METHOD WITH ABC - RANGE DETERMINATION
////////////////////////////////////////////////////////////////////////////////

void CN_solver_ABC_range(double complex R[Nx][Ny][2], double complex Psi_old[Nx][Ny], double complex Psi_new[Nx][Ny], double Vs[Nx][Ny], double dx, double dy, double m, char name[64]){

    // Parameters for the interior
    double ax = dt/(4*m*dx*dx);
    double ay = dt/(4*m*dy*dy);
    double complex bs[Nx][Ny];
    for (int i = 0; i < Nx; i++){
        for (int j = 0; j < Ny; j++){
            bs[i][j] = 1 + (2*ax + 2*ay + dt*Vs[i][j]/2.0)*I;
        }
    }


    // Parameters for the boundary
    double g1 = (sqrt(2.0*m*a2) - sqrt(2.0*m*a1))/(a2 - a1);
    double g2 = (a2*sqrt(2.0*m*a1) - a1*sqrt(2*m*a2))/(a2 - a1);

    double complex el_1s_x[Nx];
    double complex el_2s_x[Nx];
    double complex er_1s_x[Nx];
    double complex er_2s_x[Nx];

    double complex el_1s_y[Ny];
    double complex el_2s_y[Ny];
    double complex er_1s_y[Ny];
    double complex er_2s_y[Ny];

    for (int i = 0; i < Nx; i++){
        // Conditions at y = 0
        el_1s_x[i] = 1.0 - 2.0*dt/(g1*dx) + I*dt*(g2/g1 - Vs[i][0]);
        el_2s_x[i] = 1.0 + 2.0*dt/(g1/dx) + I*dt*(g2/g1 - Vs[i][0]);
        // Conditions at y = Ly
        er_1s_x[i] = 1.0 + 2.0*dt/(g1*dx) + I*dt*(g2/g1 - Vs[i][Ny-1]);
        er_2s_x[i] = 1.0 - 2.0*dt/(g1*dx) + I*dt*(g2/g1 - Vs[i][Ny-1]);
    }

    for (int j = 0; j < Ny; j++){
        // Conditions at x = 0
        el_1s_y[j] = 1.0 - 2.0*dt/(g1*dx) + I*dt*(g2/g1 - Vs[0][j]);
        el_2s_y[j] = 1.0 + 2.0*dt/(g1/dx) + I*dt*(g2/g1 - Vs[0][j]);
        // Conditions at x = Lx
        er_1s_y[j] = 1.0 + 2.0*dt/(g1*dx) + I*dt*(g2/g1 - Vs[Nx-1][j]);
        er_2s_y[j] = 1.0 - 2.0*dt/(g1*dx) + I*dt*(g2/g1 - Vs[Nx-1][j]);
    }


    // Initial conditions
    for (int i = 0; i < Nx; i++){
        for (int j = 0; j < Ny; j++){
            R[i][j][0] = Psi_old[i][j];
        }
    }


    // Computation
    FILE *arxiu;
    arxiu = fopen(name, "w");

    for (int n = 0; n < M; n++){
        for (int k = 0; k < iterations; k++){


            // Left boundary conditions
            for (int i = 0; i < Nx; i++){
                Psi_new[i][0] = 0.0 - Psi_old[i][1] + el_1s_x[i]*R[i][0][0] + el_2s_x[i]*R[i][1][0];
                Psi_old[i][0] = Psi_new[i][0];
                R[i][0][1] = Psi_new[i][0];
            }
            for (int j = 0; j < Ny; j++){
                Psi_new[0][j] = 0.0 - Psi_old[1][j] + el_1s_y[j]*R[0][j][0] + el_2s_y[j]*R[1][j][0];
                Psi_old[0][j] = Psi_new[0][j];
                R[0][j][1] = Psi_new[0][j];
            }


            // Interior
            for (int i = 1; i < Nx - 1; i++){
                for (int j = 1; j < Ny - 1; j++){
                    Psi_new[i][j] = conj(bs[i][j])/bs[i][j]*R[i][j][0];
                    Psi_new[i][j] = Psi_new[i][j] + ax*I/(bs[i][j])*(Psi_old[i+1][j] + Psi_old[i-1][j] + R[i+1][j][0] + R[i-1][j][0]);
                    Psi_new[i][j] = Psi_new[i][j] + ay*I/(bs[i][j])*(Psi_old[i][j+1] + Psi_old[i][j-1] + R[i][j+1][0] + R[i][j-1][0]);

                    Psi_old[i][j] = Psi_new[i][j];
                    R[i][j][1] = Psi_new[i][j];
                }
            }

            // Right boundary conditions
            for (int i = 0; i < Nx; i++){
                Psi_new[i][Ny-1] = 0.0 - Psi_old[i][Ny-2] + er_1s_x[i]*R[i][Ny-2][0] + er_2s_x[i]*R[i][Ny-1][0];
                Psi_old[i][Ny-1] = Psi_new[i][Ny-1];
                R[i][Ny-1][1] = Psi_new[i][Ny-1];
            }
            for (int j = 0; j < Ny; j++){
                Psi_new[Nx-1][j] = 0.0 - Psi_old[Nx-2][j] + er_1s_y[j]*R[Nx-2][j][0] + er_2s_y[j]*R[Nx-1][j][0];
                Psi_old[Nx-1][j] = Psi_new[Nx-1][j];
                R[Nx-1][j][1] = Psi_new[Nx-1][j];
            }
        }

        // Writing the result
        if (n == 440){

            for (int j = 0; j < Ny; j++){
                for (int i = 0; i < Nx; i++){
                    fprintf(arxiu, "%lf\t%lf\t", creal(R[i][j][0]), cimag(R[i][j][0]));
                }
            }
            fprintf(arxiu, "\n");
        }

        for (int i = 0; i < Nx; i++){
            for (int j = 0; j < Ny; j++){
                R[i][j][0] = R[i][j][1];
            }
        }
    }

    fclose(arxiu);

    return;
}
