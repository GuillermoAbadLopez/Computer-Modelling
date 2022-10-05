#include "2D-TDSE-CN-solver.h"
#include "2D-configurations.h"

/*INITIAL FUNCTIONS*/

double V_func(double x, double y);
double complex Psi_func(double x, double y);



/*MAIN CODE*/

int main(){
    // Initial parameters
    double dx = Lx/(Nx - 1);
    double dy = Ly/(Ny - 1);

    // Configuration of the system
    write_positions(dx, dy, "positions.txt");
    set_V(Vs, V_func, dx, dy);
    write_V(Vs, "potential.txt");
    set_psi(Psi_old, Psi_func, dx, dy);
    write_psi(Psi_old, "initial-state.txt");

    // Resolution of the system
    CN_solver(R, Psi_old, Psi_new, Vs, dx, dy, 1.0, "results-RBC.txt");
    CN_solver_ABC(R, Psi_old, Psi_new, Vs, dx, dy, 1.0, "results-ABC.txt");

    return 0;
}

/*FUNCTIONS DECLARATION*/

double V_func(double x, double y){
    double V_graphene(double x, double y);
}

double complex Psi_func(double x, double y){
    return G_packet_2D(x, y, 2.5, 7.5, 15.0, 0.0, 1);
}
