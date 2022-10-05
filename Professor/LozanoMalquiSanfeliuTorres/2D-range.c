#include "2D-TDSE-CN-solver.h"
#include "2D-configurations.h"

/*INITIAL FUNCTIONS*/

double V_func(double x, double y);
double complex Psi_func(double x, double y);

double V_d1_00(double x, double y);
double V_d1_50(double x, double y);
double V_d2_00(double x, double y);
double V_d2_50(double x, double y);
double V_d3_00(double x, double y);
double V_d3_50(double x, double y);
double V_d4_00(double x, double y);
double V_d4_50(double x, double y);
double V_d5_00(double x, double y);
double V_d5_50(double x, double y);
double V_d6_00(double x, double y);
double V_d6_50(double x, double y);
double V_d7_00(double x, double y);
double V_d7_50(double x, double y);
double V_d8_00(double x, double y);
double V_d8_50(double x, double y);
double V_d9_00(double x, double y);
double V_d9_50(double x, double y);
double V_d10_00(double x, double y);


struct potential{
    double (*func)(double, double);
    char location[64];
};

struct potential sample[19] =  {{V_d1_00, "results-d1-00.txt"},
                                {V_d1_50, "results-d1-50.txt"},
                                {V_d2_00, "results-d2-00.txt"},
                                {V_d2_50, "results-d2-50.txt"},
                                {V_d3_00, "results-d3-00.txt"},
                                {V_d3_50, "results-d3-50.txt"},
                                {V_d4_00, "results-d4-00.txt"},
                                {V_d4_50, "results-d4-50.txt"},
                                {V_d5_00, "results-d5-00.txt"},
                                {V_d5_50, "results-d5-50.txt"},
                                {V_d6_00, "results-d6-00.txt"},
                                {V_d6_50, "results-d6-50.txt"},
                                {V_d7_00, "results-d7-00.txt"},
                                {V_d7_50, "results-d7-50.txt"},
                                {V_d8_00, "results-d8-00.txt"},
                                {V_d8_50, "results-d8-50.txt"},
                                {V_d9_00, "results-d9-00.txt"},
                                {V_d9_50, "results-d9-50.txt"},
                                {V_d10_00, "results-d10-00.txt"}};
int sample_len = (int) sizeof(sample)/sizeof(struct potential);



/*CONSTANTS TO DETERMINE THE RANGE*/
double Z = 10.0;
double pos_x = 13.0;
double momentum = 15.0;



/*MAIN CODE*/

int main(){
    // Initial parameters
    double dx = Lx/(Nx - 1);
    double dy = Ly/(Ny - 1);

    // Configuration of the system
    write_positions(dx, dy, "positions.txt");
    set_psi(Psi_old, Psi_func, dx, dy);
    write_psi(Psi_old, "initial-state.txt");

    // Resolution of the system
    for (int i = 0; i < sample_len; i++){
        // Configuration of the system
        set_V(Vs, sample[i].func, dx, dy);
        set_psi(Psi_old, Psi_func, dx, dy);

        // Resolution of the system
        CN_solver_ABC_range(R, Psi_old, Psi_new, Vs, dx, dy, 1.0, sample[i].location);
        printf("Cas %i fet.\n\n", i + 1);
    }


    return 0;
}

/*FUNCTIONS DECLARATION*/

double V_func(double x, double y){
    double distance = 0.0;
    double result = 0.0;
    result = result + V_coulomb(x, y, pos_x, Ly/2.0 - distance/2.0, 2.0);
    result = result + V_coulomb(x, y, pos_x, Ly/2.0 + distance/2.0, 10.0);
    return result;
}

double complex Psi_func(double x, double y){
    return G_packet_2D(x, y, 2.5, 7.5, 15.0, 0.0, 1);
}

double V_d1_00(double x, double y){
    double distance = 1.00;
    double result = 0.0;
    result = result + V_coulomb(x, y, pos_x, Ly/2.0 - distance/2.0, Z);
    result = result + V_coulomb(x, y, pos_x, Ly/2.0 + distance/2.0, Z);
    return result;
}
double V_d1_50(double x, double y){
    double distance = 1.50;
    double result = 0.0;
    result = result + V_coulomb(x, y, pos_x, Ly/2.0 - distance/2.0, Z);
    result = result + V_coulomb(x, y, pos_x, Ly/2.0 + distance/2.0, Z);
    return result;
}
double V_d2_00(double x, double y){
    double distance = 2.00;
    double result = 0.0;
    result = result + V_coulomb(x, y, pos_x, Ly/2.0 - distance/2.0, Z);
    result = result + V_coulomb(x, y, pos_x, Ly/2.0 + distance/2.0, Z);
    return result;
}
double V_d2_50(double x, double y){
    double distance = 2.50;
    double result = 0.0;
    result = result + V_coulomb(x, y, pos_x, Ly/2.0 - distance/2.0, Z);
    result = result + V_coulomb(x, y, pos_x, Ly/2.0 + distance/2.0, Z);
    return result;
}
double V_d3_00(double x, double y){
    double distance = 3.00;
    double result = 0.0;
    result = result + V_coulomb(x, y, pos_x, Ly/2.0 - distance/2.0, Z);
    result = result + V_coulomb(x, y, pos_x, Ly/2.0 + distance/2.0, Z);
    return result;
}
double V_d3_50(double x, double y){
    double distance = 3.50;
    double result = 0.0;
    result = result + V_coulomb(x, y, pos_x, Ly/2.0 - distance/2.0, Z);
    result = result + V_coulomb(x, y, pos_x, Ly/2.0 + distance/2.0, Z);
    return result;
}
double V_d4_00(double x, double y){
    double distance = 4.00;
    double result = 0.0;
    result = result + V_coulomb(x, y, pos_x, Ly/2.0 - distance/2.0, Z);
    result = result + V_coulomb(x, y, pos_x, Ly/2.0 + distance/2.0, Z);
    return result;
}
double V_d4_50(double x, double y){
    double distance = 4.50;
    double result = 0.0;
    result = result + V_coulomb(x, y, pos_x, Ly/2.0 - distance/2.0, Z);
    result = result + V_coulomb(x, y, pos_x, Ly/2.0 + distance/2.0, Z);
    return result;
}
double V_d5_00(double x, double y){
    double distance = 5.00;
    double result = 0.0;
    result = result + V_coulomb(x, y, pos_x, Ly/2.0 - distance/2.0, Z);
    result = result + V_coulomb(x, y, pos_x, Ly/2.0 + distance/2.0, Z);
    return result;
}
double V_d5_50(double x, double y){
    double distance = 5.50;
    double result = 0.0;
    result = result + V_coulomb(x, y, pos_x, Ly/2.0 - distance/2.0, Z);
    result = result + V_coulomb(x, y, pos_x, Ly/2.0 + distance/2.0, Z);
    return result;
}
double V_d6_00(double x, double y){
    double distance = 6.00;
    double result = 0.0;
    result = result + V_coulomb(x, y, pos_x, Ly/2.0 - distance/2.0, Z);
    result = result + V_coulomb(x, y, pos_x, Ly/2.0 + distance/2.0, Z);
    return result;
}
double V_d6_50(double x, double y){
    double distance = 6.50;
    double result = 0.0;
    result = result + V_coulomb(x, y, pos_x, Ly/2.0 - distance/2.0, Z);
    result = result + V_coulomb(x, y, pos_x, Ly/2.0 + distance/2.0, Z);
    return result;
}
double V_d7_00(double x, double y){
    double distance = 7.00;
    double result = 0.0;
    result = result + V_coulomb(x, y, pos_x, Ly/2.0 - distance/2.0, Z);
    result = result + V_coulomb(x, y, pos_x, Ly/2.0 + distance/2.0, Z);
    return result;
}
double V_d7_50(double x, double y){
    double distance = 7.50;
    double result = 0.0;
    result = result + V_coulomb(x, y, pos_x, Ly/2.0 - distance/2.0, Z);
    result = result + V_coulomb(x, y, pos_x, Ly/2.0 + distance/2.0, Z);
    return result;
}
double V_d8_00(double x, double y){
    double distance = 8.00;
    double result = 0.0;
    result = result + V_coulomb(x, y, pos_x, Ly/2.0 - distance/2.0, Z);
    result = result + V_coulomb(x, y, pos_x, Ly/2.0 + distance/2.0, Z);
    return result;
}
double V_d8_50(double x, double y){
    double distance = 8.50;
    double result = 0.0;
    result = result + V_coulomb(x, y, pos_x, Ly/2.0 - distance/2.0, Z);
    result = result + V_coulomb(x, y, pos_x, Ly/2.0 + distance/2.0, Z);
    return result;
}
double V_d9_00(double x, double y){
    double distance = 9.00;
    double result = 0.0;
    result = result + V_coulomb(x, y, pos_x, Ly/2.0 - distance/2.0, Z);
    result = result + V_coulomb(x, y, pos_x, Ly/2.0 + distance/2.0, Z);
    return result;
}
double V_d9_50(double x, double y){
    double distance = 9.50;
    double result = 0.0;
    result = result + V_coulomb(x, y, pos_x, Ly/2.0 - distance/2.0, Z);
    result = result + V_coulomb(x, y, pos_x, Ly/2.0 + distance/2.0, Z);
    return result;
}
double V_d10_00(double x, double y){
    double distance = 10.00;
    double result = 0.0;
    result = result + V_coulomb(x, y, pos_x, Ly/2.0 - distance/2.0, Z);
    result = result + V_coulomb(x, y, pos_x, Ly/2.0 + distance/2.0, Z);
    return result;
}
