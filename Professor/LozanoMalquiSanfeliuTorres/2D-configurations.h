/*POTENTIAL FUNCTIONS*/

double V_doubleslit(double x, double y, double x0, double y1, double y2, double width_x, double width_y, double amplitude){
    if (fabs(x - x0) < width_x/2.0){
        if ((fabs(y - y1) > width_y/2.0) && (fabs(y - y2) > width_y/2.0)){
            return amplitude;
        }
        else {
            return 0;
        }
    }
    return 0;
}

double V_gaussian(double x, double y, double x0, double y0, double sigma, double amplitude){
    double rr = (x - x0)*(x - x0) + (y - y0)*(y - y0);
    return amplitude*exp(-rr/(2.0*sigma*sigma));
}

double V_coulomb(double x, double y, double x0, double y0, double Z){
    double rr = sqrt((x - x0)*(x - x0) + (y - y0)*(y - y0));
    return -Z/rr;
}

double V_zero(double x, double y){
    return 0.0;
}

double V_graphene(double x, double y){
    double distance = 2.70;
    double Z = 6.0;
    double result = 0.0;
    result = result + V_coulomb(x, y, 14.0 - distance*cos(M_PI/6.0), Ly/2.0 - distance/2.0, Z);
    result = result + V_coulomb(x, y, 14.0 - distance*cos(M_PI/6.0), Ly/2.0 + distance/2.0, Z);
    result = result + V_coulomb(x, y, 14.0, Ly/2.0 + distance/2.0 + distance*sin(M_PI/6.0), Z);
    result = result + V_coulomb(x, y, 14.0, Ly/2.0 - distance/2.0 - distance*sin(M_PI/6.0), Z);
    result = result + V_coulomb(x, y, 14.0 + distance*cos(M_PI/6.0), Ly/2.0 - distance/2.0, Z);
    result = result + V_coulomb(x, y, 14.0 + distance*cos(M_PI/6.0), Ly/2.0 + distance/2.0, Z);
    return result;
}

/*INITIAL STATES*/

double complex G_packet_2D(double x, double y, double x0, double y0, double p0x, double p0y, double sigma){
    double rr = (x - x0)*(x - x0) + (y - y0)*(y - y0);
    double pr = p0x*(x - x0) + p0y*(y - y0);
    double A = 1.0/pow(2.0*M_PI*sigma*sigma, 0.5);
    return A*exp(-rr/(4.0*sigma*sigma))*cexp(I*pr);
}

