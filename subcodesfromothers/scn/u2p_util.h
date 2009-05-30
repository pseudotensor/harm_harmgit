

extern void primtoU_g( double prim[], double gcov[][4], double gcon[][4],  double U[] );
extern void ucon_calc_g(double prim[],double gcov[][4],double gcon[][4],double ucon[]);
extern void raise_g(double vcov[], double gcon[][4], double vcon[]);
extern void lower_g(double vcon[], double gcov[][4], double vcov[]);
extern void ncov_calc(double gcon[][4],double ncov[]) ;
extern void bcon_calc_g(double prim[],double ucon[],double ucov[],double ncov[],double bcon[]); 
extern double pressure_rho0_u(double rho0, double u);
extern double pressure_rho0_w(double rho0, double w);
extern void dualfprintf(FILE* fileptr, char *format, ...);
