double addProtect3(double a, double b, double c);
double addProtect4(double a, double b, double c, double d);
double addProtectN(double a[], int len);
double addProtect2(double a, double b);
void ComputeP11_norm(unsigned numWin, int tk_l, double* P1, double* PP1, double** fw, double** bw, double* fw_bw_norm, double** emis);
void ComputeP22_norm(unsigned numWind, int tk_l, double** P, double* PP2, double** fw, double** bw,double* fw_bw_norm,  double** emis);
void ComputeP33_norm(unsigned numWind, int tk_l, double* P3, double* PP3, double** fw, double** bw,double* fw_bw_norm,  double** emis);
void ComputeP44_norm(unsigned numWind, int tk_l, double* P4, double* PP4, double** fw, double** bw,double* fw_bw_norm,  double** emis);
void ComputeP55_norm(unsigned numWind, int tk_l, double** P, double* PP5, double** fw, double** bw,double* fw_bw_norm,  double* stationary, double** emis);
void ComputeP66_norm(unsigned numWind, int tk_l, double** P, double* PP6, double** fw, double** bw,double* fw_bw_norm,  double* stationary, double** emis);
void ComputeP77_norm(unsigned numWind, int tk_l, double** P, double* PP7, double** fw, double** bw,double* fw_bw_norm,  double* stationary, double** emis);
void ComputeP11(unsigned numWin, int tk_l, double* P1, double* PP1, double** fw, double** bw,  double** emis);
void ComputeP22(unsigned numWind, int tk_l, double** P, double* PP2, double** fw, double** bw,  double** emis);
void ComputeP33(unsigned numWind, int tk_l, double* P3, double* PP3, double** fw, double** bw,  double** emis);
void ComputeP44(unsigned numWind, int tk_l, double* P4, double* PP4, double** fw, double** bw,  double** emis);
void ComputeP55(unsigned numWind, int tk_l, double** P, double* PP5, double** fw, double** bw,  double* stationary, double** emis);
void ComputeP66(unsigned numWind, int tk_l, double** P, double* PP6, double** fw, double** bw,  double* stationary, double** emis);
void ComputeP77(unsigned numWind, int tk_l, double** P, double* PP7, double** fw, double** bw,  double* stationary, double** emis);
void ComputeP1(double* tk, int tk_l, double* P, const double* epsize, double rho);
void ComputeP5(double* tk, int tk_l, double* P, const double* epsize);
void ComputeP6(double* tk, int tk_l, double* P, const double* epsize, double rho);
void ComputeP2(int tk_l, double* P2, double* P5);
void ComputeP3(double* tk, int tk_l, double* P3, const double* epsize, double rho);
void ComputeP4(double* tk, int tk_l, double* P4, const double* epsize, double rho);
void ComputeP7(double* tk, int tk_l, double* P7, double* P3, const double* epsize, double rho);
void ComputeP0(int tk_l, double* P0, double* P5); //probability P(T > i)
void ComputeExpectedCoalTime(double* tk, int tk_l, double* expectCoalT, const double* epsize);
double ComputeXXX_norm(int i, int j, double** P);
double calc_trans_norm(int k, int j, double** P);
double ComputeXXX(int i, int j, double** P);
double calc_trans(int k, int j, double** P);
double smartsize1(int tk_l, int numWind, double** fw, double** bw, double* tk, double* P1, double** emis, double pix, double rho, double* newEpSize);
