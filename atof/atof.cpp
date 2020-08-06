#include <iostream>

constexpr double my_atof (const char* p)
{
    while (*p == ' ' || *p == '\t') ++p;

    double sign = (*p == '-') ? -1.0 : 1.0;
    if (*p == '-' || *p == '+') ++p;

    double r = 0.;
    while (*p >= '0' && *p <= '9') {
        r = r*10. + double(*(p++) - '0');
    }

    if (*p == '.') {
        ++p;
        double r2 = 0.;
        double d = 1.;
        while (*p >= '0' && *p <= '9') {
            r2 = r2*10. + double(*(p++) - '0');
            d *= 10;
        }
        r += r2/d;
    }

    if (*p == 'e') {
        ++p;
        int sexp = (*p == '-') ? -1 : 1;
        if (*p == '-' || *p == '+') ++p;
        int iexp = 0;
        while (*p >= '0' && *p <= '9') {
            iexp = iexp*10 + (*(p++) - '0');
        }
        // need to compute 10**iexp = 10**(\sum 2**n) = \prod 10**(2**n)
        int nmax = 0;
        unsigned short int tmp = iexp;
        while (tmp >>= 1) ++nmax;
        double d = 1.0;
        constexpr double powers[] = {10.,100.,1.e4,1.e8,1.e16,1.e32,1.e64,1.e128,1.e256};
        for (int n = 0; n <= nmax; ++n) {
            if (iexp & 0x1) {
                d *= powers[n];
            }
            iexp >>= 1;
        }
        if (sexp == -1) {
            r /= d;
        } else {
            r *= d;
        }
    }

    return sign*r;
}

constexpr double
operator"" _rt( const char* x )
{
    return double( my_atof(x) );
}

constexpr double
operator"" _rt( unsigned long long int x )
{
    return double( x );
}

int main (int argc, char* argv[])
{
    std::cout << std::scientific
              << 300_rt << std::endl
              << -423.56e3_rt << std::endl
              << -.56e0_rt << std::endl
              << -7e-3_rt << std::endl
              << 7.8e+0_rt << std::endl
              << 0.314e+10_rt << std::endl
              << std::endl;
}
