#include <iostream>
#include <string>
#include <vector>
#include <array>
#include <cstdlib>

namespace {
    std::string to_symbol (int i) {
        if (i == -1) {
            return "m";
        } else if (i == 0) {
            return "0";
        } else if (i == 1) {
            return "p";
        } else {
            std::cout << "to_symbol: error\n";
            std::abort();
            return "";
        }
    }

    std::string to_symbol (int i, int j, int k) {
        return to_symbol(i)+to_symbol(j)+to_symbol(k);
    }

    std::string to_symbol (std::array<int,3> const& a) {
        return to_symbol(a[0])+to_symbol(a[1])+to_symbol(a[2]);
    }

    std::string to_indices (int i, int j, int k) {
        return "("+std::to_string(i)+","+std::to_string(j)+","+std::to_string(k)+")";
    }

    std::string to_indices (std::array<int,3> const& a) {
        return to_indices(a[0],a[1],a[2]);
    }

    std::string to_offset (int i) {
        if (i == -1) {
            return "-1";
        } else if (i == 0) {
            return "";
        } else if (i == 1) {
            return "+1";
        } else {
            std::cout << "to_offset: error\n";
            std::abort();
            return "";
        }
    }

    std::string to_A (std::array<int,3> const& base, std::array<int,3> const& off) {
        return "A"+to_symbol(off)+"(ii"+to_offset(base[0])+",jj"+to_offset(base[1])
            +",kk"+to_offset(base[2])+")";
    }

    std::string to_offset_fixed_width (int i) {
        if (i == -1) {
            return "-1";
        } else if (i == 0) {
            return " 0";
        } else if (i == 1) {
            return "+1";
        } else {
            std::cout << "to_offset_fixed_width: error\n";
            std::abort();
            return "";
        }
    }

    std::string to_p (std::array<int,3> const& ijk) {
        return "p("+to_offset_fixed_width(ijk[0])+","+to_offset_fixed_width(ijk[1])+","
            +to_offset_fixed_width(ijk[2])+")";
    }

    std::string to_ap (std::array<int,3> const& ijk) {
        return "ap("+to_offset_fixed_width(ijk[0])+","+to_offset_fixed_width(ijk[1])+","
            +to_offset_fixed_width(ijk[2])+")";
    }
}

int main (int argc, char* argv[])
{
    std::vector<std::array<int,3> > stencil_offset{{1,0,0},
                                                   {0,1,0},
                                                   {0,0,1},
                                                   {1,1,0},
                                                   {1,0,1},
                                                   {0,1,1},
                                                   {1,1,1}};

    std::string space(13,' ');

    for (auto const& sos : stencil_offset) {
        std::string ist = "ist_" + to_symbol(sos[0]) + to_symbol(sos[1]) + to_symbol(sos[2]);
        std::cout << "! csten(i,j,k," + ist << ")\n";

        int i_src = sos[0] * 2;
        int j_src = sos[1] * 2;
        int k_src = sos[2] * 2;

        std::vector<std::array<int,3> > off_ap;

        for         (int k_ap = -1; k_ap <= 1; ++k_ap) {
            for     (int j_ap = -1; j_ap <= 1; ++j_ap) {
                for (int i_ap = -1; i_ap <= 1; ++i_ap) {
                    std::vector<std::array<int,3> > aoff;
                    std::vector<std::array<int,3> > poff;
                    for         (int koff = -1; koff <= 1; ++koff) {
                        for     (int joff = -1; joff <= 1; ++joff) {
                            for (int ioff = -1; ioff <= 1; ++ioff) {
                                int ioff_src = i_ap + ioff - i_src;
                                int joff_src = j_ap + joff - j_src;
                                int koff_src = k_ap + koff - k_src;
                                if (ioff_src >= -1 and ioff_src <= 1 and
                                    joff_src >= -1 and joff_src <= 1 and
                                    koff_src >= -1 and koff_src <= 1)
                                {
                                    poff.push_back({ioff_src, joff_src, koff_src});
                                    aoff.push_back( {ioff, joff, koff});
                                }
                            }
                        }
                    }
                    if (!aoff.empty()) {
                        off_ap.push_back({i_ap,j_ap,k_ap});
                        std::cout << space << "ap"+to_indices(i_ap,j_ap,k_ap) + " = &\n";
                        for (int item = 0; item < aoff.size(); ++item) {
                            if (item == 0) { 
                                std::cout << space << "  & ";
                            } else {
                                std::cout << space << "  + ";
                            }
                            std::cout << space << to_A({i_ap,j_ap,k_ap},aoff[item]) + " * "
                                + to_p(poff[item]);
                            if (item == aoff.size()-1) {
                                std::cout << "\n";
                            } else {
                                std::cout << " &\n";
                            }
                        }
                    }
                }
            }
        }

        std::cout << "\n";

        std::cout << space << "csten(i,j,k,"+ist+") = 0.125d0 * &\n";
        for (int item = 0; item < off_ap.size(); ++item) {
            if (item == 0) {
                std::cout << space << "  ( ";
            } else {
                std::cout << space << "  + ";
            }
            std::cout << "restrict_from_" + to_symbol(off_ap[item]) + "_to(ii,jj,kk) * "
                      << to_ap(off_ap[item]);
            if (item == off_ap.size()-1) {
                std::cout << ")\n";
            } else {
                std::cout << " &\n";
            }
        }

        std::cout << "\n";
    }
}
