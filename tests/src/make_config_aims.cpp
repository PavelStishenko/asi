#include <fstream>
#include <boost/format.hpp>

using boost::format;

extern const double test_E = -75.668289;
extern const double test_forces[] = {0.004411, 0.000006, 0.021649, -0.004700, -0.014259, -0.009786, -0.004700, 0.014259, -0.009786};
extern const double test_charges[] = {0.243408, -0.121661, -0.121661};
extern const double test_calc_esp[] = {-0.001278, 0.250445};

const std::string control_in = R"--(
KS_method elsi
output   hirshfeld-I
final_forces_cleaned .false.

%1%

xc                                 pw-lda


    sc_accuracy_rho  1E-4
    sc_accuracy_eev  1E-4
    sc_accuracy_etot 1E-6
    sc_iter_limit    100 
    sc_accuracy_forces 5E-2

compensate_multipole_errors .true.

density_update_method density_matrix

spin collinear
charge_mix_param      0.05000000
packed_matrix_format index
multiplicity 3


#===============================================================================

################################################################################
#
#  FHI-aims code project
#  VB, Fritz-Haber Institut, 2009
#
#  Suggested "light" defaults for C atom (to be pasted into control.in file)
#  Be sure to double-check any results obtained with these settings for post-processing,
#  e.g., with the "tight" defaults and larger basis sets.
#
################################################################################
  species        C
#     global species definitions
    nucleus             6
    mass                12.0107
#
    l_hartree           4
#
    cut_pot             3.5  1.5  1.0
    basis_dep_cutoff    1e-4
#
    radial_base         34 5.0
    radial_multiplier   1
    angular_grids specified
      division   0.3326   50
      division   0.5710  110
      division   0.7727  194
      division   0.8772  302
#      division   0.9334  434
#      division   0.9625  590
#      division   0.9924  770
#      division   1.0230  974
#      division   1.4589 1202
#     outer_grid  974
     outer_grid 302
################################################################################
#
#  Definition of "minimal" basis
#
################################################################################
#     valence basis states
    valence      2  s   2.
    valence      2  p   2.
#     ion occupancy
    ion_occ      2  s   1.
    ion_occ      2  p   1.
################################################################################
#
#  Suggested additional basis functions. For production calculations, 
#  uncomment them one after another (the most important basis functions are
#  listed first).
#
#  Constructed for dimers: 1.0 A, 1.25 A, 1.5 A, 2.0 A, 3.0 A
#
################################################################################
#  "First tier" - improvements: -1214.57 meV to -155.61 meV
     hydro 2 p 1.7
     hydro 3 d 6
     hydro 2 s 4.9
#  "Second tier" - improvements: -67.75 meV to -5.23 meV
#     hydro 4 f 9.8
#     hydro 3 p 5.2
#     hydro 3 s 4.3
#     hydro 5 g 14.4
#     hydro 3 d 6.2
#  "Third tier" - improvements: -2.43 meV to -0.60 meV
#     hydro 2 p 5.6
#     hydro 2 s 1.4
#     hydro 3 d 4.9
#     hydro 4 f 11.2
#  "Fourth tier" - improvements: -0.39 meV to -0.18 meV
#     hydro 2 p 2.1
#     hydro 5 g 16.4
#     hydro 4 d 13.2
#     hydro 3 s 13.6
#     hydro 4 f 17.6
#  Further basis functions - improvements: -0.08 meV and below
#     hydro 3 s 2
#     hydro 3 p 6
#     hydro 4 d 20
################################################################################
#
#  FHI-aims code project
#  VB, Fritz-Haber Institut, 2009
#
#  Suggested "light" defaults for N atom (to be pasted into control.in file)
#  Be sure to double-check any results obtained with these settings for post-processing,
#  e.g., with the "tight" defaults and larger basis sets.
#
################################################################################
  species        N
#     global species definitions
    nucleus             7
    mass                14.0067
#
    l_hartree           4
#
    cut_pot             3.5  1.5  1.0
    basis_dep_cutoff    1e-4
#
    radial_base         35 5.0
    radial_multiplier   1
    angular_grids       specified
      division   0.2599   50
      division   0.4601  110
      division   0.5885  194
      division   0.6503  302
#      division   0.6939  434
#      division   0.7396  590
#      division   0.7632  770
#      division   0.8122  974
#      division   1.1604 1202
#      outer_grid  974
      outer_grid  302
################################################################################
#
#  Definition of "minimal" basis
#
################################################################################
#     valence basis states
    valence      2  s   2.
    valence      2  p   3.
#     ion occupancy
    ion_occ      2  s   1.
    ion_occ      2  p   2.
################################################################################
#
#  Suggested additional basis functions. For production calculations, 
#  uncomment them one after another (the most important basis functions are
#  listed first).
#
#  Constructed for dimers: 1.0 A, 1.1 A, 1.5 A, 2.0 A, 3.0 A
#
################################################################################
#  "First tier" - improvements: -1193.42 meV to -220.60 meV
     hydro 2 p 1.8
     hydro 3 d 6.8
     hydro 3 s 5.8
#  "Second tier" - improvements: -80.21 meV to -6.86 meV
#     hydro 4 f 10.8
#     hydro 3 p 5.8
#     hydro 1 s 0.8
#     hydro 5 g 16
#     hydro 3 d 4.9
#  "Third tier" - improvements: -4.29 meV to -0.53 meV
#     hydro 3 s 16
#     ionic 2 p auto
#     hydro 3 d 6.6
#     hydro 4 f 11.6
#  "Fourth tier" - improvements: -0.75 meV to -0.25 meV
#     hydro 2 p 4.5
#     hydro 2 s 2.4
#     hydro 5 g 14.4
#     hydro 4 d 14.4
#     hydro 4 f 16.8
# Further basis functions - -0.21 meV and below
#     hydro 3 p 14.8
#     hydro 3 s 4.4
#     hydro 3 d 19.6
#     hydro 5 g 12.8
################################################################################
#
#  FHI-aims code project
#  VB, Fritz-Haber Institut, 2009
#
#  Suggested "light" defaults for H atom (to be pasted into control.in file)
#  Be sure to double-check any results obtained with these settings for post-processing,
#  e.g., with the "tight" defaults and larger basis sets.
#
################################################################################
  species        H
#     global species definitions
    nucleus             1
    mass                1.00794
#
    l_hartree           4
#
    cut_pot             3.5  1.5  1.0
    basis_dep_cutoff    1e-4
#     
    radial_base         24 5.0
    radial_multiplier   1
    angular_grids       specified
      division   0.2421   50
      division   0.3822  110
      division   0.4799  194
      division   0.5341  302
#      division   0.5626  434
#      division   0.5922  590
#      division   0.6542  770
#      division   0.6868 1202
#      outer_grid  770
      outer_grid  302
################################################################################
#
#  Definition of "minimal" basis
#
################################################################################
#     valence basis states
    valence      1  s   1.
#     ion occupancy
    ion_occ      1  s   0.5
################################################################################
#
#  Suggested additional basis functions. For production calculations, 
#  uncomment them one after another (the most important basis functions are
#  listed first).
#
#  Basis constructed for dimers: 0.5 A, 0.7 A, 1.0 A, 1.5 A, 2.5 A
#
################################################################################
#  "First tier" - improvements: -1014.90 meV to -62.69 meV
     hydro 2 s 2.1
     hydro 2 p 3.5
#  "Second tier" - improvements: -12.89 meV to -1.83 meV
#     hydro 1 s 0.85
#     hydro 2 p 3.7
#     hydro 2 s 1.2
#     hydro 3 d 7
#  "Third tier" - improvements: -0.25 meV to -0.12 meV
#     hydro 4 f 11.2
#     hydro 3 p 4.8
#     hydro 4 d 9
#     hydro 3 s 3.2
################################################################################
#
#  FHI-aims code project
#  VB, Fritz-Haber Institut, 2009
#
#  Suggested "light" defaults for O atom (to be pasted into control.in file)
#  Be sure to double-check any results obtained with these settings for post-processing,
#  e.g., with the "tight" defaults and larger basis sets.
#
################################################################################
  species        O
#     global species definitions
    nucleus             8
    mass                15.9994
#
    l_hartree           4
#
    cut_pot             3.5  1.5  1.0
    basis_dep_cutoff    1e-4
#
    radial_base         36 5.0
    radial_multiplier   1
     angular_grids specified
      division   0.2659   50
      division   0.4451  110
      division   0.6052  194
      division   0.7543  302
#      division   0.8014  434
#      division   0.8507  590
#      division   0.8762  770
#      division   0.9023  974
#      division   1.2339 1202
#      outer_grid 974
      outer_grid 302
################################################################################
#
#  Definition of "minimal" basis
#
################################################################################
#     valence basis states
    valence      2  s   2.
    valence      2  p   4.
#     ion occupancy
    ion_occ      2  s   1.
    ion_occ      2  p   3.
################################################################################
#
#  Suggested additional basis functions. For production calculations, 
#  uncomment them one after another (the most important basis functions are
#  listed first).
#
#  Constructed for dimers: 1.0 A, 1.208 A, 1.5 A, 2.0 A, 3.0 A
#
################################################################################
#  "First tier" - improvements: -699.05 meV to -159.38 meV
     hydro 2 p 1.8
     hydro 3 d 7.6
     hydro 3 s 6.4
#  "Second tier" - improvements: -49.91 meV to -5.39 meV
#     hydro 4 f 11.6
#     hydro 3 p 6.2
#     hydro 3 d 5.6
#     hydro 5 g 17.6
#     hydro 1 s 0.75
#  "Third tier" - improvements: -2.83 meV to -0.50 meV
#     ionic 2 p auto
#     hydro 4 f 10.8
#     hydro 4 d 4.7
#     hydro 2 s 6.8
#  "Fourth tier" - improvements: -0.40 meV to -0.12 meV
#     hydro 3 p 5
#     hydro 3 s 3.3
#     hydro 5 g 15.6
#     hydro 4 f 17.6
#     hydro 4 d 14
# Further basis functions - -0.08 meV and below
#     hydro 3 s 2.1
#     hydro 4 d 11.6
#     hydro 3 p 16
#     hydro 2 s 17.2
)--";

const char *geometry_in = R"--(
atom 0.        0.        0.119262 O
initial_moment 2
atom 0.        0.763239 -0.477047 H
atom 0.       -0.763239 -0.477047 H
)--";

const char *pbc_geometry_in = R"--(
lattice_vector 10.0000000000000000 0.0000000000000000 0.0000000000000000 
lattice_vector 0.0000000000000000 10.0000000000000000 0.0000000000000000 
lattice_vector 0.0000000000000000 0.0000000000000000 10.0000000000000000 
atom_frac 0.5000000000000000 0.5000000000000000 0.5298154500000000 O
initial_moment 2
atom_frac 0.5000000000000000 0.5763239000000001 0.4701845500000000 H
atom_frac 0.5000000000000000 0.4236761000000000 0.4701845500000000 H
)--";

void make_config_files(bool pbc=false)
{
  const char *kgrid = pbc ? "k_grid                 2 1 1" : "#no k_grid";

  std::ofstream fo_control("control.in");
  fo_control << format(control_in) % kgrid;

  std::ofstream fo_geometry("geometry.in");
  fo_geometry << (pbc ? pbc_geometry_in : geometry_in);

}
