#include <complex>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>
#include <memory>

#include <mpi.h>
#include <memory>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <ctime>

#include "asi.h"
#include "utils.hpp"
#include "blacsutils.h"

#include <boost/multi_array.hpp>


int mpi_provided_threading, world_size, world_rank;


template <typename T>
void gather(int *blacs_descr, T *blacs_data, int *dest_desc, T *dest)
{
  int blacs_ctx = blacs_descr[CTXT_];
  int sys_ctx = get_system_context(blacs_ctx);

  int gatherer_blacs_ctx = make_blacs_context(sys_ctx, 1, 1);
  blacs_desc_init(blacs_descr[M_], blacs_descr[N_], gatherer_blacs_ctx, dest_desc);
  
  int nprow, npcol, myrow, mycol;
  blacs_gridinfo_(&blacs_ctx, &nprow, &npcol, &myrow, &mycol);
  assert((dest_desc[CTXT_] != -1) == ((myrow == 0) && (mycol == 0)));
  
  int ONE = 1;
  pXgemr2d(&blacs_descr[M_], &blacs_descr[N_], 
            blacs_data, &ONE, &ONE, blacs_descr, 
            dest,  &ONE, &ONE, dest_desc, &blacs_ctx);
}


struct expdata_base_t
{
  virtual ~expdata_base_t() {}
  virtual void print_all() = 0;
};

template<typename T>
T univ_conj(T x);

template<>
double univ_conj(double x)
{
  return x;
}

template<>
std::complex<double> univ_conj(std::complex<double> x)
{
  return std::conj(x);
}

double kill_negative_zero(double x, double eps)
{
    return fabs(x) < eps ? 0.0 : x;
}

std::complex<double> kill_negative_zero(std::complex<double> x, double eps)
{
  double xr = kill_negative_zero(x.real(), eps);
  double xi = kill_negative_zero(x.imag(), eps);
  return std::complex<double>(xr, xi);
}

template <typename T>
struct expdata_t: public expdata_base_t
{
  expdata_t(size_t n_kpts, size_t n_spin, size_t n_basis): n_kpts(n_kpts), n_spin(n_spin), n_basis(n_basis), 
      dm(boost::extents[n_kpts][n_spin][n_basis][n_basis]), 
      overlap(boost::extents[n_kpts][n_spin][n_basis][n_basis]), 
      hamiltonian(boost::extents[n_kpts][n_spin][n_basis][n_basis])
  {
    ASI_register_dm_callback(&expdata_t<T>::static_dm_callback, this);
    ASI_register_overlap_callback(&expdata_t<T>::static_overlap_callback, this);
    if (ASI_flavour() != ASI_FHI_AIMS) // not implemented for AIMS
    {
      ASI_register_hamiltonian_callback(&expdata_t<T>::static_hamiltonian_callback, this);
    }
  }
  size_t n_kpts, n_spin, n_basis;
  boost::multi_array<T, 4> dm, overlap, hamiltonian;
  int gathered_desc_dm[DLEN_], gathered_desc_s[DLEN_], gathered_desc_h[DLEN_];
  
  template<typename... Args>
  static void static_dm_callback(void *this_ptr, Args... args)
  {
    reinterpret_cast<expdata_t<T>*>(this_ptr)->dm_callback(args...);
  }
  void dm_callback(int iK, int iS, int *blacs_descr, void *blacs_data)
  {
    if (blacs_descr==0)
    {
      memcpy(&(dm[iK-1][iS-1][0][0]), blacs_data, sizeof(T)*n_basis*n_basis);
    }
    else
    {
      memset(&(dm[iK-1][iS-1][0][0]), 0, sizeof(T)*n_basis*n_basis);
      gather(blacs_descr, reinterpret_cast<T*>(blacs_data), gathered_desc_dm, &(dm[iK-1][iS-1][0][0]));
    }
  }

  template<typename... Args>
  static void static_overlap_callback(void *this_ptr, Args... args)
  {
    reinterpret_cast<expdata_t<T>*>(this_ptr)->overlap_callback(args...);
  }
  void overlap_callback(int iK, int iS, int *blacs_descr, void *blacs_data)
  {
    if (blacs_descr==0) // allways TRUE form AIMS
    {
      memcpy(&(overlap[iK-1][iS-1][0][0]), blacs_data, sizeof(T)*n_basis*n_basis);
    }
    else
    {
      memset(&(overlap[iK-1][iS-1][0][0]), 0, sizeof(T)*n_basis*n_basis);
      gather(blacs_descr, reinterpret_cast<T*>(blacs_data), gathered_desc_s, &(overlap[iK-1][iS-1][0][0]));
    }
  }

  template<typename... Args>
  static void static_hamiltonian_callback(void *this_ptr, Args... args)
  {
    reinterpret_cast<expdata_t<T>*>(this_ptr)->hamiltonian_callback(args...);
  }
  void hamiltonian_callback(int iK, int iS, int *blacs_descr, void *blacs_data)
  {
    if (blacs_descr==0)
    {
      memcpy(&(hamiltonian[iK-1][iS-1][0][0]), blacs_data, sizeof(T)*n_basis*n_basis);
    }
    else
    {
      memset(&(hamiltonian[iK-1][iS-1][0][0]), 0, sizeof(T)*n_basis*n_basis);
      gather(blacs_descr, reinterpret_cast<T*>(blacs_data), gathered_desc_h, &(hamiltonian[iK-1][iS-1][0][0]));
    }
  }
  
  template<typename F, typename M>
  void print_matrix(F &f, const M &m)
  {
    for (int i = 0; i < n_basis; ++i)
    {
      for (int j = 0; j < n_basis; ++j)
      {
        std::complex<double> cc = m[i][j];
        f << m[i][j] << " ";
      }
      f << std::endl;
    }
  }
  
  void print_all()
  {
    for (int iK = 0; iK < n_kpts; ++iK)
    {
      for (int iS = 0; iS < n_spin; ++iS)
      {
        std::ofstream fo(STR("dmhs_"<<iK<<"_"<<iS << "_" << world_rank <<".out").c_str());
        fo << "DM" << std::endl;
        print_matrix(fo, dm[iK][iS]);
        fo << "overlap" << std::endl;
        print_matrix(fo, overlap[iK][iS]);
        fo << "hamiltonian" << std::endl;
        print_matrix(fo, hamiltonian[iK][iS]);
      }
    }
    
    int n_local_ks = ASI_get_n_local_ks();
    int local_ks[2 * n_local_ks];
    ASI_get_local_ks(local_ks);

    std::vector<T> SDM(n_kpts * n_spin), HDM(n_kpts * n_spin);
    std::vector<T> SDM_reduced(n_kpts * n_spin), HDM_reduced(n_kpts * n_spin);
    std::fill(SDM.begin(), SDM.end(), 0);
    std::fill(HDM.begin(), HDM.end(), 0);

    for (int i_local_ks = 0; i_local_ks < n_local_ks; ++i_local_ks)
    {
      int iK = local_ks[i_local_ks*2] - 1;
      int iS = local_ks[i_local_ks*2 + 1] - 1;
      assert((iK >= 0) and (iK < n_kpts));
      assert((iS >= 0) and (iS < n_spin));
      T S = 0, H = 0;
      for (int i = 0; i < n_basis; ++i)
      {
        for (int j = 0; j < n_basis; ++j)
        {
          if (ASI_flavour() != ASI_DFTBP)
          {
            assert(fabs(abs(dm[iK][iS][i][j]) - abs(dm[iK][iS][j][i])) < 1e-9); // Hermitian, excluding DFTB+
          }
          T dm_ij = (i < j) ? dm[iK][iS][i][j] : univ_conj(dm[iK][iS][j][i]); // necessary for DFTB+ only
          T  s_ij = overlap[iK][iS][i][j];
          S += s_ij * dm_ij;
          H += hamiltonian[iK][iS][i][j] * dm_ij;
        }
      }
      SDM.at(iK*n_spin + iS) = S;
      HDM.at(iK*n_spin + iS) = H;
    }// i_local_ks
    
    MPI_Reduce(&(SDM[0]), &(SDM_reduced[0]), SDM.size(), MPI_DOUBLE_COMPLEX,  MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&(HDM[0]), &(HDM_reduced[0]), HDM.size(), MPI_DOUBLE_COMPLEX,  MPI_SUM, 0, MPI_COMM_WORLD);
    if (world_rank == 0)
    {
      for (int iK = 0; iK < n_kpts; ++iK)
      {
        for (int iS = 0; iS < n_spin; ++iS)
        {
          std::cout << "S*DM " << iK << " " << iS << " "  << kill_negative_zero(SDM_reduced[iK*n_spin + iS], 1e-6) << "  ";
          std::cout << "H*DM " << iK << " " << iS << " "  << kill_negative_zero(HDM_reduced[iK*n_spin + iS], 1e-6) << std::endl;
        }
      }
      std::cout << "N_el = "  << kill_negative_zero(std::accumulate(SDM_reduced.begin(), SDM_reduced.end(), T(0.0)), 1e-6) << std::endl;
    }
  }// print_all
};



int main(int argc, char *argv[])
{
  MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &mpi_provided_threading); // instead of MPI_Init
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
 
  assert(argc == 2);

  const MPI_Fint f_mpi_comm = MPI_Comm_c2f(MPI_COMM_WORLD);
  if (world_rank == 0) std::cout << "world_size == " << world_size << std::endl;

  ASI_init(argv[1], "asi.log", f_mpi_comm);
  if ((ASI_flavour() == ASI_DFTBP) &&  // not implemented properly for AIMS
    (ASI_n_atoms() == 2))              //  only for Si2
  {
  
    /* Coordinates in row major format, atomic units */
    double coords_si2[] = {
      0.0000000000000000, 0.0000000000000000, 0.0000000000000000,
      2.2639291987021915, 2.4639291987021915, 2.5639291987021915
    };

    /* Lattice vectors in row major format, atomic units */
    double latvecs_si2[] = {
      5.2278583974043830, 5.1278583974043830, 0.0000000000000000,
      0.0000000000000000, 5.3278583974043830, 5.1278583974043830,
      5.1278583974043830, 0.0000000000000000, 5.4278583974043830
    };
    ASI_set_geometry(coords_si2, 2, latvecs_si2);
  }

  if (world_rank == 0) std::cout << "ASI_init complete" << std::endl;
  bool is_real = ASI_is_hamiltonian_real();
  int   n_kpts = ASI_get_nkpts();
  int   n_spin = ASI_get_nspin();
  int  n_basis = ASI_get_basis_size();
  int n_local_ks = ASI_get_n_local_ks();
  
  if (world_rank == 0)
  {
    std::cout << "n_kpts == " << n_kpts << std::endl;
    std::cout << "n_spin == " << n_spin << std::endl;
    std::cout << "n_basis == " << n_basis << std::endl;
    std::cout << "is_real == " << is_real << std::endl;
    std::cout << "n_local_ks == " << n_local_ks << std::endl;
  }
  
  
  std::unique_ptr<expdata_base_t> export_data;
  if (is_real)
  {
    export_data = std::make_unique<expdata_t<double> >(n_kpts, n_spin, n_basis);
  }
  else
  {
    export_data = std::make_unique<expdata_t<std::complex<double> > >(n_kpts, n_spin, n_basis);
  }

  ASI_run();
  if (world_rank == 0) std::cout << "ASI_run complete" << std::endl;
  assert(ASI_is_hamiltonian_real() == is_real);
  assert(ASI_get_nkpts() == n_kpts);
  assert(ASI_get_nspin() == n_spin);
  assert(ASI_get_basis_size() == n_basis);

  auto E = ASI_energy();
  std::cout << std::setprecision(6) << std::fixed;
  if (world_rank == 0) std::cout << "Energy == " << E  << " Ha = " << E * 27.2113845 << " eV" << std::endl;
  const double *forces = ASI_forces();
  if ((world_rank == 0) && (forces))
  {
    std::cout << "Forces == ";
    for (size_t i = 0; i < 3 * ASI_n_atoms(); ++i )
      std::cout << kill_negative_zero(forces[i], 1e-6) << " ";
    std::cout << std::endl;
  }
  
  export_data->print_all();

  ASI_finalize(); 
  MPI_Finalize();
  return 0;
}
