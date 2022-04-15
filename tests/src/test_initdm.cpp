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
void distrib(int *src_descr, T *src_data, int *dest_desc, T *dest_data)
{
  int ONE = 1;
  int M = dest_desc[M_];
  int N = dest_desc[N_];
  int blacs_ctx = dest_desc[CTXT_];
  int sys_ctx = get_system_context(blacs_ctx);
  int src_blacs_ctx = make_blacs_context(sys_ctx, 1, 1);
  blacs_desc_init(M, N, src_blacs_ctx, src_descr);
  
  pXgemr2d(&M, &N, 
            src_data, &ONE, &ONE, src_descr, 
            dest_data,  &ONE, &ONE, dest_desc, &blacs_ctx);
}


struct expdata_base_t
{
  virtual ~expdata_base_t() {}
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


std::vector<double> load_txt(const char *fn)
{
  std::ifstream fi(fn);
  std::vector<double> v;
  double tmp;
  while(true)
  {
    fi >>  tmp;
    if (fi.eof()) break;
    v.push_back(tmp);
  }
  fi.close();
  return v;
}

template<typename V, typename S>
void copy_vector2square(const V &src, S dst)
{
  size_t n = sqrt(src.size());
  assert(n*n == src.size());
  for (size_t i = 0; i < n; ++i)
  for (size_t j = 0; j < n; ++j)
    dst[i][j] = src[i*n + j];
}



template <typename T>
struct expdata_t: public expdata_base_t
{
  expdata_t(size_t n_kpts, size_t n_spin, size_t n_basis): n_kpts(n_kpts), n_spin(n_spin), n_basis(n_basis), 
      dm(boost::extents[n_kpts][n_spin][n_basis][n_basis])
  {
    for (int iK = 1; iK <= n_kpts; ++iK)
    for (int iS = 1; iS <= n_spin; ++iS)
    {
      auto dm_loaded = load_txt(STR("dm_"<<iK<<"_"<<iS<<".init").c_str());
      assert(dm_loaded.size() == n_basis * n_basis);
      copy_vector2square(dm_loaded, dm[iK-1][iS-1]);
    }
    ASI_register_dm_init_callback(&expdata_t<T>::static_dm_init_callback, this);
  }
  size_t n_kpts, n_spin, n_basis;
  boost::multi_array<T, 4> dm;
  
  template<typename... Args>
  static void static_dm_init_callback(void *this_ptr, Args... args)
  {
    reinterpret_cast<expdata_t<T>*>(this_ptr)->dm_init_callback(args...);
  }
  void dm_init_callback(int iK, int iS, int *blacs_descr, void *blacs_data)
  {
    if (world_rank == 0) std::cout << "dm_init_callback " << iK <<  " " << iS << " blacs_descr=" << (blacs_descr==0 ? "NULL":"True") << std::endl;
    if ((blacs_descr==0) or (world_size == 1))
    {
      memcpy(blacs_data, &(dm[iK-1][iS-1][0][0]), sizeof(T)*n_basis*n_basis);
    }
    else
    {
      //assert(false);
      //memset(&(dm[iK-1][iS-1][0][0]), 0, sizeof(T)*n_basis*n_basis);
      
      int src_desc[DLEN_];
      distrib(src_desc, &(dm[iK-1][iS-1][0][0]), blacs_descr, reinterpret_cast<T*>(blacs_data));

    }
  }
};

int main(int argc, char *argv[])
{
  MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &mpi_provided_threading); // instead of MPI_Init
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
 
  assert(argc == 1);

  const MPI_Fint f_mpi_comm = MPI_Comm_c2f(MPI_COMM_WORLD);

  ASI_init("dummy", "asi.log", f_mpi_comm);
  if (world_rank == 0) std::cout << "ASI_init complete" << std::endl;
  bool is_real = ASI_is_hamiltonian_real();
  int   n_kpts = ASI_get_nkpts();
  int   n_spin = ASI_get_nspin();
  int  n_basis = ASI_get_basis_size();
  
  if (world_rank == 0)
  {
    std::cout << "n_kpts == " << n_kpts << std::endl;
    std::cout << "n_spin == " << n_spin << std::endl;
    std::cout << "n_basis == " << n_basis << std::endl;
    std::cout << "is_real == " << is_real << std::endl;
  }
  
  
  std::unique_ptr<expdata_base_t> export_data;
  if (is_real)
  {
    export_data = std::make_unique<expdata_t<double> >(n_kpts, n_spin, n_basis);
  }
  else
  {
    //export_data = std::make_unique<expdata_t<std::complex<double> > >(n_kpts, n_spin, n_basis);
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
      std::cout << forces[i] << " ";
    std::cout << std::endl;
  }

  ASI_finalize(); 
  MPI_Finalize();
  return 0;
}
