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


int mpi_provided_threading, world_size, world_rank;


template<typename A>
void loadtxt(const char *fn, int n, A &a)
{
  std::ifstream f(fn);
  for (size_t i = 0; i < n; ++i)
    f >> a[i];
}

template<typename A>
void savetxt(const char *fn, int n, const A &a)
{
  std::ofstream f(fn);
  for (size_t i = 0; i < n; ++i)
    f << a[i] << std::endl;
}


struct DM_callbacks_t
{
  DM_callbacks_t(size_t n_basis):  n_basis(n_basis), dm(n_basis*n_basis)
  {
    ASI_register_dm_init_callback(&DM_callbacks_t::static_dm_init_callback, this);
    ASI_register_dm_callback(&DM_callbacks_t::static_dm_callback, this);
  }
  size_t n_basis; 
  std::vector<double> dm; // density matrix n_basis*n_basis
  
  
  template<typename... Args>
  static void static_dm_init_callback(void *this_ptr, Args... args) // static callback wrapper
  {
    reinterpret_cast<DM_callbacks_t*>(this_ptr)->dm_init_callback(args...);
  }
  void dm_init_callback(int iK, int iS, int *blacs_descr, void *data)
  {
    std::cout << "dm_init_callback" << std::endl;
    assert(("Only 1 k-point supported", iK==1));
    assert(("Only spin-paired case supported", iS==1));
    if ((blacs_descr==0) or (world_size == 1))
    {
      memcpy(data, &(dm[0]), sizeof(double)*n_basis*n_basis);
    }
    else
    {
      assert(("Scalapack not implemented. Use pXgemr2d for redistribution or pick smaller number of MPI processes.", false));
    }
  }
  
  template<typename... Args>
  static void static_dm_callback(void *this_ptr, Args... args)  // static callback wrapper
  {
    reinterpret_cast<DM_callbacks_t*>(this_ptr)->dm_callback(args...);
  }
  void dm_callback(int iK, int iS, int *blacs_descr, void *data)
  {
    std::cout << "dm_callback" << std::endl;
    assert(("Only 1 k-point supported", iK==1));
    assert(("Only spin-paired case supported", iS==1));
    if ((blacs_descr==0) or (world_size == 1))
    {
      memcpy(&(dm[0]), data, sizeof(double)*n_basis*n_basis);
    }
    else
    {
      assert(("Scalapack not implemented. Use pXgemr2d for redistribution or pick smaller number of MPI processes.", false));
    }
  }
};

int main(int argc, char *argv[])
{
  MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &mpi_provided_threading); // instead of MPI_Init
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  const MPI_Fint f_mpi_comm = MPI_Comm_c2f(MPI_COMM_WORLD);

  ASI_init("dummy", "aims.out", f_mpi_comm); // read geometry.in and control.in
  if (world_rank == 0) std::cout << "ASI_init complete" << std::endl;

  int  n_basis = ASI_get_basis_size();
  if (world_rank == 0) std::cout << "n_basis == " << n_basis << std::endl;
  
  DM_callbacks_t dm_callbacks(n_basis);  // register DM init callbacks here!
  loadtxt("dm.txt", n_basis*n_basis, dm_callbacks.dm); // load initial DM

  ASI_run();    // DO CALCULATIONS!
  if (world_rank == 0) std::cout << "ASI_run complete" << std::endl;
  
  savetxt("out_dm.txt", n_basis*n_basis, dm_callbacks.dm); // save final DM

  auto E = ASI_energy();
  if (world_rank == 0) std::cout << std::setprecision(6) << "Energy == " << E  << " Ha = " << E * 27.2113845 << " eV" << std::endl;

  ASI_finalize(); 
  MPI_Finalize();
  return 0;
}
