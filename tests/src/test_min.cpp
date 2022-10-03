#include <iostream>
#include <iomanip>
#include <mpi.h>
#include <cassert>

#include "asi.h"

int mpi_provided_threading, world_size, world_rank;


int main(int argc, char *argv[])
{
  MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &mpi_provided_threading); // instead of MPI_Init
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  const MPI_Fint f_mpi_comm = MPI_Comm_c2f(MPI_COMM_WORLD);

  ASI_init("dftb_in.hsd", "asi.log", f_mpi_comm); // read geometry.in and control.in OR dftb_in.hsd
  if (world_rank == 0) std::cout << "ASI_init complete" << std::endl;

  int  n_basis = ASI_get_basis_size();
  if (world_rank == 0) std::cout << "n_basis == " << n_basis << std::endl;
  
  ASI_run();    // DO CALCULATIONS!
  if (world_rank == 0) std::cout << "ASI_run complete" << std::endl;

  auto E = ASI_energy();
  if (world_rank == 0) std::cout << std::setprecision(6) << "Energy == " << E  << " Ha = " << E * 27.2113845 << " eV" << std::endl;

  ASI_finalize(); 
  MPI_Finalize();
  return 0;
}
