#include <complex>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>

#include <asi.h>

#include <mpi.h>

#include <cassert>
#include <cmath>
#include <cstdint>
#include <ctime>


#include "codespec.hpp"
#include "utils.hpp"

int mpi_provided_threading, world_size, world_rank;



double homog_pot(void *aux_ptr, double x, double y, double z)
{
  double *hom_field = reinterpret_cast<double*>(aux_ptr);
  double P = -(x * hom_field[0] + y * hom_field[1] + z * hom_field[2]);
  //std::cout << "homog_pot @ "  << x << " " << y << " " << z << " = " << P << std::endl;
  return P;
}

double pnt_pot_1(void *aux_ptr, double x, double y, double z)
{
  double *charge_and_coords = reinterpret_cast<double*>(aux_ptr);
  
  //std::cout << charge_and_coords[0] << " " << charge_and_coords[1] << " " << charge_and_coords[2] << " " << charge_and_coords[3] << std::endl;
  const double chg = charge_and_coords[0];
  const double *coords = &(charge_and_coords[1]);
  const double r = sqrt((x - coords[0])*(x - coords[0]) + (y - coords[1])*(y - coords[1]) + (z - coords[2])*(z - coords[2]));
  const double P = chg / r;
  //std::cout << "pnt_pot @ "  << x << " " << y << " " << z << " = " << P << std::endl;
  return P;
}

void pnt_pot(void *aux_ptr, int n, const double *coords, double *potential, double *potential_grad)
{
  const double d = 0.0001;
  for(int i  = 0; i < n; ++i)
  {
    double p = pnt_pot_1(aux_ptr, coords[i*3 + 0], coords[i*3 + 1], coords[i*3 + 2]);
    if (potential)
    {
      potential[i] = p;
    }
    
    if (potential_grad)
    {
      potential_grad[i*3 + 0] = (pnt_pot_1(aux_ptr, coords[i*3 + 0] + d, coords[i*3 + 1], coords[i*3 + 2]) - p) / d;
      potential_grad[i*3 + 1] = (pnt_pot_1(aux_ptr, coords[i*3 + 0], coords[i*3 + 1] + d, coords[i*3 + 2]) - p) / d;
      potential_grad[i*3 + 2] = (pnt_pot_1(aux_ptr, coords[i*3 + 0], coords[i*3 + 1], coords[i*3 + 2] + d) - p) / d;
    }
  }
}

std::vector<double> h2o_coords()
{
  std::vector<double> coords {
      1, 0, -1,
      1,  0.78306400, 0,
      1, -0.78306400, 0};
  for (auto &x : coords)
  {
    x /= 0.52917721; // convert to Ang -> Bohr
  }
  return coords;
}

int main(int argc, char *argv[])
{
  MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &mpi_provided_threading); // instead of MPI_Init
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
 
  assert(argc == 1);

  
  const MPI_Fint f_mpi_comm = MPI_Comm_c2f(MPI_COMM_WORLD);

  make_config_files(false);

  //ASI_set_atom_coords(h2o_coords().data(), 3);
  ASI_init("dftb_in.hsd", "asi.log", f_mpi_comm);  // input file ignored for AIMS
  ASI_set_atom_coords(h2o_coords().data(), 3);

  double esp_field[3] = {0, 0.01, 0};
  double pnt_charge_and_coords[4] = {1.0, 1, 0, 2.0};


  {
    //ASI_set_external_potential(homog_pot, esp_field);
    //ASI_register_external_potential(homog_pot, esp_field);
    ASI_register_external_potential(pnt_pot, pnt_charge_and_coords);
    //ASI_set_external_potential(pnt_pot, pnt_charge_and_coords);
  }

  if (world_rank == 0) std::cout << "ASI_run() ..." << std::endl;
  auto t0 = rdtsc();
  ASI_run();
  if (world_rank == 0) std::cout << "ASI_run() DONE!" << std::endl;

  t0 = rdtsc();
  auto E = ASI_energy();
  std::cout << std::setprecision(6) << std::fixed;
  if (world_rank == 0) std::cout << "Energy == " << E  << " Ha = " << E * 27.2113845 << " eV" << std::endl;
  
  const double *forces = ASI_forces();
  if (world_rank == 0) {
    std::cout << "Forces == ";
    for (size_t i = 0; i < 3 * ASI_n_atoms(); ++i )
      std::cout << forces[i] << " ";
    std::cout << std::endl;
  }

  const double *atomic_charges = ASI_atomic_charges();
  assert(atomic_charges != 0);
  if (world_rank == 0) {
    std::cout << "atomic_charges == ";
    for (size_t i = 0; i < ASI_n_atoms(); ++i )
      std::cout << atomic_charges[i] << " ";
    std::cout << std::endl;
  }

  double esp_coords[2*3] = { 10, 10, 10,    1,  2,  0}; 
  double esp[2];
  
  ASI_calc_esp(2, esp_coords, esp, 0);
  
  if (world_rank == 0) std::cout << "EP1=" << esp[0] << " = " << esp[0] * 27.2113845 << " V" << std::endl;
  if (world_rank == 0) std::cout << "EP2=" << esp[1] << " = " << esp[1] * 27.2113845 << " V"  << std::endl;
  
  assert_all(E, forces, atomic_charges, esp);
  
  ASI_finalize(); 
  MPI_Finalize();
  return 0;
}
