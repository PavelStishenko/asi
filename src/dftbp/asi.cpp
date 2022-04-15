//#include <cstdbool>
#define _Bool bool
#include <dftbplus.h>
#include <cstdlib>
#include <cassert>
#include <cmath>
#include <vector>
#include <array>
#include <numeric>
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <complex>

#include "asi.h"

using namespace std;

const double bohr = 0.52917721; // bohrs in 1 Ang.
const double hartree = 1./27.2113845; // Hartree in 1 eV.

bool ASI_initialized = false;
DftbPlus calculator;
DftbPlusInput input;
int major, minor, patch;
std::vector<std::array<double, 3> > atom_coords;
std::vector<std::array<double, 3> > saved_total_forces;
std::vector<std::array<double, 3> > saved_atomic_charges;

ASI_ext_pot_func_t ext_pot_func = 0;
void *ext_pot_func_aux_ptr = 0;

void ext_pot_func_dftbp(void *refptr, double * /*dqatom unused*/, double *extpotatom);
void ext_pot_grad_func_dftbp(void *refptr, double */*dqatom unused*/, double *extpotatomgrad);

int ASI_flavour()
{
  return ASI_DFTBP;
}

void ASI_init(const char *inputpath, const char *outputfilename, int mpiComm)
{
  dftbp_api(&major, &minor, &patch);
  dftbp_init_mpi(&calculator, outputfilename, mpiComm);
  dftbp_get_input_from_file(&calculator, inputpath, &input);
  dftbp_process_input(&calculator, &input);
  
  ASI_initialized = true;
}

int ASI_get_basis_size()
{
  return dftbp_get_basis_size(&calculator);
}

bool ASI_is_hamiltonian_real()
{
  return dftbp_is_hs_real(&calculator);
}

void ASI_register_dm_callback(ASI_dmhs_callback_t dm_callback, void *aux_ptr)
{
  dftbp_register_dm_callback(&calculator, dm_callback, aux_ptr);
}

void ASI_register_overlap_callback(ASI_dmhs_callback_t hs_callback, void *aux_ptr)
{
  dftbp_register_s_callback(&calculator, hs_callback, aux_ptr);
}

void ASI_register_hamiltonian_callback(ASI_dmhs_callback_t hs_callback, void *aux_ptr)
{
  dftbp_register_h_callback(&calculator, hs_callback, aux_ptr);
}

/*
  coords assumed in Bohrs
*/
void ASI_set_atom_coords(const double *coords, int /*n_atoms unused*/)
{
  const int n = ASI_n_atoms();
  atom_coords.resize(n);
 
  // copy and convert in bohrs // NOT CONVERT
  std::transform(coords, coords + 3*n, atom_coords[0].data(), [](double x){return x/1.0;} );
  
  dftbp_set_coords(&calculator, atom_coords[0].data());
  
  // save in Angstroms
  std::copy(coords, coords + 3*n, atom_coords[0].data());

}

/*
  coords assumed in Bohrs
*/
void ASI_set_geometry(const double *coords, int n_atoms, const double *lattice)
{
  const int n = ASI_n_atoms();
  assert(n == n_atoms);
  
  atom_coords.resize(n);
 
  // copy and convert in bohrs // NOT CONVERT
  std::transform(coords, coords + 3*n, atom_coords[0].data(), [](double x){return x/1.0;} );
  
  //dftbp_set_coords(&calculator, atom_coords[0].data());
  dftbp_set_coords_and_lattice_vecs(&calculator, atom_coords[0].data(), lattice);
  
  // save in Angstroms
  std::copy(coords, coords + 3*n, atom_coords[0].data());

}


void ASI_set_external_potential(ASI_ext_pot_func_t ext_pot_func, void * ext_pot_func_aux_ptr)
{
  const int n = ASI_n_atoms();
  assert(atom_coords.size() == n);
  std::vector<double> extpot(n);
  std::vector<std::array<double, 3> > extpotgrad(n);
  
  ext_pot_func(ext_pot_func_aux_ptr, n, atom_coords[0].data(), extpot.data(), extpotgrad[0].data());
  
  for (int i = 0; i < n; ++i)
  {
    extpot[i] *= -1;
    extpotgrad[i][0] *= -1;
    extpotgrad[i][1] *= -1;
    extpotgrad[i][2] *= -1;
  }
  
  dftbp_set_external_potential(&calculator, extpot.data(), extpotgrad[0].data());
}

void ASI_register_external_potential(ASI_ext_pot_func_t _ext_pot_func, void * _ext_pot_func_aux_ptr)
{
  ext_pot_func = _ext_pot_func;
  ext_pot_func_aux_ptr = _ext_pot_func_aux_ptr;
  
  dftbp_register_ext_pot_generator(&calculator, ext_pot_func_aux_ptr, ext_pot_func_dftbp, ext_pot_grad_func_dftbp);
}

void ASI_run()
{
  assert(ASI_initialized && "ASI not initialized");

  double mermin_energy;
  dftbp_get_energy(&calculator, &mermin_energy); 
}

double ASI_energy()
{
  double mermin_energy;
  dftbp_get_energy(&calculator, &mermin_energy);
  return mermin_energy;
}

int ASI_n_atoms()
{
  return dftbp_get_nr_atoms(&calculator);
}

int ASI_get_nspin()
{
  return dftbp_get_nr_spin(&calculator);
}

int ASI_get_nkpts()
{
  return dftbp_nr_kpoints(&calculator);
}

int ASI_get_n_local_ks()
{
  return dftbp_get_nr_local_ks(&calculator);
}

int ASI_get_local_ks(int *local_ks)
{
  return dftbp_get_local_ks(&calculator, local_ks);
}

void ASI_get_dm(int i_spin, int i_kpnt, int **dm_desc, double **dm)
{
}

void ASI_set_dm(int **dm_descs, double **dms, int n_kpts) //  TODO: get rid of n_kpts here
{
}


const double * ASI_forces()
{
  saved_total_forces.resize(ASI_n_atoms());
  dftbp_get_gradients(&calculator, saved_total_forces[0].data());
  for (size_t i = 0; i < saved_total_forces.size(); ++i)
  {
    saved_total_forces.at(i)[0] *= -1;
    saved_total_forces.at(i)[1] *= -1;
    saved_total_forces.at(i)[2] *= -1;
  }
  return saved_total_forces[0].data();
}

const double* ASI_atomic_charges(int scheme)
{
  saved_atomic_charges.resize(ASI_n_atoms());
  dftbp_get_gross_charges(&calculator, saved_atomic_charges[0].data());
  return saved_atomic_charges[0].data();
}

void ASI_calc_esp(int n, const double *coords, double *potential, double *potential_grad)
{
  std::vector<double> potential_buf; // empty for a while
  
  if (not potential)
  {
    potential_buf.resize(n);
    potential = &(potential_buf[0]);
  }

  dftbp_get_elstat_potential(&calculator, n, potential, coords); //need it in any way
  for(int i = 0; i < n; ++i)
  {
      potential[i] *= -1; // DFTB+ operates by electronic potential, so invertion is necessary to get conventional protonic potential
  }

  if (potential_grad) // TODO: replace with dftbp_get_potential_gradient
  {
    const double d = 0.001;
    std::vector<std::array<double, 3> > grad_coords(n*3);
    for(int i = 0; i < n; ++i)
    {
      for(int j = 0; j < 3; ++j)
      {
        grad_coords[i*3 + j] = {coords[i*3 + 0], coords[i*3 + 1], coords[i*3 + 2]};
        grad_coords[i*3 + j][j] += d;
      }
    }
    dftbp_get_elstat_potential(&calculator, n*3, potential_grad, grad_coords[0].data()); 
    for(int i = 0; i < n; ++i)
    {
      for(int j = 0; j < 3; ++j)
      {
        potential_grad[i*3 + j] *=  -1;
        potential_grad[i*3 + j] = (potential_grad[i*3 + j] - potential[i]) / d;
      }
    }
  }
}

void ASI_get_esp(int *n, double **coords, double **potential, double **potential_grad)
{
  *n = 0;
}

void ASI_finalize()
{
  dftbp_final(&calculator);
}

// Local functions


void ext_pot_func_dftbp(void *refptr, double * /*dqatom unused*/, double *extpotatom)
{
  int n = dftbp_get_nr_atoms(&calculator);
  assert(refptr == ext_pot_func_aux_ptr); // just in case
  ext_pot_func(ext_pot_func_aux_ptr, n, atom_coords[0].data(), extpotatom, 0);
  for (int i = 0; i < n; ++i)
  {
    extpotatom[i] *= -1;
  }
}

void ext_pot_grad_func_dftbp(void *refptr, double */*dqatom unused*/, double *extpotatomgrad)
{
  int n = dftbp_get_nr_atoms(&calculator);
  assert(refptr == ext_pot_func_aux_ptr); // just in case
  ext_pot_func(ext_pot_func_aux_ptr, n, atom_coords[0].data(), 0, extpotatomgrad);
  for (int i = 0; i < n*3; ++i)
  {
    extpotatomgrad[i] *= -1;
  }  
}


