#ifndef ASI_H
#define ASI_H

#define ASI_FHI_AIMS	1
#define ASI_DFTBP 	2

typedef void (*ASI_ext_pot_func_t)(void *aux_ptr, int n, const double *coords, double *potential, double *potential_grad);
typedef void (*ASI_dmhs_callback_t)(void *aux_ptr, int iK, int iS, int *blacs_descr, void *blacs_data); // iK and iS are 1-based indices

extern "C" int ASI_flavour();

// generate all input files before ASI_init; no calls before ASI_init
extern "C" void ASI_init(const char *inputpath, const char *outputfilename, int mpiComm); // mpiComm must be obtained from MPI_Comm_c2f(); it is recommended to read input files here

extern "C" void ASI_set_atom_coords(const double *coords, int n_atoms = -1); // AIMS: cannot be used after ASI_init
extern "C" void ASI_set_geometry(const double *coords, int n_atoms = -1, const double *lattice=0); // AIMS: cannot be used after ASI_init
extern "C" void ASI_set_external_potential(ASI_ext_pot_func_t, void *); // potential of positive charges; is called *before* ASI_run(), calls callback on its own.
extern "C" void ASI_register_external_potential(ASI_ext_pot_func_t, void *);  // potential of positive charges;  is called *during* ASI_run()

extern "C" void ASI_run(); // do long-time calculations; input files also can be read here but it is discouraged

extern "C" int ASI_n_atoms();
extern "C" double ASI_energy();
extern "C" const double* ASI_forces(); // may return NULL pointer is forces are not available AIMS: use `compute_forces .true.` DFTB+: use `Analysis { CalculateForces = Yes }`

extern "C" const double* ASI_atomic_charges(int scheme=-1); // AIMS:  Hirschfeld by default; DFTB+ Mulliken by default
extern "C" void ASI_calc_esp(int n, const double *coords, double *potential, double *potential_grad); // calculate electrostatic potential (ESP) for positive charges in specified coordinates

extern "C" int ASI_get_nspin();
extern "C" int ASI_get_nkpts();
extern "C" int ASI_get_n_local_ks();
extern "C" int ASI_get_local_ks(int *local_ks); // local_ks[0, i] == i_kpnt, local_ks[1, i] == i_spin; i < ASI_get_n_local_ks();  iK and iS are 1-based indices
extern "C" int ASI_get_basis_size();
extern "C" bool ASI_is_hamiltonian_real();
extern "C" void ASI_register_dm_callback(ASI_dmhs_callback_t , void *aux_ptr); // AIMS: `density_update_method density_matrix` option is necessary; cannot be used if forces are computed (returns garbage)
extern "C" void ASI_register_overlap_callback(ASI_dmhs_callback_t , void *aux_ptr);
extern "C" void ASI_register_hamiltonian_callback(ASI_dmhs_callback_t , void *aux_ptr);
extern "C" void ASI_register_dm_init_callback(ASI_dmhs_callback_t , void *aux_ptr); // AIMS: elsi_restart read is necessary

extern "C" void ASI_finalize();
// no calls after ASI_finalize

#endif // ASI_H
