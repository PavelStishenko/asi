#ifndef ASI_H
#define ASI_H

/*! Constant returned by ASI_flavour() for FHI-aims implementaion */
#define ASI_FHI_AIMS	1
/*! Constant returned by ASI_flavour() for DFTB+ implementaion */
#define ASI_DFTBP 	2
/// @file

/*!
  Signature of callback functions for electrostatic potential (ESP) import.
  
  Signature of callback functions for electrostatic potential (ESP) import.
  Meant to be used with ASI_set_external_potential() or ASI_register_external_potential() functions.
  ESP sign is supposed to be as filled by positive charge.
  
  \param aux_ptr pointer to auxilary object specified on callback registration.
  \param n number of points where the ESP is requested.
  \param coords pointer to array in row major format with shape `[n, 3]`. Contains coordinates of points where ESP is requested.
  \param potential pointer to array of size `n` or `NULL`. Points to output array, where ESP potential values should be placed.
                    If NULL, then ESP is not requested, that may be if ESP gradient is requested.
  \param potential_grad pointer to array in row major format with shape `[n, 3]` or `NULL`. 
                    Points to output array, where ESP potential gradient values should be placed. 
                    If NULL, then ESP gradient is not requested, that may be if ESP only is requested.
*/
typedef void (*ASI_ext_pot_func_t)(void *aux_ptr, int n, const double *coords, double *potential, double *potential_grad);

/*!
  Signature of callback functions for export of density, Hamiltonian or overlap matrices.
  
  Signature of callback functions for export of density (DM), Hamiltonian (H), or overlap (S) matrices.
  
  \param aux_ptr pointer to auxilary object specified on callback registration.
  \param iK k-point index, 1-based.
  \param iS spin channel index, 1-based.
  \param blacs_descr pointer to BLACS array descriptor or `NULL`.
  \param blacs_data pointer to the exported matrix. Data type of the matrix depends on ASI_is_hamiltonian_real() value.
               Size of the matrix  depends on `blacs_descr` content or ASI_get_basis_size() value.
*/
typedef void (*ASI_dmhs_callback_t)(void *aux_ptr, int iK, int iS, int *blacs_descr, void *blacs_data); // iK and iS are 1-based indices

/*!
  Returns integer constant to unique for each API implementaion
  \return Unique implementaion identifier
*/
extern "C" int ASI_flavour();

/*!
  Initialize calculations
  
  User must generate all input files before ASI_init().
  Implementations must read input files here.
  No ASI API calls permited before this function.
  
  \param inputpath path to file or directory with input files (configs, initial geometry, etc)
  \param outputfilename path to output file
  \param mpiComm MPI communicator in Fortran format. Can be obtained from MPI_Comm_c2f().
*/
extern "C" void ASI_init(const char *inputpath, const char *outputfilename, int mpiComm);

/*!
  Set coordinates of atoms
  
  Set coordinates of atoms.
  \param coords pointer to array in row major format with shape [n_atoms, 3]. Contains new atomic coordinates.
  \param n_atoms number of atoms. Implementations may ignore this parameter if number if atoms is known at invocation moment.
*/
extern "C" void ASI_set_atom_coords(const double *coords, int n_atoms = -1);

/*!
  Set coordinates of atoms and lattice vectors
  
  Set coordinates of atoms and lattice vectors. Same as ASI_set_atom_coords() but with additional
  parameter for lattice vectors.
  
  \param coords pointer to array in row major format with shape [n_atoms, 3]. Contains new atomic coordinates.
  \param n_atoms number of atoms. Implementations may ignore this parameter if number if atoms is known at invocation moment.
  \param lattice pointer to array in row major format with shape [3, 3]. Contains new lattice vectors with vectors as rows.
*/
extern "C" void ASI_set_geometry(const double *coords, int n_atoms = -1, const double *lattice=0);

/*!
  Set external electrostatic potential (ESP).
  
  Set external electrostatic potential (ESP). The `callback` is not stored and will be invoked only before the function returns.
  
  \param callback callback function of \ref ASI_ext_pot_func_t type.
  \param aux_ptr auxilary pointer that will be passed in every `callback` call.
*/
extern "C" void ASI_set_external_potential(ASI_ext_pot_func_t callback, void *aux_ptr);

/*!
  Register callback for external electrostatic potential (ESP) calculation.
  
  Register callback for external electrostatic potential (ESP) calculation. The `callback` will be stored and will be invoked during ASI_run() call.
  
  \param callback callback function of \ref ASI_ext_pot_func_t type.
  \param aux_ptr auxilary pointer that will be passed in every `callback` call.
*/
extern "C" void ASI_register_external_potential(ASI_ext_pot_func_t callback, void *aux_ptr);

/*!
  Do actual calculations.
  
  Do actual calculations. Long-time calculations supposed to be done here. 
  Functions such as ASI_energy(), ASI_forces(), ASI_atomic_charges() must be invoked only after ASI_run()
*/
extern "C" void ASI_run();

/*!
  Return number of atoms
  
  Return number of atoms
  \return number of atoms
*/
extern "C" int ASI_n_atoms();

/*!
  Return total energy
  
  Return total energy
  \return total energy
*/
extern "C" double ASI_energy();

/*!
  Return forces acting on each atom
  
  Return forces acting on each atom. May return `NULL` if forces calculation was not configured on ASI_init() call.
  For FHI-aims use `compute_forces .true.`. For DFTB+ use `Analysis { CalculateForces = Yes }`.
  
  \return pointer to array in row major format of shape `[ASI_n_atoms(), 3]` or `NULL`.
*/
extern "C" const double* ASI_forces();

/*!
  Return atomic charges
  
  Return atomic charges calculated by requested algorithm (charge partitioning scheme). 
  May return `NULL` if charge.partitioning was not configured on ASI_init() call.
  
  \return pointer to array of size `ASI_n_atoms()` or `NULL`.
*/
extern "C" const double* ASI_atomic_charges(int scheme=-1);

/*!
  Calculate electrostatic potential (ESP) for positive charges in specified coordinates.
  
  Calculate electrostatic potential (ESP) and/or ESP gradient for positive charges in specified coordinates.
  Parameters are similar to \ref ASI_ext_pot_func_t callback type.
  
  \param n number of points for which ESP calculation is requested
  \param coords pointer to array in row major format with shape `[n, 3]`. Contains coordinates of points where ESP is requested.
  \param potential pointer to array of size `n` or `NULL`. Points to output array, where ESP potential values should be placed.
                    If NULL, then ESP is not requested, that may be if ESP gradient is requested.
  \param potential_grad pointer to array in row major format with shape `[n, 3]` or `NULL`. 
                    Points to output array, where ESP potential gradient values should be placed. 
                    If NULL, then ESP gradient is not requested, that may be if ESP only is requested.
*/
extern "C" void ASI_calc_esp(int n, const double *coords, double *potential, double *potential_grad);

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
