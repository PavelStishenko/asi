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
  
  Set coordinates of atoms. Coordinates are copied into interal storage, so the pointer do not have to be valid after return.
  \param coords pointer to array in row major format with shape [n_atoms, 3]. Contains new atomic coordinates.
  \param n_atoms number of atoms. Implementations may ignore this parameter if number if atoms is known at invocation moment.
*/
extern "C" void ASI_set_atom_coords(const double *coords, int n_atoms = -1);

/*!
  Set coordinates of atoms and lattice vectors
  
  Set coordinates of atoms and lattice vectors. Same as ASI_set_atom_coords() but with additional
  parameter for lattice vectors. Coordinates and lattive vectors are copied into interal storage, so the pointer do not have to be valid after return.
  
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
  Returns number of atoms
  
  Returns number of atoms
  \return number of atoms
*/
extern "C" int ASI_n_atoms();

/*!
  Returns total energy
  
  Returns total energy
  \return total energy
*/
extern "C" double ASI_energy();

/*!
  Returns forces acting on each atom
  
  Returns forces acting on each atom. May return `NULL` if forces calculation was not configured on ASI_init() call.
  For FHI-aims use `compute_forces .true.`. For DFTB+ use `Analysis { CalculateForces = Yes }`.
  
  \return pointer to array in row major format of shape `[ASI_n_atoms(), 3]` or `NULL`.
*/
extern "C" const double* ASI_forces();

/*!
  Returns stress tensor
  
  Returns stress tensor. May return `NULL` for non-periodic system or if stress calculation was not 
  configured on ASI_init() call.
  
  \return pointer to array in row major format of shape `[3, 3]` or `NULL`.
*/
extern "C" const double* ASI_stress();

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

/*!
  Returns number of spin channels.
  
  Returns number of spin channels. Returns 1 for spin-paired calculations.
  
  \return number of spin channels. Returns 1 for spin-paired calculations.
*/
extern "C" int ASI_get_nspin();

/*!
  Returns number of k-points.
  
  Returns number of k-points. Returns 1 for non-periodic or Gamma-point calculations.
  
  \return number of k-points. Returns 1 for non-periodic or Gamma-point calculations.
*/
extern "C" int ASI_get_nkpts();

/*!
  Returns number of (k-point, spin channel) pairs processed by current process.
  
  Meant to be used in case of distributed (MPI) calculation, when parallelization is performed along k-points 
  and spin channels. Returns number of (k-point, spin channel) pairs processed by current process.
  
  \return number of (k-point, spin channel) pairs processed by current process.
*/
extern "C" int ASI_get_n_local_ks();

/*!
  Returns list of (k-point, spin channel) pairs processed by current process.
  
  Meant to be used in case of distributed (MPI) calculation, when parallelization is performed along k-points 
  and spin channels. Returns list and number of (k-point, spin channel) pairs processed by current process.
  
  \param pointer to output array 
  \return number of (k-point, spin channel) pairs processed by current process.
*/
extern "C" int ASI_get_local_ks(int *local_ks); // local_ks[0, i] == i_kpnt, local_ks[1, i] == i_spin; i < ASI_get_n_local_ks();  iK and iS are 1-based indices

/*!
  Returns basis set size
  
  Returns basis set size. Defines total size of DM, H, and S matrices passed into callbacks of \ref ASI_dmhs_callback_t type.
  Note, that in case of distributed calculation local matrix size is defined by BLACS descriptor.
  \return basis set size. 
*/
extern "C" int ASI_get_basis_size();

/*!
  If hamiltonias, overlap, and density matrices are real.
  
  Returns True if hamiltonias, overlap, and density matrices are real. Returns False otherwise. 
  It defines type of matrices passed into callbacks of \ref ASI_dmhs_callback_t type.

  \return if hamiltonias, overlap, and density matrices are real.
*/
extern "C" bool ASI_is_hamiltonian_real();

/*!
  Register callback for density matrix evaluation.
  
  Register callback for density matrix evaluation. 
  Callback will be invoked on each density matrix evaluation, including evaluations within SCF loop.
  Note, that FHI-aims evaulates density matrix in SCF loop only if option `density_update_method density_matrix` is in config file.
  Note, that if FHI-aims evaulates forces then the last evaluations of the density matrix do not correspond to the calculated geometry.
  
  \param callback callback of type \ref ASI_dmhs_callback_t
  \param aux_ptr auxilary pointer that will be passed in every `callback` call.
*/
extern "C" void ASI_register_dm_callback(ASI_dmhs_callback_t callback, void *aux_ptr);

/*!
  Register callback for overlap matrix evaluation.
  
  Register callback for overlap matrix evaluation. 
  Callback will be invoked on each overlap matrix evaluation. Although overlap matrix is the same for all k-point and spin channels, 
  the callback will be invoked for each (k-point, spin channel) pair. It is up to client to decide when is the best time to work with the matrix. 
  Computational cost of callback invocation is negligible if no work is done in callback itself.
  
  \param callback callback of type \ref ASI_dmhs_callback_t
  \param aux_ptr auxilary pointer that will be passed in every `callback` call.
*/
extern "C" void ASI_register_overlap_callback(ASI_dmhs_callback_t , void *aux_ptr);

/*!
  Register callback for hamiltonian matrix evaluation.
  
  Register callback for hamiltonian matrix evaluation. 
  Callback will be invoked on each hamiltonian matrix evaluation.
  
  \param callback callback of type \ref ASI_dmhs_callback_t
  \param aux_ptr auxilary pointer that will be passed in every `callback` call.
*/
extern "C" void ASI_register_hamiltonian_callback(ASI_dmhs_callback_t , void *aux_ptr);

/*!
  Register callback for density matrix initial guess.
  
  Register callback for density matrix initial guess. Callback is expected to fill the matrix argument with initial guess of the density matrix.
  It is a good place for calculation restart or density matrix prediction.
  
  \param callback callback of type \ref ASI_dmhs_callback_t
  \param aux_ptr auxilary pointer that will be passed in every `callback` call.
*/
extern "C" void ASI_register_dm_init_callback(ASI_dmhs_callback_t , void *aux_ptr);

/*!
  Finalize calculation.
  
  Free memory, close openned files. Library may be unloaded after that. No API calls are permited after that call.
*/
extern "C" void ASI_finalize();
// no calls after ASI_finalize

#endif // ASI_H
