// AMReX Interface for NGA2
// C++ bridge between NGA2 Fortran and AMReX C++
// Provides amrcore_* and amrmfab_* functions callable from Fortran

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_FAmrCore.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_PhysBCFunct.H>
#include <AMReX_PlotFileUtil.H>
#ifdef AMREX_USE_HDF5
#include <AMReX_PlotFileUtilHDF5.H>
#endif

namespace nga2 {

//=============================================================================
// Callback function types - these are bind(c) Fortran procedures
// They receive the owner pointer and handle dispatching in Fortran
//=============================================================================
extern "C" {
// Level callbacks: owner, level, time, boxarray ptr, distromap ptr
using LevelDispatcher = void (*)(void *owner, int lev, double time,
                                 const amrex::BoxArray *,
                                 const amrex::DistributionMapping *);
// Clear callback: owner, level
using ClearDispatcher = void (*)(void *owner, int lev);
// Tag callback: owner, level, tagboxarray ptr, time, tagval, clearval
using TagDispatcher = void (*)(void *owner, int lev, amrex::TagBoxArray *,
                               double time, char tagval, char clearval);
// Post-regrid callback: owner only
using PostRegridDispatcher = void (*)(void *owner);
}

//=============================================================================
// NGA2AmrCore - extends FAmrCore with Fortran dispatcher callbacks
// Each amrcore has ONE owner and ONE dispatcher per event type
//=============================================================================
class NGA2AmrCore : public amrex::FAmrCore {
public:
  void *owner = nullptr;

  LevelDispatcher on_init_dispatch = nullptr;
  LevelDispatcher on_coarse_dispatch = nullptr;
  LevelDispatcher on_remake_dispatch = nullptr;
  ClearDispatcher on_clear_dispatch = nullptr;
  TagDispatcher on_tag_dispatch = nullptr;
  PostRegridDispatcher on_postregrid_dispatch = nullptr;

  // Public override of regrid to call post-regrid dispatch after base class
  void regrid(int lbase, amrex::Real time, bool initial = false) override {
    amrex::FAmrCore::regrid(lbase, time, initial);
    if (on_postregrid_dispatch && owner) {
      on_postregrid_dispatch(owner);
    }
  }

protected:
  void MakeNewLevelFromScratch(int lev, amrex::Real time,
                               const amrex::BoxArray &ba,
                               const amrex::DistributionMapping &dm) override {
    if (on_init_dispatch && owner) {
      on_init_dispatch(owner, lev, time, &ba, &dm);
    }
  }

  void MakeNewLevelFromCoarse(int lev, amrex::Real time,
                              const amrex::BoxArray &ba,
                              const amrex::DistributionMapping &dm) override {
    if (on_coarse_dispatch && owner) {
      on_coarse_dispatch(owner, lev, time, &ba, &dm);
    }
  }

  void RemakeLevel(int lev, amrex::Real time, const amrex::BoxArray &ba,
                   const amrex::DistributionMapping &dm) override {
    if (on_remake_dispatch && owner) {
      on_remake_dispatch(owner, lev, time, &ba, &dm);
    }
  }

  void ClearLevel(int lev) override {
    if (on_clear_dispatch && owner) {
      on_clear_dispatch(owner, lev);
    }
  }

  void ErrorEst(int lev, amrex::TagBoxArray &tags, amrex::Real time,
                int /*ngrow*/) override {
    if (on_tag_dispatch && owner) {
      on_tag_dispatch(owner, lev, &tags, time, amrex::TagBox::SET,
                      amrex::TagBox::CLEAR);
    }
  }
};

//=============================================================================
// FillPatch BC Callback - receives solver context, can access all fields
//=============================================================================

// Fortran BC callback signature: receives solver context, multifab, geometry,
// time
using FillPatchBCDispatcher = void (*)(void *solver_ctx, amrex::MultiFab *mf,
                                       const amrex::Geometry *geom,
                                       amrex::Real time, int scomp, int ncomp);

// Functor that AMReX FillPatch calls; forwards to Fortran with context
class NGA2BCFunctor {
public:
  void *solver_ctx = nullptr;
  FillPatchBCDispatcher bc_dispatch = nullptr;
  const amrex::Geometry *geom = nullptr;

  NGA2BCFunctor() = default;
  NGA2BCFunctor(void *ctx, FillPatchBCDispatcher f, const amrex::Geometry *g)
      : solver_ctx(ctx), bc_dispatch(f), geom(g) {}

  // Called by AMReX FillPatch when physical BCs needed
  void operator()(amrex::MultiFab &mf, int icomp, int ncomp,
                  amrex::IntVect const & /*nghost*/, amrex::Real time,
                  int /*bccomp*/) {
    if (bc_dispatch && solver_ctx) {
      bc_dispatch(solver_ctx, &mf, geom, time, icomp, ncomp);
    }
  }
};

} // namespace nga2

//=============================================================================
// C Interface for Fortran - all functions use amrcore_* prefix
//=============================================================================
extern "C" {

//-----------------------------------------------------------------------------
// Lifecycle
//-----------------------------------------------------------------------------
void *amrcore_create() { return new nga2::NGA2AmrCore(); }

void amrcore_destroy(void *core) {
  delete static_cast<nga2::NGA2AmrCore *>(core);
}

//-----------------------------------------------------------------------------
// Set Owner
//-----------------------------------------------------------------------------
void amrcore_set_owner(void *core, void *owner) {
  static_cast<nga2::NGA2AmrCore *>(core)->owner = owner;
}

//-----------------------------------------------------------------------------
// Set Dispatchers
//-----------------------------------------------------------------------------
void amrcore_set_on_init_dispatch(void *core, nga2::LevelDispatcher f) {
  static_cast<nga2::NGA2AmrCore *>(core)->on_init_dispatch = f;
}

void amrcore_set_on_coarse_dispatch(void *core, nga2::LevelDispatcher f) {
  static_cast<nga2::NGA2AmrCore *>(core)->on_coarse_dispatch = f;
}

void amrcore_set_on_remake_dispatch(void *core, nga2::LevelDispatcher f) {
  static_cast<nga2::NGA2AmrCore *>(core)->on_remake_dispatch = f;
}

void amrcore_set_on_clear_dispatch(void *core, nga2::ClearDispatcher f) {
  static_cast<nga2::NGA2AmrCore *>(core)->on_clear_dispatch = f;
}

void amrcore_set_on_tag_dispatch(void *core, nga2::TagDispatcher f) {
  static_cast<nga2::NGA2AmrCore *>(core)->on_tag_dispatch = f;
}

void amrcore_set_on_postregrid_dispatch(void *core,
                                        nga2::PostRegridDispatcher f) {
  static_cast<nga2::NGA2AmrCore *>(core)->on_postregrid_dispatch = f;
}

//-----------------------------------------------------------------------------
// Grid Operations
//-----------------------------------------------------------------------------
void amrcore_init_from_scratch(void *core, double time) {
  static_cast<nga2::NGA2AmrCore *>(core)->InitFromScratch(time);
}

void amrcore_regrid(void *core, int base_level, double time) {
  static_cast<nga2::NGA2AmrCore *>(core)->regrid(base_level, time);
}

//-----------------------------------------------------------------------------
// Level Queries
//-----------------------------------------------------------------------------
int amrcore_finest_level(void *core) {
  return static_cast<nga2::NGA2AmrCore *>(core)->finestLevel();
}

int amrcore_max_level(void *core) {
  return static_cast<nga2::NGA2AmrCore *>(core)->maxLevel();
}

void amrcore_set_finest_level(void *core, int lev) {
  static_cast<nga2::NGA2AmrCore *>(core)->SetFinestLevel(lev);
}

//-----------------------------------------------------------------------------
// Geometry Access
//-----------------------------------------------------------------------------
void amrcore_get_geometry(void **geom_ptr, int lev, void *core) {
  auto *amr = static_cast<nga2::NGA2AmrCore *>(core);
  *geom_ptr = const_cast<amrex::Geometry *>(&(amr->Geom(lev)));
}

void amrcore_get_ref_ratio(int *ref_ratio, void *core) {
  auto *amr = static_cast<nga2::NGA2AmrCore *>(core);
  for (int lev = 0; lev < amr->maxLevel(); ++lev) {
    ref_ratio[lev] = amr->refRatio(lev)[0];
  }
}

//-----------------------------------------------------------------------------
// BoxArray and DistributionMapping Access
//-----------------------------------------------------------------------------
void amrcore_get_boxarray(void **ba_ptr, int lev, void *core) {
  auto *amr = static_cast<nga2::NGA2AmrCore *>(core);
  *ba_ptr = const_cast<amrex::BoxArray *>(&(amr->boxArray(lev)));
}

void amrcore_get_distromap(void **dm_ptr, int lev, void *core) {
  auto *amr = static_cast<nga2::NGA2AmrCore *>(core);
  *dm_ptr =
      const_cast<amrex::DistributionMapping *>(&(amr->DistributionMap(lev)));
}

//=============================================================================
// MultiFab Operations - amrmfab_* prefix
//=============================================================================

//-----------------------------------------------------------------------------
// FillPatch for Level 0 (single level, just physical BCs)
// Note: scomp/dcomp are 1-indexed (Fortran convention), converted to 0-indexed
// here
//-----------------------------------------------------------------------------
void amrmfab_fillpatch_single(void *mf_ptr, double time_old, void *mf_old_ptr,
                              double time_new, void *mf_new_ptr, void *geom_ptr,
                              void *solver_ctx,
                              nga2::FillPatchBCDispatcher bc_dispatch,
                              double time, int scomp, int dcomp, int ncomp) {
  auto *mf = static_cast<amrex::MultiFab *>(mf_ptr);
  auto *mf_old = static_cast<amrex::MultiFab *>(mf_old_ptr);
  auto *mf_new = static_cast<amrex::MultiFab *>(mf_new_ptr);
  auto *geom = static_cast<amrex::Geometry *>(geom_ptr);

  amrex::Vector<amrex::MultiFab *> smf = {mf_old, mf_new};
  amrex::Vector<amrex::Real> stime = {time_old, time_new};

  nga2::NGA2BCFunctor bc_functor(solver_ctx, bc_dispatch, geom);

  // Convert from 1-indexed (Fortran) to 0-indexed (C++)
  amrex::FillPatchSingleLevel(*mf, time, smf, stime, scomp - 1, dcomp - 1,
                              ncomp, *geom, bc_functor, 0);
}
//-----------------------------------------------------------------------------
// FillPatch for Fine Levels (two levels: coarse + fine)
// Note: scomp/dcomp are 1-indexed (Fortran convention), converted to 0-indexed
// here
//-----------------------------------------------------------------------------
void amrmfab_fillpatch_two(void *mf_ptr, double time_old_c, void *mf_old_c_ptr,
                           double time_new_c, void *mf_new_c_ptr,
                           void *geom_c_ptr, double time_old_f,
                           void *mf_old_f_ptr, double time_new_f,
                           void *mf_new_f_ptr, void *geom_f_ptr,
                           void *solver_ctx,
                           nga2::FillPatchBCDispatcher bc_dispatch, double time,
                           int scomp, int dcomp, int ncomp, int ref_ratio,
                           int interp_type, int *lo_bc, int *hi_bc, int nbc) {
  auto *mf = static_cast<amrex::MultiFab *>(mf_ptr);
  auto *mf_old_c = static_cast<amrex::MultiFab *>(mf_old_c_ptr);
  auto *mf_new_c = static_cast<amrex::MultiFab *>(mf_new_c_ptr);
  auto *mf_old_f = static_cast<amrex::MultiFab *>(mf_old_f_ptr);
  auto *mf_new_f = static_cast<amrex::MultiFab *>(mf_new_f_ptr);
  auto *geom_c = static_cast<amrex::Geometry *>(geom_c_ptr);
  auto *geom_f = static_cast<amrex::Geometry *>(geom_f_ptr);

  // Look up interpolator from type constant (matches AMReX Fortran constants)
  amrex::Interpolater *interp = nullptr;
  switch (interp_type) {
  case 0: // amrex_interp_pc (piecewise constant)
    interp = &amrex::pc_interp;
    break;
  case 1: // amrex_interp_cell_cons (cell conservative)
  default:
    interp = &amrex::cell_cons_interp;
    break;
  }

  amrex::Vector<amrex::MultiFab *> cmf = {mf_old_c, mf_new_c};
  amrex::Vector<amrex::Real> ctime = {time_old_c, time_new_c};
  amrex::Vector<amrex::MultiFab *> fmf = {mf_old_f, mf_new_f};
  amrex::Vector<amrex::Real> ftime = {time_old_f, time_new_f};

  // Build BCRec from lo_bc/hi_bc arrays
  amrex::Vector<amrex::BCRec> bcs(ncomp);
  for (int i = 0; i < ncomp; ++i) {
    bcs[i].set(amrex::Orientation(0, amrex::Orientation::low),
               lo_bc[i * 3 + 0]);
    bcs[i].set(amrex::Orientation(0, amrex::Orientation::high),
               hi_bc[i * 3 + 0]);
    bcs[i].set(amrex::Orientation(1, amrex::Orientation::low),
               lo_bc[i * 3 + 1]);
    bcs[i].set(amrex::Orientation(1, amrex::Orientation::high),
               hi_bc[i * 3 + 1]);
    bcs[i].set(amrex::Orientation(2, amrex::Orientation::low),
               lo_bc[i * 3 + 2]);
    bcs[i].set(amrex::Orientation(2, amrex::Orientation::high),
               hi_bc[i * 3 + 2]);
  }

  nga2::NGA2BCFunctor bc_functor_c(solver_ctx, bc_dispatch, geom_c);
  nga2::NGA2BCFunctor bc_functor_f(solver_ctx, bc_dispatch, geom_f);

  amrex::IntVect ratio(AMREX_D_DECL(ref_ratio, ref_ratio, ref_ratio));

  // Convert from 1-indexed (Fortran) to 0-indexed (C++)
  amrex::FillPatchTwoLevels(*mf, time, cmf, ctime, fmf, ftime, scomp - 1,
                            dcomp - 1, ncomp, *geom_c, *geom_f, bc_functor_c, 0,
                            bc_functor_f, 0, ratio, interp, bcs, 0);
}

//-----------------------------------------------------------------------------
// FillCoarsePatch - fill fine level from coarse only (no fine data needed)
// Note: scomp/dcomp are 1-indexed (Fortran convention), converted to 0-indexed
// here
//-----------------------------------------------------------------------------
void amrmfab_fillcoarsepatch(void *mf_f_ptr, double time, void *mf_c_ptr,
                             void *geom_c_ptr, void *geom_f_ptr,
                             void *solver_ctx,
                             nga2::FillPatchBCDispatcher bc_dispatch, int scomp,
                             int dcomp, int ncomp, int ref_ratio,
                             int interp_type, int *lo_bc, int *hi_bc, int nbc) {
  auto *mf_f = static_cast<amrex::MultiFab *>(mf_f_ptr);
  auto *mf_c = static_cast<amrex::MultiFab *>(mf_c_ptr);
  auto *geom_c = static_cast<amrex::Geometry *>(geom_c_ptr);
  auto *geom_f = static_cast<amrex::Geometry *>(geom_f_ptr);

  // Look up interpolator
  amrex::Interpolater *interp = nullptr;
  switch (interp_type) {
  case 0:
    interp = &amrex::pc_interp;
    break;
  case 1:
  default:
    interp = &amrex::cell_cons_interp;
    break;
  }

  // Build BCRec
  amrex::Vector<amrex::BCRec> bcs(ncomp);
  for (int i = 0; i < ncomp; ++i) {
    bcs[i].set(amrex::Orientation(0, amrex::Orientation::low),
               lo_bc[i * 3 + 0]);
    bcs[i].set(amrex::Orientation(0, amrex::Orientation::high),
               hi_bc[i * 3 + 0]);
    bcs[i].set(amrex::Orientation(1, amrex::Orientation::low),
               lo_bc[i * 3 + 1]);
    bcs[i].set(amrex::Orientation(1, amrex::Orientation::high),
               hi_bc[i * 3 + 1]);
    bcs[i].set(amrex::Orientation(2, amrex::Orientation::low),
               lo_bc[i * 3 + 2]);
    bcs[i].set(amrex::Orientation(2, amrex::Orientation::high),
               hi_bc[i * 3 + 2]);
  }

  nga2::NGA2BCFunctor bc_functor_c(solver_ctx, bc_dispatch, geom_c);
  nga2::NGA2BCFunctor bc_functor_f(solver_ctx, bc_dispatch, geom_f);

  amrex::IntVect ratio(AMREX_D_DECL(ref_ratio, ref_ratio, ref_ratio));

  // Convert from 1-indexed (Fortran) to 0-indexed (C++)
  // Use InterpFromCoarseLevel which doesn't require fine-level source data
  amrex::InterpFromCoarseLevel(*mf_f, time, *mf_c, scomp - 1, dcomp - 1, ncomp,
                               *geom_c, *geom_f, bc_functor_c, 0, bc_functor_f,
                               0, ratio, interp, bcs, 0);
}

//-----------------------------------------------------------------------------
// Plotfile Output (Native and HDF5)
//-----------------------------------------------------------------------------

void amrplotfile_write_native(const char *name, int nlevels, void **mf_ptrs,
                              const char **varnames, int ncomp,
                              void **geom_ptrs, double time, int *level_steps,
                              int *ref_ratios) {
  // Convert void* arrays to AMReX types
  amrex::Vector<const amrex::MultiFab *> mf_vec(nlevels);
  amrex::Vector<amrex::Geometry> geom_vec(nlevels);
  amrex::Vector<int> steps_vec(nlevels);
  amrex::Vector<amrex::IntVect> rr_vec(nlevels - 1);
  amrex::Vector<std::string> varname_vec(ncomp);

  for (int lev = 0; lev < nlevels; ++lev) {
    mf_vec[lev] = static_cast<amrex::MultiFab *>(mf_ptrs[lev]);
    geom_vec[lev] = *static_cast<amrex::Geometry *>(geom_ptrs[lev]);
    steps_vec[lev] = level_steps[lev];
  }
  for (int lev = 0; lev < nlevels - 1; ++lev) {
    rr_vec[lev] = amrex::IntVect(
        AMREX_D_DECL(ref_ratios[lev], ref_ratios[lev], ref_ratios[lev]));
  }
  for (int i = 0; i < ncomp; ++i) {
    varname_vec[i] = std::string(varnames[i]);
  }

  amrex::WriteMultiLevelPlotfile(std::string(name), nlevels, mf_vec,
                                 varname_vec, geom_vec, time, steps_vec,
                                 rr_vec);
}

#ifdef AMREX_USE_HDF5
void amrplotfile_write_hdf5(const char *name, int nlevels, void **mf_ptrs,
                            const char **varnames, int ncomp, void **geom_ptrs,
                            double time, int *level_steps, int *ref_ratios,
                            const char *compression) {
  // Convert void* arrays to AMReX types
  amrex::Vector<const amrex::MultiFab *> mf_vec(nlevels);
  amrex::Vector<amrex::Geometry> geom_vec(nlevels);
  amrex::Vector<int> steps_vec(nlevels);
  amrex::Vector<amrex::IntVect> rr_vec(nlevels - 1);
  amrex::Vector<std::string> varname_vec(ncomp);

  for (int lev = 0; lev < nlevels; ++lev) {
    mf_vec[lev] = static_cast<amrex::MultiFab *>(mf_ptrs[lev]);
    geom_vec[lev] = *static_cast<amrex::Geometry *>(geom_ptrs[lev]);
    steps_vec[lev] = level_steps[lev];
  }
  for (int lev = 0; lev < nlevels - 1; ++lev) {
    rr_vec[lev] = amrex::IntVect(
        AMREX_D_DECL(ref_ratios[lev], ref_ratios[lev], ref_ratios[lev]));
  }
  for (int i = 0; i < ncomp; ++i) {
    varname_vec[i] = std::string(varnames[i]);
  }

  std::string comp_str = compression ? std::string(compression) : "";
  amrex::WriteMultiLevelPlotfileHDF5(std::string(name), nlevels, mf_vec,
                                     varname_vec, geom_vec, time, steps_vec,
                                     rr_vec, comp_str);
}
#endif

//-----------------------------------------------------------------------------
// Checkpoint I/O (VisMF)
//-----------------------------------------------------------------------------
#include <AMReX_Utility.H>
#include <AMReX_VisMF.H>

// Write a MultiFab to disk using VisMF
void amrmfab_vismf_write(void *mf, const char *path) {
  amrex::VisMF::Write(*static_cast<amrex::MultiFab *>(mf), std::string(path));
}

// Read a MultiFab from disk using VisMF
// NOTE: The MultiFab must already be defined with correct BoxArray/DistroMap
void amrmfab_vismf_read(void *mf, const char *path) {
  amrex::VisMF::Read(*static_cast<amrex::MultiFab *>(mf), std::string(path));
}

// Create directory hierarchy for checkpoints
void amrcheckpoint_prebuild_dirs(const char *dirname, const char *subdir_prefix,
                                 int nlevels) {
  amrex::PreBuildDirectorHierarchy(std::string(dirname),
                                   std::string(subdir_prefix), nlevels, true);
}

// Get MultiFab file prefix for a level (e.g., "chk00100/Level_0/phi")
// Returns the full path that can be used with VisMF::Write/Read
void amrcheckpoint_mfab_prefix(char *result, int result_len, int lev,
                               const char *dirname, const char *level_prefix,
                               const char *mfab_name) {
  std::string prefix = amrex::MultiFabFileFullPrefix(lev, std::string(dirname),
                                                     std::string(level_prefix),
                                                     std::string(mfab_name));
  // Copy to result buffer
  size_t len = std::min(prefix.size(), static_cast<size_t>(result_len - 1));
  std::copy(prefix.begin(), prefix.begin() + len, result);
  result[len] = '\0';
}

// Set number of output files for VisMF (controls I/O aggregation)
// nfiles=1 gives single file output; nfiles=nranks gives file-per-rank
void amrvismf_set_noutfiles(int nfiles) { amrex::VisMF::SetNOutFiles(nfiles); }

// Get current number of output files setting
int amrvismf_get_noutfiles() { return amrex::VisMF::GetNOutFiles(); }

//=============================================================================
// HDF5 Plotfile Utilities
//=============================================================================
#ifdef AMREX_USE_HDF5
#include <hdf5.h>

// Read the time attribute from an AMReX HDF5 plotfile
// Returns the time value, or -1.0 if file doesn't exist or can't be read
double amrplotfile_read_time(const char *filename) {
  double time = -1.0;

  // Check if file exists
  hid_t file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file_id < 0) {
    return -1.0; // File doesn't exist or can't be opened
  }

  // Read the "time" attribute from root group
  if (H5Aexists(file_id, "time") > 0) {
    hid_t attr_id = H5Aopen(file_id, "time", H5P_DEFAULT);
    if (attr_id >= 0) {
      // Read as double (time is stored as array of size 1)
      double time_arr[1];
      H5Aread(attr_id, H5T_NATIVE_DOUBLE, time_arr);
      time = time_arr[0];
      H5Aclose(attr_id);
    }
  }

  H5Fclose(file_id);
  return time;
}

#else
// Stub when HDF5 not available
double amrplotfile_read_time(const char *filename) { return -1.0; }
#endif

} // extern "C"
