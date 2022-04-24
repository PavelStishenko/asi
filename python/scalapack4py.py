# dl_iterate_phdr
# LoadLibrary
# 

from ctypes import cdll, CDLL, RTLD_LOCAL, Structure
from ctypes import POINTER, byref, c_int, c_int64, c_char, c_bool, c_char_p, c_double, c_void_p, CFUNCTYPE, py_object, cast, byref
import numpy as np
from numpy.ctypeslib import ndpointer

from mpiprint import *
#libdl = CDLL("libdl.so")
#print(libdl.dl_iterate_phdr

#libdl = CDLL("libdl.so")
#blacs = CDLL("/home/mazay/prg/aims/build-so-3/libaims.220309.scalapack.mpi.so")
#default_lib = None


def nullable_ndpointer(*args, **kwargs):
  # solution from here: https://stackoverflow.com/a/37664693/3213940
  # TODO look for this isssue: https://github.com/numpy/numpy/issues/6239
  base = ndpointer(*args, **kwargs)
  def from_param(cls, obj):
    if obj is None:
      return obj
    return base.from_param(obj)
  return type(base.__name__, (base,), {'from_param': classmethod(from_param)})

DLEN_ = 9
class blacs_desc(Structure): # https://www.netlib.org/scalapack/slug/node68.html
  _fields_ = [("dtype", c_int),
              ("ctxt", c_int),
              ("m", c_int),
              ("n", c_int),
              ("mb", c_int),
              ("nb", c_int),
              ("rsrc", c_int),
              ("csrc", c_int),
              ("lld", c_int)]
  def __init__(self, lib, blacs_ctx=-1, m=0, n=0, mb=1, nb=1, rsrc=0, csrc=0, lld=None, buf=None):
    super().__init__()
    self.lib = lib
    
    if buf is not None:
      if not hasattr(buf, '__get_index__'):
        buf = np.ctypeslib.as_array(buf, shape=(DLEN_,))
      self.dtype, self.ctxt, self.m, self.n, self.mb, self.nb, self.rsrc, self.csrc, self.lld = (*buf,)
      return # wrapped copy

    if blacs_ctx == -1: # stub for non-participating process
      self.dtype=1
      self.ctxt=blacs_ctx
      self.m = 0
      self.n = 0
      self.mb = 0
      self.nb = 0
      self.rsrc = 0
      self.csrc = 0
      self.lld = 0
      return # default, zero-intitalized

    assert m > 0, m
    assert n > 0, n
    if lld is None: # calc default lld
      nprow, _, myrow, _ = self.lib.blacs_gridinfo(blacs_ctx)
      lld = self.lib.numroc(m, mb, myrow, rsrc, nprow)
      lld = max(1, lld) # ldd is not less than 1

    self.lib.descinit(self, m, n, mb, nb, rsrc, csrc, blacs_ctx, lld)
    
  @property
  def myrow(self):
    _, _, myrow, _ = self.lib.blacs_gridinfo(self.ctxt)
    return myrow

  @property
  def mycol(self):
    _, _, _, mycol = self.lib.blacs_gridinfo(self.ctxt)
    return mycol

  @property
  def nprow(self):
    nprow, _, _, _ = self.lib.blacs_gridinfo(self.ctxt)
    return nprow

  @property
  def npcol(self):
    _, npcol, _, _ = self.lib.blacs_gridinfo(self.ctxt)
    return npcol

  @property
  def is_distributed(self):
    nprow, npcol, _, _ = self.lib.blacs_gridinfo(self.ctxt)
    return (nprow * npcol) > 1

  @property
  def locrow(self):
    nprow, _, myrow, _ = self.lib.blacs_gridinfo(self.ctxt)
    return self.lib.numroc(self.m, self.mb, myrow, self.rsrc, nprow)

  @property
  def loccol(self):
    _, npcol, _, mycol = self.lib.blacs_gridinfo(self.ctxt)
    return self.lib.numroc(self.n, self.nb, mycol, self.csrc, npcol)

  def alloc_zeros(self, dtype):
    res = self.alloc(dtype)
    res[:,:] = 0
    return res
  
  def alloc(self, dtype):
    return np.ndarray(shape=(self.locrow, self.loccol), order='F', dtype=dtype)

  def __repr__(self):
    return f"{self.dtype} {self.ctxt} {self.m} {self.n} {self.mb} {self.nb} {self.rsrc} {self.csrc} {self.lld}"


class ScaLAPACK4py:
  def __init__(self, blacs):
    """
      blacs: Library loaded via CTypes with BLACS and ScaLAPACK functions
    """
    self.blacs = blacs
    blacs.blacs_get_.restype = None
    blacs.blacs_get_.argtypes = [POINTER(c_int), POINTER(c_int), POINTER(c_int)]
    blacs.blacs_gridinit_.restype = None
    blacs.blacs_gridinit_.argtypes = [POINTER(c_int), POINTER(c_char), POINTER(c_int), POINTER(c_int)]
    blacs.blacs_gridexit_.restype = None
    blacs.blacs_gridexit_.argtypes = [POINTER(c_int)]
    blacs.blacs_gridinfo_.restype = None
    blacs.blacs_gridinfo_.argtypes = [POINTER(c_int), POINTER(c_int), POINTER(c_int), POINTER(c_int), POINTER(c_int)]
    blacs.numroc_.restype = c_int
    blacs.numroc_.argtypes = [POINTER(c_int), POINTER(c_int), POINTER(c_int), POINTER(c_int), POINTER(c_int)]

    blacs.descinit_.restype = None
    #descinit_:   DESC, M, N, MB, NB, IRSRC, ICSRC, ICTXT, LLD, INFO
    blacs.descinit_.argtypes = [POINTER(blacs_desc), POINTER(c_int), POINTER(c_int), POINTER(c_int), POINTER(c_int), POINTER(c_int), POINTER(c_int), POINTER(c_int), POINTER(c_int), POINTER(c_int)] 

    #EXTERN void pdgemr2d_(int *m , int *n , double *a , int *ia , int *ja , int *desca , double *b , int *ib , int *jb , int *descb , int *ictxt);
    blacs.pdgemr2d_.restype = None
    blacs.pdgemr2d_.argtypes = [POINTER(c_int), POINTER(c_int), 
                                nullable_ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'), POINTER(c_int), POINTER(c_int), POINTER(blacs_desc), 
                                nullable_ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'), POINTER(c_int), POINTER(c_int), POINTER(blacs_desc), 
                                POINTER(c_int)]
    #EXTERN void pzgemr2d_(int *m , int *n , double_complex_t *a , int *ia , int *ja , int *desca , double_complex_t *b , int *ib , int *jb , int *descb , int *ictxt);
  
  def descinit(self, desc, m, n, mb, nb, rsrc, csrc, blacs_ctx, lld):
    m = c_int(m)
    n = c_int(n)
    mb = c_int(mb)
    nb = c_int(nb)
    irsrc = c_int(rsrc)
    icsrc = c_int(csrc)
    blacs_ctx = c_int(blacs_ctx)
    lld = c_int(lld)
    info = c_int()
    self.blacs.descinit_(desc, m, n, mb, nb, irsrc, icsrc, blacs_ctx, lld, info)
    assert info.value == 0, info.value
  
  def pdgemr2d(self, m, n, a, ia, ja, desca, b, ib, jb, descb, blacs_ctx):
    m = c_int(m)
    n = c_int(n)
    ia = c_int(ia)
    ja = c_int(ja)
    ib = c_int(ib)
    jb = c_int(jb)
    blacs_ctx = c_int(blacs_ctx)
    self.blacs.pdgemr2d_(m, n, a, ia, ja, desca, b, ib, jb, descb, blacs_ctx)


  def get_default_system_context(self):
    zero = c_int(0)
    system_context = c_int()
    self.blacs.blacs_get_(zero, zero, system_context)
    return system_context.value

  def get_system_context(self, blacs_ctx): # get system context of the blacs_context
    blacs_ctx = c_int(blacs_ctx)
    what = c_int(10)
    val = c_int()
    self.blacs.blacs_get_(blacs_ctx, what, val)
    return val.value

  def make_blacs_context(self, sys_context, MP, NP):
    order = c_char(b'R')
    context_inout = c_int(sys_context)
    c_MP = c_int(MP)
    c_NP = c_int(NP)
    self.blacs.blacs_gridinit_(context_inout, byref(order), c_MP, c_NP)
    return context_inout.value
  
  def close_blacs_context(self, blacs_ctx):
    if blacs_ctx == -1:
      return
    blacs_ctx = c_int(blacs_ctx)
    self.blacs.blacs_gridexit_(blacs_ctx)

  def make_blacs_desc(self, blacs_ctx, m, n, mb=1, nb=1, rsrc=0, csrc=0, lld=None):
    return blacs_desc(self, blacs_ctx, m, n, mb, nb, rsrc, csrc, lld)

  def wrap_blacs_desc(self, buf):
    return blacs_desc(self, buf=buf)

  def blacs_gridinfo(self, blacs_ctx):
    c_blacs_ctx = c_int(blacs_ctx)
    nprow = c_int()
    npcol = c_int()
    myrow = c_int()
    mycol = c_int()
    self.blacs.blacs_gridinfo_(c_blacs_ctx, nprow, npcol, myrow, mycol)
    return nprow.value, npcol.value, myrow.value, mycol.value

  def numroc(self, n, nb, iproc, isrcproc, nprocs):
    n = c_int(n)
    nb = c_int(nb)
    iproc = c_int(iproc)
    isrcproc = c_int(isrcproc)
    nprocs = c_int(nprocs)
    return self.blacs.numroc_(n, nb, iproc, isrcproc, nprocs)

  def scatter(self, src_data, dest_desc, dest_data):
    m = dest_desc.m
    n = dest_desc.m
    common_blacs_ctx = dest_desc.ctxt
    src_blacs_ctx = self.make_blacs_context(self.get_system_context(common_blacs_ctx), 1, 1)
    if src_blacs_ctx != -1:
      assert m == src_data.shape[0]
      assert n == src_data.shape[1]
    else:
      assert src_data is None, f"src_data must be None if src_blacs_ctx == -1"
    src_desc = self.make_blacs_desc(src_blacs_ctx, m, n)
    self.pdgemr2d(m, n, src_data, 1, 1, src_desc, dest_data, 1, 1, dest_desc, common_blacs_ctx)
    self.close_blacs_context(src_blacs_ctx)

  def gather(self, src_desc, src_data, dest_data):
    m = src_desc.m
    n = src_desc.m
    common_blacs_ctx = src_desc.ctxt
    gatherer_blacs_ctx = self.make_blacs_context(self.get_system_context(common_blacs_ctx), 1, 1)
    if gatherer_blacs_ctx != -1:
      assert m == dest_data.shape[0]
      assert n == dest_data.shape[1]
    else:
      assert dest_data is None, dest_data
    dest_desc = self.make_blacs_desc(gatherer_blacs_ctx, m, n)
    self.pdgemr2d(m, n, src_data, 1, 1, src_desc, dest_data, 1, 1, dest_desc, common_blacs_ctx)
    self.close_blacs_context(gatherer_blacs_ctx)

