from mpi4py import MPI
import io
def print_to_string(*args, **kwargs):
  # inspired by https://stackoverflow.com/a/39823534
  with io.StringIO() as output:
    print(*args, file=output, **kwargs)
    return output.getvalue()

def parprint(*args, **kwargs):
  if MPI.COMM_WORLD.Get_rank() == 0:
    print (*args, **kwargs)

def ordprint(*args):
  """
    Print from all processes in rank order
  """
  s = print_to_string(*args)
  strings = MPI.COMM_WORLD.gather(s)
  if MPI.COMM_WORLD.Get_rank() == 0:
    for s in strings:
      print(s, end="")
