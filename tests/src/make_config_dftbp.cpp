#include <fstream>
#include <boost/format.hpp>

using boost::format;


extern const double test_E = -3.943259;
extern const double test_forces[] = {-0.001499, -0.000000, 0.168818, 0.002580, -0.021285, -0.086934, 0.002580, 0.021285, -0.086934};
extern const double test_charges[] = {-0.106997, 0.053498, 0.053498};
extern const double test_calc_esp[] = {0.000461 , 0.029803}; // TODO fill out


const char *dftb_in_hsd = R"--(
Geometry = GenFormat {
3 %1%
  O H

  1 1  1.00000000000E+00 -0.10000000000E+01  0.00000000000E+00
  2 2  1.00000000000E+00  0.00000000000E+00  0.78306400000E+00
  3 2  1.00000000000E+00  0.00000000000E+00 -0.78306400000E+00
%2%
}

Hamiltonian = DFTB {
  Scc = Yes
  SlaterKosterFiles {
    O-O = "./O-O.skf"
    O-H = "./O-H.skf"
    H-O = "./H-O.skf"
    H-H = "./H-H.skf"
  }
  MaxAngularMomentum {
    O = "p"
    H = "s"
  }

  Solver = DivideAndConquer {}

  SpinPolarisation = Colinear {
    RelaxTotalSpin = Yes
  }
  
  SpinConstants = {
    O={
    #Wpp
    -0.028
    }
    H={
    # Wss
    -0.072
    }
  }

%3%
}

Options {
  #WriteHS = Yes
  #WriteRealHS = Yes
}

Analysis {
  CalculateForces = Yes
  WriteBandOut = Yes
  ElectrostaticPotential {
    Points {
      10 10 10
      5  5  5
      1.0 2.0 0.1
      1.0 2.0 0.0
    }
  }
}

Parallel {
Groups = 2
Blacs {
BlockSize = 4
}
}

ParserOptions {
  ParserVersion = 9
}
)--";

/*
Add in Hamiltonian section for testing:
	ElectricField = {
		PointCharges = {
			CoordsAndCharges [Bohr] = {
			 1.0 2.0 0.0  0.1
			}
		}
	}
*/

void make_config_files(bool pbc)
{
  std::ofstream in_hsd("dftb_in.hsd");
  const char *c_s, *cell, *kgrid;
  if (pbc)
  {
    c_s = "S";
    cell = R"--(0 0 0 
30 0 0 
0 30 0 
0 0 30)--";
    kgrid = R"--(  KPointsAndWeights = SupercellFolding { 
  1 0 0 
  0 1 0 
  0 0 2 
  0.5 0.5 0.0 
  })--";
  }
  else
  {
    c_s = "C";
    cell = "";
    kgrid = "";
  }
  
  auto s = format(dftb_in_hsd) % c_s % cell % kgrid;
  in_hsd << s;
}
