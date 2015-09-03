
// Initialization

#include <iostream>
#include "getfem/getfem_model_solvers.h"
#include "getfem/getfem_export.h"
#include "gmm/gmm.h"
#include "getfem/getfem_mesher.h"
#include "getfem/getfem_generic_assembly.h"

using namespace std;
using bgeot::size_type;
using bgeot::base_node;
using bgeot::base_small_vector;
typedef getfem::model_real_plain_vector plain_vector;

int main(void) {

// Parameters of the model

  double epsilon = 1.; // Thickness of the plate (cm)
  double E = 21E6;     // Young Modulus (N/cm^2)
  double nu = 0.3;     // Poisson ratio
  double clambda = E*nu/((1+nu)*(1-2*nu));
  double cmu = E/(2*(1+nu));
  double clambdastar = 2*clambda*cmu/(clambda+2*cmu);
  double F = 100E2;    // Force density at the right boundary (N/cm^2)
  double kappa = 4.;   // Thermal conductivity (W/(cm K))
  double D = 10.;      // Heat transfert coefficient (W/(K cm^2))
  double air_temp = 20;// Temperature of the air in oC.
  double alpha_th = 16.6E-6; // Thermal expansion coefficient (/K)
  double T0 = 20.;     // Reference temperature in oC
  double rho_0 = 1.754E-8; // Resistance temperature coefficient at T0
  double alpha = 0.0039; // Second resistance temperature coefficient

  double h = 2;        // Approximate mesh size
  bgeot::dim_type elements_degree = 2; // Degree of the finite element methods

// Mesh generation

  getfem::mesh mesh;
  getfem::mesher_rectangle mo1(base_node(0., 0.), base_node(100., 25.));
  getfem::mesher_ball mo2(base_node(25., 12.5), 8.);
  getfem::mesher_ball mo3(base_node(50., 12.5), 8.);
  getfem::mesher_ball mo4(base_node(75., 12.5), 8.);
  getfem::mesher_union mo5(mo2, mo3, mo4);
  getfem::mesher_setminus mo(mo1, mo5);

  std::vector<getfem::base_node> fixed;
  getfem::build_mesh(mesh, mo, h, fixed, 2, -2);

// Mesh draw

  getfem::vtk_export exp("mesh.vtk", false);
  exp.exporting(mesh);
  exp.write_mesh();
  // You can view the mesh for instance with
  // mayavi2 -d mesh.vtk -f ExtractEdges -m Surface

// Definition of finite element methods and integration method

  getfem::mesh_fem mfu(mesh, 2);
  mfu.set_classical_finite_element(elements_degree);
  getfem::mesh_fem mft(mesh, 1);
  mft.set_classical_finite_element(elements_degree);
  getfem::mesh_fem mfvm(mesh, 1);
  mfvm.set_classical_discontinuous_finite_element(elements_degree);

  getfem::mesh_im  mim(mesh);
  mim.set_integration_method(bgeot::dim_type(gmm::sqr(elements_degree)));

}
