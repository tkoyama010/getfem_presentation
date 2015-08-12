
// Initialization

#include "getfem/getfem_model_solvers.h"
#include "getfem/getfem_export.h"
#include "gmm/gmm.h"
#include "getfem/getfem_mesher.h"
#include "getfem/getfem_generic_assembly.h"

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

// Boundary selection

  getfem::mesh_region border_faces;
  getfem::outer_faces_of_mesh(mesh, border_faces);
  getfem::mesh_region fb1
    = getfem::select_faces_in_box(mesh, border_faces, base_node(1., 1.),
                   base_node(99., 24.));
  getfem::mesh_region fb2
    = getfem::select_faces_of_normal(mesh, border_faces,
                       base_small_vector( 1., 0.), 0.01);
  getfem::mesh_region fb3
    = getfem::select_faces_of_normal(mesh, border_faces,
                       base_small_vector(-1., 0.), 0.01);
  getfem::mesh_region fb4
    = getfem::select_faces_of_normal(mesh, border_faces,
                       base_small_vector(0.,  1.), 0.01);
  getfem::mesh_region fb5
    = getfem::select_faces_of_normal(mesh, border_faces,
                       base_small_vector(0., -1.), 0.01);

  size_type RIGHT_BOUND=1, LEFT_BOUND=2, TOP_BOUND=3, BOTTOM_BOUND=4;
  mesh.region( RIGHT_BOUND) = getfem::mesh_region::subtract(fb2, fb1);
  mesh.region(  LEFT_BOUND) = getfem::mesh_region::subtract(fb3, fb1);
  mesh.region(   TOP_BOUND) = getfem::mesh_region::subtract(fb4, fb1);
  mesh.region(BOTTOM_BOUND) = getfem::mesh_region::subtract(fb5, fb1);

// Mesh draw

  getfem::vtk_export exp("mesh.vtk", false);
  exp.exporting(mesh);
  exp.write_mesh();
  // You can view the mesh for instance with
  // mayavi2 -d mesh.vtk -f ExtractEdges -m Surface
}
