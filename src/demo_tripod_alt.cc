#include "getfem/getfem_model_solvers.h"
#include "getfem/getfem_import.h"
#include "getfem/getfem_export.h"
#include "gmm/gmm.h"
#include "getfem/getfem_mesher.h"
#include "getfem/getfem_generic_assembly.h"

using bgeot::size_type;
using bgeot::base_node;
using bgeot::base_small_vector;
typedef getfem::model_real_plain_vector plain_vector;

int main(void) {
	char *file_msh = "tripod.GiD.msh";
	int degree = 2;
	double E = 1e3;
	double Nu = 0.3;
	double Lambda = E*Nu/((1+Nu)*(1-2*Nu));
	double Mu = E/(2*(1+Nu));
	getfem::mesh m;
	getfem::import_mesh("gid:../mesh/tripod.GiD.msh",m);
	m.optimize_structure();
	getfem::vtk_export exp("m.vtk", false);
	exp.exporting(m);
	exp.write_mesh();
	getfem::mesh_fem mfu(m,3);
	getfem::mesh_fem mfd(m,1);
	mfu.set_finite_element(getfem::PK_fem(getfem::dim_type(3), getfem::short_type(1)));
	mfd.set_finite_element(getfem::PK_fem(getfem::dim_type(3), getfem::short_type(0)));
	getfem::mesh_im mim;
	mim.set_integration_method(m.convex_index(getfem::dim_type(5)), getfem::int_method_descriptor("IM_TETRAHEDRON(5)"));
	int P = m.nb_points();
	return 0;
}

