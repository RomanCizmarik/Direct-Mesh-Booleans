

#ifdef _MSC_VER // Workaround for known bugs and issues on MSVC
#define _HAS_STD_BYTE 0  // https://developercommunity.visualstudio.com/t/error-c2872-byte-ambiguous-symbol/93889
#define NOMINMAX // https://stackoverflow.com/questions/1825904/error-c2589-on-stdnumeric-limitsdoublemin
#endif

#include <iostream>

#include "SimpleMesh.h"

#include "Booleans.h"

int main(int argc, char** argv)
{
	//std::cout << "Hello world" << std::endl;

    //std::vector<double> values = { 3.4, 5.6, 1.2, 7.8, 2.0, 4.5 };
    //double X = 4.0;

    //// Partition: elements <= X will go to the front, > X will go to the back
    //auto it = std::partition(values.begin(), values.end(),
    //    [X](double val) { return val <= X; });

    //// Print partitioned vector
    //std::cout << "Partitioned vector: ";
    //for (double val : values)
    //{
    //    std::cout << val << ", ";
    //}
    //std::cout << "\n";

    //// Print split index
    //std::cout << "Split index: " << std::distance(values.begin(), it) << "\n";

    //std::vector<double> left(values.begin(), it);
    //std::vector<double> right(it, values.end());

    //std::cout << "left: ";
    //for (auto v : left)
    //{
    //    std::cout << v << ", ";
    //}
    //std::cout<<std::endl;

    //std::cout << "right: ";
    //for (auto v : right)
    //{
    //    std::cout << v << ", ";
    //}
    //std::cout << std::endl;


	//DMB::SimpleMesh m;
	////OpenMesh::IO::read_mesh(m, "C:/skola/PhD/Samples/booleans/result19.obj");
	//OpenMesh::IO::read_mesh(m, "C:/skola/PhD/Samples/booleans/sphere_9x6.obj");

	//DMB::AABBTree< DMB::SimpleMesh > tree;
	//tree.build(m);

 //   tree.debug_outputAsMesh();
    DMB::SimpleMesh m;
    DMB::meshUnion(m, "C:/skola/PhD/VUT/booleans_paper/evaluation/Validity/input_models/101954.stl", "C:/skola/PhD/VUT/booleans_paper/evaluation/Validity/input_models/894835.stl");
    
    //std::string filename("C:/skola/PhD/VUT/booleans_paper/evaluation/Validity/input_models/101954.stl");

    //std::vector<double> in_coords, out_coords;
    //std::vector<uint> in_tris, out_tris;
    //std::vector<genericPoint*> gen_points;
    //point_arena arena;

    //load(filename, in_coords, in_tris);

    ///*-------------------------------------------------------------------
    // * There are 4 versions of the solveIntersections function. Please
    // * refer to the solve_intersections.h file to see how to use them. */

    //solveIntersections(in_coords, in_tris, arena, gen_points, out_tris);

    //DMB::TriangleSoup soup{};
    //DMB::MeshArrangement<MeshType> ma{};
    //solveIntersections(soup.coordinates, soup.triangles, soup.labels, arena, ma.m_coordinatesImplicit, ma.m_triangles, ma.m_labels/*, true*/);

    //        ! Triplets (x,y,z) of world-space coordinates
    //std::vector<double> coordinates;
    //! Triplets of indices to coordinates
    //std::vector<uint> triangles;
    //! Per-triangle model label (model index)
    //std::vector<uint> labels;

    //std::vector<double> out_coordinates;
    //std::vector<genericPoint*> out_coordinatesImplicit;
    //std::vector<uint> out_triangles;
    //std::vector<std::bitset<NBIT>> out_labels;

    //solveIntersections(coordinates, triangles, labels, arena, out_coordinatesImplicit, out_triangles, out_labels);
}