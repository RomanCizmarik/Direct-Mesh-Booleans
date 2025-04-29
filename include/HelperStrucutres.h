//TODO: License

#pragma once

#include <bitset>
#include <Eigen/SparseCore>
#include <OpenMesh/Core/Mesh/PolyConnectivity.hh>

namespace DMB
{
    using tSparseIntMatrix = Eigen::SparseMatrix<int>;
    using tTriplet = Eigen::Triplet<int>;

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    //! \brief   Encapsulation of input data for FaRMA algorithms.
    //! The library uses vectors for everything. Consult the library API for more info on these.
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    struct TriangleSoup
    {
        //! Triplets (x,y,z) of world-space coordinates
        std::vector<double> coordinates;
        //! Triplets of indices to coordinates
        std::vector<uint> triangles;
        //! Per-triangle model label (model index)
        std::vector<uint> labels;
    };

    struct InputTriangleMesh
    {
        //! Triplets (x,y,z) of world-space coordinates
        std::vector<double> coordinates;
        //! Triplets of indices to coordinates
        std::vector<uint> triangles;
        //! Per-triangle model label (model index)
        uint label;
    };

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    //! \brief   Structure for measuring time of the operations. For debugging only.
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    struct MeshBooleansTimings
    {
        float brewingSoup;                 //! How long it takes to add meshes to the soup.
        float resolvingIntersections;      //! How long the library call is.
        float performingBooleanOperation;  //! How long it takes to figure out what is inside and what is outside.
    };


    struct SIntersectionCurve
    {
        std::vector<OpenMesh::SmartHalfedgeHandle> loopHalfedges;
        bool closed = false;
        bool closedWithBoundary = false; // the curve is considered to be closed, but it is closed using boundary edges, which are not part of intersection curve
        bool closedWithCoplanarRegions = false; // the curve is considered to be closed, but it is closed using coplanar regions
    };
}
