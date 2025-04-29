//TODO: License



#pragma once

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#include <array>

#pragma warning (push, 0)
#pragma warning(disable: 2220)
#pragma warning(disable: 4701)
#pragma warning(disable: 4267)
#pragma warning(disable: 4715)
#pragma warning(disable: 4172)
#include <solve_intersections.h>
#include <common.h> //for NBIT define
#include <io_functions.h>
#pragma warning (pop)

#include "HelperStrucutres.h"
#include "HelperFunctions.h"
#include "MeshArrangement.h"
#include "MatrixMesh.h"

//#include "common.h"
#include "io_functions.h"

namespace DMB
{
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    //! \brief   Performs generic two mesh boolean operation. 
    //! The predicate is used to indicate which triangles shall be part of the output mesh and in which orientation.
    //! Utilizes paralelization with omp.
    //!
    //! \param   output      Reference to output mesh. It will be cleaned at the beginning of the operations!
    //! \param   soup        Input triangle soup.
    //! \param   predicate   Predicate in form of "int predicate(const std::array<typename MeshType::Point, 3> &triangle, int label)"
    //!                      which for an output triangle and a its label determines if it should be discarded (0), 
    //!                      part of the output mesh (>0) or part of the output mesh but flipped (<0).
    //! \param   progress    Progress callback functor, if this functor returns false, the computation will be interrupted.
    //! \param   timings     Optional MeshBooleansTimings for measuring time of the operations.
    //! \return              SUCCESS - everything ok, FAILED - operation failed, INTERRUPTED - operation was interrupted.
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    template <typename MeshType, typename F>
    bool meshBoolean(MeshType& output,
        const TriangleSoup& soup,
        const F& predicate, const std::function<bool()>& progress, MeshBooleansTimings* timings = nullptr)
    {
        std::chrono::steady_clock::time_point point;
        using tTriangle = std::array<typename MeshType::Point, 3>;
        using tTriangleVertexHandles = std::array<typename MeshType::VertexHandle, 3>;
        using tTriangleIndices = std::array<uint, 3>;
        using tPoint = typename MeshType::Point;
        using tVertexHandle = typename MeshType::VertexHandle;
        using tFaceHandle = typename MeshType::FaceHandle;
        using tEdgeHandle = typename MeshType::EdgeHandle;

        point = std::chrono::steady_clock::now();

        MeshArrangement<MeshType> ma{};

        point_arena arena;

        //try to solve mesh arrangement, it can throw if it fails internally
        try
        {

            //TODO: make this part of MeshArrangement class
            solveIntersections(soup.coordinates, soup.triangles, soup.labels, arena, ma.m_coordinatesImplicit, ma.m_triangles, ma.m_labels, true);

            //TODO: there is memory corruption in FaRMA
            //solveIntersections(soup.coordinates, soup.triangles, soup.labels, arena, ma.m_coordinates, ma.m_triangles, ma.m_labels, true);
        }
        catch (const std::runtime_error& e)
        {
            // do stuff with exception...
            return false;
        }

        ma.m_multiplier = ma.m_coordinatesImplicit.back()->toExplicit3D().X();
        computeApproximateCoordinates(ma.m_coordinatesImplicit, ma.m_coordinates);

        //leave only the sign from m_multiplier (1 or -1) and add it to the orientation result,
        ma.m_multiplier = ma.m_multiplier / std::abs(ma.m_multiplier);

        //ma.buildDebugMesh();


        if (timings)
        {
            timings->resolvingIntersections = (float)std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - point).count();
        }

        point = std::chrono::steady_clock::now();

        auto copyFunctor = [&ma](uint tId, uint newTId, MatrixMesh<MeshType>& matrixMesh)
            {
                matrixMesh.m_intersectionFaceProp[newTId] = ma.m_intersectionEdgeFaceProp[tId];
                matrixMesh.m_coplanarFaceProp[newTId] = ma.m_coplanarFace[tId];
                matrixMesh.m_faceIOLabeling[newTId] = ma.m_faceIOLabeling[tId];

                for (int i = 0; i < 3; ++i)
                {
                    matrixMesh.m_intersectionEdgeProp[matrixMesh.getFaceEdgeId(newTId, i)] = ma.m_intersectionEdgeProp[ma.getFaceEdgeId(tId, i)];
                }
            };


        std::bitset<NBIT> leftLabel;
        leftLabel[0] = 1;

        std::bitset<NBIT> rightLabel;
        rightLabel[1] = 1;

        MatrixMesh<MeshType> left(0, ma,
            [&ma, &leftLabel](uint tId)
            {
                auto faceLabel = ma.m_labels[tId];

                bool belongsToThisOperand = (faceLabel & leftLabel) != 0;

                return belongsToThisOperand;
            });

        if (!left.buildManifoldMesh())
        {
            return false;
        }

        MatrixMesh<MeshType> right(1, ma,
            [&ma, &rightLabel](uint tId)
            {
                auto faceLabel = ma.m_labels[tId];

                bool belongsToThisOperand = (faceLabel & rightLabel) != 0;

                return belongsToThisOperand;
            });


        if (!right.buildManifoldMesh())
        {
            return false;
        }

        ma.updateMatrices();
        ma.detectCoplanarFaces();
        ma.detectIntersectionEdges();

        ma.initComponentsVolume();
        left.computeComponentsVolume(ma);
        right.computeComponentsVolume(ma);

        //TODO: resolve classification for multiple faces connected to an edge (although it shouldn't happen)
        try
        {
            ma.classifyFaces();
        }
        catch (const std::exception& e)
        {
            std::cerr << "meshBoolean: " << e.what() << std::endl;
            return false;
        }

        left.copyMAProperties(copyFunctor);
        right.copyMAProperties(copyFunctor);

        if (!left.disconnectComponents(ma))
        {
            return false;
        }

        if (!right.disconnectComponents(ma))
        {
            return false;
        }

        left.classifyIsolatedComponents(right);
        right.classifyIsolatedComponents(left);


        left.classifyMeshArrangement(ma, 0, predicate);
        right.classifyMeshArrangement(ma, 1, predicate);

        MatrixMesh<MeshType> result(2, ma,
            [&ma, &predicate](uint tId) -> int
            {
                if (ma.m_coplanarFace[tId])
                {
                    return predicate(0/*TODO: make this the lowest set bit? Or even better, adjust the predicate*/, ma.m_labels[tId]);
                }

                return ma.m_boPredicates[tId];
            });

        if (!result.buildManifoldMesh(true))
        {
            return false;
        }

        if (!result.resolveNmfVerticesInHalfedgeMesh())
        {
            return false;
        }



        //TODO: enable after testing
#if 0
            //TODO: This is for mesh cut, REDO THIS!!
            //make the functor an optional input
        OpenMesh::FProp<std::bitset<NBIT>> faceOriginLabel(result.getMeshNonCost(), "bo_faceOriginLabel");
        result.copyMeshProperties([&faceOriginLabel](MatrixMesh<MeshType>& matrixMesh, MeshType& mesh, tFaceHandle fh, uint tId)
            {
                faceOriginLabel[fh] = matrixMesh.m_faceOriginLabel[tId];
            });
#endif


        output = result.getClearMesh();

        if (timings)
        {
            timings->performingBooleanOperation = (float)std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - point).count();
            std::cout << "Performing Boolean operation: " << timings->performingBooleanOperation << " ms" << std::endl;
        }

        //output.update_normals();
        return true;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////
//! \brief   Performs mesh union using new connected components analysis approach.
//!
//! \param [in/out]  output         Reference to output mesh.
//! \param           lhs            Path to first mesh.
//! \param           rhs            Path to second mesh.
//! \param           timings        Optional MeshBooleansTimings for measuring time of the operations.
//! 
//! \return                         True if successful, false otherwise. 
    template <typename MeshType>
    bool meshUnion(MeshType& output,
        std::string lhs, std::string rhs,
        MeshBooleansTimings* timings = nullptr)
    {
        auto predicateAdd = [&](int meshOperandLabel, const std::bitset<NBIT>& faceLabel) -> int
            {
                if (meshOperandLabel == 0 && faceLabel[1] == 0 && faceLabel[0] == 0)
                {
                    return 1;
                }

                if (meshOperandLabel == 1 && faceLabel[0] == 0 && faceLabel[1] == 0)
                {
                    return 1;
                }

                //and we want the coplanars as well
                if (faceLabel[0] == 1 && faceLabel[1] == 1)
                {
                    return 1;
                }

                return 0;
        };

        TriangleSoup soup{};

        uint label = 0;
        for (std::string filename : {lhs, rhs})
        {
            InputTriangleMesh m;
            m.label = label;
            load(filename, m.coordinates, m.triangles);
            addMesh(soup, m);

            ++label;
        }

        auto progress = []() { return true; };

        return meshBoolean(output, soup, predicateAdd, progress, timings);
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////
//! \brief   Performs mesh intersection using new connected components analysis approach.
//!
//! \param [in/out]  output         Reference to output mesh.
//! \param           lhs            Path to first mesh.
//! \param           rhs            Path to second mesh.
//! \param           timings        Optional MeshBooleansTimings for measuring time of the operations.
//! 
//! \return                         True if successful, false otherwise. 
    template <typename MeshType>
    bool meshIntersection(MeshType& output,
        std::string lhs, std::string rhs,
        MeshBooleansTimings* timings = nullptr)
    {
        auto predicateInt = [](int meshOperandLabel, const std::bitset<NBIT>& faceLabel) -> int
            {
                if (meshOperandLabel == 0 && faceLabel[1] == 1 && faceLabel[0] == 0)
                {
                    return 1;
                }

                if (meshOperandLabel == 1 && faceLabel[0] == 1 && faceLabel[1] == 0)
                {
                    return 1;
                }

                //and we want the coplanars as well
                if (faceLabel[0] == 1 && faceLabel[1] == 1)
                {
                    return 1;
                }

                return 0;
            };

        TriangleSoup soup{};

        uint label = 0;
        for (std::string filename : {lhs, rhs})
        {
            InputTriangleMesh m;
            m.label = label;
            load(filename, m.coordinates, m.triangles);
            addMesh(soup, m);
            //save("C:/skola/PhD/Samples/booleans/input"+std::to_string(label) + ".obj", m.coordinates, m.triangles);

            ++label;
        }

        auto progress = []() { return true; };

        return meshBoolean(output, soup, predicateInt, progress, timings);
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////
//! \brief   Performs mesh subtraction using new connected components analysis approach.
//!
//! \param [in/out]  output         Reference to output mesh.
//! \param           lhs            Path to first mesh.
//! \param           rhs            Path to second mesh.
//! \param           timings        Optional MeshBooleansTimings for measuring time of the operations.
//! 
//! \return                         True if successful, false otherwise. 
    template <typename MeshType>
    bool meshSubtraction(MeshType& output,
        std::string lhs, std::string rhs,
        MeshBooleansTimings* timings = nullptr)
    {
        auto predicateSub = [](int meshOperandLabel, const std::bitset<NBIT>& faceLabel) -> int
            {
                if (meshOperandLabel == 0 && faceLabel[1] == 0 && faceLabel[0] == 0)
                {
                    return 1;
                }

                if (meshOperandLabel == 1 && faceLabel[0] == 1 && faceLabel[1] == 0)
                {
                    return -1;
                }

                return 0;
            };

        TriangleSoup soup{};

        uint label = 0;
        for (std::string filename : {lhs, rhs})
        {
            InputTriangleMesh m;
            m.label = label;
            load(filename, m.coordinates, m.triangles);
            addMesh(soup, m);
            ++label;
        }

        auto progress = []() { return true; };

        return meshBoolean(output, soup, predicateSub, progress, timings);
    }
}//namespace