//TODO: License

#pragma once

#include "MeshArrangement.h"
//#include "HelperDefines.h"

#include <Eigen/SparseCore>

namespace DMB
{
    template<typename MeshType>
    class MeshOperand
    {
    public: //types
        using tSparseIntMatrix = Eigen::SparseMatrix<int>;
        using tPoint = typename MeshType::Point;
        using tVertexHandle = typename MeshType::VertexHandle;
        using tFaceHandle = typename MeshType::FaceHandle;
        using tEdgeHandle = typename MeshType::EdgeHandle;
        using tTriangle = std::array<tPoint, 3>;
        using tTriangleIndices = std::array<uint, 3>;

    public: //methods
        //TODO: unify names of properties
        //TODO: unify result/output usage
        //TODO: make label bitset
        MeshOperand(int label);

        //! Careful, this can change the mesh arrangement!!!
        bool createMesh(MeshArrangement& ma);

        void processComponent(const MeshArrangement& ma, const tSparseIntMatrix& edgeAdjacencyMatrix);

        void classifyComponents(const MeshArrangement& ma,const MeshType& otherMesh, const std::vector<tVertexHandle>& maToVh, int otherMeshLabel);

        static bool combineToResult(MeshArrangement& ma, MeshType& result);

        template<typename F>
        void classifyMeshArrangement(MeshArrangement& ma, const F& predicate);

        const MeshType& getMesh() const
        {
            return m_mesh;
        }

        MeshType& getMesh()
        {
            return m_mesh;
        }

        const std::vector<tVertexHandle>& meshArrangementToVhs() const
        {
            return m_meshArrangementToVhs;
        }

        bool disconnectComponents();

        int getLabel() const
        {
            return m_meshOperandLabel;
        }

    private: //methods
        void handleCoplanarFaces();
        void detectBoundaries();
        void detectIntersectionEdges(const MeshArrangement& ma);
        void discardUnusableComponents();
        uint addMeshArrangementVertex(MeshArrangement& ma, uint vertexHandle, uint triangleVertexHandle);

        //utils
        int findHighestSetBit(const std::bitset<NBIT> bitset);
        int findLowestSetBit(const std::bitset<NBIT> bitset);

    private: //attributes
        //TODO: make label bitset
        int m_meshOperandLabel;

        tSparseIntMatrix m_edgeValanceMatrix;
        std::vector<tVertexHandle> m_meshArrangementToVhs;

        //! Actual mesh of this operand
        MeshType m_mesh;

        std::vector<SIntersectionCurve> m_intersectionCurves;
        //pair <triangle Id, vertex ids>
        std::vector<std::pair<uint,tTriangleIndices>> m_coplanarTriangles;

        OpenMesh::VPropHandleT<int> m_pboundaryId;
        OpenMesh::VPropHandleT<uint> m_pOperandToMeshArrangementIndex;
        OpenMesh::VPropHandleT<bool> m_pCoplanarVertex;
        OpenMesh::EPropHandleT<bool> m_pIntersectionEdge;
        OpenMesh::EPropHandleT<bool> m_pCoplanarEdge;
        OpenMesh::FPropHandleT < std::bitset<NBIT> > m_pLabeling;
        OpenMesh::FPropHandleT<bool> m_pCoplanarFace;
        OpenMesh::FPropHandleT<bool> m_pIntersectionFace;
        OpenMesh::FPropHandleT<uint> m_pFhToMaFh;

    };

}

#include "MeshOperand.hxx"
