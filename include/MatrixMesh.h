//TODO: License

#pragma once

#include "HelperStrucutres.h"
//#include "HelperDefines.h"
#include "MeshArrangement.h"
#include "FastWindingNumber.h"

#include <Eigen/SparseCore>

namespace DMB
{
    //TODO debug: remove this template and all references to a OpenMesh...
    template <typename MeshType>
    class MatrixMesh
    {
    public: //types
        using tSparseIntMatrix = Eigen::SparseMatrix<int>;
        using tPoint = typename MeshType::Point;
        using tExactPoint = std::array<bigfloat, 3>;
        using tVertexHandle = typename MeshType::VertexHandle;
        using tFaceHandle = typename MeshType::FaceHandle;
        using tEdgeHandle = typename MeshType::EdgeHandle;
        using tHalfedgeHandle = typename MeshType::HalfedgeHandle;
        using tScalar = typename MeshType::Scalar;
        using tTriangle = std::array<tPoint, 3>;
        using tTriangleIndices = std::array<uint, 3>;
        using tFWN = typename FastWindingNumber<MeshType>;

    public:
        MatrixMesh(int label, const MeshArrangement<MeshType>& ma);

        //TODO: make a version without copyFunction
        template<typename F>
        MatrixMesh(int intLabel, const MeshArrangement<MeshType>& ma, const F& filterFunction):
            m_intLabel(intLabel),
            m_multiplier(ma.m_multiplier)
        {
            std::bitset<NBIT> bitsetLabel;
            bitsetLabel[intLabel] = 1;
            m_label = bitsetLabel;

            buildOperand(ma, filterFunction);
        }
        uint getFaceEdgeId(uint tId, int i);

        const MeshType& getMesh() const
        {
            return m_mesh;
        }

        MeshType& getMeshNonCost()
        {
            return m_mesh;
        }

        //CAREFULL!! This will screw up all mappings, use only when done with this operand!
        const MeshType& getClearMesh()
        {
            m_mesh.garbage_collection();
            return m_mesh;
        }

        void saveMatrixMesh(std::string outputFileName);

        bool buildManifoldMesh(bool considerOrigin = false);
        bool resolveNmfVerticesInHalfedgeMesh();
        void buildDebugMesh();

        void flipToConsistentOrientation(MeshArrangement<MeshType>& ma);
        void adjustOrientation(MeshArrangement<MeshType>& ma);

        template<typename C>
        void copyMAProperties(const C& copyFunctor)
        {
            //TODO: MAKE SOME PROPERTY MANAGEMENT SYSTEM!!
            m_coplanarFaceProp.clear();
            m_coplanarFaceProp.resize(m_triangles.size() / 3, false);
            m_intersectionFaceProp.clear();
            m_intersectionFaceProp.resize(m_triangles.size() / 3, false);
            m_faceIOLabeling.clear();
            m_faceIOLabeling.resize(m_triangles.size() / 3, {});
            m_intersectionEdgeProp.clear();
            m_intersectionEdgeProp.resize(m_nEdges + 1, false);
            //m_faceOriginLabel.clear();
            //m_faceOriginLabel.resize(m_triangles.size() / 3);

            for (int i = 0; i < std::size(m_triangles); i += 3)
            {
                uint tId = i / 3;
                copyFunctor(m_tIdToOriginalTId[tId], tId, *this);
            }
        }

        template<typename C>
        void copyMeshProperties(const C& copyFunctor)
        {
            OpenMesh::FProp<uint> pFhToMaFh(m_mesh, m_pFhToMaFh);

            for (auto fh : m_mesh.faces())
            {
                copyFunctor(*this, m_mesh, fh, pFhToMaFh[fh]);
            }
        }

        bool disconnectComponents(MeshArrangement<MeshType>& ma);

        void classifyIsolatedComponents(MatrixMesh<MeshType>& other);
        std::shared_ptr<tFWN> getFWN();

        template<typename F>
        void classifyMeshArrangement(MeshArrangement<MeshType>& ma, int operandLabel, const F& predicate);

        template<typename F>
        void classifyMeshArrangementWithMAFaceHandle(MeshArrangement<MeshType>& ma, int operandLabel, const F& predicate);

        void markSelfIntersectingFaces();

        void resolveConflictingComponents(MeshArrangement<MeshType>& ma);
        void computeComponentsVolume(MeshArrangement<MeshType>& ma);

        //TODO: MAKE ACCESSORS FOR ThIS!!!!
        std::vector<bool> m_coplanarFaceProp;
        std::vector<bool> m_intersectionFaceProp;
        std::vector<std::bitset<NBIT>> m_faceIOLabeling;
        std::vector<bool> m_intersectionEdgeProp;
        std::vector<std::bitset<NBIT>> m_faceOriginLabel;
        std::bitset < NBIT > m_label;

    private: //methods

        template<typename F>
        void buildOperand(const MeshArrangement<MeshType>& ma, const F& filterFunction)
        {
            std::vector<int> maToOperandVertices;
            maToOperandVertices.resize(ma.m_coordinatesImplicit.size(), -1);

            m_tIdToOriginalTId.reserve(ma.m_triangles.size());

            uint newVhId = 0;

            for (int i = 0; i < std::size(ma.m_triangles); i += 3)
            {
                uint tId = i / 3;

                /*
                 * filterResult ==  0 -> not part of mesh
                 * filterResult ==  1 -> part of mesh as is
                 * filterResult == -1 -> part of mesh but flipped
                */
                int filterResult = filterFunction(tId);

                if (filterResult == 0)
                {
                    continue;
                }

                uint maVh0 = ma.m_triangles[i + 0];
                uint maVh1 = ma.m_triangles[i + 1];
                uint maVh2 = ma.m_triangles[i + 2];

                std::array<uint, 3> vertices;
                vertices[0] = maVh0;
                vertices[1] = maVh1;
                vertices[2] = maVh2;

                if (filterResult == -1)
                {
                    //flip triangle
                    vertices[1] = maVh2;
                    vertices[2] = maVh1;
                }

                for(auto maVh : vertices)
                {
                    //uint maVh = ma.m_triangles[i + j];

                    if (maToOperandVertices[maVh] < 0)
                    {
                        maToOperandVertices[maVh] = newVhId;
                        ++newVhId;

                        auto* gp = ma.m_coordinatesImplicit[maVh];

                        m_coordinatesImplicit.push_back(gp);
                        m_coordinates.push_back(ma.m_coordinates[maVh * 3 + 0]);
                        m_coordinates.push_back(ma.m_coordinates[maVh * 3 + 1]);
                        m_coordinates.push_back(ma.m_coordinates[maVh * 3 + 2]);
                    }

                    uint vh = maToOperandVertices[maVh];

                    m_triangles.push_back(vh);
                }

                //TODO: this is a mess...
                m_faceOriginLabel.push_back(ma.m_labels[tId]);
                m_tIdToOriginalTId.push_back(tId);
            }

            updateMatrices();
        }

        //TODO: this needs to be done in a more inteligent way...
        void updateMatrices();
        void buildVertexEdgeMatrix();
        void buildEdgeFaceAdjacency();
        void buildVertexFaceAdjacency();
        void buildEdgeVertices();

        void buildOperand(const MeshArrangement<MeshType>& ma);

        void detectNonManifolds();
        void detectNonManifoldEdges();
        void detectNonManifoldVertices();
        void detectNonManifoldFaces();
        void detectNonOrientableEdges();

        void addManifoldFaces();
        void addNonManifoldEdges(bool considerOrigin = false);
        void addNonManifoldVertices();
        void resolveNonOrientableEdges();

        uint duplicateVertex(uint vId);

        //mesh functions
        void transferArrangementPropertiesToMesh();
        void detectBoundaries();
        void handleCoplanarFaces();

        double calcSignedVolume(const std::vector<OpenMesh::SmartFaceHandle>& component);
        bigfloat calcSignedVolumeExact(const std::vector<OpenMesh::SmartFaceHandle>& component);
        int calcVolumeSignExact(const std::vector<OpenMesh::SmartFaceHandle>& component);

        void debug_CouldNotAddFace(tVertexHandle v0, tVertexHandle v1, tVertexHandle v2);
        void debug_CurrentMesh();

        template<typename F>
        void debug_showMesh(const F& func);

    private:
        //TODO: debug - delete
        int m_intLabel;


        double m_multiplier;

        //TODO: debug?
        std::vector<double> m_coordinates;
        std::vector<genericPoint*> m_coordinatesImplicit;
        std::vector<uint> m_triangles;

        //connectivity
        tSparseIntMatrix m_vertexEdgeMatrix;
        std::vector<std::vector<uint>> m_edgeFaceAdjacency;
        std::vector<std::vector<uint>> m_vertexFaceAdjacency;
        std::vector<std::pair<uint, uint>> m_edgeVertices;

        size_t m_nEdges;
        size_t m_nFaces;
        size_t m_nVertices;

        int m_nNonManifoldVertices;

        //properties - TODO: rework OpenMesh style
        std::vector<bool> m_nonManifoldVertexProp;
        std::vector<bool> m_nonManifoldEdgeProp;
        std::vector<bool> m_nonManifoldEdgeVertexProp;
        std::vector<bool> m_nonManifoldEdgeFaceProp;
        std::vector<bool> m_nonManifoldFaceProp;
        std::vector<bool> m_nonOrientableEdgeProp;
        std::vector<bool> m_nonOrientableEdgeFaceProp;
        std::vector<bool> m_faceAdded;
        std::vector<bool> m_faceFlipped;
        std::vector<bool> m_selfIntersectingFace;

        std::vector<uint> m_tIdToOriginalTId;

        //TODO delete debug
        std::vector<tVertexHandle> m_matrixVhToOMVh;

        MeshType m_mesh;

        //mesh properties
        OpenMesh::FPropHandleT<uint> m_pFhToMaFh;
        OpenMesh::FPropHandleT<bool> m_pCoplanarFace;
        OpenMesh::FPropHandleT < std::bitset<NBIT> > m_pFaceIOLabeling;
        OpenMesh::FPropHandleT<bool> m_pIntersectionFace;
        OpenMesh::VPropHandleT<int> m_pboundaryId;
        OpenMesh::EPropHandleT<bool> m_pIntersectionEdge;
        OpenMesh::VPropHandleT<uint> m_pVhToMaVId;
        OpenMesh::VPropHandleT<bool> m_pCoplanarVertex;
        OpenMesh::EPropHandleT<bool> m_pCoplanarEdge;

        //can this be optimized somehow? this can take up a lot of memory, maybe it's an unnecessary copy
        std::vector<std::vector<tFaceHandle>> m_isolatedComponents;
    };


}

#include "MatrixMesh.hxx"
