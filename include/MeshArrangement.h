//TODO: License

#pragma once

//#include "HelperDefines.h"
#include "HelperStrucutres.h"

#include <bitset>


#include <Eigen/SparseCore>

////////////////////////////////////////////////////////////////////////////////////////////////////
//! \brief   Encapsulation of output data for FaRMA algorithms.
//! The library uses vectors for everything. Consult the library API for more info on these.
////////////////////////////////////////////////////////////////////////////////////////////////////

namespace DMB
{
    //TODO: remove template in the future
    template<typename MeshType>
    class MeshArrangement
    {


    public:
        MeshArrangement();

        void updateMatrices();
        void detectCoplanarFaces();
        void detectIntersectionEdges();

        void classifyFaces();

        void buildDebugMesh();

        uint getFaceEdgeId(uint tId, int i);
        uint getVertexEdgeId(uint v0, uint v1) const;
        bool isEdgeManifold(uint eId) const;
        const std::vector<uint>& getFaceAdjacentEdges(uint eId) const;
        //TODO: make this more intelligent...
    public:
        std::vector<double> m_coordinates;
        std::vector<genericPoint*> m_coordinatesImplicit;
        std::vector<uint> m_triangles;
        std::vector<std::bitset<NBIT>> m_labels;
        std::vector<int> m_boPredicates;
        double m_multiplier;

        //todo: make accessors
        //properties - TODO: rework OpenMesh style
        std::vector<bool> m_intersectionEdgeProp;
        std::vector<bool> m_intersectionEdgeFaceProp;
        std::vector<bool> m_intersectionVertexProp;
        std::vector<bool> m_coplanarFace;
        std::vector<std::bitset<NBIT>> m_faceIOLabeling;

        std::vector<bool> m_conflictingFace;
        std::vector<bool> m_checkOrientation;
        std::vector<double> m_componentsVolume;
        std::vector<int> m_faceComponentId;
        std::vector<int> m_faceClassifiedByComponent;
        void initComponentsVolume();


    private: //methods
        void buildVertexEdgeMatrix();
        void buildEdgeFaceAdjacency();
        void buildEdgeVertices();

        template<typename F>
        void debug_showMesh(const F& func);

    private: //properties

        //connectivity
        tSparseIntMatrix m_vertexEdgeMatrix;
        std::vector<std::vector<uint>> m_edgeFaceAdjacency;
        //std::vector<std::vector<uint>> m_vertexFaceAdjacency;
        std::vector<std::pair<uint, uint>> m_edgeVertices;

        size_t m_nEdges;
        size_t m_nFaces;
        size_t m_nVertices;



    };
}

#include "MeshArrangement.hxx"
