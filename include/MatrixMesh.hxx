//TODO: License

#include "MatrixMesh.h"
#include <unordered_set>
#include <vector>
#include <stack>
#include <array>
#include <atomic>

/*
* self intersections https://diglib.eg.org/bitstream/handle/10.2312/conf.EG2012.tutorials.t4/t4.pdf?sequence=1
*     [CK10] use a intermediate BSP
*/

template<typename MeshType>
inline DMB::MatrixMesh<MeshType>::MatrixMesh(int label, const MeshArrangement<MeshType>& ma):
    m_nEdges(0),
    m_nFaces(0),
    m_nVertices(0),
    m_nNonManifoldVertices(0),
    m_multiplier(ma.m_multiplier)
{
    std::bitset<NBIT> bitsetLabel;
    bitsetLabel[label] = 1;
    m_label = bitsetLabel;

    //buildOperand(ma);
    //buildManifoldMesh();
}




template<typename MeshType>
inline void DMB::MatrixMesh<MeshType>::updateMatrices()
{
    m_nFaces = m_triangles.size();
    m_nVertices = m_coordinatesImplicit.size();

    buildVertexEdgeMatrix();
    buildEdgeFaceAdjacency();
    buildVertexFaceAdjacency();
    buildEdgeVertices();
}

template<typename MeshType>
inline void DMB::MatrixMesh<MeshType>::buildVertexEdgeMatrix()
{
    std::vector<tTriplet> tripletList;
    tripletList.reserve(std::size(m_coordinates));

    for (int i = 0; i < std::size(m_triangles); i += 3)
    {
        const uint v0 = m_triangles[i + 0];
        const uint v1 = m_triangles[i + 1];
        const uint v2 = m_triangles[i + 2];

        tripletList.emplace_back(v0, v1, -1);
        tripletList.emplace_back(v1, v2, -1);
        tripletList.emplace_back(v2, v0, -1);
    }

    //TODO: maybe the copying of matrices is not neccessary - but I am not sure about the Eigen resize, if it clears all previously set values
    //m_vertexEdgeMatrix.resize(std::size(m_coordinatesImplicit), std::size(m_coordinatesImplicit));
    //work with m_vertexEdgeMatrix only
    //....

    tSparseIntMatrix edges(std::size(m_coordinatesImplicit), std::size(m_coordinatesImplicit));
    edges.setFromTriplets(tripletList.begin(), tripletList.end());

    int edgeId = 1;
    m_nEdges = 0;

    for (int k = 0; k < edges.outerSize(); ++k)
    {
        for (tSparseIntMatrix::InnerIterator it(edges, k); it; ++it)
        {
            auto r = it.row();
            auto c = it.col();

            //existing edge without assigned Id -> assign new Id
            if (edges.coeff(r, c) < 0)
            {
                edges.coeffRef(r, c) = edgeId;

                //try to assign the same Id for twin edge
                if (edges.coeff(c, r) != 0)
                {
                    edges.coeffRef(c, r) = edgeId;
                }

                ++edgeId;
                ++m_nEdges;
            }
        }
    }

    m_vertexEdgeMatrix = edges;
}

template<typename MeshType>
inline void DMB::MatrixMesh<MeshType>::buildEdgeFaceAdjacency()
{
    m_edgeFaceAdjacency.clear();
    m_edgeFaceAdjacency.resize(m_nEdges + 1); // edge Ids starts from 1

    // preallocate to minize reallocations
    for (auto& adjacency : m_edgeFaceAdjacency)
    {
        adjacency.reserve(16);
    }

    for (uint i = 0; i < std::size(m_triangles); i += 3)
    {
        uint tId = i / 3;

        uint v0 = m_triangles[i + 0];
        uint v1 = m_triangles[i + 1];
        uint v2 = m_triangles[i + 2];

        int e0 = m_vertexEdgeMatrix.coeff(v0, v1);
        int e1 = m_vertexEdgeMatrix.coeff(v1, v2);
        int e2 = m_vertexEdgeMatrix.coeff(v2, v0);

        m_edgeFaceAdjacency[e0].push_back(tId);
        m_edgeFaceAdjacency[e1].push_back(tId);
        m_edgeFaceAdjacency[e2].push_back(tId);
    }
}

template<typename MeshType>
inline void DMB::MatrixMesh<MeshType>::buildVertexFaceAdjacency()
{
    m_vertexFaceAdjacency.clear();
    m_vertexFaceAdjacency.resize(std::size(m_coordinatesImplicit));

    // preallocate to minize reallocations
    for (auto& adjacency : m_vertexFaceAdjacency)
    {
        adjacency.reserve(16);
    }

    for (uint i = 0; i < std::size(m_triangles); i += 3)
    {
        uint tId = i / 3;

        uint v0 = m_triangles[i + 0];
        uint v1 = m_triangles[i + 1];
        uint v2 = m_triangles[i + 2];

        m_vertexFaceAdjacency[v0].push_back(tId);
        m_vertexFaceAdjacency[v1].push_back(tId);
        m_vertexFaceAdjacency[v2].push_back(tId);
    }
}

template<typename MeshType>
inline void DMB::MatrixMesh<MeshType>::buildEdgeVertices()
{
    m_edgeVertices.clear();
    m_edgeVertices.resize(m_nEdges + 1, {0,0});

    for (int k = 0; k < m_vertexEdgeMatrix.outerSize(); ++k)
    {
        for (tSparseIntMatrix::InnerIterator it(m_vertexEdgeMatrix, k); it; ++it)
        {
            Eigen::Index r = it.row();
            Eigen::Index c = it.col();
            auto eId = m_vertexEdgeMatrix.coeff(r, c);

            if (eId > 0) //This edge id is twice in the m_vertexEdgeMatrix, so it will be overwritten, is it ok?
            {
                m_edgeVertices[eId] = std::make_pair((uint)r, (uint)c);
            }
        }
    }
}

template<typename MeshType>
inline void DMB::MatrixMesh<MeshType>::detectNonManifoldEdges()
{
    int nNonManifoldEdges = 0;
    for(uint eId = 0; eId < m_edgeFaceAdjacency.size(); ++eId)
    {
        if (m_edgeFaceAdjacency[eId].size() > 2)
        {
            m_nonManifoldEdgeProp[eId] = true;
            ++nNonManifoldEdges;

            for (auto tId : m_edgeFaceAdjacency[eId])
            {
                m_nonManifoldFaceProp[tId] = true;
                m_nonManifoldEdgeFaceProp[tId] = true;

                //TODO: make edge vertex matrix for this, this is stupid...
                uint v0 = m_triangles[tId * 3 + 0];
                uint v1 = m_triangles[tId * 3 + 1];
                uint v2 = m_triangles[tId * 3 + 2];

                uint e0 = m_vertexEdgeMatrix.coeff(v0, v1);
                uint e1 = m_vertexEdgeMatrix.coeff(v1, v2);
                uint e2 = m_vertexEdgeMatrix.coeff(v2, v0);

                //TODO: use vertex edge adjacency for this
                if (eId == e0)
                {
                    m_nonManifoldEdgeVertexProp[v0] = true;
                    m_nonManifoldEdgeVertexProp[v1] = true;
                }

                if (eId == e1)
                {
                    m_nonManifoldEdgeVertexProp[v1] = true;
                    m_nonManifoldEdgeVertexProp[v2] = true;
                }

                if (eId == e2)
                {
                    m_nonManifoldEdgeVertexProp[v2] = true;
                    m_nonManifoldEdgeVertexProp[v0] = true;
                }
            }
        }
    }

    //std::cout << "nmf edges: " << nNonManifoldEdges << std::endl;
}

template<typename MeshType>
inline void DMB::MatrixMesh<MeshType>::detectNonManifoldVertices()
{
    m_nNonManifoldVertices = 0;

    //tmp prop
    // Note: initialized to bool{} (false) due to C++20 constructor of std::atomic<bool>
    // see: https://en.cppreference.com/w/cpp/atomic/atomic/atomic
    //std::vector<std::atomic<bool>> faceVisited(m_nFaces);
    std::vector<bool> faceVisited(m_nFaces, false);

    {
        std::stack<uint, std::vector<uint>> facesToVisit;

        for (int ivId = 0; ivId < m_vertexFaceAdjacency.size(); ++ivId)
        {
            uint vId = static_cast<uint>(ivId);

            assert(facesToVisit.empty());
            //TODO: do i need the components?
            //std::vector<std::vector<uint>> components;
            unsigned componentsCount = 0;
            bool isVertexClosed = false;

            for (auto vertexFace : m_vertexFaceAdjacency[vId])
            {
                if (faceVisited[vertexFace])
                {
                    continue;
                }

                facesToVisit.push(vertexFace);

                //std::vector<uint> component;
                bool componentClosed = true;

                while (!facesToVisit.empty())
                {
                    auto topFace = facesToVisit.top();
                    facesToVisit.pop();

                    if (faceVisited[topFace])
                    {
                        continue;
                    }

                    //component.push_back(topFace);

                    std::array<uint, 3> faceVertices;

                    faceVertices[0] = (m_triangles[topFace * 3 + 0]);
                    faceVertices[1] = (m_triangles[topFace * 3 + 1]);
                    faceVertices[2] = (m_triangles[topFace * 3 + 2]);

                    auto mid = std::find_if(faceVertices.begin(), faceVertices.end(),
                        [&](const auto& faceVertex)
                        {
                            return faceVertex == vId;
                        });

                    //vector sorted in such a way, that vId vertex is first
                    std::rotate(faceVertices.begin(), mid, faceVertices.end());

                    int e0 = m_vertexEdgeMatrix.coeff(vId, faceVertices[1]);
                    int e1 = m_vertexEdgeMatrix.coeff(vId, faceVertices[2]);

                    //component is closed if all edges are closed
                    componentClosed = componentClosed && (m_edgeFaceAdjacency[e0].size() == 2 && m_edgeFaceAdjacency[e1].size() == 2);

                    //TODO: skip non-manifold edges?
                    if (!m_nonManifoldEdgeProp[e0])
                    {
                        for (auto fId : m_edgeFaceAdjacency[e0])
                        {
                            facesToVisit.push(fId);
                        }
                    }

                    if (!m_nonManifoldEdgeProp[e1])
                    {
                        for (auto fId : m_edgeFaceAdjacency[e1])
                        {
                            facesToVisit.push(fId);
                        }
                    }


                    faceVisited[topFace] = true;
                }

                isVertexClosed = isVertexClosed || componentClosed; // at least one component is closed

                componentsCount++;
            }

            //reset tag
            for (auto vertexFace : m_vertexFaceAdjacency[vId])
            {
                faceVisited[vertexFace] = false;
            }

            //this is non-manifold vertex
            if (componentsCount > 1 && isVertexClosed)
            {
                {
                    m_nonManifoldVertexProp[vId] = true;
                    ++m_nNonManifoldVertices;
                }
            }
        }
    }

    //std::cout << "nmf vertices: " << m_nNonManifoldVertices << std::endl;

}

template<typename MeshType>
inline void DMB::MatrixMesh<MeshType>::detectNonManifoldFaces()
{
    //std::vector<bool> faceVisited(m_nFaces, false);

    //TODO: parallel
#pragma omp parallel for
    for (int i = 0; i < std::size(m_triangles); i += 3)
    {
        uint tId = i / 3;

        for (int j = 0; j < 3; ++j)
        {
            uint v = m_triangles[tId * 3 + j];

            if (m_nonManifoldVertexProp[v] || m_nonManifoldEdgeVertexProp[v])
            {
                m_nonManifoldFaceProp[tId] = true;

            }
        }
    }
}

template<typename MeshType>
inline void DMB::MatrixMesh<MeshType>::detectNonOrientableEdges()
{
    int nNonOrientableEdges = 0;

    std::vector<tTriplet> tripletList;
    tripletList.reserve(std::size(m_triangles));

    for (int i = 0; i < std::size(m_triangles); i += 3)
    {
        uint v0 = m_triangles[i + 0];
        uint v1 = m_triangles[i + 1];
        uint v2 = m_triangles[i + 2];

        tripletList.emplace_back(v0, v1, 1);
        tripletList.emplace_back(v1, v2, 1);
        tripletList.emplace_back(v2, v0, 1);
    }

    //Edge adjacency matrx (row major sparse matrix) - rows and columns are indices of vertices, the value indicates how many edges are between this vertices
    tSparseIntMatrix edgeAdjacency(std::size(m_coordinatesImplicit), std::size(m_coordinatesImplicit));
    edgeAdjacency.setFromTriplets(tripletList.begin(), tripletList.end());

    tSparseIntMatrix edgeValanceMatrix = tSparseIntMatrix(edgeAdjacency.transpose()) + edgeAdjacency;

    m_nonOrientableEdgeProp.clear();
    m_nonOrientableEdgeProp.resize(m_nEdges + 1, false);
    m_nonOrientableEdgeFaceProp.clear();
    m_nonOrientableEdgeFaceProp.resize(m_nFaces, false);

    for (int k = 0; k < edgeAdjacency.outerSize(); ++k)
    {
        for (tSparseIntMatrix::InnerIterator it(edgeAdjacency, k); it; ++it)
        {
            auto r = it.row();
            auto c = it.col();

            if (edgeAdjacency.coeff(r, c) > 1)
            {
                //std::cout << "edge valence [" << r << ", " << c <<"]: " << edgeAdjacency.coeff(r, c) << std::endl;
                if (edgeAdjacency.coeff(r, c) != edgeAdjacency.coeff(c, r))
                {
                    m_nonOrientableEdgeProp[m_vertexEdgeMatrix.coeff(r, c)] = true;
                    ++nNonOrientableEdges;
                }
            }
        }
    }

    for (auto eId = 0; eId < m_edgeFaceAdjacency.size(); ++eId)
    {
        if (!m_nonOrientableEdgeProp[eId])
        {
            continue;
        }

        //std::cout << "This non orientable edge has: " << m_edgeFaceAdjacency[eId].size() << " faces" << std::endl;

        for (auto tId : m_edgeFaceAdjacency[eId])
        {
            m_nonOrientableEdgeFaceProp[tId] = true;
        }
    }
}

template<typename MeshType>
inline void DMB::MatrixMesh<MeshType>::addManifoldFaces()
{
    OpenMesh::FProp<uint> pFhToMaFh(m_mesh, m_pFhToMaFh);
    OpenMesh::VProp<uint> pVhToMaVId(m_mesh, m_pVhToMaVId);
    
    auto addFace = [this, &pFhToMaFh, &pVhToMaVId](uint tId) -> OpenMesh::SmartFaceHandle
    {

        uint v0 = m_triangles[tId * 3 + 0];
        uint v1 = m_triangles[tId * 3 + 1];
        uint v2 = m_triangles[tId * 3 + 2];

        tTriangle triangle =
        {
            tPoint { m_coordinates[v0 * 3 + 0], m_coordinates[v0 * 3 + 1], m_coordinates[v0 * 3 + 2] },
            tPoint { m_coordinates[v1 * 3 + 0], m_coordinates[v1 * 3 + 1], m_coordinates[v1 * 3 + 2] },
            tPoint { m_coordinates[v2 * 3 + 0], m_coordinates[v2 * 3 + 1], m_coordinates[v2 * 3 + 2] },
        };

        //add vertices, if not already present in the mesh
        int j = 0;
        for (auto v : { v0, v1, v2 })
        {
            if (!m_matrixVhToOMVh[v].is_valid())
            {
                auto newVh = m_mesh.add_vertex(triangle[j]);
                m_matrixVhToOMVh[v] = newVh;
                pVhToMaVId[newVh] = v;
            }
            ++j;
        }

        auto vh0 = m_matrixVhToOMVh[v0];
        auto vh1 = m_matrixVhToOMVh[v1];
        auto vh2 = m_matrixVhToOMVh[v2];


        //vhToMatrixId[vh0] = v0;
        //vhToMatrixId[vh1] = v1;
        //vhToMatrixId[vh2] = v2;

        auto fh = m_mesh.add_face(vh0, vh1, vh2);

        if (!fh.is_valid())
        {
            //std::cout << "Could not add manifold face" << std::endl;
            return fh;
        }

        pFhToMaFh[fh] = tId;
        m_faceAdded[tId] = true;

        return fh;
    };



    for (uint i = 0; i < std::size(m_triangles); i += 3)
    {
        uint tId = i / 3;


        if (m_nonManifoldFaceProp[tId])
        {
            continue;
        }

        uint v0 = m_triangles[i + 0];
        uint v1 = m_triangles[i + 1];
        uint v2 = m_triangles[i + 2];

        int e0 = m_vertexEdgeMatrix.coeff(v0, v1);
        int e1 = m_vertexEdgeMatrix.coeff(v1, v2);
        int e2 = m_vertexEdgeMatrix.coeff(v2, v0);

        addFace(tId);

    }

    //DMB::saveMesh(m_mesh, "C:/T3D/Samples/Booleans/meshWithoutNonManifolds.obj");

}

template<typename MeshType>
inline void DMB::MatrixMesh<MeshType>::addNonManifoldEdges(bool considerOrigin /*= false*/)
{
    auto isTwoManifoldInMesh = [this](uint eId, std::bitset<NBIT> origin)
        {
            if (!m_nonManifoldEdgeProp[eId])
            {
                return true;
            }

            int nFaces = 0;

            for (auto tId : m_edgeFaceAdjacency[eId])
            {
                //std::cout << "current label: " << m_faceOriginLabel[tId] << std::endl;
                //std::cout << "label to check: " << (m_faceOriginLabel[tId] & origin) << std::endl;

                if ((m_faceOriginLabel[tId] & origin).count() == 1)
                {
                    ++nFaces;
                }
            }

            //std::cout << "faces around edge with the same origin: " << nFaces << " faces total: " << m_edgeFaceAdjacency[eId].size() << std::endl;

            return nFaces == 2;
        };


    //split faces around non-manifold edges into manifold connected components
    std::vector<std::vector<uint>> components;

    //which component references this face?
    std::vector<std::vector<size_t>> faceComponentsMap;
    faceComponentsMap.resize(m_nFaces);

    //create components
    {
        std::vector<bool> componenetFaceVisited(m_nFaces, false);

        std::vector<uint> componentIdToEh;


        for (auto eId = 0; eId < m_edgeFaceAdjacency.size(); ++eId)
        {
            if (!m_nonManifoldEdgeProp[eId])
            {
                continue;
            }

            //this is non manifold edge

            //get edge vertices of this edge
            auto edgeVertices = m_edgeVertices[eId];

            auto v0 = edgeVertices.first;
            auto v1 = edgeVertices.second;

            //TODO: edge adjacent faces would be enough
            //gather all faces around end vertices
            std::vector<uint> allFaces;

            for (auto tId : m_edgeFaceAdjacency[eId])
            {
                if (!m_nonManifoldFaceProp[tId])
                {
                    continue;
                }

                allFaces.push_back(tId);
            }

            std::vector<uint> allVisitedFaces;

            //sort faces into separated connected components

            for (auto seedFh : allFaces)
            {
                if (componenetFaceVisited[seedFh])
                {
                    continue;
                }

                std::vector<uint> component;
                std::queue<uint> facesToVisit;
                facesToVisit.push(seedFh);

                auto seedOriginLabel = m_faceOriginLabel[seedFh];

                while (!facesToVisit.empty())
                {
                    auto topFace = facesToVisit.front();
                    facesToVisit.pop();

                    if (componenetFaceVisited[topFace] /*|| !m_nonManifoldFaceProp[topFace]*/) //Not neccessary, all faces should be non manifold
                    {
                        continue;
                    }

                    component.push_back(topFace);
                    allVisitedFaces.push_back(topFace);
                    componenetFaceVisited[topFace] = true;


                    std::vector<uint> faceVertices;
                    faceVertices.reserve(3);
                    faceVertices.push_back(m_triangles[topFace * 3 + 0]);
                    faceVertices.push_back(m_triangles[topFace * 3 + 1]);
                    faceVertices.push_back(m_triangles[topFace * 3 + 2]);

                    int e0 = m_vertexEdgeMatrix.coeff(faceVertices[0], faceVertices[1]);
                    int e1 = m_vertexEdgeMatrix.coeff(faceVertices[1], faceVertices[2]);
                    int e2 = m_vertexEdgeMatrix.coeff(faceVertices[2], faceVertices[0]);


                    for (auto e : { e0,e1,e2 })
                    {
                        if (considerOrigin && m_nonManifoldEdgeProp[e] && isTwoManifoldInMesh(e, seedOriginLabel))
                        {
                            for (auto adjFace : m_edgeFaceAdjacency[e])
                            {
                                if ((m_faceOriginLabel[adjFace] & seedOriginLabel).count() >= 1) //this takes coplanars into account, at least one bit must be set
                                {
                                    facesToVisit.push(adjFace);
                                }
                            }

                            continue;
                        }

                        //skip nmf edges
                        if (m_nonManifoldEdgeProp[e])
                        {
                            continue;
                        }

                        //TODO: skip non orientable edges as well
                        if (m_nonOrientableEdgeProp[e])
                        {
                            continue;
                        }

                        //skip edges non connected to v0 nor v1
                        bool connectedToV0 = m_edgeVertices[e].first == v0 || m_edgeVertices[e].second == v0;
                        bool connectedToV1 = m_edgeVertices[e].first == v1 || m_edgeVertices[e].second == v1;

                        if (!(connectedToV0 || connectedToV1))
                        {
                            continue;
                        }

                        for (auto adjFace : m_edgeFaceAdjacency[e])
                        {
                            
                            facesToVisit.push(adjFace);
                        }

                    }
                }

                if (!component.empty())
                {
                    auto componentId = components.size();

                    componentIdToEh.push_back(eId);
                    components.push_back(component);

                    //mark that all the faces are referenced by this component
                    for (auto componentFace : component)
                    {
                        faceComponentsMap[componentFace].push_back(componentId);
                    }

                }
            }

            //reset tagged status
            for (auto visitedFace : allVisitedFaces)
            {
                componenetFaceVisited[visitedFace] = false;
            }
        }
    }

    std::vector<std::vector<uint>> mergedComponents;
    //merge components
    {
        std::vector<bool> componentProcessed;
        componentProcessed.resize(components.size(), false);

        std::vector<bool> componentFaceAdded;
        componentFaceAdded.resize(m_nFaces, false);

        for (int i = 0; i < components.size(); ++i)
        {
            if (componentProcessed[i])
            {
                continue;
            }

            std::queue<size_t> componentsToVisit;
            std::vector<uint> mergedComponent;

            std::map<uint, int> referencedNonManifoldEdges;
            //std::unordered_multiset<uint> referencedNonManfioldEdges;

            auto addToMergedComponent = [this, &referencedNonManifoldEdges, &mergedComponent](uint tId)
            {
                std::vector<uint> faceVertices;
                faceVertices.reserve(3);
                faceVertices.push_back(m_triangles[tId * 3 + 0]);
                faceVertices.push_back(m_triangles[tId * 3 + 1]);
                faceVertices.push_back(m_triangles[tId * 3 + 2]);

                int e0 = m_vertexEdgeMatrix.coeff(faceVertices[0], faceVertices[1]);
                int e1 = m_vertexEdgeMatrix.coeff(faceVertices[1], faceVertices[2]);
                int e2 = m_vertexEdgeMatrix.coeff(faceVertices[2], faceVertices[0]);


                for (auto e : { e0,e1,e2 })
                {
                    if (m_nonManifoldEdgeProp[e])
                    {
                        if (referencedNonManifoldEdges.find(e) == referencedNonManifoldEdges.end())
                        {
                            referencedNonManifoldEdges[e] = 0;
                        }

                        referencedNonManifoldEdges[e] += 1;
                    }
                }

                mergedComponent.push_back(tId);
            };


            //let's take a look at this component
            for (auto tId : components[i])
            {
                if (componentFaceAdded[tId])
                {
                    continue;
                }

                addToMergedComponent(tId);
                componentFaceAdded[tId] = true;

                //which components reference this face?
                for (auto referencingComponentId : faceComponentsMap[tId])
                {
                    componentsToVisit.push(referencingComponentId);
                }
            }

            //this component was processed
            componentProcessed[i] = true;

            //now go through all 
            while (!componentsToVisit.empty())
            {
                auto topComponent = componentsToVisit.front();
                componentsToVisit.pop();

                if (componentProcessed[topComponent])
                {
                    continue;
                }

                //add all faces from adjacent component (that were not added yet)
                for (auto tId : components[topComponent])
                {
                    if (componentFaceAdded[tId])
                    {
                        continue;
                    }

                    //mergedComponent.push_back(tId);
                    addToMergedComponent(tId);
                    componentFaceAdded[tId] = true;

                    //which components reference this face?
                    for (auto referencingComponentId : faceComponentsMap[tId])
                    {
                        componentsToVisit.push(referencingComponentId);
                    }
                }

                //this component was processed as well
                componentProcessed[topComponent] = true;
            }

            std::unordered_set<uint> referencedVertices;

            for (auto [key, count] : referencedNonManifoldEdges)
            {
                referencedVertices.insert(m_edgeVertices[key].first);
                referencedVertices.insert(m_edgeVertices[key].second);
               
            }

            std::map<uint, uint> vertexDuplicateMap;

            //duplicate vertices of edges that became manifold
            for (uint vId : referencedVertices)
            {
                uint newVId = duplicateVertex(vId);
                vertexDuplicateMap[vId] = newVId;
            }

            for (auto tId : mergedComponent)
            {
                for (int j = 0; j < 3; ++j)
                {
                    uint vId = m_triangles[tId * 3 + j];
                    if (referencedVertices.find(vId) != referencedVertices.end())
                    {
                        uint newVId = vertexDuplicateMap[vId];
                        m_triangles[tId * 3 + j] = newVId;
                    }
                }
            }


            mergedComponents.push_back(mergedComponent);
        }
    }

    updateMatrices();
    detectNonManifolds();

    //now only non manifold edges that could not be merged remain
    //duplicate edges
    for (auto eId = 0; eId < m_edgeFaceAdjacency.size(); ++eId)
    {
        if (!m_nonManifoldEdgeProp[eId] || m_edgeFaceAdjacency[eId].empty())
        {
            continue;
        }

        //this is non manifold edge

        //get edge vertices of this edge
        auto edgeVertices = m_edgeVertices[eId];

        auto v0 = edgeVertices.first;
        auto v1 = edgeVertices.second;

        //get faces connected to this edge and duplicate the vertices
        for (int i = 1 /*skip first face?*/; i < m_edgeFaceAdjacency[eId].size(); ++i)
        {
            uint tId = m_edgeFaceAdjacency[eId][i];

            for (int j = 0; j < 3; ++j)
            {
                uint newVId = duplicateVertex(m_triangles[tId * 3 + j]);
                m_triangles[tId * 3 + j] = newVId;
            }
        }

    }

    updateMatrices();
    detectNonManifolds();

    //m_mesh.delete_isolated_vertices();
    //m_mesh.garbage_collection();
    //DMB::saveMesh(m_mesh, "C:/skola/PhD/Samples/booleans/meshHalfWay"+ std::to_string(m_intLabel) +".obj");
}

template<typename MeshType>
inline void DMB::MatrixMesh<MeshType>::addNonManifoldVertices()
{

    //TODO: do this only if there are any non manifold vertices!!

    using tComponent = std::vector<uint>;
    using tComponents = std::vector<tComponent>;
    //tmp prop
    std::vector<bool> faceVisited(m_nFaces, false);
    std::unordered_map<uint, tComponents> vhComponentsMap;

    for (uint vId = 0; vId < m_vertexFaceAdjacency.size(); ++vId)
    {
        if (!m_nonManifoldVertexProp[vId])
        {
            continue;
        }

        std::queue<uint> facesToVisit;
        tComponents components;

        for (auto vertexFace : m_vertexFaceAdjacency[vId])
        {
            if (faceVisited[vertexFace])
            {
                continue;
            }

            facesToVisit.push(vertexFace);

            tComponent component;

            while (!facesToVisit.empty())
            {
                auto topFace = facesToVisit.front();
                facesToVisit.pop();

                if (faceVisited[topFace])
                {
                    continue;
                }

                component.push_back(topFace);

                std::vector<uint> faceVertices;
                faceVertices.reserve(3);

                faceVertices.push_back(m_triangles[topFace * 3 + 0]);
                faceVertices.push_back(m_triangles[topFace * 3 + 1]);
                faceVertices.push_back(m_triangles[topFace * 3 + 2]);

                auto mid = std::find_if(faceVertices.begin(), faceVertices.end(),
                    [&](const auto& faceVertex)
                    {
                        return faceVertex == vId;
                    });

                //vector sorted in such a way, that vId vertex is first
                std::rotate(faceVertices.begin(), mid, faceVertices.end());

                int e0 = m_vertexEdgeMatrix.coeff(vId, faceVertices[1]);
                int e1 = m_vertexEdgeMatrix.coeff(vId, faceVertices[2]);

                //TODO: skip non-manifold edges?
                for (auto fId : m_edgeFaceAdjacency[e0])
                {
                    facesToVisit.push(fId);
                }

                for (auto fId : m_edgeFaceAdjacency[e1])
                {
                    facesToVisit.push(fId);
                }


                faceVisited[topFace] = true;
            }

            components.push_back(component);
        }

        //reset tag
        for (auto vertexFace : m_vertexFaceAdjacency[vId])
        {
            faceVisited[vertexFace] = false;
        }

        //sanity check if it is non-manifold vertex
        if (components.size() > 1)
        {
            vhComponentsMap[vId] = components;
        }
    }

    for (auto [vId, components] : vhComponentsMap)
    {
        //skip first component and replace the vertex in remaining components?
        for (int i = 1; i < components.size(); ++i)
        {
            uint newVId = duplicateVertex(vId);

            for (auto tId : components[i])
            {
                for (int j = 0; j < 3; ++j)
                {
                    if (m_triangles[tId * 3 + j] == vId)
                    {
                        //replace this vertex
                        m_triangles[tId * 3 + j] = newVId;
                    }
                }
            }
        }
    }

    updateMatrices();
    //TODO: no need to detect non manifolds?
    detectNonManifolds();
}

template<typename MeshType>
inline void DMB::MatrixMesh<MeshType>::resolveNonOrientableEdges()
{
#if 0
    for (auto eId = 0; eId < m_edgeFaceAdjacency.size(); ++eId)
    {
        if (!m_nonOrientableEdgeProp[eId])
        {
            continue;
        }

        //std::cout << "This non orientable edge has: " << m_edgeFaceAdjacency[eId].size() << " faces" << std::endl;

        for (auto tId : m_edgeFaceAdjacency[eId])
        {

            uint v0 = m_triangles[tId * 3 + 0];
            uint v1 = m_triangles[tId * 3 + 1];
            uint v2 = m_triangles[tId * 3 + 2];

            int e0 = m_vertexEdgeMatrix.coeff(v0, v1);
            int e1 = m_vertexEdgeMatrix.coeff(v1, v2);
            int e2 = m_vertexEdgeMatrix.coeff(v2, v0);

            int additionalNonOrientableEdges = 0;
            int additionalBoundaryEdges = 0;

            for (auto e : { e0,e1,e2 })
            {
                if (e == eId)
                {
                    continue;
                }

                if (m_nonOrientableEdgeProp[e])
                {
                    ++additionalNonOrientableEdges;
                }

                if (m_edgeFaceAdjacency[e].size() == 1)
                {
                    ++additionalBoundaryEdges;
                }

            }

            //std::cout << "This non orientable face has: " << additionalNonOrientableEdges << " additional non orientable edges" << std::endl;

            if (additionalNonOrientableEdges == 2 || (additionalNonOrientableEdges == 1 && additionalBoundaryEdges == 1))
            {
                //flipping this face will resolve the orientation
                m_triangles[tId * 3 + 1] = v2;
                m_triangles[tId * 3 + 2] = v1;

                m_faceFlipped[tId] = !m_faceFlipped[tId];
            }
        }
    }

    updateMatrices();

    std::cout << "after first round of non orientable edge resolution" << std::endl;
    detectNonOrientableEdges();
#endif

    //nothing more we can do - disconnect the faces
    for (auto eId = 0; eId < m_edgeFaceAdjacency.size(); ++eId)
    {
        if (!m_nonOrientableEdgeProp[eId])
        {
            continue;
        }

        //std::cout << "This non orientable edge has: " << m_edgeFaceAdjacency[eId].size() << " faces" << std::endl;

        //get edge vertices of this edge
        auto edgeVertices = m_edgeVertices[eId];

        auto v0 = edgeVertices.first;
        auto v1 = edgeVertices.second;

        //for (auto tId : m_edgeFaceAdjacency[eId])
        for (int i = 1 /*skip first face*/; i < m_edgeFaceAdjacency[eId].size(); ++i)
        {
            uint tId = m_edgeFaceAdjacency[eId][i];

            uint v0New = duplicateVertex(v0);
            uint v1New = duplicateVertex(v1);

            for (int j = 0; j < 3; ++j)
            {
                uint v = m_triangles[tId * 3 + j];

                if (v == v0)
                {
                    m_triangles[tId * 3 + j] = v0New;
                }

                if (v == v1)
                {
                    m_triangles[tId * 3 + j] = v1New;
                }

            }
        }
    }

    updateMatrices();
    detectNonOrientableEdges();
}


template<typename MeshType>
inline void DMB::MatrixMesh<MeshType>::markSelfIntersectingFaces()
{
    //TODO: remove this call...
    detectNonManifolds();
    m_selfIntersectingFace.clear();
    m_selfIntersectingFace.resize(m_nFaces, false);

    for (uint i = 0; i < std::size(m_triangles); i += 3)
    {
        uint tId = i / 3;

        m_selfIntersectingFace[tId] = m_nonManifoldEdgeFaceProp[tId];
    }
}

template<typename MeshType>
inline void DMB::MatrixMesh<MeshType>::resolveConflictingComponents(MeshArrangement<MeshType>& ma)
{
    //TODO: make this a method of the class
    auto getConnectedComponents = [](MeshType& mesh)
        {
            using tComponent = std::vector<OpenMesh::SmartFaceHandle>;
            using tComponents = std::vector<tComponent>;

            tComponents components;

            OpenMesh::FProp<bool> visited(false, mesh);

            for (auto fh : mesh.faces())
            {
                if (!fh.is_valid() || mesh.status(fh).deleted() || visited[fh])
                {
                    continue;
                }

                tComponent component;

                //It is a face, which is not visited => new component
                std::queue<OpenMesh::SmartFaceHandle> facesToVisit;
                facesToVisit.push(fh);

                while (!facesToVisit.empty())
                {
                    auto topFace = facesToVisit.front();
                    facesToVisit.pop();

                    if (!topFace.is_valid() || mesh.status(topFace).deleted() || visited[topFace])
                    {
                        continue;
                    }

                    component.push_back(topFace);
                    visited[topFace] = true;

                    for (auto ffh : topFace.faces())
                    {
                        facesToVisit.push(ffh);
                    }
                }

                components.push_back(component);

            }

            return components;
        };

    OpenMesh::FProp<uint> pFhToMaFh(m_mesh, m_pFhToMaFh);

    auto flipComponent = [&](MeshType& mesh, const std::vector<OpenMesh::SmartFaceHandle>& component)
        {

            //flip this component in matrix rep., in halfedge rep. and in mesh arrangement as well
            std::vector<std::vector<OpenMesh::SmartVertexHandle>> newFaces;
            //std::map<tFaceHandle, uint> fhToTIdMap;
            vector<uint> fhToTIdMap;

            for (auto fh : component)
            {

                //flip in matrix rep.
                uint tId = pFhToMaFh[fh];
                std::array<uint, 3> triangleVertices;
                triangleVertices[0] = m_triangles[tId * 3 + 0];
                triangleVertices[1] = m_triangles[tId * 3 + 1];
                triangleVertices[2] = m_triangles[tId * 3 + 2];
                m_triangles[tId * 3 + 1] = triangleVertices[2];
                m_triangles[tId * 3 + 2] = triangleVertices[1];

                //flip in ma
                uint maTId = m_tIdToOriginalTId[tId];
                std::array<uint, 3> maTriangleVertices;
                maTriangleVertices[0] = ma.m_triangles[maTId * 3 + 0];
                maTriangleVertices[1] = ma.m_triangles[maTId * 3 + 1];
                maTriangleVertices[2] = ma.m_triangles[maTId * 3 + 2];
                ma.m_triangles[maTId * 3 + 1] = maTriangleVertices[2];
                ma.m_triangles[maTId * 3 + 2] = maTriangleVertices[1];


                //prepare for halfedge flip
                //fhToTIdMap[fh] = tId;
                auto faceVertices = fh.vertices().to_vector();
                newFaces.push_back(faceVertices);
                fhToTIdMap.push_back(tId);
                mesh.delete_face(fh, false);
            }

            //flip in halfedge mesh
            int j = 0;
            for (auto faceVertices : newFaces)
            {
                auto newFh = mesh.add_face(faceVertices[0], faceVertices[2], faceVertices[1]);
                pFhToMaFh[newFh] = fhToTIdMap[j];

                ++j;
            }
        };

    auto components = getConnectedComponents(m_mesh);

    std::unordered_set<int> componentsToCheck;

    for (int i = 0; i < components.size(); ++i)
    {
        const auto& component = components[i];
        for (auto fh : component)
        {
            if (ma.m_checkOrientation[m_tIdToOriginalTId[pFhToMaFh[fh]]])
            {
                componentsToCheck.insert(i);
            }
        }
    }

    for (int componentId : componentsToCheck)
    {
        const auto& component = components[componentId];
        if (calcVolumeSignExact(component) < 0)
        {
            flipComponent(m_mesh, component);
        }
    }

    updateMatrices();
}

template<typename MeshType>
inline void DMB::MatrixMesh<MeshType>::computeComponentsVolume(MeshArrangement<MeshType>& ma)
{
    //TODO: make this a method of the class
    auto getConnectedComponents = [](MeshType& mesh)
        {
            using tComponent = std::vector<OpenMesh::SmartFaceHandle>;
            using tComponents = std::vector<tComponent>;

            tComponents components;

            OpenMesh::FProp<bool> visited(false, mesh);

            for (auto fh : mesh.faces())
            {
                if (!fh.is_valid() || mesh.status(fh).deleted() || visited[fh])
                {
                    continue;
                }

                tComponent component;

                //It is a face, which is not visited => new component
                std::queue<OpenMesh::SmartFaceHandle> facesToVisit;
                facesToVisit.push(fh);

                while (!facesToVisit.empty())
                {
                    auto topFace = facesToVisit.front();
                    facesToVisit.pop();

                    if (!topFace.is_valid() || mesh.status(topFace).deleted() || visited[topFace])
                    {
                        continue;
                    }

                    component.push_back(topFace);
                    visited[topFace] = true;

                    for (auto ffh : topFace.faces())
                    {
                        facesToVisit.push(ffh);
                    }
                }

                components.push_back(component);

            }

            return components;
        };

    auto components = getConnectedComponents(m_mesh);

    OpenMesh::FProp<uint> pFhToMaFh(m_mesh, m_pFhToMaFh);
    
    for (auto component : components)
    {
        for (auto fh : component)
        {
            ma.m_faceComponentId[m_tIdToOriginalTId[pFhToMaFh[fh]]] = ma.m_componentsVolume.size();
        }

        int sign = calcVolumeSignExact(component);
        double volume = calcSignedVolume(component);
        double volumeSize = std::abs(volume);
        ma.m_componentsVolume.push_back(volumeSize * sign);
    }
}

template<typename MeshType>
inline uint DMB::MatrixMesh<MeshType>::duplicateVertex(uint vId)
{
    uint newVId = (uint)m_coordinatesImplicit.size();

    auto* gp = m_coordinatesImplicit[vId];

    m_coordinatesImplicit.push_back(gp);
    m_coordinates.push_back(m_coordinates[vId * 3 + 0]);
    m_coordinates.push_back(m_coordinates[vId * 3 + 1]);
    m_coordinates.push_back(m_coordinates[vId * 3 + 2]);

    tVertexHandle invalidVh;
    invalidVh.invalidate();
    m_matrixVhToOMVh.push_back(invalidVh);
    //TODO: make some property management!
    m_nonManifoldEdgeVertexProp.push_back(false);
    m_nonManifoldVertexProp.push_back(false);

    ++m_nVertices;

    return newVId;
}

template<typename MeshType>
inline void DMB::MatrixMesh<MeshType>::detectNonManifolds()
{
    m_nonManifoldVertexProp.clear();
    m_nonManifoldEdgeProp.clear();
    m_nonManifoldEdgeVertexProp.clear();
    m_nonManifoldFaceProp.clear();
    m_nonManifoldEdgeFaceProp.clear();

    //TODO: THIS IS MADNESS!! Make some reasonable property manager
    m_nonManifoldVertexProp.resize(m_nVertices, false);
    m_nonManifoldEdgeProp.resize(m_nEdges + 1, false);
    m_nonManifoldEdgeVertexProp.resize(m_nVertices, false);
    m_nonManifoldFaceProp.resize(m_nFaces, false);
    m_nonManifoldEdgeFaceProp.resize(m_nFaces, false);

    detectNonManifoldEdges();
    detectNonManifoldVertices();
    detectNonManifoldFaces();
}

template<typename MeshType>
inline bool DMB::MatrixMesh<MeshType>::buildManifoldMesh(bool considerOrigin/* = false*/)
{
    detectNonManifolds();

    detectNonManifolds();

    //TODO: THIS IS MADNESS!! Make some reasonable property manager
    tVertexHandle invalidVh;
    invalidVh.invalidate();
    m_matrixVhToOMVh.clear();
    m_matrixVhToOMVh.resize(m_nVertices, invalidVh);
    m_faceAdded.clear();
    m_faceAdded.resize(m_nFaces, false);
    //m_faceFlipped.clear();
    //m_faceFlipped.resize(m_nFaces, false);

    OpenMesh::FProp<uint> pFhToMaFh(0, m_mesh, "bo_fhToMaFh");
    m_pFhToMaFh = pFhToMaFh.getRawProperty();
    OpenMesh::VProp<uint> pVhToMaVId(0, m_mesh, "bo_operandToMeshArrangementIndex"); //TODO: rename, use header with defines for the names of properties
    m_pVhToMaVId = pVhToMaVId.getRawProperty();

    //addManifoldFaces();

    detectNonOrientableEdges();

    addNonManifoldEdges(considerOrigin);
    

    //TODO: WHY ARE NOT NON MANIFOLD VVERTICES RESOLVED???
    while (m_nNonManifoldVertices != 0)
    {
        addNonManifoldVertices();
    }

    detectNonOrientableEdges();
    resolveNonOrientableEdges();

    detectNonOrientableEdges();

    //TODO: MOVE THIS!!!
    bool success = true;
    auto addFace = [this, &pFhToMaFh, &pVhToMaVId, &success](uint tId) -> OpenMesh::SmartFaceHandle
    {

        uint v0 = m_triangles[tId * 3 + 0];
        uint v1 = m_triangles[tId * 3 + 1];
        uint v2 = m_triangles[tId * 3 + 2];

        tTriangle triangle =
        {
            tPoint { m_coordinates[v0 * 3 + 0], m_coordinates[v0 * 3 + 1], m_coordinates[v0 * 3 + 2] },
            tPoint { m_coordinates[v1 * 3 + 0], m_coordinates[v1 * 3 + 1], m_coordinates[v1 * 3 + 2] },
            tPoint { m_coordinates[v2 * 3 + 0], m_coordinates[v2 * 3 + 1], m_coordinates[v2 * 3 + 2] },
        };

        //add vertices, if not already present in the mesh
        int j = 0;
        for (auto v : { v0, v1, v2 })
        {
            if (!m_matrixVhToOMVh[v].is_valid())
            {
                auto newVh = m_mesh.add_vertex(triangle[j]);
                m_matrixVhToOMVh[v] = newVh;
                pVhToMaVId[newVh] = v;
            }
            ++j;
        }

        auto vh0 = m_matrixVhToOMVh[v0];
        auto vh1 = m_matrixVhToOMVh[v1];
        auto vh2 = m_matrixVhToOMVh[v2];

        auto fh = m_mesh.add_face(vh0, vh1, vh2);

        if (!fh.is_valid())
        {
            //std::cout << "Could not add face" << std::endl;
            success = false;
            return fh;
        }

        pFhToMaFh[fh] = tId;
        m_faceAdded[tId] = true;

        return fh;
    };

    for (uint i = 0; i < std::size(m_triangles); i += 3)
    {
        uint tId = i / 3;

        if (m_faceAdded[tId])
        {
            continue;
        }

        uint v0 = m_triangles[i + 0];
        uint v1 = m_triangles[i + 1];
        uint v2 = m_triangles[i + 2];

        addFace(tId);
    }

    return success;
}

template<typename MeshType>
inline bool DMB::MatrixMesh<MeshType>::resolveNmfVerticesInHalfedgeMesh()
{
    std::vector<OpenMesh::SmartVertexHandle>verticesToResolve;
    for (auto vh : m_mesh.vertices())
    {
        if (!m_mesh.is_manifold(vh))
        {
            verticesToResolve.push_back(vh);
        }
    }

    if (verticesToResolve.empty())
    {
        return true;
    }

    for (auto vh : verticesToResolve)
    {
        std::vector<std::vector<OpenMesh::SmartHalfedgeHandle>> outHeGroups;

        std::vector<OpenMesh::SmartHalfedgeHandle> group;
        for (auto he : vh.outgoing_halfedges())
        {
            if (he.is_boundary())
            {
                if (!group.empty())
                {
                    outHeGroups.push_back(group);
                }
                group.clear();
            }

            group.push_back(he);
        }
        outHeGroups.push_back(group);


        if (outHeGroups.size() < 2)
        {
            continue;
        }

        for (int i = 1; i < outHeGroups.size(); ++i)
        {
            auto newVh = m_mesh.add_vertex(m_mesh.point(vh));

            std::vector<OpenMesh::SmartFaceHandle> facesToResolve;

            for (auto he : outHeGroups[i])
            {
                if (!he.is_boundary())
                {
                    facesToResolve.push_back(he.face());
                }

            }

            std::vector<std::vector<OpenMesh::SmartVertexHandle>> newFaces;

            for (auto fh : facesToResolve)
            {
                std::vector<OpenMesh::SmartVertexHandle> newVertices;
                for (auto vvh : fh.vertices())
                {
                    if (vvh == vh)
                    {
                        newVertices.push_back(newVh);
                    }
                    else
                    {
                        newVertices.push_back(vvh);
                    }
                }


                newFaces.push_back(newVertices);
            }

            for (auto fh : facesToResolve)
            {
                m_mesh.delete_face(fh, false);
            }

            for (auto newFace : newFaces)
            {
                if (!m_mesh.add_face(newFace).is_valid())
                {
                    return false;
                }
            }
        }
    }

    return true;
}

template<typename MeshType>
inline uint DMB::MatrixMesh<MeshType>::getFaceEdgeId(uint tId, int i)
{
    Eigen::Index v0 = m_triangles[tId * 3 + 0];
    Eigen::Index v1 = m_triangles[tId * 3 + 1];
    Eigen::Index v2 = m_triangles[tId * 3 + 2];

    std::array<uint, 3> edges;

    edges[0] = m_vertexEdgeMatrix.coeff(v0, v1);
    edges[1] = m_vertexEdgeMatrix.coeff(v1, v2);
    edges[2] = m_vertexEdgeMatrix.coeff(v2, v0);

    return edges[i];
}

template<typename MeshType>
inline void DMB::MatrixMesh<MeshType>::buildOperand(const MeshArrangement<MeshType>& ma)
{
    std::vector<int> maToOperandVertices;
    maToOperandVertices.resize(ma.m_coordinatesImplicit.size(), -1);

    uint newVhId = 0;

    for (int i = 0; i < std::size(ma.m_triangles); i += 3)
    {
        uint tId = i / 3;
        auto faceLabel = ma.m_labels[tId];

        bool belongsToThisOperand = (faceLabel & m_label) != 0;

        if (!belongsToThisOperand)
        {
            continue;
        }

        std::vector<tVertexHandle> vertices;

        for (int j = 0; j < 3; ++j)
        {
            uint maVh = ma.m_triangles[i + j];

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


        //TODO: THIS IS A MESS...
        m_faceOriginLabel.push_back( ma.m_labels[tId]);
    }

}

template<typename MeshType>
inline void DMB::MatrixMesh<MeshType>::buildDebugMesh()
{

    transferArrangementPropertiesToMesh();
    detectBoundaries();
    handleCoplanarFaces();

    //std::string path = "C:/T3D/Samples/Booleans/";
    std::string path = "C:/skola/PhD/Samples/booleans/components/";

#if 0

    MeshType debugMesh;

    for (int i = 0; i < std::size(m_triangles); i += 3)
    {
        uint tId = i / 3;
        uint maVh0 = m_triangles[i + 0];
        uint maVh1 = m_triangles[i + 1];
        uint maVh2 = m_triangles[i + 2];

        tTriangle triangle =
        {
            tPoint { m_coordinates[maVh0 * 3 + 0], m_coordinates[maVh0 * 3 + 1], m_coordinates[maVh0 * 3 + 2] },
            tPoint { m_coordinates[maVh1 * 3 + 0], m_coordinates[maVh1 * 3 + 1], m_coordinates[maVh1 * 3 + 2] },
            tPoint { m_coordinates[maVh2 * 3 + 0], m_coordinates[maVh2 * 3 + 1], m_coordinates[maVh2 * 3 + 2] },
        };

        auto vh0 = debugMesh.add_vertex(triangle[0]);
        auto vh1 = debugMesh.add_vertex(triangle[1]);
        auto vh2 = debugMesh.add_vertex(triangle[2]);

        auto fh = debugMesh.add_face(vh0, vh1, vh2);
        debugMesh.set_color(fh, { 192,192,192 });

        
        
    }

    OpenMesh::IO::Options opt = OpenMesh::IO::Options::Default;
    opt += OpenMesh::IO::Options::FaceColor;
    DMB::saveMesh(debugMesh, path + "debugMAMesh"+std::to_string(m_intLabel)+".stl");
#endif


#if 1

    {
        OpenMesh::FProp< bool > intersectionFace(m_mesh, m_pIntersectionFace);

        for (auto fh : m_mesh.faces())
        {
            m_mesh.set_color(fh, { 162,162,162 });

            if (intersectionFace[fh])
            {
                m_mesh.set_color(fh, { 0,0,250 });
            }
        }

        OpenMesh::IO::Options opt = OpenMesh::IO::Options::Default;
        opt += OpenMesh::IO::Options::FaceColor;
        //DMB::saveMesh(m_mesh, path + "intersecionFaces_operand_" + std::to_string(m_intLabel) + ".ply", opt);
    }

    {
        OpenMesh::FProp<bool> coplanarFace(m_mesh, m_pCoplanarFace);

        for (auto fh : m_mesh.faces())
        {
            m_mesh.set_color(fh, { 162,162,162 });

            if (coplanarFace[fh])
            {
                m_mesh.set_color(fh, { 255,0,250 });
            }
        }

        OpenMesh::IO::Options opt = OpenMesh::IO::Options::Default;
        opt += OpenMesh::IO::Options::FaceColor;
        //DMB::saveMesh(m_mesh, path + "coplanarFaces_operand_" + std::to_string(m_intLabel) + ".ply", opt);
    }

    {
        OpenMesh::FProp< bool > intersectionFace(m_mesh, m_pIntersectionFace);
        OpenMesh::FProp< std::bitset<NBIT>> labeling(m_mesh, m_pFaceIOLabeling);

        for (auto fh : m_mesh.faces())
        {
            m_mesh.set_color(fh, { 162,162,162 });

            if (intersectionFace[fh])
            {
                typename MeshType::Color c = labeling[fh].count() > 0 ? MeshType::Color(255,0,0) : MeshType::Color(0, 255, 0);
                m_mesh.set_color(fh, c);

            }
        }

        OpenMesh::IO::Options opt = OpenMesh::IO::Options::Default;
        opt += OpenMesh::IO::Options::FaceColor;
       // DMB::saveMesh(m_mesh, path + "IOFaces_operand_" + std::to_string(m_intLabel) + ".ply", opt);
    }

    {
        auto getConnectedComponents = [](MeshType& mesh)
            {
                using tComponent = std::vector<OpenMesh::FaceHandle>;
                using tComponents = std::vector<tComponent>;

                tComponents components;

                OpenMesh::FProp<bool> visited(false, mesh);

                for (auto fh : mesh.faces())
                {
                    if (!fh.is_valid() || mesh.status(fh).deleted() || visited[fh])
                    {
                        continue;
                    }

                    tComponent component;

                    //It is a face, which is not visited => new component
                    std::queue<OpenMesh::SmartFaceHandle> facesToVisit;
                    facesToVisit.push(fh);

                    while (!facesToVisit.empty())
                    {
                        auto topFace = facesToVisit.front();
                        facesToVisit.pop();

                        if (!topFace.is_valid() || mesh.status(topFace).deleted() || visited[topFace])
                        {
                            continue;
                        }

                        component.push_back(topFace);
                        visited[topFace] = true;

                        for (auto ffh : topFace.faces())
                        {
                            facesToVisit.push(ffh);
                        }
                    }

                    components.push_back(component);

                }

                return components;
            };

        auto components = getConnectedComponents(m_mesh);

        static int cc = 0;

        for (const auto& component : components)
        {
            MeshType meshPart;
            //DMB::copyMeshPart<MeshType>(m_mesh, meshPart, component);

            OpenMesh::IO::Options opt = OpenMesh::IO::Options::Default;
            //DMB::saveMesh(meshPart, "C:/skola/PhD/Samples/booleans/components/cmp_" + std::to_string(m_intLabel) + "_" + std::to_string(cc) + ".obj", opt);
            ++cc;
        }
    }

#endif

#if 0
    //std::string path = "D:/skola/PhD/VUT/tezy/images and models/bool principle/";
    //std::string path = "C:/skola/PhD/Samples/booleans/components/";

    //intersection lines
    {
        OpenMesh::EProp<bool> pIntersectionEdge(m_mesh, m_pIntersectionEdge);

        OpenMesh::VProp<bool> pIntersectionVertex(false, m_mesh);

        for (auto eh : m_mesh.edges())
        {
            if (pIntersectionEdge[eh])
            {
                auto v0 = eh.v0();
                auto v1 = eh.v1();
                pIntersectionVertex[v0] = true;
                pIntersectionVertex[v1] = true;
            }
        }

        std::map<int, int> meshToObjVertexMap;
        int objVertexId = 1;
        std::ofstream lineObj(path + "intEdges" + std::to_string(m_intLabel) + ".obj");

        for (auto vh : m_mesh.vertices())
        {
            if (pIntersectionVertex[vh])
            {
                meshToObjVertexMap[vh.idx()] = objVertexId;

                auto p = m_mesh.point(vh);
                lineObj << "v " << std::to_string(p[0]) << " " << std::to_string(p[1]) << " " << std::to_string(p[2]) << std::endl;

                ++objVertexId;
            }
        }

        for (auto eh : m_mesh.edges())
        {
            if (pIntersectionEdge[eh])
            {
                auto v0 = eh.v0();
                auto v1 = eh.v1();
                lineObj << "l " << std::to_string(meshToObjVertexMap[v0.idx()]) << " " << std::to_string(meshToObjVertexMap[v1.idx()]) << std::endl;
            }
        }
    }

    //intersection faces
    {
        OpenMesh::FProp< std::bitset<NBIT>> labeling(m_mesh, m_pFaceIOLabeling);
        OpenMesh::FProp< bool > intersectionFace(m_mesh, m_pIntersectionFace);

        MeshType insideMesh;
        MeshType outsideMesh;

        for(auto fh : m_mesh.faces())
        {
            if (!intersectionFace[fh])
            {
                continue;
            }

            auto l = labeling[fh];
            l[NBIT - 1] = 0;
            MeshType& mesh = (l.count() > 0 ? insideMesh : outsideMesh);

            auto vertices = fh.vertices().to_vector();
            auto vh0 = mesh.add_vertex(m_mesh.point(vertices[0]));
            auto vh1 = mesh.add_vertex(m_mesh.point(vertices[1]));
            auto vh2 = mesh.add_vertex(m_mesh.point(vertices[2]));

            mesh.add_face(vh0, vh1, vh2);
        }
        insideMesh.update_normals();
        outsideMesh.update_normals();

        DMB::saveMesh(insideMesh, path +"meshOperand" + std::to_string(m_intLabel) + "IntFacesInside.stl");
        DMB::saveMesh(outsideMesh, path +"meshOperand" + std::to_string(m_intLabel) + "IntFacesOutside.stl");
    }

    //inside outside components
    {
        OpenMesh::FProp< std::bitset<NBIT>> labeling(m_mesh, m_pFaceIOLabeling);

        MeshType insideMesh;
        MeshType outsideMesh;

        for (auto fh : m_mesh.faces())
        {
            

            MeshType& mesh = (labeling[fh].count() > 0 ? insideMesh : outsideMesh);

            auto vertices = fh.vertices().to_vector();
            auto vh0 = mesh.add_vertex(m_mesh.point(vertices[0]));
            auto vh1 = mesh.add_vertex(m_mesh.point(vertices[1]));
            auto vh2 = mesh.add_vertex(m_mesh.point(vertices[2]));

            mesh.add_face(vh0, vh1, vh2);
        }

        insideMesh.update_normals();
        outsideMesh.update_normals();
        DMB::saveMesh(insideMesh, path + "insideMesh" + std::to_string(m_intLabel) + ".stl");
        DMB::saveMesh(outsideMesh, path + "outsideMesh" + std::to_string(m_intLabel) + ".stl");
    }

    //whole mesh
    m_mesh.update_normals();
    DMB::saveMesh(m_mesh, path + "meshOperand" + std::to_string(m_intLabel) + ".obj");
#endif
}

template<typename MeshType>
inline void DMB::MatrixMesh<MeshType>::flipToConsistentOrientation(MeshArrangement<MeshType>& ma)
{
    
    OpenMesh::FProp<uint> pFhToMaFh(m_mesh, m_pFhToMaFh);

    auto getConnectedComponents = [](MeshType& mesh)
    {
        using tComponent = std::vector<OpenMesh::SmartFaceHandle>;
        using tComponents = std::vector<tComponent>;

        tComponents components;

        OpenMesh::FProp<bool> visited(false, mesh);

        for (auto fh : mesh.faces())
        {
            if (!fh.is_valid() || mesh.status(fh).deleted() || visited[fh])
            {
                continue;
            }

            tComponent component;

            //It is a face, which is not visited => new component
            std::queue<OpenMesh::SmartFaceHandle> facesToVisit;
            facesToVisit.push(fh);

            while (!facesToVisit.empty())
            {
                auto topFace = facesToVisit.front();
                facesToVisit.pop();

                if (!topFace.is_valid() || mesh.status(topFace).deleted() || visited[topFace])
                {
                    continue;
                }

                component.push_back(topFace);
                visited[topFace] = true;

                for (auto ffh : topFace.faces())
                {
                    facesToVisit.push(ffh);
                }
            }

            components.push_back(component);

        }

        return components;
    };

    auto areNormalsFlipped = [this](std::vector<OpenMesh::SmartFaceHandle> component)
    {
        std::vector<tFaceHandle> tmpComponent;
        tmpComponent.resize(component.size());

        for (int i = 0; i < component.size(); ++i)
        {
            tmpComponent[i] = component[i];
        }

        return calcSignedVolume(tmpComponent) <= 0;
    };

    auto components = getConnectedComponents(m_mesh);

    //for (auto fh : m_mesh.faces())
    //{
    //    m_mesh.set_color(fh, { 192,192,192 });
    //}

    //needed for areNormalsFlipped check - it uses mean scalar curvature
    m_mesh.update_normals();

    for (auto component : components)
    {
        //is this self intersecting component?
        bool isSelfIntersecting = false;
        for(auto fh : component)
        {
            if (m_selfIntersectingFace[pFhToMaFh[fh]])
            {
                isSelfIntersecting = true;
            }
        }

        if (!isSelfIntersecting)
        {
            continue;
        }

        //get component vertices
        std::vector<tVertexHandle> componentVertices;
        for (auto fh : component)
        {
            for (auto vh : fh.vertices())
            {
                componentVertices.push_back(vh);
            }
        }

        DMB::makeUniqueVector(componentVertices);

        if (areNormalsFlipped(component))
        {
            //flip this component in matrix rep., in halfedge rep. and in mesh arrangement as well
            std::vector<std::vector<OpenMesh::SmartVertexHandle>> newFaces;
            //std::map<tFaceHandle, uint> fhToTIdMap;
            vector<uint> fhToTIdMap;

            for (auto fh : component)
            {

                //flip in matrix rep.
                uint tId = pFhToMaFh[fh];
                std::array<uint, 3> triangleVertices;
                triangleVertices[0] = m_triangles[tId * 3 + 0];
                triangleVertices[1] = m_triangles[tId * 3 + 1];
                triangleVertices[2] = m_triangles[tId * 3 + 2];
                m_triangles[tId * 3 + 1] = triangleVertices[2];
                m_triangles[tId * 3 + 2] = triangleVertices[1];

                //flip in ma
                uint maTId = m_tIdToOriginalTId[tId];
                std::array<uint, 3> maTriangleVertices;
                maTriangleVertices[0] = ma.m_triangles[maTId * 3 + 0];
                maTriangleVertices[1] = ma.m_triangles[maTId * 3 + 1];
                maTriangleVertices[2] = ma.m_triangles[maTId * 3 + 2];
                ma.m_triangles[maTId * 3 + 1] = maTriangleVertices[2];
                ma.m_triangles[maTId * 3 + 2] = maTriangleVertices[1];


                //prepare for halfedge flip
                //fhToTIdMap[fh] = tId;
                auto faceVertices = fh.vertices().to_vector();
                newFaces.push_back(faceVertices);
                fhToTIdMap.push_back(tId);
                m_mesh.delete_face(fh, false);
            }

            //flip in halfedge mesh
            int j = 0;
            for (auto faceVertices : newFaces)
            {
                auto newFh = m_mesh.add_face(faceVertices[0], faceVertices[2], faceVertices[1]);
                pFhToMaFh[newFh] = fhToTIdMap[j];

                ++j;
            }
        }

    }

    updateMatrices();
}

template<typename MeshType>
inline void DMB::MatrixMesh<MeshType>::adjustOrientation(MeshArrangement<MeshType>& ma)
{
    for (int i = 0; i < std::size(m_triangles); i += 3)
    {
        uint tId = i / 3;
        if (m_faceFlipped[tId])
        {
            //flip in ma
            uint maTId = m_tIdToOriginalTId[tId];
            std::array<uint, 3> maTriangleVertices;
            maTriangleVertices[0] = ma.m_triangles[maTId * 3 + 0];
            maTriangleVertices[1] = ma.m_triangles[maTId * 3 + 1];
            maTriangleVertices[2] = ma.m_triangles[maTId * 3 + 2];
            ma.m_triangles[maTId * 3 + 1] = maTriangleVertices[2];
            ma.m_triangles[maTId * 3 + 2] = maTriangleVertices[1];
        }
    }
    
}

template<typename MeshType>
inline bool DMB::MatrixMesh<MeshType>::disconnectComponents(MeshArrangement<MeshType>& ma)
{
    transferArrangementPropertiesToMesh();
    detectBoundaries();
    handleCoplanarFaces();

    OpenMesh::EProp<bool> pIntersectionEdge(m_mesh, m_pIntersectionEdge);
    //OpenMesh::VProp<uint> pOperandToMeshArrangementIndex(m_mesh, m_pOperandToMeshArrangementIndex);
    OpenMesh::EProp<bool> pCoplanarEdge(m_mesh, m_pCoplanarEdge);
    OpenMesh::VProp<bool> coplanarVertex(m_mesh, m_pCoplanarVertex);

    OpenMesh::FProp< std::bitset<NBIT>> labeling(m_mesh, m_pFaceIOLabeling);
    OpenMesh::FProp< bool > intersectionFace(m_mesh, m_pIntersectionFace);

    OpenMesh::FProp<uint> pFhToMaFh(m_mesh, m_pFhToMaFh);
    OpenMesh::FProp<bool> coplanarFace(m_mesh, m_pCoplanarFace);


    OpenMesh::VProp<bool> intersectionVertex(false, m_mesh);

    //mark all intersection vertices
    for (auto eh : m_mesh.edges())
    {
        if (pIntersectionEdge[eh] && !pCoplanarEdge[eh])
        {
            intersectionVertex[eh.v0()] = true;
            intersectionVertex[eh.v1()] = true;
        }
    }

    OpenMesh::VProp<int> intersectionValance(0, m_mesh);

    for (auto vh : m_mesh.vertices())
    {
        if (intersectionVertex[vh])
        {
            int val = 0;

            for (auto outHe : vh.outgoing_halfedges())
            {
                if (pIntersectionEdge[outHe.edge()] || coplanarVertex[outHe.to()] || outHe.to().is_boundary())
                {
                    ++val;
                }
            }

            intersectionValance[vh] = val;
        }
    }

    //reset tagged - just to be sure
    for (auto fh : m_mesh.faces())
    {
        m_mesh.status(fh).set_tagged(false);
    }

    std::vector<std::vector<tFaceHandle>> components;

    //find components
    for (auto fh : m_mesh.faces())
    {
        if (fh.deleted())
        {
            continue;
        }

        if (coplanarFace[fh])
        {
            continue;
        }

        //found a new component seed face
        if (!fh.tagged())
        {
            std::queue<tFaceHandle> facesToVisit;
            std::vector<tFaceHandle> componentFaces;

            facesToVisit.push(fh);

            while (!facesToVisit.empty())
            {
                //get top face
                auto topFace = OpenMesh::make_smart(facesToVisit.front(), m_mesh);
                facesToVisit.pop();

                if (topFace.tagged() || topFace.deleted() || coplanarFace[topFace])
                {
                    continue;
                }


                //store new component face
                componentFaces.push_back(topFace);

                //continue with adjacent face, but adjacencies through intersection edges are not allowed
                for (auto he : topFace.halfedges())
                {

                    if (pIntersectionEdge[he.edge()] || pCoplanarEdge[he.edge()])
                    {
                        continue;
                    }

                    auto adjFace = he.opp().face();

                    if (adjFace.is_valid() && !adjFace.tagged())
                    {
                        facesToVisit.push(adjFace);
                    }

                }

                //mark this face as processed
                m_mesh.status(topFace).set_tagged(true);
            }

            //store this component
            components.push_back(componentFaces);
        }
    }

#if 0
    {
        int cc = 0;

        for (const auto& component : components)
        {
            MeshType meshPart;
            DMB::copyMeshPart<MeshType>(m_mesh, meshPart, component);

            OpenMesh::IO::Options opt = OpenMesh::IO::Options::Default;
            DMB::saveMesh(meshPart, "C:/skola/PhD/Samples/booleans/components/cut_cmp_" + std::to_string(m_intLabel) + "_" + std::to_string(cc) + ".obj", opt);
            ++cc;
        }
    }
#endif

    std::vector<uint> notLabeledComponents;
    std::vector<std::bitset<NBIT>> componentsLabels;

    //now spread the labeling information over the connected components
    //for (const auto& component : components)
    for (uint  componentId = 0; componentId < components.size(); ++componentId)
    {
        const auto& component = components[componentId];

        //find seed label
        std::bitset<NBIT> seedLabel;

        bool conflictingLabeling = false;
        bool unclosedCurve = false;
        bool componentLabeled = false;

        double largestComponentVolume = std::numeric_limits<double>::lowest();
        double largestComponentId = -1;


        OpenMesh::FaceHandle seedFh;
        for (auto fh : component)
        {
            if (intersectionFace[fh])
            {
                if (componentLabeled && seedLabel != labeling[fh])
                {
                    //std::cout << "conflicting labeling: " << seedLabel << " | " << labeling[fh] << std::endl;
                    //std::cout << "coplanar faces: " << (coplanarFace[seedFh] ? "true" : "false")  << " | " << (coplanarFace[fh] ? "true" : "false") << std::endl;
                    conflictingLabeling = true;
                }

                if (!componentLabeled)
                {
                    seedFh = fh;
                    seedLabel = labeling[fh];
                    componentLabeled = true;
                }

                for (auto vh : OpenMesh::make_smart(fh, m_mesh).vertices())
                {
                    if (intersectionValance[vh] == 1)
                    {
                        unclosedCurve = true;
                    }
                }
            }
        }


        //TODO: label somehow that it wasn't labeled? DOES NOT CONTAIN INTERSECTIONS!
        if (!componentLabeled)
        {
            std::bitset<NBIT> notLabeledLabel;
            notLabeledLabel[NBIT - 2] = 1;

            for (auto cfh : component)
            {
                labeling[cfh] = labeling[cfh] | notLabeledLabel;
            }

            notLabeledComponents.push_back(componentId);
            componentsLabels.push_back(notLabeledLabel);

            continue;
        }


        if (conflictingLabeling)
        {
            //volume based voting
            std::unordered_map<std::bitset<NBIT>, std::unordered_set<int>> labelToLabelingComponentMap;

            for (auto fh : component)
            {
                if (intersectionFace[fh])
                {
                    int componentId = ma.m_faceClassifiedByComponent[m_tIdToOriginalTId[pFhToMaFh[fh]]];
                    auto l = labeling[fh];

                    labelToLabelingComponentMap[l].insert(componentId);
                }
            }

            std::unordered_map < std::bitset<NBIT>, double > labelToLabelingVolumeMap;

            for (auto [label, componentsId] : labelToLabelingComponentMap)
            {
                //TODO: should be bigfloat
                double v = 0;
                for (auto id : componentsId)
                {
                    v += ma.m_componentsVolume[id];
                }

                labelToLabelingVolumeMap[label] = v;
            }

            double maxVolume = std::numeric_limits<double>::lowest();

            //assume conflicting labeling will be resolved
            conflictingLabeling = false;

            for (auto [key, val] : labelToLabelingVolumeMap)
            {
                if (val > maxVolume)
                {
                    maxVolume = val;
                    seedLabel = key;
                    conflictingLabeling = false;
                }

                if (val == maxVolume && seedLabel != key)
                {
                    //still conflicting => e-surface
                    //conflictingLabeling = true;
                    //TODO: I should put 0 in place of the input mesh, that labeled this face/component
                    seedLabel = 0;
                }
            }

        }

        //set labeling
        if (!conflictingLabeling)
        {
            for (auto fh : component)
            {
                labeling[fh] = seedLabel;
            }

            componentsLabels.push_back(seedLabel);
        }
    }

    //reset tagged
    for (auto fh : m_mesh.faces())
    {
        m_mesh.status(fh).set_tagged(false);
    }

    //deal with dangling not labeled components
    if (!notLabeledComponents.empty())
    {
        std::map<uint, uint> maTriangleToComponentMap;
        std::vector<double> componentsVolume;

        for (uint componentId = 0; componentId < components.size(); ++componentId)
        {
            const auto& component = components[componentId];
            std::vector<OpenMesh::SmartFaceHandle> smartComponent;

            for (auto fh : component)
            {
                uint maTriangle = m_tIdToOriginalTId[pFhToMaFh[fh]];
                maTriangleToComponentMap[maTriangle] = componentId;
                smartComponent.push_back(OpenMesh::make_smart(fh, m_mesh));
            }

            int sign = calcVolumeSignExact(smartComponent);
            double volume = calcSignedVolume(smartComponent);
            double volumeSize = std::abs(volume);
            componentsVolume.push_back(volumeSize* sign);
        }

        //first, merge not labeled components
        std::vector<std::vector<tFaceHandle>> mergedUnlabeledComponents;
        std::unordered_set<int> visitedComponents;

        for (auto componentId : notLabeledComponents)
        {
            if (visitedComponents.find(componentId) != visitedComponents.end())
            {
                continue;
            }

            const auto& currentComponent = components[componentId];
            visitedComponents.insert(componentId);

            std::vector<tFaceHandle> mergedComponent;
            std::queue<tFaceHandle> facesToVisit;

            for (auto fh : currentComponent)
            {
                facesToVisit.push(fh);
            }

            while(!facesToVisit.empty())
            {
                auto topFace = OpenMesh::make_smart(facesToVisit.front(), m_mesh);
                facesToVisit.pop();

                if (topFace.tagged() || topFace.deleted())
                {
                    continue;
                }

                //add to this merged component
                mergedComponent.push_back(topFace);

                uint maTriangle = m_tIdToOriginalTId[pFhToMaFh[topFace]];
                auto seedOrigin = ma.m_labels[maTriangle];

                std::vector<uint> faceVertices;
                faceVertices.reserve(3);
                faceVertices.push_back(ma.m_triangles[maTriangle * 3 + 0]);
                faceVertices.push_back(ma.m_triangles[maTriangle * 3 + 1]);
                faceVertices.push_back(ma.m_triangles[maTriangle * 3 + 2]);

                int e0 = ma.getVertexEdgeId(faceVertices[0], faceVertices[1]);
                int e1 = ma.getVertexEdgeId(faceVertices[1], faceVertices[2]);
                int e2 = ma.getVertexEdgeId(faceVertices[2], faceVertices[0]);

                for (auto e : { e0,e1,e2 })
                {
                    //TODO: should it be only non manifold edges? maybe not
                    if (!ma.m_intersectionEdgeProp[e]/* && !ma.isEdgeManifold(e)*/) //this is self-intersecting edge
                    {
                        //look at adjacent faces and take their labeling with component volume
                        for (auto adjFace : ma.getFaceAdjacentEdges(e))
                        {
                            //skip different meshes
                            if (ma.m_labels[adjFace] != seedOrigin)
                            {
                                continue;
                            }

                            auto otherComponentId = maTriangleToComponentMap[adjFace];

                            if (otherComponentId == componentId)
                            {
                                continue;
                            }

                            //skip not labeled component as well
                            if (componentsLabels[otherComponentId][NBIT - 2] == 1 && visitedComponents.find(otherComponentId) == visitedComponents.end())
                            {
                                const auto& otherUnlabeledComponent = components[otherComponentId];

                                for (auto otherCmpFh : otherUnlabeledComponent)
                                {
                                    facesToVisit.push(otherCmpFh);
                                }

                                visitedComponents.insert(otherComponentId);
                            }
                        }
                    }
                }


                //mark this face as processed
                m_mesh.status(topFace).set_tagged(true);
            }

            mergedUnlabeledComponents.push_back(mergedComponent);
        }

        for(const auto& mergedComponent : mergedUnlabeledComponents)
        {
            double largestVolume = std::numeric_limits<double>::lowest();
            int largestNeighbouringComponent = -1;

            for (auto fh : mergedComponent)
            {
                uint maTriangle = m_tIdToOriginalTId[pFhToMaFh[fh]];
                auto seedOrigin = ma.m_labels[maTriangle];

                std::vector<uint> faceVertices;
                faceVertices.reserve(3);
                faceVertices.push_back(ma.m_triangles[maTriangle * 3 + 0]);
                faceVertices.push_back(ma.m_triangles[maTriangle * 3 + 1]);
                faceVertices.push_back(ma.m_triangles[maTriangle * 3 + 2]);

                int e0 = ma.getVertexEdgeId(faceVertices[0], faceVertices[1]);
                int e1 = ma.getVertexEdgeId(faceVertices[1], faceVertices[2]);
                int e2 = ma.getVertexEdgeId(faceVertices[2], faceVertices[0]);

                for (auto e : { e0,e1,e2 })
                {
                    if (!ma.m_intersectionEdgeProp[e] /*&& !ma.isEdgeManifold(e)*/) //this is self-intersecting edge
                    {
                        //look at adjacent faces and take their labeling with component volume
                        for (auto adjFace : ma.getFaceAdjacentEdges(e))
                        {
                            //skip different meshes
                            if (ma.m_labels[adjFace] != seedOrigin)
                            {
                                continue;
                            }

                            auto otherComponentId = maTriangleToComponentMap[adjFace];

                            //if (otherComponentId == componentId)
                            //{
                            //    continue;
                            //}

                            //skip not labeled component as well
                            if (componentsLabels[otherComponentId][NBIT - 2] == 1)
                            {
                                continue;
                            }

                            auto otherComponentVolume = componentsVolume[otherComponentId];

                            if (otherComponentVolume > largestVolume)
                            {
                                largestVolume = otherComponentVolume;
                                largestNeighbouringComponent = otherComponentId;
                            }
                        }
                    }
                }
            }

            if (largestNeighbouringComponent != -1)
            {
                auto bestLabel = componentsLabels[largestNeighbouringComponent];

                for (auto fh : mergedComponent)
                {
                    labeling[fh] = bestLabel;
                }
            }
            else
            {
                //This is completely isolated component - this needs to be classified via GWN
                m_isolatedComponents.push_back({ mergedComponent.begin(), mergedComponent.end() });
            }
        }

        //reset tagged
        for (auto fh : m_mesh.faces())
        {
            m_mesh.status(fh).set_tagged(false);
        }

    }

    return true;
}

template<typename MeshType>
inline void DMB::MatrixMesh<MeshType>::classifyIsolatedComponents(MatrixMesh<MeshType>& other)
{
    //no isolated components, nothing to do
    if (m_isolatedComponents.empty())
    {
        return;
    }

    OpenMesh::FProp< std::bitset<NBIT>> labeling(m_mesh, m_pFaceIOLabeling);

    //TODO: This is UNUSABLE for variadic booleans, as it always creates and destroys the FWN, but for now
    // it is better, as it frees the memory after the use...
    // TODO: REWORK for variadic booleans

    //we need to build FWN of the other mesh
    auto acc = other.getFWN();

    static const int nSamples = 10;

    //classifiy the isolated components
    for (auto componentFaces : m_isolatedComponents)
    {
        //pick 5 random sample vertices from this component
        std::random_device rd{};
        std::default_random_engine rng(rd());
        std::uniform_int_distribution<int> randomSamples(0, static_cast<int>(componentFaces.size() - 1));

        std::array<bool, nSamples> votes;

        for (int i = 0; i < nSamples; ++i)
        {
            auto fh = OpenMesh::make_smart(componentFaces[randomSamples(rng)], m_mesh);

            // Get a property manager of the points property of the mesh to use as functor
            const auto& points = OpenMesh::getPointsProperty(m_mesh);
            tPoint centroid = fh.vertices().avg(points);

            votes[i] = acc->windingNumber(centroid) > 0.5;
        }

        int insideCount = std::count_if(votes.begin(), votes.end(), [](bool b) {return b; });

        std::bitset<NBIT> IOLabel;

        //is inside the other mesh
        if (insideCount > nSamples / 2)
        {
            IOLabel |= other.m_label;
        }

        //label this component
        for (auto fh : componentFaces)
        {
            labeling[fh][NBIT - 2] = 0;
            labeling[fh] |= IOLabel;
        }
    }
}

template<typename MeshType>
inline std::shared_ptr <typename DMB::MatrixMesh<MeshType>::tFWN> DMB::MatrixMesh<MeshType>::getFWN()
{
    // TODO: REWORK for variadic booleans

    m_mesh.update_normals();
    auto acc = std::make_shared<tFWN>();
    acc->build(m_mesh);

    return acc;
}

template<typename MeshType>
inline void DMB::MatrixMesh<MeshType>::detectBoundaries()
{
    //reset tagged flag for edges
    for (auto eh : m_mesh.edges())
    {
        m_mesh.status(eh).set_tagged(false);
    }

    auto pBoundaryId = OpenMesh::VProp<int>(-1, m_mesh, "bo_pBoundaryId");
    m_pboundaryId = pBoundaryId.getRawProperty();
    int boundaryId = 0;

    //set boundary curves ids
    for (auto eh : m_mesh.edges())
    {
        if (eh.is_boundary() && !eh.tagged())
        {
            //start new boundary
            auto seedHe = eh.h0().is_boundary() ? eh.h0() : eh.h1();
            auto nextHe = seedHe.next();

            m_mesh.status(eh).set_tagged(true);
            pBoundaryId[seedHe.from()] = boundaryId;
            pBoundaryId[seedHe.to()] = boundaryId;

            while (nextHe != seedHe)
            {
                m_mesh.status(nextHe.edge()).set_tagged(true);
                pBoundaryId[nextHe.to()] = boundaryId;

                nextHe = nextHe.next();
            }

        }

        ++boundaryId;
    }

    //reset tagged flag for edges
    for (auto eh : m_mesh.edges())
    {
        m_mesh.status(eh).set_tagged(false);
    }
}

template<typename MeshType>
inline void DMB::MatrixMesh<MeshType>::handleCoplanarFaces()
{
    OpenMesh::FProp<bool> coplanarFace(m_mesh, m_pCoplanarFace);

    OpenMesh::EProp<bool> pCoplanarEdge(false, m_mesh, "bo_coplanarEdge");
    m_pCoplanarEdge = pCoplanarEdge.getRawProperty();

    OpenMesh::VProp<bool> coplanarVertex(false, m_mesh, "bo_coplanarVertex");
    m_pCoplanarVertex = coplanarVertex.getRawProperty();


    for (auto fh : m_mesh.faces())
    {
        if (coplanarFace[fh])
        {
            for (auto vh : fh.vertices())
            {
                coplanarVertex[vh] = true;
            }

            for (auto eh : fh.edges())
            {
                pCoplanarEdge[eh] = true;
            }

            //TODO: tmp
            //m_mesh.delete_face(fh);
        }
    }

}

template<typename MeshType>
inline double DMB::MatrixMesh<MeshType>::calcSignedVolume(const std::vector<OpenMesh::SmartFaceHandle>& component)
{

    //http://chenlab.ece.cornell.edu/Publication/Cha/icip01_Cha.pdf
    //https://www.ams.org/journals/mcom/1986-46-173/S0025-5718-1986-0815838-7/S0025-5718-1986-0815838-7.pdf
    //https://dsp.stackexchange.com/questions/7856/calculating-the-volume-of-a-triangular-mesh

    double volume = 0;

    tPoint centroid(0, 0, 0);

    {
        std::unordered_set<tVertexHandle> vertices;

        for (auto fh : component)
        {
            for (auto vh : OpenMesh::make_smart(fh, m_mesh).vertices())
            {
                vertices.insert(vh);
            }
        }

        for (auto vh : vertices)
        {
            centroid += m_mesh.point(vh);
        }

        centroid /= vertices.size();
    }

    for (auto fh : component)
    {
        std::vector<tPoint> pts{};

        for (auto vh : OpenMesh::make_smart(fh, m_mesh).vertices())
        {
            pts.push_back(m_mesh.point(vh) - centroid);
        }
        
        typename MeshType::Scalar v =  (
            -pts[2][0] * pts[1][1] * pts[0][2]
            + pts[1][0] * pts[2][1] * pts[0][2]
            + pts[2][0] * pts[0][1] * pts[1][2]
            - pts[0][0] * pts[2][1] * pts[1][2]
            - pts[1][0] * pts[0][1] * pts[2][2]
            + pts[0][0] * pts[1][1] * pts[2][2]
            );

        volume += v;
    }

    return volume; // should be multiplied by (1.0 / 6.0), but it does not make a difference for my use case
}

template<typename MeshType>
inline bigfloat DMB::MatrixMesh<MeshType>::calcSignedVolumeExact(const std::vector<OpenMesh::SmartFaceHandle>& component)
{
    bigfloat volume = 0;

    tExactPoint centroid;
    centroid[0] = bigfloat(0);
    centroid[1] = bigfloat(0);
    centroid[2] = bigfloat(0);


    {
        std::unordered_set<tVertexHandle> vertices;

        for (auto fh : component)
        {
            for (auto vh : OpenMesh::make_smart(fh, m_mesh).vertices())
            {
                vertices.insert(vh);
            }
        }

        for (auto vh : vertices)
        {
            centroid[0] = centroid[0] + bigfloat(m_mesh.point(vh)[0]);
            centroid[1] = centroid[1] + bigfloat(m_mesh.point(vh)[1]);
            centroid[2] = centroid[2] + bigfloat(m_mesh.point(vh)[2]);
        }

        centroid[0] = centroid[0] * bigfloat((1.0 / (double)vertices.size()));
        centroid[1] = centroid[1] * bigfloat((1.0 / (double)vertices.size()));
        centroid[2] = centroid[2] * bigfloat((1.0 / (double)vertices.size()));
    }


    for (auto fh : component)
    {
        std::vector<tExactPoint> pts;

        for (auto vh : OpenMesh::make_smart(fh, m_mesh).vertices())
        {
            tExactPoint p;
            p[0] = bigfloat(m_mesh.point(vh)[0]) - centroid[0];
            p[1] = bigfloat(m_mesh.point(vh)[1]) - centroid[1];
            p[2] = bigfloat(m_mesh.point(vh)[2]) - centroid[2];

            pts.push_back(p);
        }

        bigfloat v = (
            - pts[2][0] * pts[1][1] * pts[0][2]
            + pts[1][0] * pts[2][1] * pts[0][2]
            + pts[2][0] * pts[0][1] * pts[1][2]
            - pts[0][0] * pts[2][1] * pts[1][2]
            - pts[1][0] * pts[0][1] * pts[2][2]
            + pts[0][0] * pts[1][1] * pts[2][2]
            );

        volume = volume + v;
    }

    return volume; // should be multiplied by (1.0 / 6.0), but it does not make a difference for my use case
}

template<typename MeshType>
inline int DMB::MatrixMesh<MeshType>::calcVolumeSignExact(const std::vector<OpenMesh::SmartFaceHandle>& component)
{
    return calcSignedVolumeExact(component).sgn();
}

template<typename MeshType>
inline void DMB::MatrixMesh<MeshType>::debug_CouldNotAddFace(tVertexHandle v0, tVertexHandle v1, tVertexHandle v2)
{
    for (auto fh : m_mesh.faces())
    {
        m_mesh.set_color(fh, { 192,192,192 });
    }

    m_mesh.update_normals();
    T3D_VDT_STORE_NAMED_MESH_AND_WAIT("", m_mesh);

    MeshType failedFace;
    auto vh0 = failedFace.add_vertex(m_mesh.point(v0));
    auto vh1 = failedFace.add_vertex(m_mesh.point(v1));
    auto vh2 = failedFace.add_vertex(m_mesh.point(v2));

    auto fh = failedFace.add_face(vh0, vh1, vh2);

    failedFace.set_color(fh, { 255,0,0 });

    failedFace.update_normals();
    T3D_VDT_STORE_NAMED_MESH_AND_WAIT("", failedFace);

}

template<typename MeshType>
inline void DMB::MatrixMesh<MeshType>::debug_CurrentMesh()
{
    for (auto fh : m_mesh.faces())
    {
        m_mesh.set_color(fh, { 192,192,192 });
    }

    m_mesh.update_normals();
    T3D_VDT_STORE_NAMED_MESH_AND_WAIT("", m_mesh);
}

template<typename MeshType>
template<typename F>
inline void DMB::MatrixMesh<MeshType>::classifyMeshArrangement(MeshArrangement<MeshType>& ma, int operandLabel, const F& predicate)
{
    OpenMesh::FProp<std::bitset<NBIT>> labeling(m_mesh, m_pFaceIOLabeling);
    ////OpenMesh::VProp<uint> operandToMeshArrangementIndex(m_mesh, m_pOperandToMeshArrangementIndex);
    OpenMesh::FProp<uint> pFhToMaFh(m_mesh, m_pFhToMaFh);

    for (auto fh : m_mesh.faces())
    {
        if (!fh.is_valid() || fh.deleted())
        {
            continue;
        }

        /*
         * predicateResult ==  0 -> not part of result
         * predicateResult ==  1 -> part of result as is
         * predicateResult == -1 -> part of result but flipped
        */
        auto faceLabel = labeling[fh];
        auto predicateResult = predicate(operandLabel, faceLabel);


        ma.m_boPredicates[m_tIdToOriginalTId[pFhToMaFh[fh]]] = predicateResult;
    }
}

template<typename MeshType>
template<typename F>
inline void DMB::MatrixMesh<MeshType>::classifyMeshArrangementWithMAFaceHandle(MeshArrangement<MeshType>& ma, int operandLabel, const F& predicate)
{
    OpenMesh::FProp<uint> pFhToMaFh(m_mesh, m_pFhToMaFh);

    for (auto fh : m_mesh.faces())
    {
        if (!fh.is_valid() || fh.deleted())
        {
            continue;
        }

        /*
         * predicateResult ==  0 -> not part of result
         * predicateResult ==  1 -> part of result as is
         * predicateResult == -1 -> part of result but flipped
        */
        auto predicateResult = predicate(operandLabel, m_tIdToOriginalTId[pFhToMaFh[fh]]);

        ma.m_boPredicates[m_tIdToOriginalTId[pFhToMaFh[fh]]] = predicateResult;
    }
}

template<typename MeshType>
template<typename F>
inline void DMB::MatrixMesh<MeshType>::debug_showMesh(const F& func)
{
    using tPoint = typename MeshType::Point;
    using tVertexHandle = typename MeshType::VertexHandle;
    using tFaceHandle = typename MeshType::FaceHandle;
    using tEdgeHandle = typename MeshType::EdgeHandle;
    using tHalfedgeHandle = typename MeshType::HalfedgeHandle;
    using tTriangle = std::array<tPoint, 3>;
    using tTriangleIndices = std::array<uint, 3>;

    MeshType debugMesh;


    for (int i = 0; i < std::size(m_triangles); i += 3)
    {
        uint tId = i / 3;
        uint maVh0 = m_triangles[i + 0];
        uint maVh1 = m_triangles[i + 1];
        uint maVh2 = m_triangles[i + 2];

        tTriangle triangle =
        {
            tPoint{ m_coordinates[maVh0 * 3 + 0], m_coordinates[maVh0 * 3 + 1], m_coordinates[maVh0 * 3 + 2] },
            tPoint{ m_coordinates[maVh1 * 3 + 0], m_coordinates[maVh1 * 3 + 1], m_coordinates[maVh1 * 3 + 2] },
            tPoint{ m_coordinates[maVh2 * 3 + 0], m_coordinates[maVh2 * 3 + 1], m_coordinates[maVh2 * 3 + 2] },
        };

        auto vh0 = debugMesh.add_vertex(triangle[0]);
        auto vh1 = debugMesh.add_vertex(triangle[1]);
        auto vh2 = debugMesh.add_vertex(triangle[2]);

        auto fh = debugMesh.add_face(vh0, vh1, vh2);

        debugMesh.set_color(fh, { 192,192,192 });

        func(debugMesh, fh, tId);


    }

    debugMesh.update_normals();
    T3D_VDT_STORE_NAMED_MESH_AND_WAIT("", debugMesh);
}


template<typename MeshType>
inline void DMB::MatrixMesh<MeshType>::transferArrangementPropertiesToMesh()
{
    OpenMesh::FProp<uint> pFhToMaFh(m_mesh, m_pFhToMaFh);
    OpenMesh::VProp<uint> pVhToMaVId(m_mesh, m_pVhToMaVId);

    OpenMesh::FProp<bool> coplanarFace(false, m_mesh, "bo_coplanarFace");
    m_pCoplanarFace = coplanarFace.getRawProperty();
    OpenMesh::FProp< std::bitset<NBIT>> faceIOLabeling(0, m_mesh, "bo_labeling");
    m_pFaceIOLabeling = faceIOLabeling.getRawProperty();
    OpenMesh::FProp< bool > intersectionFace(false, m_mesh, "bo_intersectionFace");
    m_pIntersectionFace = intersectionFace.getRawProperty();

    OpenMesh::EProp<bool> pIntersectionEdge(false, m_mesh, "bo_pIntersectionEdge");
    m_pIntersectionEdge = pIntersectionEdge.getRawProperty();

    for (auto fh : m_mesh.faces())
    {
        uint tId = pFhToMaFh[fh];
        coplanarFace[fh] = m_coplanarFaceProp[tId];
        faceIOLabeling[fh] = m_faceIOLabeling[tId];
        intersectionFace[fh] = m_intersectionFaceProp[tId];

        for (auto eh : fh.edges())
        {
            auto vh0 = eh.v0();
            auto vh1 = eh.v1();

            uint v0 = pVhToMaVId[vh0];
            uint v1 = pVhToMaVId[vh1];

            uint eId = m_vertexEdgeMatrix.coeff(v0, v1);

            pIntersectionEdge[eh] = m_intersectionEdgeProp[eId];
        }
    }
}

template<typename MeshType>
inline void DMB::MatrixMesh<MeshType>::saveMatrixMesh(std::string outputFileName)
{
    save(outputFileName, m_coordinates, m_triangles);
}
