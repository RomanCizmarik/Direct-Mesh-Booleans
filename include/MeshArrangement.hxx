//TODO: License

#include "MeshArrangement.h"

template<typename MeshType>
DMB::MeshArrangement<MeshType>::MeshArrangement() :
    m_multiplier(0),
    m_nEdges(0),
    m_nFaces(0),
    m_nVertices(0)
{
}

template<typename MeshType>
inline void DMB::MeshArrangement<MeshType>::updateMatrices()
{
    m_nFaces = m_triangles.size();
    m_nVertices = m_coordinatesImplicit.size();

    //TODO: make some property management system
    m_boPredicates.resize(std::size(m_labels), 0);


    buildVertexEdgeMatrix();
    buildEdgeFaceAdjacency();
    buildEdgeVertices();

}

template<typename MeshType>
inline void DMB::MeshArrangement<MeshType>::detectCoplanarFaces()
{
    m_coplanarFace.clear();
    m_coplanarFace.resize(m_nFaces, false);

    for (uint i = 0; i < m_triangles.size();  i += 3)
    {
        uint tId = i / 3;

        auto faceLabel = m_labels[tId];

        if (faceLabel.count() > 1)
        {
            m_coplanarFace[tId] = true;
        }
    }
}

template<typename MeshType>
inline void DMB::MeshArrangement<MeshType>::detectIntersectionEdges()
{
    m_intersectionEdgeProp.clear();
    m_intersectionEdgeProp.resize(m_nEdges +1, false);

    m_intersectionEdgeFaceProp.clear();
    m_intersectionEdgeFaceProp.resize(m_nFaces, false);

    m_intersectionVertexProp.clear();
    m_intersectionVertexProp.resize(m_nVertices, false);

    for (uint eId = 0; eId < m_edgeFaceAdjacency.size(); ++eId)
    {
        if (m_edgeFaceAdjacency[eId].size() <= 2)
        {
            continue;
        }

        std::unordered_set<uint> adjacentFacesLabels;

        std::bitset<NBIT> totalBitset;

        for (auto tId : m_edgeFaceAdjacency[eId])
        {
            auto faceLabel = m_labels[tId];
            totalBitset = (faceLabel | totalBitset);
        }

        if (totalBitset.count() > 1)
        {
            m_intersectionEdgeProp[eId] = true;

            for (auto tId : m_edgeFaceAdjacency[eId])
            {
                m_intersectionEdgeFaceProp[tId] = true;
            }

            auto edgeVertices = m_edgeVertices[eId];
            m_intersectionVertexProp[edgeVertices.first] = true;
            m_intersectionVertexProp[edgeVertices.second] = true;
            
        }
    }
}

template<typename MeshType>
inline void DMB::MeshArrangement<MeshType>::classifyFaces()
{
    m_faceIOLabeling.clear();
    m_faceIOLabeling.resize(m_nFaces);

    m_conflictingFace.clear();
    m_conflictingFace.resize(m_nFaces, false);

    m_faceClassifiedByComponent.clear();
    m_faceClassifiedByComponent.resize(m_nFaces, -1);

    
    auto remainingVId = [this](uint tId, uint v0, uint v1) -> int
    {
        std::vector<uint> faceVertices;
        faceVertices.reserve(3);
        faceVertices.push_back(m_triangles[tId * 3 + 0]);
        faceVertices.push_back(m_triangles[tId * 3 + 1]);
        faceVertices.push_back(m_triangles[tId * 3 + 2]);

        auto it = std::find_if(faceVertices.begin(), faceVertices.end(), [v0, v1](uint tVId) { return tVId != v0 && tVId != v1; });

        if (it == faceVertices.end())
        {
            return -1;
        }

        return (*it);
    };

    auto sortedPoints = [this](uint tId, uint v0, uint v1) 
    {
       

        std::array<uint, 3> faceVertices;
        faceVertices[0] = m_triangles[tId * 3 + 0];
        faceVertices[1] = m_triangles[tId * 3 + 1];
        faceVertices[2] = m_triangles[tId * 3 + 2];

        //std::cout << "before rotate:\t\t";
        //for (auto v : faceVertices)
        //    std::cout << v << ' ';
        //std::cout << '\n';


        auto it = std::find_if(faceVertices.begin(), faceVertices.end(), [v0, v1](uint tVId) { return tVId != v0 && tVId != v1; });
        //std::rotate(faceVertices.rbegin(), std::make_reverse_iterator(it), faceVertices.rend());
        std::rotate(faceVertices.begin(), it, faceVertices.end());

        //std::cout << "after rotate:\t\t";
        //for (auto v : faceVertices)
        //    std::cout << v << ' ';
        //std::cout << '\n';

        std::array<genericPoint*, 3> triangle;

        triangle[0] = m_coordinatesImplicit[faceVertices[0]];
        triangle[1] = m_coordinatesImplicit[faceVertices[1]];
        triangle[2] = m_coordinatesImplicit[faceVertices[2]];

        return triangle;
    };

    auto trianglePoints = [this](uint tId)
    {
        std::array<genericPoint*, 3> triangle;
        for (int j = 0; j < 3; ++j)
        {
            triangle[j] = m_coordinatesImplicit[m_triangles[tId * 3 + j]];
        }

        return triangle;
    };

    auto insideConvex = [this, &trianglePoints](uint leftTriId, uint rightTriId, genericPoint* p) ->bool
    {
        std::array<genericPoint*, 3 > leftTri = trianglePoints(leftTriId);
        std::array<genericPoint*, 3 > rightTri = trianglePoints(rightTriId);

        //std::cout << (genericPoint::orient3D(*leftTri[0], *leftTri[1], *leftTri[2], *p) * m_multiplier) << std::endl;
        //std::cout << (genericPoint::orient3D(*rightTri[0], *rightTri[1], *rightTri[2], *p) * m_multiplier) << std::endl;
        //std::cout << (genericPoint::orient3D(*leftTri[0], *leftTri[1], *leftTri[2], *p)) << std::endl;
        //std::cout << (genericPoint::orient3D(*rightTri[0], *rightTri[1], *rightTri[2], *p)) << std::endl;

        int orientationL = (genericPoint::orient3D(*leftTri[0], *leftTri[1], *leftTri[2], *p) * (int)m_multiplier) ;
        int orientationR = (genericPoint::orient3D(*rightTri[0], *rightTri[1], *rightTri[2], *p) * (int)m_multiplier) ;

        //inside, if left of both, outside otherwise
        //orientation sum < 0 => inside
        //return orientationL && orientationR;
        return (orientationL + orientationR) < 0;
    };

    auto insideConcave = [this, &trianglePoints](uint leftTriId, uint rightTriId, genericPoint* p) ->bool
    {
        std::array<genericPoint*, 3 > leftTri = trianglePoints(leftTriId);
        std::array<genericPoint*, 3 > rightTri = trianglePoints(rightTriId);


        int orientationL = (genericPoint::orient3D(*leftTri[0], *leftTri[1], *leftTri[2], *p) * (int)m_multiplier);
        int orientationR = (genericPoint::orient3D(*rightTri[0], *rightTri[1], *rightTri[2], *p) * (int)m_multiplier);

        //std::cout << "insideConcave" << std::endl;

        //std::cout << (genericPoint::orient3D(*leftTri[0], *leftTri[1], *leftTri[2], *p) * m_multiplier) << std::endl;
        //std::cout << (genericPoint::orient3D(*rightTri[0], *rightTri[1], *rightTri[2], *p) * m_multiplier) << std::endl;
        //std::cout << (genericPoint::orient3D(*leftTri[0], *leftTri[1], *leftTri[2], *p)) << std::endl;
        //std::cout << (genericPoint::orient3D(*rightTri[0], *rightTri[1], *rightTri[2], *p)) << std::endl;
        //std::cout << "----------------" << std::endl;


        //outside, if right of both, inside otherwise
        //return !(!orientationL && !orientationR);
        //orientation sum > 0 => outside
        return !((orientationL + orientationR) > 0);
    };


    std::vector<bool> faceClassified;
    faceClassified.resize(m_nFaces, false);

    for (int eId = 0; eId < m_edgeFaceAdjacency.size(); ++eId)
    {
        if (!m_intersectionEdgeProp[eId])
        {
            continue;
        }

        const auto& adjacentFaces = m_edgeFaceAdjacency[eId];

        //if (adjacentFaces.size() != 4)
        //{
        //    std::cout << "not possible" << std::endl;

        //    debug_showMesh<MeshType>([this, &adjacentFaces](MeshType& mesh, OpenMesh::SmartFaceHandle fh, uint triangleId)
        //    {
        //        mesh.set_color(fh, { 192,192,192 });

        //        std::bitset<NBIT> left;
        //        std::bitset<NBIT> right;
        //        std::bitset<NBIT> coplanar;

        //        left[0] = 1;
        //        right[1] = 1;
        //        coplanar[0] = 1;
        //        coplanar[1] = 1;

        //        if (std::find(adjacentFaces.begin(), adjacentFaces.end(), triangleId) != adjacentFaces.end())
        //        {
        //            if(m_labels[triangleId] == left)
        //                mesh.set_color(fh, { 255,0,0 });
        //            if (m_labels[triangleId] == right)
        //                mesh.set_color(fh, { 0,255,0 });
        //            if (m_labels[triangleId] == coplanar)
        //                mesh.set_color(fh, { 0,0,255 });
        //        }

        //    });
        //}

        auto edgeVertices = m_edgeVertices[eId];

        auto inside = [this, &insideConvex, &insideConcave, &trianglePoints, &remainingVId, &edgeVertices](uint leftTriId, uint rightTriId, genericPoint* p) ->bool
            {
                std::array<genericPoint*, 3 > leftTri = trianglePoints(leftTriId);
                std::array<genericPoint*, 3 > rightTri = trianglePoints(rightTriId);


                bool orientationRL = (genericPoint::orient3D(*leftTri[0], *leftTri[1], *leftTri[2], *m_coordinatesImplicit[remainingVId(rightTriId, edgeVertices.first, edgeVertices.second)]) * m_multiplier) <= 0;
                bool orientationLR = (genericPoint::orient3D(*rightTri[0], *rightTri[1], *rightTri[2], *m_coordinatesImplicit[remainingVId(leftTriId, edgeVertices.first, edgeVertices.second)]) * m_multiplier) <= 0;

                if (orientationRL != orientationLR)
                {
                    std::cout << "unorientable edge, should not happen!" << std::endl;
                    assert(false);
                }

                //convex
                if (orientationRL && orientationLR)
                {
                    //std::cout << "convex" << std::endl;
                    //std::cout << (genericPoint::orient3D(*leftTri[0], *leftTri[1], *leftTri[2], *m_coordinatesImplicit[remainingVId(rightTriId, edgeVertices.first, edgeVertices.second)]) * m_multiplier) << std::endl;
                    //std::cout << (genericPoint::orient3D(*leftTri[0], *leftTri[1], *leftTri[2], *m_coordinatesImplicit[remainingVId(rightTriId, edgeVertices.first, edgeVertices.second)]) * m_multiplier) << std::endl;
                    //std::cout << (genericPoint::orient3D(*rightTri[0], *rightTri[1], *rightTri[2], *m_coordinatesImplicit[remainingVId(leftTriId, edgeVertices.first, edgeVertices.second)])) << std::endl;
                    //std::cout << (genericPoint::orient3D(*rightTri[0], *rightTri[1], *rightTri[2], *m_coordinatesImplicit[remainingVId(leftTriId, edgeVertices.first, edgeVertices.second)])) << std::endl;
                    //std::cout << "----------------" << std::endl;


                    return insideConvex(leftTriId, rightTriId, p);
                }
                //concave
                else
                {
                    //std::cout << "concave" << std::endl;

                    //std::cout << (genericPoint::orient3D(*leftTri[0], *leftTri[1], *leftTri[2], *m_coordinatesImplicit[remainingVId(rightTriId, edgeVertices.first, edgeVertices.second)]) * m_multiplier)   << std::endl;
                    //std::cout << (genericPoint::orient3D(*leftTri[0], *leftTri[1], *leftTri[2], *m_coordinatesImplicit[remainingVId(rightTriId, edgeVertices.first, edgeVertices.second)]) * m_multiplier)   << std::endl; 
                    //std::cout << (genericPoint::orient3D(*rightTri[0], *rightTri[1], *rightTri[2], *m_coordinatesImplicit[remainingVId(leftTriId, edgeVertices.first, edgeVertices.second)])) << std::endl; 
                    //std::cout << (genericPoint::orient3D(*rightTri[0], *rightTri[1], *rightTri[2], *m_coordinatesImplicit[remainingVId(leftTriId, edgeVertices.first, edgeVertices.second)])) << std::endl; 

                    //std::cout << "----------------" << std::endl;

                    //debug_showMesh<MeshType>([this, &leftTriId, &rightTriId](MeshType& mesh, OpenMesh::SmartFaceHandle fh, uint triangleId)
                    //{
                    //    mesh.set_color(fh, { 192,192,192 });

                    //    if (triangleId == leftTriId)
                    //    {
                    //        mesh.set_color(fh, { 255,0,0 });
                    //    }
                    //    if (triangleId == rightTriId)
                    //    {
                    //        mesh.set_color(fh, { 0,255,0 });
                    //    }

                    //});


                    return insideConcave(leftTriId, rightTriId, p);
                }

                return false;
            };

        auto insideSimplified = [this, &trianglePoints](uint leftTriId, uint rightTriId, genericPoint* p)
            {
                std::array<genericPoint*, 3 > leftTri = trianglePoints(leftTriId);
                std::array<genericPoint*, 3 > rightTri = trianglePoints(rightTriId);

                int orientationL = (genericPoint::orient3D(*leftTri[0], *leftTri[1], *leftTri[2], *p) * (int)m_multiplier);
                int orientationR = (genericPoint::orient3D(*rightTri[0], *rightTri[1], *rightTri[2], *p) * (int)m_multiplier);

                return (orientationL + orientationR) <= 0;
            };

        auto radialSort = [this, &edgeVertices, &trianglePoints, &remainingVId](uint tId, const std::vector<uint>& trianglesToSort)
            {
                std::array<genericPoint*, 3 > triPoints = trianglePoints(tId);

                std::vector<uint> sortedTris = trianglesToSort;

                std::sort(sortedTris.begin(), sortedTris.end(), [this, &edgeVertices, &triPoints, &remainingVId](uint leftTId, uint rightTId)
                    {
                        uint leftVh = remainingVId(leftTId, edgeVertices.first, edgeVertices.second);
                        genericPoint* pLeft = m_coordinatesImplicit[leftVh];

                        uint rightVh = remainingVId(rightTId, edgeVertices.first, edgeVertices.second);
                        genericPoint* pRight = m_coordinatesImplicit[rightVh];

                        int orientationL = (genericPoint::orient3D(*triPoints[0], *triPoints[1], *triPoints[2], *pLeft) * (int)m_multiplier);
                        int orientationR = (genericPoint::orient3D(*triPoints[0], *triPoints[1], *triPoints[2], *pRight) * (int)m_multiplier);

                        //left is under the triangle, right is above
                        if (orientationL < 0 && orientationR > 0)
                        {
                            //left is "smaller" than right
                            return true;
                        }

                        //left is above the triangle, right is under
                        if (orientationL > 0 && orientationR < 0)
                        {
                            //left is "greater" than right
                            return false;
                        }

                        int orientationLR = (genericPoint::orient3D(*triPoints[0], *triPoints[1], *pLeft, *pRight) * (int)m_multiplier);
                        return orientationLR < 0;
                        /*
                        //both are above
                        if (orientationL > 0 && orientationR > 0)
                        {
                            //the one closer to the original tri is greater
                            int orientationLR = (genericPoint::orient3D(*triPoints[0], *triPoints[1], *pLeft, *pRight) * (int)m_multiplier);
                            return orientationLR < 0;
                        }

                        //both are under or coplanar
                        if (orientationL > 0 && orientationR < 0)
                        {
                            //left is "greater" than right
                            int orientationLR = (genericPoint::orient3D(*triPoints[0], *triPoints[1], *pLeft, *pRight) * (int)m_multiplier);
                            return orientationLR < 0;
                        }

                        //default, should never happen
                        return false;
                        */
                    });

                return sortedTris;

            };

        auto debugSaveTriangle = [this](uint tId, std::string name)
            {
                using tPoint = typename MeshType::Point;
                using tTriangle = std::array<tPoint, 3>;

                MeshType debugMesh;

                //uint tId = i / 3;
                uint maVh0 = m_triangles[tId * 3 + 0];
                uint maVh1 = m_triangles[tId * 3 + 1];
                uint maVh2 = m_triangles[tId * 3 + 2];

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

               // DMB::saveMesh(debugMesh, "C:/skola/PhD/Samples/booleans/radial/"+name+".obj");

            };

        std::unordered_map<std::bitset<NBIT>, std::vector<uint>> triangleGroups;

        //TODO: now I assume binary input only
        for (int i = 0; i < adjacentFaces.size(); ++i)
        {

            //TODO: check if face was already processed
            uint tId = adjacentFaces[i];

            if (m_coplanarFace[tId])
            {
                continue;
            }

            auto triangleLabel = m_labels[tId];

            triangleGroups[triangleLabel].push_back(tId);
        }

        for (auto [triangleLabel, triangles] : triangleGroups)
        {

            for (auto tId : triangles)
            {
                uint vIdToInvestigate = remainingVId(tId, edgeVertices.first, edgeVertices.second);
                genericPoint* p = m_coordinatesImplicit[vIdToInvestigate];

                std::bitset<NBIT> IOLabel;
                int currentComponentId = -1;

                for (auto const& pair : triangleGroups)
                {
                    auto otherTriangleLabel = pair.first;

                    if (otherTriangleLabel == triangleLabel)
                    {
                        continue;
                    }

                    const auto& otherTriangles = pair.second;

                    bool isInside = false;

                    if (otherTriangles.size() == 2)
                    {
                        //isInside = inside(otherTriangles[0], otherTriangles[1], p);
                        isInside = inside(otherTriangles[0], otherTriangles[1], p);

                        //pick component with max volume - but generally it should be the same component
                        currentComponentId = m_componentsVolume[m_faceComponentId[otherTriangles[0]]] > m_componentsVolume[m_faceComponentId[otherTriangles[1]]] ? m_faceComponentId[otherTriangles[0]] : m_faceComponentId[otherTriangles[1]];

                    }
                    else if (otherTriangles.size() == 1)
                    {
                        std::array<genericPoint*, 3 > leftTri = trianglePoints(otherTriangles[0]);

                        bool orientation = (genericPoint::orient3D(*leftTri[0], *leftTri[1], *leftTri[2], *p) * (int)m_multiplier) < 0;
                        isInside = orientation;

                        currentComponentId = m_faceComponentId[otherTriangles[0]];
                    }
                    else
                    {
                        auto sortedTris = radialSort(tId, otherTriangles);

                        std::vector<uint> enclosingTris = { sortedTris[0], sortedTris[sortedTris.size() - 1] };

                        isInside = inside(enclosingTris[0], enclosingTris[1], p);

                        //pick component with max volume - but generally it should be the same component
                        currentComponentId = m_componentsVolume[m_faceComponentId[enclosingTris[0]]] > m_componentsVolume[m_faceComponentId[enclosingTris[1]]] ? m_faceComponentId[enclosingTris[0]] : m_faceComponentId[enclosingTris[1]];

                        //std::cout << "unclassified for now" << std::endl;
                        //assert(false);
                        //throw std::runtime_error(std::string("Failed: Unclassified MA triangle."));

                        //debugSaveTriangle(tId, "baseTri");

                        //int i = 0;
                        //for (auto stId : sortedTris)
                        //{
                        //    debugSaveTriangle(stId, std::to_string(i));
                        //    ++i;
                        //}

                    }

                    if (isInside)
                    {
                        IOLabel = IOLabel | otherTriangleLabel;
                    }
                }

                if (faceClassified[tId])
                {
                    //this triangles will be classified by a component with the larges volume
                    if (m_componentsVolume[currentComponentId] < m_componentsVolume[m_faceClassifiedByComponent[tId]])
                    {
                        currentComponentId = m_faceClassifiedByComponent[tId];
                        IOLabel = m_faceIOLabeling[tId];
                    }
                }

                //TODO: lets assume that each classification came from another component -> copy component ID from MatrixMesh to MA
                // face is outside, only if it is outside of all components
                //so that should be logical OR
                //m_faceIOLabeling[tId] = (m_faceIOLabeling[tId] | IOLabel);

                m_faceIOLabeling[tId] = IOLabel;
                faceClassified[tId] = true;

                m_faceClassifiedByComponent[tId] = currentComponentId;
            }

        }

    }
}


template<typename MeshType>
inline uint DMB::MeshArrangement<MeshType>::getFaceEdgeId(uint tId, int i)
{
    uint v0 = m_triangles[tId * 3 + 0];
    uint v1 = m_triangles[tId * 3 + 1];
    uint v2 = m_triangles[tId * 3 + 2];

    std::array<uint, 3> edges;

    edges[0] = m_vertexEdgeMatrix.coeff(v0, v1);
    edges[1] = m_vertexEdgeMatrix.coeff(v1, v2);
    edges[2] = m_vertexEdgeMatrix.coeff(v2, v0);

    return edges[i];
}

template<typename MeshType>
inline uint DMB::MeshArrangement<MeshType>::getVertexEdgeId(uint v0, uint v1) const
{
    return m_vertexEdgeMatrix.coeff(v0, v1);
}

template<typename MeshType>
inline void DMB::MeshArrangement<MeshType>::initComponentsVolume()
{
    m_faceComponentId.clear();
    m_faceComponentId.resize(m_nFaces, -1);
    m_componentsVolume.clear();
}

template<typename MeshType>
inline bool DMB::MeshArrangement<MeshType>::isEdgeManifold(uint eId) const
{
    assert(eId < m_edgeFaceAdjacency.size());

    return m_edgeFaceAdjacency[eId].size() < 3;
}

template<typename MeshType>
inline const std::vector<uint>& DMB::MeshArrangement<MeshType>::getFaceAdjacentEdges(uint eId) const
{
    assert(eId < m_edgeFaceAdjacency.size());

    return m_edgeFaceAdjacency[eId];
}

template<typename MeshType>
inline void DMB::MeshArrangement<MeshType>::buildVertexEdgeMatrix()
{
    std::vector<tTriplet> tripletList;
    tripletList.reserve(std::size(m_coordinates));

    for (int i = 0; i < std::size(m_triangles); i += 3)
    {
        uint v0 = m_triangles[i + 0];
        uint v1 = m_triangles[i + 1];
        uint v2 = m_triangles[i + 2];

        tripletList.push_back(tTriplet(v0, v1, -1));
        tripletList.push_back(tTriplet(v1, v2, -1));
        tripletList.push_back(tTriplet(v2, v0, -1));
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
inline void DMB::MeshArrangement<MeshType>::buildEdgeFaceAdjacency()
{
    m_edgeFaceAdjacency.clear();
    m_edgeFaceAdjacency.resize(m_nEdges + 1); // edge Ids starts from 1

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
inline void DMB::MeshArrangement<MeshType>::buildEdgeVertices()
{
    m_edgeVertices.clear();
    m_edgeVertices.resize(m_nEdges + 1, { 0,0 });

    for (int k = 0; k < m_vertexEdgeMatrix.outerSize(); ++k)
    {
        for (tSparseIntMatrix::InnerIterator it(m_vertexEdgeMatrix, k); it; ++it)
        {
            auto r = it.row();
            auto c = it.col();
            auto eId = m_vertexEdgeMatrix.coeff(r, c);

            if (eId > 0) //This edge id is twice in the m_vertexEdgeMatrix, so it will be overwritten, is it ok?
            {
                m_edgeVertices[eId] = std::make_pair((uint)r, (uint)c);
            }
        }
    }
}

template<typename MeshType>
inline void DMB::MeshArrangement<MeshType>::buildDebugMesh()
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
            tPoint { m_coordinates[maVh0 * 3 + 0], m_coordinates[maVh0 * 3 + 1], m_coordinates[maVh0 * 3 + 2] },
            tPoint { m_coordinates[maVh1 * 3 + 0], m_coordinates[maVh1 * 3 + 1], m_coordinates[maVh1 * 3 + 2] },
            tPoint { m_coordinates[maVh2 * 3 + 0], m_coordinates[maVh2 * 3 + 1], m_coordinates[maVh2 * 3 + 2] },
        };

        auto vh0 = debugMesh.add_vertex(triangle[0]);
        auto vh1 = debugMesh.add_vertex(triangle[1]);
        auto vh2 = debugMesh.add_vertex(triangle[2]);

        auto fh = debugMesh.add_face(vh0, vh1, vh2);

        debugMesh.set_color(fh, { 192,192,192 });



        //if (m_intersectionEdgeFaceProp[tId])
        //{
        //    auto IOLabel = m_faceIOLabeling[tId];
        //    auto faceLabel = m_labels[tId];


        //    if (IOLabel.count() > 0)
        //    {
        //        debugMesh.set_color(fh, { 0,255,0 });
        //    }
        //    else
        //    {
        //        debugMesh.set_color(fh, { 255,0,0 });

        //    }


        ////    //if (faceLabel[0] == 1 && IOLabel[1] == 0 && IOLabel[0] == 0)
        ////    //{
        ////    //    debugMesh.set_color(fh, { 0,255,0 });
        ////    //}

        ////    //if (faceLabel[1] == 1 && IOLabel[0] == 1 && IOLabel[1] == 0)
        ////    //{
        ////    //    debugMesh.set_color(fh, { 0,255,0 });
        ////    //}

        //}

        //if (m_coplanarFace[tId])
        //{
        //    debugMesh.set_color(fh, { 0,0,255 });
        //}


        
    }

    //debugMesh.update_normals();
    //T3D_VDT_STORE_NAMED_MESH_AND_WAIT("", debugMesh);

    OpenMesh::IO::Options opt = OpenMesh::IO::Options::Default;
    //opt += OpenMesh::IO::Options::FaceColor;

    debugMesh.update_normals();
   // DMB::saveMesh(debugMesh, "C:/skola/PhD/Samples/booleans/debugMA.obj", opt);
}

template<typename MeshType>
template<typename F>
inline void DMB::MeshArrangement<MeshType>::debug_showMesh(const F& func)
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
            tPoint { m_coordinates[maVh0 * 3 + 0], m_coordinates[maVh0 * 3 + 1], m_coordinates[maVh0 * 3 + 2] },
            tPoint { m_coordinates[maVh1 * 3 + 0], m_coordinates[maVh1 * 3 + 1], m_coordinates[maVh1 * 3 + 2] },
            tPoint { m_coordinates[maVh2 * 3 + 0], m_coordinates[maVh2 * 3 + 1], m_coordinates[maVh2 * 3 + 2] },
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
