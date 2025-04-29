//TODO: License

#include "MeshOperand.h"

#include <OpenMesh/Core/Utils/PropertyManager.hh>


//TODO: for property names use HelperDefines.h, maybe rename it
template<typename MeshType>
inline DMB::MeshOperand<MeshType>::MeshOperand(int label) :
    m_meshOperandLabel(label)
{
    assert(label < NBIT);
}

template<typename MeshType>
inline void DMB::MeshOperand<MeshType>::processComponent(const MeshArrangement& ma, const tSparseIntMatrix& edgeAdjacencyMatrix)
{
    m_edgeValanceMatrix = tSparseIntMatrix(edgeAdjacencyMatrix.transpose()) + edgeAdjacencyMatrix;

    detectBoundaries();
    detectIntersectionEdges(ma);
    handleCoplanarFaces();
}

template<typename MeshType>
inline void DMB::MeshOperand<MeshType>::classifyComponents(const MeshArrangement& ma,const MeshType& otherMesh, const std::vector<tVertexHandle>& maToVh, int otherMeshLabel)
{
    OpenMesh::EProp<bool> pIntersectionEdge(m_mesh, m_pIntersectionEdge);
    OpenMesh::VProp<uint> pOperandToMeshArrangementIndex(m_mesh, m_pOperandToMeshArrangementIndex);
    OpenMesh::EProp<bool> pCoplanarEdge(m_mesh, m_pCoplanarEdge);
    OpenMesh::FProp<bool> pCoplanarFace(m_mesh, m_pCoplanarFace);

    OpenMesh::FProp< std::bitset<NBIT>> labeling(0, m_mesh, "bo_labeling");
    m_pLabeling = labeling.getRawProperty();
    OpenMesh::FProp< bool > intersectionFace(false, m_mesh, "bo_intersectionFace");
    m_pIntersectionFace = intersectionFace.getRawProperty();

    OpenMesh::VPropHandleT<uint> pOtherOperandToMAHandle;
    otherMesh.get_property_handle(pOtherOperandToMAHandle, "bo_operandToMeshArrangementIndex");


    auto isInside = [this, &otherMesh, &maToVh, &pOperandToMeshArrangementIndex, &pOtherOperandToMAHandle, &ma](OpenMesh::SmartHalfedgeHandle he) -> bool
    {
        auto vhToInvestigate = he.next().to();
        auto pToInvestigate = m_mesh.point(vhToInvestigate);

        auto heOtherMesh = otherMesh.find_halfedge(maToVh[pOperandToMeshArrangementIndex[he.from()]], maToVh[pOperandToMeshArrangementIndex[he.to()]]);

        if (heOtherMesh.is_boundary())
        {
            heOtherMesh = heOtherMesh.opp();
        }

        if (!heOtherMesh.is_valid())
        {
            std::cout << "neighbour he not found" << std::endl;
            assert(false);
            return false;
        }

        auto vh0 = heOtherMesh.from();
        auto vh1 = heOtherMesh.to();
        auto vh2 = heOtherMesh.next().to();

        auto p0 = otherMesh.point(vh0);
        auto p1 = otherMesh.point(vh1);
        auto p2 = otherMesh.point(vh2);

        genericPoint* p0Generic = ma.m_coordinatesImplicit[otherMesh.property(pOtherOperandToMAHandle, vh0)];
        genericPoint* p1Generic = ma.m_coordinatesImplicit[otherMesh.property(pOtherOperandToMAHandle, vh1)];
        genericPoint* p2Generic = ma.m_coordinatesImplicit[otherMesh.property(pOtherOperandToMAHandle, vh2)];
        genericPoint* pToInvestigateGeneric = ma.m_coordinatesImplicit[pOperandToMeshArrangementIndex[vhToInvestigate]];

        auto orientation = genericPoint::orient3D(*p0Generic, *p1Generic, *p2Generic, *pToInvestigateGeneric) * ma.m_multiplier;

        assert(orientation != 0); //there shouldn't be any coplanar faces

        if (orientation == 0)
        {
            std::cout << "not handled coplanar face" << std::endl;
            assert(false);
        }

        //sanity check
        if (!heOtherMesh.is_boundary())
        {
            auto vh3 = heOtherMesh.opp().next().to();
            auto p = otherMesh.point(vh3);

            genericPoint* p3Generic = ma.m_coordinatesImplicit[otherMesh.property(pOtherOperandToMAHandle, vh3)];
            auto orientationOpp = genericPoint::orient3D(*p1Generic, *p0Generic, *p3Generic, *pToInvestigateGeneric) * ma.m_multiplier;

            //TODO: this is not true - and edge can be lying in a triangle face, then the orientation is the same...
            /*
            *   \  / 
            * ___\/___
            * 
            */
            if (orientation != orientationOpp)
            {
                std::cout << "this is unresolved intersection" << std::endl;
                assert(false);
            }
        }

        return orientation < 0;
    };


    //TODO: parallel
    //for (auto eh : m_mesh.edges())
#pragma omp parallel for
    for (int i = 0; i < m_mesh.n_edges(); ++i)
    {
        auto eh = OpenMesh::make_smart(m_mesh.edge_handle(i), m_mesh);

        if (!eh.is_valid() || m_mesh.status(eh).deleted())
        {
            continue;
        }

        if (pIntersectionEdge[eh]/* && !pCoplanarEdge[eh]*/) //TODO: coplanar edges should not be skipped...
        {
            auto he = eh.h0();

            if (he.is_boundary())
            {
                he = eh.h1();
            }

            bool inside = isInside(he);

            std::bitset<NBIT> label;
            label[otherMeshLabel] = inside;

            labeling[he.face()] = label;
            intersectionFace[he.face()] = true;

            if (!eh.is_boundary())
            {
                bool oppInside = isInside(he.opp());

                label[otherMeshLabel] = oppInside;
                labeling[he.opp().face()] = label;
                intersectionFace[he.opp().face()] = true;

                if (oppInside == inside)
                {
                    std::cout << "this is not an intersection edge, it is actually non-manifold!" << std::endl;

                    pIntersectionEdge[eh] = false;
                    intersectionFace[he.face()] = false;
                    intersectionFace[he.opp().face()] = false;
                }

            }
        }
    }

    for (auto fh : m_mesh.faces())
    {
        m_mesh.set_color(fh, { 192,192,192 });

        if (intersectionFace[fh])
        {
            m_mesh.set_color(fh, { 255,0,0 });
        }
    }

    m_mesh.update_normals();
    T3D_VDT_STORE_NAMED_MESH_AND_WAIT("intersection faces before disconnect", m_mesh);

}

template<typename MeshType>
inline bool DMB::MeshOperand<MeshType>::disconnectComponents()
{
    OpenMesh::EProp<bool> pIntersectionEdge(m_mesh, m_pIntersectionEdge);
    OpenMesh::VProp<uint> pOperandToMeshArrangementIndex(m_mesh, m_pOperandToMeshArrangementIndex);
    OpenMesh::EProp<bool> pCoplanarEdge(m_mesh, m_pCoplanarEdge);
    OpenMesh::VProp<bool> coplanarVertex(m_mesh, m_pCoplanarVertex);

    OpenMesh::FProp< std::bitset<NBIT>> labeling(m_mesh, m_pLabeling);
    OpenMesh::FProp< bool > intersectionFace(m_mesh, m_pIntersectionFace);

    OpenMesh::FProp<uint> pFhToMaFh(m_mesh, m_pFhToMaFh);

    //reset tagged
    for (auto fh : m_mesh.faces())
    {
        m_mesh.status(fh).set_tagged(false);
    }

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

    //spread labeling
    for (auto vh : m_mesh.vertices())
    {
        if (intersectionVertex[vh] /*&& intersectionValance[vh] > 1*/)
        {
            auto edges = vh.edges_ccw().to_vector();

            //find first intersection face
            auto mid = std::find_if(edges.begin(), edges.end(),
                [&](const auto& eh)
                {
                    return pIntersectionEdge[eh];
                });

            if (mid == edges.end())
            {
                continue;
            }

            //vector sorted in such a way, that an intersection face is first
            std::rotate(edges.begin(), mid, edges.end());

            std::bitset<NBIT> currentLabel;
            bool labelValid = true;
            std::vector<OpenMesh::SmartFaceHandle> unlabeledFaces;

            for (auto eh : edges)
            {
                if (pIntersectionEdge[eh])
                {
                    continue;
                }

                OpenMesh::FaceHandle f0, f1;

                if (!eh.h0().is_boundary())
                {
                    f0 = eh.h0().face();
                }

                if (!eh.h1().is_boundary())
                {
                    f1 = eh.h1().face();
                }

                labelValid = false;

                if (f0.is_valid() && intersectionFace[f0])
                {
                    currentLabel = labeling[f0];
                    labelValid = true;
                }

                if (f1.is_valid() && intersectionFace[f1])
                {
                    currentLabel = labeling[f1];
                    labelValid = true;
                }

                if (labelValid)
                {
                    unlabeledFaces.push_back(OpenMesh::make_smart(f0, m_mesh));
                    unlabeledFaces.push_back(OpenMesh::make_smart(f1, m_mesh));

                    for (auto unlabeledFace : unlabeledFaces)
                    {
                        if (!unlabeledFace.is_valid())
                        {
                            continue;
                        }

                        if (unlabeledFace.tagged())
                        {
                            continue;
                        }

                        labeling[unlabeledFace] = currentLabel;
                        intersectionFace[unlabeledFace] = true;
                        m_mesh.status(unlabeledFace).set_tagged(true);
                    }

                    unlabeledFaces.clear();
                }
                else
                {
                    unlabeledFaces.push_back(OpenMesh::make_smart(f0, m_mesh));
                    unlabeledFaces.push_back(OpenMesh::make_smart(f1, m_mesh));
                }
            }

            auto eh = *edges.begin();
            if (eh.h0().face().is_valid() && intersectionFace[eh.h0().face()])
            {
                currentLabel = labeling[eh.h0().face()];
            }

            if (eh.h1().face().is_valid() && intersectionFace[eh.h1().face()])
            {
                currentLabel = labeling[eh.h1().face()];
            }

            for (auto unlabeledFace : unlabeledFaces)
            {
                if (!unlabeledFace.is_valid())
                {
                    continue;
                }

                labeling[unlabeledFace] = currentLabel;
                intersectionFace[unlabeledFace] = true;
                m_mesh.status(unlabeledFace).set_tagged(true);
            }
        }
    }

    OpenMesh::VProp<int> intersectionValance(0, m_mesh);

    for (auto vh : m_mesh.vertices())
    {
        if (intersectionVertex[vh])
        {
            int val = 0;
            for (auto vvh : vh.vertices())
            {
                if (intersectionVertex[vvh] || coplanarVertex[vvh] || vvh.is_boundary())
                {
                    ++val;
                }
            }

            intersectionValance[vh] = val;
        }
    }

    //reset tagged
    for (auto fh : m_mesh.faces())
    {
        m_mesh.status(fh).set_tagged(false);
    }

    std::vector<std::vector<tFaceHandle>> components;

    //find components
    for (auto fh : m_mesh.faces())
    {
        //found a new component seed face - get all connected intersection faces with the same labeling
        if (intersectionFace[fh] && !fh.tagged())
        {
            //seed label
            auto seedLabel = labeling[fh];

            std::queue<tFaceHandle> facesToVisit;
            std::vector<tFaceHandle> componentFaces;

            facesToVisit.push(fh);

            while (!facesToVisit.empty())
            {
                //get top face
                auto topFace = OpenMesh::make_smart(facesToVisit.front(), m_mesh);
                facesToVisit.pop();

                if (topFace.tagged())
                {
                    continue;
                }

                //store new component face
                componentFaces.push_back(topFace);

                //continue with adjacent face with the same labeling
                for (auto ffh : topFace.faces())
                {
                    if (!ffh.tagged() && labeling[ffh] == seedLabel)
                    {
                        facesToVisit.push(ffh);
                    }
                }

                //mark this face as processed
                m_mesh.status(topFace).set_tagged(true);
            }

            //store this component
            components.push_back(componentFaces);
        }
    }

    //reset tagged
    for (auto fh : m_mesh.faces())
    {
        m_mesh.status(fh).set_tagged(false);
    }

    OpenMesh::VProp<tVertexHandle> duplicatedVertex(m_mesh);

    //disconnect components
    for (auto component : components)
    {
        std::vector<tVertexHandle> verticesToDuplicate;
        //gather component vertices to duplicate
        for (auto fh : component)
        {
            for (auto vh : OpenMesh::make_smart(fh, m_mesh).vertices())
            {
                if (intersectionVertex[vh])
                {
                    verticesToDuplicate.push_back(vh);
                }
            }
        }

        DMB::makeUniqueVector(verticesToDuplicate);

        //duplicate vertices
        for (auto vh : verticesToDuplicate)
        {
            auto p = m_mesh.point(vh);
            auto newVh = m_mesh.add_vertex(p);
            duplicatedVertex[vh] = newVh;
            pOperandToMeshArrangementIndex[newVh] = pOperandToMeshArrangementIndex[vh];
            intersectionValance[newVh] = intersectionValance[vh];

        }

        //disconnect faces
        for (auto fh : component)
        {
            //face properties to copy
            auto label = labeling[fh];
            bool isIntersectionFace = intersectionFace[fh];
            uint maFh = pFhToMaFh[fh];

            std::vector<tVertexHandle> fVertices;

            for (auto vh : OpenMesh::make_smart(fh, m_mesh).vertices())
            {
                if (intersectionVertex[vh])
                {
                    auto newVh = duplicatedVertex[vh];

                    if (!newVh.is_valid())
                    {
                        std::cout << "could not add new vertex, not possible" << std::endl;
                        assert(false);
                    }

                    fVertices.push_back(newVh);
                }
                else
                {
                    fVertices.push_back(vh);
                }
            }

            m_mesh.delete_face(fh, false);
            auto newFh = m_mesh.add_face(fVertices);

            if (!newFh.is_valid())
            {
                std::cout << "disconnect failed" << std::endl;
                return false;
            }

            intersectionFace[newFh] = isIntersectionFace;
            labeling[newFh] = label;
            pFhToMaFh[newFh] = maFh;

        }

        for (auto vh : verticesToDuplicate)
        {
            intersectionVertex[duplicatedVertex[vh]] = intersectionVertex[vh];
        }

    }

    for (auto vh : m_mesh.vertices())
    {
        m_mesh.set_color(vh, { 0,0,0 });
    }
    for (auto fh : m_mesh.faces())
    {
        m_mesh.set_color(fh, { 192,192,192 });
    }

    for (auto fh : m_mesh.faces())
    {
        if (intersectionFace[fh])
        {
            m_mesh.set_color(fh, { 255,0,0 });
        }
    }

    m_mesh.update_normals();
    T3D_VDT_STORE_NAMED_MESH_AND_WAIT("intersection faces after disconnect", m_mesh);

    std::bitset<NBIT> inA;
    inA[0] = 1;

    std::bitset<NBIT> inB;
    inB[1] = 1;

    for (auto fh : m_mesh.faces())
    {
        if (intersectionFace[fh])
        {
            auto label = labeling[fh];

            if ((label & inA) != 0)
            {
                m_mesh.set_color(fh, { 0,255,0 });
            }

            if ((label & inB) != 0)
            {
                m_mesh.set_color(fh, { 0,0,255 });
            }

        }
    }

    m_mesh.update_normals();
    T3D_VDT_STORE_NAMED_MESH_AND_WAIT("", m_mesh);

    std::vector < std::vector<tFaceHandle> > componentsToDelete;

    //now spread the labeling information over the connected components
    for (auto fh : m_mesh.faces())
    {
        //found a new component seed face - get all connected faces
        if (intersectionFace[fh] && !fh.tagged())
        {
            //seed label
            auto seedLabel = labeling[fh];

            std::queue<tFaceHandle> facesToVisit;
            std::vector<tFaceHandle> componentFaces;

            facesToVisit.push(fh);

            bool conflictingLabeling = false;
            bool unclosedCurve = false;

            while (!facesToVisit.empty())
            {
                //get top face
                auto topFace = OpenMesh::make_smart(facesToVisit.front(), m_mesh);
                facesToVisit.pop();

                if (topFace.tagged())
                {
                    continue;
                }

                componentFaces.push_back(topFace);

                if (intersectionFace[topFace])
                {
                    if (labeling[topFace] != seedLabel)
                    {
                        conflictingLabeling = true;

                        std::cout << "seedLabel: " << seedLabel << std::endl;
                        std::cout << "labeling[topFace]: " << labeling[topFace] << std::endl;

                        m_mesh.set_color(fh, { 255,0,0 });
                        m_mesh.set_color(topFace, { 0,255,0 });
                        m_mesh.update_normals();
                        T3D_VDT_STORE_NAMED_MESH_AND_WAIT("", m_mesh);

                        for (auto dbfh : m_mesh.faces())
                        {
                            if (intersectionFace[dbfh])
                            {
                                auto label = labeling[dbfh];

                                if ((label & inA) != 0)
                                {
                                    m_mesh.set_color(dbfh, { 0,255,0 });
                                }

                                if ((label & inB) != 0)
                                {
                                    m_mesh.set_color(dbfh, { 0,0,255 });
                                }

                            }
                        }

                        m_mesh.update_normals();
                        T3D_VDT_STORE_NAMED_MESH_AND_WAIT("", m_mesh);

                    }
                }
                else
                {
                    labeling[topFace] = seedLabel;
                }

                for (auto vh : topFace.vertices())
                {
                    if (intersectionValance[vh] == 1)
                    {
                        unclosedCurve = true;

                        m_mesh.set_color(vh, { 255,0,0 });
                        m_mesh.set_color(topFace, { 255,0,0 });

                        m_mesh.update_normals();
                        T3D_VDT_STORE_NAMED_MESH_AND_WAIT("", m_mesh);

                        MeshType debugMesh;

                        std::vector<tVertexHandle> vertices;
                        for (auto fvh : topFace.vertices())
                        {
                            m_mesh.set_color(fvh, { 255,0,0 });

                            vertices.push_back(debugMesh.add_vertex(m_mesh.point(fvh)));
                        }
                        debugMesh.add_face(vertices);
                        T3D_VDT_STORE_NAMED_MESH_AND_WAIT("", debugMesh);

                    }
                }

                //continue with adjacent face with the same labeling
                for (auto ffh : topFace.faces())
                {
                    if (!ffh.tagged())
                    {
                        facesToVisit.push(ffh);
                    }
                }

                //mark this face as processed
                m_mesh.status(topFace).set_tagged(true);
            }

            if (conflictingLabeling || unclosedCurve)
            {
                componentsToDelete.push_back(componentFaces);

                for (auto _fh : m_mesh.faces())
                {
                    m_mesh.set_color(_fh, { 192,192,192 });
                }

                for (auto _fh : componentFaces)
                {
                    m_mesh.set_color(_fh, { 255,0,0 });
                }
                m_mesh.update_normals();
                T3D_VDT_STORE_NAMED_MESH_AND_WAIT("", m_mesh);

            }
        }
    }

    //reset tagged
    for (auto fh : m_mesh.faces())
    {
        m_mesh.status(fh).set_tagged(false);
    }

    //delete unusable components
    for (auto component : componentsToDelete)
    {
        for (auto fh : component)
        {
            m_mesh.delete_face(fh);
        }
    }

    return true;
}

template<typename MeshType>
template<typename F>
inline void DMB::MeshOperand<MeshType>::classifyMeshArrangement(MeshArrangement& ma, const F& predicate)
{
    OpenMesh::FProp<std::bitset<NBIT>> labeling(m_mesh, m_pLabeling);
    OpenMesh::VProp<uint> operandToMeshArrangementIndex(m_mesh, m_pOperandToMeshArrangementIndex);
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
        auto predicateResult = predicate(m_meshOperandLabel, faceLabel);

        ma.m_boPredicates[pFhToMaFh[fh]] = predicateResult;
    }

    for (auto coplanarTri : m_coplanarTriangles)
    {
        //we support binary booleans only, that's why I can do it here, but it should be a list of faces to which this face belong in variadic booleans
        std::bitset<NBIT> faceLabel;
        faceLabel[0] = 1;
        faceLabel[1] = 1;

        auto predicateResult = predicate(m_meshOperandLabel, faceLabel);
        ma.m_boPredicates[coplanarTri.first] = predicateResult;
    }

}

template<typename MeshType>
inline bool DMB::MeshOperand<MeshType>::combineToResult(MeshArrangement& ma, MeshType& result)
{
    std::vector<tTriplet> tripletList;
    tripletList.reserve(std::size(ma.m_coordinates));

    for (int i = 0; i < std::size(ma.m_triangles); i += 3)
    {
        /*
         * predicateResult ==  0 -> not part of result
         * predicateResult ==  1 -> part of result as is
         * predicateResult == -1 -> part of result but flipped
        */
        uint tId = i / 3;

        auto predicateResult = ma.m_boPredicates[tId];

        if (predicateResult == 0)
        {
            continue;
        }

        uint v0 = ma.m_triangles[i + 0];
        uint v1 = ma.m_triangles[i + 1];
        uint v2 = ma.m_triangles[i + 2];

        tripletList.push_back(tTriplet(v0, v1, -1));
        tripletList.push_back(tTriplet(v1, v2, -1));
        tripletList.push_back(tTriplet(v2, v0, -1));
    }

    tSparseIntMatrix edges(std::size(ma.m_coordinates) / 3, std::size(ma.m_coordinates) / 3);
    edges.setFromTriplets(tripletList.begin(), tripletList.end());

    int edgeId = 1;
    int cnt = 0;

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
                ++cnt;
            }
        }
    }

    std::vector<std::vector<uint>> edgeFaceAdjacency;
    edgeFaceAdjacency.resize(cnt + 1);


    for (uint i = 0; i < std::size(ma.m_triangles); i += 3)
    {
        /*
         * predicateResult ==  0 -> not part of result
         * predicateResult ==  1 -> part of result as is
         * predicateResult == -1 -> part of result but flipped
        */
        uint tId = i / 3;

        auto predicateResult = ma.m_boPredicates[tId];

        if (predicateResult == 0)
        {
            continue;
        }

        uint v0 = ma.m_triangles[i + 0];
        uint v1 = ma.m_triangles[i + 1];
        uint v2 = ma.m_triangles[i + 2];

        int e0 = edges.coeff(v0, v1);
        int e1 = edges.coeff(v1, v2);
        int e2 = edges.coeff(v2, v0);

        edgeFaceAdjacency[e0].push_back(tId);
        edgeFaceAdjacency[e1].push_back(tId);
        edgeFaceAdjacency[e2].push_back(tId);
    }


    OpenMesh::FProp<std::bitset<NBIT>> faceOriginLabel({}, result, "bo_faceOriginLabel");

    std::vector<bool> faceVisited(std::size(ma.m_labels), false);
    std::vector<typename MeshType::VertexHandle> meshArrangementToOutputIndex(std::size(ma.m_coordinates) / 3);
    std::vector<int> vertexValance;
    vertexValance.resize(std::size(ma.m_coordinates) / 3, -1);

    std::vector<bool> taggedVh;
    taggedVh.resize(std::size(ma.m_coordinates) / 3, false);

    std::vector<bool> taggedEh;
    taggedEh.resize(std::size(edgeFaceAdjacency), false);

    //using map, don't expect that there will be many non-manifold vertices in the result
    std::map < uint, std::vector<tVertexHandle>> duplicatedVertices;

    bool success = true;

    //TODO: this is bullsh*t...
    auto getConnectedFaces = [&edgeFaceAdjacency, &edges, &taggedEh, &taggedVh, &vertexValance, &ma](uint vh, uint fh, std::queue<uint>& trianglesToProcess)
    {
        if (taggedVh[vh])
        {
            return;
        }

        uint maVh0 = ma.m_triangles[fh * 3 + 0];
        uint maVh1 = ma.m_triangles[fh * 3 + 1];
        uint maVh2 = ma.m_triangles[fh * 3 + 2];

        std::queue<int> edgesToCheck;

        if(maVh0 == vh || maVh1 == vh)
            edgesToCheck.push(edges.coeff(maVh0, maVh1));

        if (maVh1 == vh || maVh2 == vh)
            edgesToCheck.push(edges.coeff(maVh1, maVh2));

        if (maVh0 == vh || maVh2 == vh)
            edgesToCheck.push(edges.coeff(maVh2, maVh0));

        std::vector<int> visitedEdges;

        //1, because this triangle is definitely connected to this vertex...
        int currenctVertexValance = 1;
        std::vector<uint> connectedFaces;

        while (!edgesToCheck.empty())
        {
            int eh = edgesToCheck.front();
            edgesToCheck.pop();

            connectedFaces.insert(connectedFaces.end(), edgeFaceAdjacency[eh].begin(), edgeFaceAdjacency[eh].end());

            if (taggedEh[eh])
            {
                continue;
            }

            for (auto adjacentTriangle : edgeFaceAdjacency[eh])
            {
                auto adjacentTrianglePredicateResult = ma.boPredicates[adjacentTriangle];

                if (adjacentTriangle != fh && adjacentTrianglePredicateResult != 0)
                {
                    trianglesToProcess.push(adjacentTriangle);
                    ++currenctVertexValance;

                    maVh0 = ma.m_triangles[adjacentTriangle * 3 + 0];
                    maVh1 = ma.m_triangles[adjacentTriangle * 3 + 1];
                    maVh2 = ma.m_triangles[adjacentTriangle * 3 + 2];

                    if (maVh0 == vh || maVh1 == vh)
                        edgesToCheck.push(edges.coeff(maVh0, maVh1));

                    if (maVh1 == vh || maVh2 == vh)
                        edgesToCheck.push(edges.coeff(maVh1, maVh2));

                    if (maVh0 == vh || maVh2 == vh)
                        edgesToCheck.push(edges.coeff(maVh2, maVh0));

                }
            }

            taggedEh[eh] = true;
            visitedEdges.push_back(eh);
        }

        for (auto eh : visitedEdges)
        {
            taggedEh[eh] = false;
        }

        taggedVh[vh] = true;

        DMB::makeUniqueVector(connectedFaces);
        vertexValance[vh] = connectedFaces.size(); 
    };


    for (uint i = 0; i < std::size(ma.m_triangles); i += 3)
    {
        if (!success)
        {
            break;
        }

        uint tId = i / 3;
        auto predicateResult = ma.m_boPredicates[tId];

        if (predicateResult == 0)
        {
            continue;
        }

        
        if (faceVisited[tId])
        {
            continue;
        }

        std::queue<uint> trianglesToProcess;
        trianglesToProcess.push(tId);

        //todo: invalidate all vertices?
        //tVertexHandle invalidVh;
        //invalidVh.invalidate();
        //meshArrangementToOutputIndex.resize(std::size(ma.m_coordinates) / 3, invalidVh);

        while (!trianglesToProcess.empty())
        {
            uint currentTriangleId = trianglesToProcess.front();
            trianglesToProcess.pop();

            auto currentPredicateResult = ma.m_boPredicates[currentTriangleId];
            auto faceLabel = ma.m_labels[currentTriangleId];

            if (currentPredicateResult == 0 || faceVisited[currentTriangleId])
            {
                faceVisited[currentTriangleId] = true;
                continue;
            }

            uint maVh0 = ma.m_triangles[currentTriangleId * 3 + 0];
            uint maVh1 = ma.m_triangles[currentTriangleId * 3 + 1];
            uint maVh2 = ma.m_triangles[currentTriangleId * 3 + 2];

            tTriangle triangle =
            {
                tPoint { ma.m_coordinates[maVh0 * 3 + 0], ma.m_coordinates[maVh0 * 3 + 1], ma.m_coordinates[maVh0 * 3 + 2] },
                tPoint { ma.m_coordinates[maVh1 * 3 + 0], ma.m_coordinates[maVh1 * 3 + 1], ma.m_coordinates[maVh1 * 3 + 2] },
                tPoint { ma.m_coordinates[maVh2 * 3 + 0], ma.m_coordinates[maVh2 * 3 + 1], ma.m_coordinates[maVh2 * 3 + 2] },
            };

            int j = 0;
            for (auto maVh : { maVh0, maVh1, maVh2 })
            {
                getConnectedFaces(maVh, currentTriangleId, trianglesToProcess);

                auto vh = meshArrangementToOutputIndex[maVh];
                if (!vh.is_valid())
                {
                    vh = result.add_vertex(triangle[j]);
                    meshArrangementToOutputIndex[maVh] = vh;

                }
                else //vertex was already referenced
                {
                    //duplicate non-manifold vertex
                    if (vertexValance[maVh] == 0)
                    {
                        duplicatedVertices[maVh].push_back(vh);

                        //std::cout << maVh << " " << vh.idx() << " actually found non manifold vertex, valance: " << vertexValance[maVh] << " is boundary " << result.is_boundary(vh) << std::endl;
                        vh = result.add_vertex(triangle[j]);
                        meshArrangementToOutputIndex[maVh] = vh;
                        taggedVh[maVh] = false;
                    }
                }


                //vertex is referenced by this triangle
                vertexValance[maVh] = vertexValance[maVh] - 1;
                
                ++j;

            }

            auto vh0 = meshArrangementToOutputIndex[maVh0];
            auto vh1 = meshArrangementToOutputIndex[maVh1];
            auto vh2 = meshArrangementToOutputIndex[maVh2];

            tFaceHandle resultFh;


            if (currentPredicateResult > 0)
            {
                resultFh = result.add_face(vh0, vh1, vh2);

                if (!resultFh.is_valid())
                {
                    std::cout << "Could not add face to result" << std::endl;
                    success = false;
                    break;
                }

            }
            else if (currentPredicateResult < 0)
            {
                resultFh = result.add_face(vh0, vh2, vh1);

                if (!resultFh.is_valid())
                {
                    std::cout << "Could not add face to result" << std::endl;
                    success = false;
                    break;
                }
            }

            faceOriginLabel[resultFh] = faceLabel;

            //this triangle was processed
            faceVisited[currentTriangleId] = true;
        }
    }

    if (!success)
    {
        return success;
    }

    //TODO: this is not tested!
    //try to merge duplicated vertices, if possible
    for(auto const& [key, vertices] : duplicatedVertices)
    {
        uint maVh = key;

        std::vector<tVertexHandle> verticesToMerge;

        auto currentVh = meshArrangementToOutputIndex[maVh];
        if (result.is_boundary(currentVh) && !result.status(currentVh).deleted())
        {
            verticesToMerge.push_back(currentVh);
        }

        for (auto vh : vertices)
        {
            if (result.is_boundary(vh) && !result.status(vh).deleted())
            {
                verticesToMerge.push_back(vh);
            }
        }

        if (verticesToMerge.empty())
        {
            continue;
        }

        auto p = result.point(verticesToMerge.back());
        auto newVh = result.add_vertex(p);


        for (auto vh : verticesToMerge)
        {
            std::vector<OpenMesh::SmartFaceHandle> facesToDelete = OpenMesh::make_smart(vh, result).faces().to_vector();

            for (auto fh : facesToDelete)
            {
                std::vector<tVertexHandle> fVertices;

                for (auto fvh : fh.vertices())
                {
                    if (fvh == vh)
                    {
                        fVertices.push_back(newVh);
                    }
                    else
                    {
                        fVertices.push_back(fvh);
                    }
                }

                result.delete_face(fh, false);

                auto newFh = result.add_face(fVertices);

                if (!newFh.is_valid())
                {
                    std::cout << "merging vertices in result failed" << std::endl;
                    assert(false);
                }

            }
        }

    }

    return success;
}

template<typename MeshType>
inline bool DMB::MeshOperand<MeshType>::createMesh(MeshArrangement& ma)
{
    //TODO: handle coplanar faces as static for this class, not in one of the operands

    OpenMesh::VProp<uint> operandToMeshArrangementIndex(0, m_mesh, "bo_operandToMeshArrangementIndex");
    m_pOperandToMeshArrangementIndex = operandToMeshArrangementIndex.getRawProperty();

    OpenMesh::VProp<tVertexHandle> nonManifoldVertexMap(m_mesh);

    OpenMesh::FProp<bool> coplanarFace(false, m_mesh, "bo_coplanarFace");
    m_pCoplanarFace = coplanarFace.getRawProperty();
    OpenMesh::FProp<uint> pFhToMaFh(0, m_mesh, "bo_fhToMaFh");
    m_pFhToMaFh = pFhToMaFh.getRawProperty();

    m_meshArrangementToVhs.resize(std::size(ma.m_coordinates) / 3);

    std::bitset<NBIT> meshArrangementLabel;
    meshArrangementLabel[m_meshOperandLabel] = 1;
    std::bitset<NBIT> meshArrangementLabelComplement = meshArrangementLabel;
    meshArrangementLabelComplement.flip();

    std::vector<bool> faceVisited(std::size(ma.m_labels), false);

    tSparseIntMatrix edges(std::size(ma.m_coordinates) / 3, std::size(ma.m_coordinates) / 3);
    std::vector<std::vector<uint>> edgeFaceAdjacency;
    {
        std::vector<tTriplet> tripletList;
        tripletList.reserve(std::size(ma.m_coordinates));


        for (int i = 0; i < std::size(ma.m_triangles); i += 3)
        {
            auto faceLabel = ma.m_labels[i / 3];

            bool belongsToThisOperand = (faceLabel & meshArrangementLabel) != 0;

            if (!belongsToThisOperand)
            {
                continue;
            }

            uint v0 = ma.m_triangles[i + 0];
            uint v1 = ma.m_triangles[i + 1];
            uint v2 = ma.m_triangles[i + 2];

            tripletList.push_back(tTriplet(v0, v1, -1));
            tripletList.push_back(tTriplet(v1, v2, -1));
            tripletList.push_back(tTriplet(v2, v0, -1));
        }

        edges.setFromTriplets(tripletList.begin(), tripletList.end());

        int edgeId = 1;
        int cnt = 0;

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
                    ++cnt;
                }
            }
        }

        edgeFaceAdjacency.resize(cnt + 1);

        for (uint i = 0; i < std::size(ma.m_triangles); i += 3)
        {
            uint tId = i / 3;
            auto faceLabel = ma.m_labels[tId];

            bool belongsToThisOperand = (faceLabel & meshArrangementLabel) != 0;

            if (!belongsToThisOperand)
            {
                continue;
            }

            uint v0 = ma.m_triangles[i + 0];
            uint v1 = ma.m_triangles[i + 1];
            uint v2 = ma.m_triangles[i + 2];

            int e0 = edges.coeff(v0, v1);
            int e1 = edges.coeff(v1, v2);
            int e2 = edges.coeff(v2, v0);

            edgeFaceAdjacency[e0].push_back(tId);
            edgeFaceAdjacency[e1].push_back(tId);
            edgeFaceAdjacency[e2].push_back(tId);
        }
    }

#if 1
    {
        MeshType wholeOperand;

        for (int i = 0; i < std::size(ma.m_triangles); i += 3)
        {
            uint tId = i / 3;
            auto faceLabel = ma.m_labels[tId];

            bool belongsToThisOperand = (faceLabel & meshArrangementLabel) != 0;

            if (!belongsToThisOperand)
            {
                continue;
            }


            uint maVh0 = ma.m_triangles[i + 0];
            uint maVh1 = ma.m_triangles[i + 1];
            uint maVh2 = ma.m_triangles[i + 2];

            tTriangle triangle =
            {
                tPoint { ma.m_coordinates[maVh0 * 3 + 0], ma.m_coordinates[maVh0 * 3 + 1], ma.m_coordinates[maVh0 * 3 + 2] },
                tPoint { ma.m_coordinates[maVh1 * 3 + 0], ma.m_coordinates[maVh1 * 3 + 1], ma.m_coordinates[maVh1 * 3 + 2] },
                tPoint { ma.m_coordinates[maVh2 * 3 + 0], ma.m_coordinates[maVh2 * 3 + 1], ma.m_coordinates[maVh2 * 3 + 2] },
            };

            auto vh0 = wholeOperand.add_vertex(triangle[0]);
            auto vh1 = wholeOperand.add_vertex(triangle[1]);
            auto vh2 = wholeOperand.add_vertex(triangle[2]);
            auto fh = wholeOperand.add_face(vh0, vh1, vh2);
            //wholeOperand.set_color(fh, { 192,192,192 });
        }

        wholeOperand.update_normals();
        DMB::saveMesh(wholeOperand, "C:/T3D/Samples/Booleans/wholeOperandBeforeEdges.stl");
    }

    std::vector<uint> nonManifoldFaces;
    std::unordered_set<uint> nonManifoldEdgeVertices;
    {
        //get faces and vertices adjacent to non-manifold edges
        uint eId = 0;
        for (auto efAdj : edgeFaceAdjacency)
        {
            if (efAdj.size() > 2)
            {
                nonManifoldFaces.insert(nonManifoldFaces.end(), efAdj.begin(), efAdj.end());
                //std::cout << "this is non manifold edge in ajd " << std::endl;

                for (auto tId : efAdj)
                {
                    //std::cout << "non manifold face: " << tId << std::endl;

                    uint v0 = ma.m_triangles[tId * 3 + 0];
                    uint v1 = ma.m_triangles[tId * 3 + 1];
                    uint v2 = ma.m_triangles[tId * 3 + 2];

                    uint e0 = edges.coeff(v0, v1);
                    uint e1 = edges.coeff(v1, v2);
                    uint e2 = edges.coeff(v2, v0);
                    //std::cout << e0 << " " << e1 << " " << e2 << " " << std::endl;
                    if (eId == e0)
                    {
                        nonManifoldEdgeVertices.insert(v0);
                        nonManifoldEdgeVertices.insert(v1);
                    }

                    if (eId == e1)
                    {
                        nonManifoldEdgeVertices.insert(v1);
                        nonManifoldEdgeVertices.insert(v2);
                    }

                    if (eId == e2)
                    {
                        nonManifoldEdgeVertices.insert(v2);
                        nonManifoldEdgeVertices.insert(v0);
                    }
                }

            }

            ++eId;
        }

        //get all faces adjacent to non-manifold vertices from previous step
        //TODO: need vertex face connectivity for this
        for (uint i = 0; i < std::size(ma.m_triangles); i += 3)
        {
            uint tId = i / 3;
            auto faceLabel = ma.m_labels[tId];

            bool belongsToThisOperand = (faceLabel & meshArrangementLabel) != 0;

            if (!belongsToThisOperand)
            {
                continue;
            }

            uint maVh0 = ma.m_triangles[tId * 3 + 0];
            uint maVh1 = ma.m_triangles[tId * 3 + 1];
            uint maVh2 = ma.m_triangles[tId * 3 + 2];

            for (auto maVh : { maVh0, maVh1, maVh2 })
            {
                if (std::find(nonManifoldEdgeVertices.begin(), nonManifoldEdgeVertices.end(), maVh) != nonManifoldEdgeVertices.end())
                {
                    nonManifoldFaces.push_back(tId);
                }
            }
        }

        DMB::makeUniqueVector(nonManifoldFaces);
    }


    {
        OpenMesh::VProp<bool> nonManifoldVertexProp(false, m_mesh);
        OpenMesh::FProp<uint> operandFhToMaTId(0, m_mesh);

        for (auto tId : nonManifoldFaces)
        {
            uint maVh0 = ma.m_triangles[tId * 3 + 0];
            uint maVh1 = ma.m_triangles[tId * 3 + 1];
            uint maVh2 = ma.m_triangles[tId * 3 + 2];

            tTriangle triangle =
            {
                tPoint { ma.m_coordinates[maVh0 * 3 + 0], ma.m_coordinates[maVh0 * 3 + 1], ma.m_coordinates[maVh0 * 3 + 2] },
                tPoint { ma.m_coordinates[maVh1 * 3 + 0], ma.m_coordinates[maVh1 * 3 + 1], ma.m_coordinates[maVh1 * 3 + 2] },
                tPoint { ma.m_coordinates[maVh2 * 3 + 0], ma.m_coordinates[maVh2 * 3 + 1], ma.m_coordinates[maVh2 * 3 + 2] },
            };

            int j = 0;
            for (auto maVh : { maVh0, maVh1, maVh2 })
            {
                //we don't have this vertex yet
                if (!m_meshArrangementToVhs[maVh].is_valid())
                {
                    auto newVh = m_mesh.add_vertex(triangle[j]);
                    m_meshArrangementToVhs[maVh] = newVh;

                    nonManifoldVertexProp[newVh] = false;

                    //this is non manifold vh
                    if (std::find(nonManifoldEdgeVertices.begin(), nonManifoldEdgeVertices.end(), maVh) != nonManifoldEdgeVertices.end())
                    {
                        nonManifoldVertexProp[newVh] = true;
                    }
                }
                //we have this vertex, but it is non-manifold
                else if (std::find(nonManifoldEdgeVertices.begin(), nonManifoldEdgeVertices.end(), maVh) != nonManifoldEdgeVertices.end())
                {
                    auto newVh = m_mesh.add_vertex(triangle[j]);
                    m_meshArrangementToVhs[maVh] = newVh;
                    nonManifoldVertexProp[newVh] = true;
                }

                ++j;
            }

            auto vh0 = m_meshArrangementToVhs[maVh0];
            auto vh1 = m_meshArrangementToVhs[maVh1];
            auto vh2 = m_meshArrangementToVhs[maVh2];

            operandToMeshArrangementIndex[vh0] = maVh0;
            operandToMeshArrangementIndex[vh1] = maVh1;
            operandToMeshArrangementIndex[vh2] = maVh2;

            tFaceHandle fh = m_mesh.add_face(vh0, vh1, vh2);
            //m_mesh.set_color(fh, { 192, 192, 192 });

            operandFhToMaTId[fh] = tId;
        }

        //DMB::saveMesh(m_mesh, "C:/T3D/Samples/Booleans/explodedMesh.obj");

        OpenMesh::EProp<bool> edgesToMerge(false, m_mesh);

        tVertexHandle invalidVh;
        invalidVh.invalidate();
        OpenMesh::VProp<tVertexHandle> replaceVertex(invalidVh, m_mesh);
        OpenMesh::VProp<tVertexHandle> targetVertex(invalidVh, m_mesh);
        
        std::vector<std::vector<tVertexHandle>> facesToAdd;

        std::vector< std::vector<tFaceHandle> > components;
        for (auto vh : m_mesh.vertices())
        {
            if (!nonManifoldVertexProp[vh] || vh.tagged())
            {
                continue;
            }

            std::vector<tFaceHandle> componentFaces;
            std::queue<OpenMesh::SmartVertexHandle> verticesToVisit;
            verticesToVisit.push(vh);

            while (!verticesToVisit.empty())
            {
                auto currentVh = verticesToVisit.front();
                verticesToVisit.pop();

                if (currentVh.tagged())
                {
                    continue;
                }

                for (auto vfh : currentVh.faces())
                {
                    componentFaces.push_back(vfh);

                    for (auto vvh : vfh.vertices())
                    {
                        verticesToVisit.push(vvh);
                    }
                }

                //this vertex was processed
                m_mesh.status(currentVh).set_tagged(true);

            }

            DMB::makeUniqueVector(componentFaces);

            components.push_back(componentFaces);
        }

        //reconnect components
        for (auto componentFaces : components)
        {
            //find originaly non-manifold vertices
            std::unordered_set<uint> nonManifoldMAVertices;

            for (auto fh : componentFaces)
            {
                for (auto fvh : OpenMesh::make_smart(fh, m_mesh).vertices())
                {
                    if (nonManifoldVertexProp[fvh])
                    {
                        nonManifoldMAVertices.insert(operandToMeshArrangementIndex[fvh]);
                    }
                }
            }

            std::map<uint, tVertexHandle> maToReplaceVertex;

            //add new replace vertices
            for (uint maVh : nonManifoldMAVertices)
            {
                tPoint p(ma.m_coordinates[maVh * 3 + 0], ma.m_coordinates[maVh * 3 + 1], ma.m_coordinates[maVh * 3 + 2]);
                auto newVh = m_mesh.add_vertex(p);
                maToReplaceVertex[maVh] = newVh;

                //add this new vertex to MA as well
                uint newMaVh = ma.m_coordinatesImplicit.size();

                auto* gp = ma.m_coordinatesImplicit[maVh];
                ma.m_coordinatesImplicit.push_back(gp);

                ma.m_coordinates.push_back(ma.m_coordinates[maVh * 3 + 0]);
                ma.m_coordinates.push_back(ma.m_coordinates[maVh * 3 + 1]);
                ma.m_coordinates.push_back(ma.m_coordinates[maVh * 3 + 2]);

                m_meshArrangementToVhs.push_back(newVh);
                operandToMeshArrangementIndex[newVh] = newMaVh;
            }

            std::map<uint, std::vector<tVertexHandle>> newFaces;
            std::vector<tFaceHandle> facesToDelete;

            for (auto fh : componentFaces)
            {
                auto tId = operandFhToMaTId[fh];

                //std::cout << "fh: " << vfh.idx() << std::endl;
                std::vector<tVertexHandle> faceVertices;
                std::vector<uint> maTIds;

                for (auto fvh : OpenMesh::make_smart(fh, m_mesh).vertices())
                {
                    if (nonManifoldVertexProp[fvh])
                    {
                        //std::cout << "replacing: " << vvh.idx() << " with " << replaceVertex[vvh].idx() << std::endl;
                        auto replaceVh = maToReplaceVertex[operandToMeshArrangementIndex[fvh]];
                        faceVertices.push_back(replaceVh);
                        newFaces[tId].push_back(replaceVh);

                        maTIds.push_back(operandToMeshArrangementIndex[replaceVh]);
                    }
                    else
                    {
                        //std::cout << "leaving: " << vvh.idx() << std::endl;
                        newFaces[tId].push_back(fvh);
                        faceVertices.push_back(fvh);
                        maTIds.push_back(operandToMeshArrangementIndex[fvh]);

                    }
                }

                //adjust indices in MA
                for (int i = 0; i < 3; ++i)
                {
                    ma.m_triangles[tId * 3 + i] = maTIds[i];
                }

            }

            //merge vertices by reinserting faces
            for (auto fh : componentFaces)
            {
                m_mesh.delete_face(fh, false);
            }

            for (auto& [tId, vertices] : newFaces)
            {
                //faceVisited[tId] = true;
                m_mesh.add_face(vertices);
            }

        }

        DMB::saveMesh(m_mesh, "C:/T3D/Samples/Booleans/mergedMesh.obj");

        //reset tagged status
        for (auto vh : m_mesh.vertices())
        {
            m_mesh.status(vh).set_tagged(false);
        }
    }


    //TODO: move this to MA class/strucutre
    //update matrices

    {
        std::vector<tTriplet> tripletList;
        tripletList.reserve(std::size(ma.m_coordinates));
        for (int i = 0; i < std::size(ma.m_triangles); i += 3)
        {
            auto faceLabel = ma.m_labels[i / 3];

            bool belongsToThisOperand = (faceLabel & meshArrangementLabel) != 0;

            if (!belongsToThisOperand)
            {
                continue;
            }

            uint v0 = ma.m_triangles[i + 0];
            uint v1 = ma.m_triangles[i + 1];
            uint v2 = ma.m_triangles[i + 2];

            tripletList.push_back(tTriplet(v0, v1, -1));
            tripletList.push_back(tTriplet(v1, v2, -1));
            tripletList.push_back(tTriplet(v2, v0, -1));
        }

        for (int k = 0; k < edges.outerSize(); ++k)
        {
            for (tSparseIntMatrix::InnerIterator it(edges, k); it; ++it)
            {
                auto r = it.row();
                auto c = it.col();
                edges.coeffRef(r, c) = 0;
            }
        }

        //tSparseIntMatrix edges2(std::size(ma.m_coordinatesImplicit), std::size(ma.m_coordinatesImplicit));

        //edges.resize(std::size(ma.m_coordinatesImplicit), std::size(ma.m_coordinatesImplicit));
        edges.resize(std::size(ma.m_coordinates) / 3, std::size(ma.m_coordinates) / 3);
        edges.setFromTriplets(tripletList.begin(), tripletList.end());

        int edgeId = 1;
        int cnt = 0;

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
                    ++cnt;
                }
            }
        }

        edgeFaceAdjacency.clear();
        edgeFaceAdjacency.resize(cnt + 1);

        for (uint i = 0; i < std::size(ma.m_triangles); i += 3)
        {
            uint tId = i / 3;
            auto faceLabel = ma.m_labels[tId];

            bool belongsToThisOperand = (faceLabel & meshArrangementLabel) != 0;

            if (!belongsToThisOperand)
            {
                continue;
            }

            uint v0 = ma.m_triangles[i + 0];
            uint v1 = ma.m_triangles[i + 1];
            uint v2 = ma.m_triangles[i + 2];

            int e0 = edges.coeff(v0, v1);
            int e1 = edges.coeff(v1, v2);
            int e2 = edges.coeff(v2, v0);

            edgeFaceAdjacency[e0].push_back(tId);
            edgeFaceAdjacency[e1].push_back(tId);
            edgeFaceAdjacency[e2].push_back(tId);
        }
    }

    //m_meshArrangementToVhs.clear();
    //m_meshArrangementToVhs.resize(std::size(ma.m_coordinates) / 3);
    //m_mesh.clear();

#endif

    //vertex-face adjacency
#if 1
    {
        MeshType wholeOperand;
        std::vector<tVertexHandle> operandToMaVh;
        tVertexHandle invalidVh;
        invalidVh.invalidate();
        operandToMaVh.resize(ma.m_coordinatesImplicit.size(), invalidVh);

        for (int i = 0; i < std::size(ma.m_triangles); i += 3)
        {
            uint tId = i / 3;
            auto faceLabel = ma.m_labels[tId];

            bool belongsToThisOperand = (faceLabel & meshArrangementLabel) != 0;

            if (!belongsToThisOperand)
            {
                continue;
            }


            uint maVh0 = ma.m_triangles[i + 0];
            uint maVh1 = ma.m_triangles[i + 1];
            uint maVh2 = ma.m_triangles[i + 2];

            tTriangle triangle =
            {
                tPoint { ma.m_coordinates[maVh0 * 3 + 0], ma.m_coordinates[maVh0 * 3 + 1], ma.m_coordinates[maVh0 * 3 + 2] },
                tPoint { ma.m_coordinates[maVh1 * 3 + 0], ma.m_coordinates[maVh1 * 3 + 1], ma.m_coordinates[maVh1 * 3 + 2] },
                tPoint { ma.m_coordinates[maVh2 * 3 + 0], ma.m_coordinates[maVh2 * 3 + 1], ma.m_coordinates[maVh2 * 3 + 2] },
            };

            int j = 0;
            for (auto maVh : { maVh0, maVh1, maVh2 })
            {
                if (!operandToMaVh[maVh].is_valid())
                {
                    auto newVh = wholeOperand.add_vertex(triangle[j]);
                    operandToMaVh[maVh] = newVh;
                }
                ++j;
            }

            auto vh0 = operandToMaVh[maVh0];
            auto vh1 = operandToMaVh[maVh1];
            auto vh2 = operandToMaVh[maVh2];

            auto fh = wholeOperand.add_face(vh0, vh1, vh2);
            wholeOperand.set_color(fh, { 192,192,192 });

            if (!fh.is_valid())
            {
                vh0 = wholeOperand.add_vertex(triangle[0]);
                vh1 = wholeOperand.add_vertex(triangle[1]);
                vh2 = wholeOperand.add_vertex(triangle[2]);

                fh = wholeOperand.add_face(vh0, vh1, vh2);
                wholeOperand.set_color(fh, { 255,0,0 });

            }

        }

        wholeOperand.update_normals();
        //T3D_VDT_STORE_NAMED_MESH_AND_WAIT("", wholeMA);
        OpenMesh::IO::Options opt = OpenMesh::IO::Options::Default;
        opt += OpenMesh::IO::Options::FaceColor;
        DMB::saveMesh(wholeOperand, "C:/T3D/Samples/Booleans/wholeOperandAfterEdges.obj", opt);


    }

    {
        std::vector<std::vector<uint>> vertexFaceAdjacency;
        vertexFaceAdjacency.resize(std::size(ma.m_coordinates) / 3);

        for (uint i = 0; i < std::size(ma.m_triangles); i += 3)
        {
            uint tId = i / 3;
            auto faceLabel = ma.m_labels[tId];

            bool belongsToThisOperand = (faceLabel & meshArrangementLabel) != 0;

            if (!belongsToThisOperand)
            {
                continue;
            }

            uint v0 = ma.m_triangles[i + 0];
            uint v1 = ma.m_triangles[i + 1];
            uint v2 = ma.m_triangles[i + 2];

            vertexFaceAdjacency[v0].push_back(tId);
            vertexFaceAdjacency[v1].push_back(tId);
            vertexFaceAdjacency[v2].push_back(tId);
        }

        std::vector<tVertexHandle> dbgMaToVh;
        tVertexHandle invalidVh;
        invalidVh.invalidate();
        dbgMaToVh.resize(ma.m_coordinatesImplicit.size(), invalidVh);
        MeshType dbgMesh;

        for (uint vId = 0; vId < vertexFaceAdjacency.size(); ++vId)
        {
            std::queue<uint> facesToVisit;
            std::vector<std::vector<uint>> components;



            for (auto vertexFace : vertexFaceAdjacency[vId])
            {
                if (faceVisited[vertexFace])
                {
                    continue;
                }

                facesToVisit.push(vertexFace);

                std::vector<uint> component;

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

                    faceVertices.push_back(ma.m_triangles[topFace * 3 + 0]);
                    faceVertices.push_back(ma.m_triangles[topFace * 3 + 1]);
                    faceVertices.push_back(ma.m_triangles[topFace * 3 + 2]);

                    auto mid = std::find_if(faceVertices.begin(), faceVertices.end(),
                        [&](const auto& faceVertex)
                        {
                            return faceVertex == vId;
                        });

                    //vector sorted in such a way, that vId vertex is first
                    std::rotate(faceVertices.begin(), mid, faceVertices.end());

                    int e0 = edges.coeff(vId, faceVertices[1]);
                    int e1 = edges.coeff(vId, faceVertices[2]);

                    for (auto fId : edgeFaceAdjacency[e0])
                    {
                        facesToVisit.push(fId);
                    }

                    for (auto fId : edgeFaceAdjacency[e1])
                    {
                        facesToVisit.push(fId);
                    }


                    faceVisited[topFace] = true;
                }

                components.push_back(component);
            }

            //reset tag
            for (auto vertexFace : vertexFaceAdjacency[vId])
            {
                faceVisited[vertexFace] = false;
            }

            auto edgeId = edgeFaceAdjacency.size();

            if (components.size() > 1)
            {
                std::cout << "non manifold vertex (will be problematic): " << vId << std::endl;

                //for (int i = 0; i < components.size(); ++i)
                //{
                //    for (auto currentTriangleId : components[i])
                //    {
                //        uint maVh0 = ma.m_triangles[currentTriangleId * 3 + 0];
                //        uint maVh1 = ma.m_triangles[currentTriangleId * 3 + 1];
                //        uint maVh2 = ma.m_triangles[currentTriangleId * 3 + 2];

                //        tTriangle triangle =
                //        {
                //            tPoint { ma.m_coordinates[maVh0 * 3 + 0], ma.m_coordinates[maVh0 * 3 + 1], ma.m_coordinates[maVh0 * 3 + 2] },
                //            tPoint { ma.m_coordinates[maVh1 * 3 + 0], ma.m_coordinates[maVh1 * 3 + 1], ma.m_coordinates[maVh1 * 3 + 2] },
                //            tPoint { ma.m_coordinates[maVh2 * 3 + 0], ma.m_coordinates[maVh2 * 3 + 1], ma.m_coordinates[maVh2 * 3 + 2] },
                //        };

                //        int j = 0;
                //        for (auto maVh : { maVh0, maVh1, maVh2 })
                //        {
                //            if (!dbgMaToVh[maVh].is_valid())
                //            {
                //                auto newVh = dbgMesh.add_vertex(triangle[j]);
                //                dbgMaToVh[maVh] = newVh;
                //            }
                //            ++j;
                //        }

                //        auto vh0 = dbgMaToVh[maVh0];
                //        auto vh1 = dbgMaToVh[maVh1];
                //        auto vh2 = dbgMaToVh[maVh2];

                //        dbgMesh.add_face(vh0, vh1, vh2);
                //    }

                    //dbgMesh.update_normals();
                    //T3D_VDT_STORE_NAMED_MESH_AND_WAIT("", dbgMesh);

                //}



                //TODO: COPY NEW VERTEX ONLY ONCE PER COMPONENT!!!
                for (int i = 1 /*skip first component*/; i < components.size(); ++i)
                {
                    //auto newMaVhId = addMeshArrangementVertex(ma, vId, currentTriangleId * 3 + j);

                    uint newMaVhId = ma.m_coordinatesImplicit.size();

                    auto* gp = ma.m_coordinatesImplicit[vId];
                    //auto* gpCopy = copyGenericPoint(gp);
                    ma.m_coordinatesImplicit.push_back(gp);

                    ma.m_coordinates.push_back(ma.m_coordinates[vId * 3 + 0]);
                    ma.m_coordinates.push_back(ma.m_coordinates[vId * 3 + 1]);
                    ma.m_coordinates.push_back(ma.m_coordinates[vId * 3 + 2]);

                    for (auto currentTriangleId : components[i])
                    {
                        std::vector<uint> faceVertices;
                        faceVertices.reserve(3);
                        faceVertices.push_back(ma.m_triangles[currentTriangleId * 3 + 0]);
                        faceVertices.push_back(ma.m_triangles[currentTriangleId * 3 + 1]);
                        faceVertices.push_back(ma.m_triangles[currentTriangleId * 3 + 2]);

                        auto mid = std::find_if(faceVertices.begin(), faceVertices.end(),
                            [&](const auto& faceVertex)
                            {
                                return faceVertex == vId;
                            });

                        //replace old ma vh with new ma vh
                        ma.m_triangles[currentTriangleId * 3 + std::distance(faceVertices.begin(), mid)] = newMaVhId;

                        //vector sorted in such a way, that vId vertex is first
                        std::rotate(faceVertices.begin(), mid, faceVertices.end());

                        //remove face from edge face adjacency
                        int e0 = edges.coeff(vId, faceVertices[1]);
                        int e1 = edges.coeff(vId, faceVertices[2]);

                        edgeFaceAdjacency[e0].erase(std::remove_if(edgeFaceAdjacency[e0].begin(), edgeFaceAdjacency[e0].end(),
                            [&](const uint& fhId)
                            {
                                return fhId == currentTriangleId;
                            }), edgeFaceAdjacency[e0].end());
                        edgeFaceAdjacency[e1].erase(std::remove_if(edgeFaceAdjacency[e1].begin(), edgeFaceAdjacency[e1].end(),
                            [&](const uint& fhId)
                            {
                                return fhId == currentTriangleId;
                            }), edgeFaceAdjacency[e1].end());



                        //TODO: edge face adjacency should be updated...
                        edges.conservativeResize(edges.rows() + 2, edges.cols() + 2);
                        ++edgeId;
                        edges.coeffRef(newMaVhId, faceVertices[1]) = edgeId;
                        edges.coeffRef(faceVertices[1], newMaVhId) = edgeId;

                        ++edgeId;
                        edges.coeffRef(newMaVhId, faceVertices[2]) = edgeId;
                        edges.coeffRef(faceVertices[2], newMaVhId) = edgeId;
                    }
                }

            }

        }
    }
#endif

    {
        std::vector<tTriplet> tripletList;
        tripletList.reserve(std::size(ma.m_coordinates));
        for (int i = 0; i < std::size(ma.m_triangles); i += 3)
        {
            auto faceLabel = ma.m_labels[i / 3];

            bool belongsToThisOperand = (faceLabel & meshArrangementLabel) != 0;

            if (!belongsToThisOperand)
            {
                continue;
            }

            uint v0 = ma.m_triangles[i + 0];
            uint v1 = ma.m_triangles[i + 1];
            uint v2 = ma.m_triangles[i + 2];

            tripletList.push_back(tTriplet(v0, v1, -1));
            tripletList.push_back(tTriplet(v1, v2, -1));
            tripletList.push_back(tTriplet(v2, v0, -1));
        }

        for (int k = 0; k < edges.outerSize(); ++k)
        {
            for (tSparseIntMatrix::InnerIterator it(edges, k); it; ++it)
            {
                auto r = it.row();
                auto c = it.col();
                edges.coeffRef(r, c) = 0;
            }
        }

        //tSparseIntMatrix edges2(std::size(ma.m_coordinatesImplicit), std::size(ma.m_coordinatesImplicit));

        //edges.resize(std::size(ma.m_coordinatesImplicit), std::size(ma.m_coordinatesImplicit));
        edges.resize(std::size(ma.m_coordinates) / 3, std::size(ma.m_coordinates) / 3);
        edges.setFromTriplets(tripletList.begin(), tripletList.end());

        int edgeId = 1;
        int cnt = 0;

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
                    ++cnt;
                }
            }
        }

        edgeFaceAdjacency.clear();
        edgeFaceAdjacency.resize(cnt + 1);

        for (uint i = 0; i < std::size(ma.m_triangles); i += 3)
        {
            uint tId = i / 3;
            auto faceLabel = ma.m_labels[tId];

            bool belongsToThisOperand = (faceLabel & meshArrangementLabel) != 0;

            if (!belongsToThisOperand)
            {
                continue;
            }

            uint v0 = ma.m_triangles[i + 0];
            uint v1 = ma.m_triangles[i + 1];
            uint v2 = ma.m_triangles[i + 2];

            int e0 = edges.coeff(v0, v1);
            int e1 = edges.coeff(v1, v2);
            int e2 = edges.coeff(v2, v0);

            edgeFaceAdjacency[e0].push_back(tId);
            edgeFaceAdjacency[e1].push_back(tId);
            edgeFaceAdjacency[e2].push_back(tId);
        }
    }

    m_meshArrangementToVhs.clear();
    m_meshArrangementToVhs.resize(std::size(ma.m_coordinates) / 3);
    m_mesh.clear();

    std::queue<uint> trianglesToProcess;

    bool success = true;

    for (uint i = 0; i < std::size(ma.m_triangles); i += 3)
    {
        uint tId = i / 3;
        auto faceLabel = ma.m_labels[tId];

        bool belongsToThisOperand = (faceLabel & meshArrangementLabel) != 0;

        if (!belongsToThisOperand || faceVisited[tId])
        {
            continue;
        }

        //if (std::find(nonManifoldFaces.begin(), nonManifoldFaces.end(), tId) != nonManifoldFaces.end())
        //{
        //    continue;
        //}

        trianglesToProcess.push(tId);

        while (!trianglesToProcess.empty())
        {
            uint currentTriangleId = trianglesToProcess.front();
            trianglesToProcess.pop();

            auto currentFaceLabel = ma.m_labels[currentTriangleId];

            bool currentFhBelongsToThisOperand = (currentFaceLabel & meshArrangementLabel) != 0;
            bool isCoplanar = (currentFaceLabel & meshArrangementLabelComplement) != 0;

            if (!currentFhBelongsToThisOperand || faceVisited[currentTriangleId])
            {
                continue;
            }

            //if (std::find(nonManifoldFaces.begin(), nonManifoldFaces.end(), currentTriangleId) != nonManifoldFaces.end())
            //{
            //    continue;
            //}
            uint maVh0 = ma.m_triangles[currentTriangleId * 3 + 0];
            uint maVh1 = ma.m_triangles[currentTriangleId * 3 + 1];
            uint maVh2 = ma.m_triangles[currentTriangleId * 3 + 2];


            if (isCoplanar && findLowestSetBit(currentFaceLabel) == m_meshOperandLabel)
            {
                std::array<uint, 3> coplanarVertices = { maVh0, maVh1, maVh2 };

                m_coplanarTriangles.push_back(std::make_pair(currentTriangleId, coplanarVertices));
            }


            tTriangle triangle =
            {
                tPoint { ma.m_coordinates[maVh0 * 3 + 0], ma.m_coordinates[maVh0 * 3 + 1], ma.m_coordinates[maVh0 * 3 + 2] },
                tPoint { ma.m_coordinates[maVh1 * 3 + 0], ma.m_coordinates[maVh1 * 3 + 1], ma.m_coordinates[maVh1 * 3 + 2] },
                tPoint { ma.m_coordinates[maVh2 * 3 + 0], ma.m_coordinates[maVh2 * 3 + 1], ma.m_coordinates[maVh2 * 3 + 2] },
            };

            int j = 0;
            for (auto maVh : { maVh0, maVh1, maVh2 })
            {
                if (!m_meshArrangementToVhs[maVh].is_valid())
                {
                    auto newVh = m_mesh.add_vertex(triangle[j]);
                    m_meshArrangementToVhs[maVh] = newVh;
                }
                ++j;
            }

            auto vh0 = m_meshArrangementToVhs[maVh0];
            auto vh1 = m_meshArrangementToVhs[maVh1];
            auto vh2 = m_meshArrangementToVhs[maVh2];

            //get adjacent triangles, before updating vertices
            int e0 = edges.coeff(maVh0, maVh1);
            int e1 = edges.coeff(maVh1, maVh2);
            int e2 = edges.coeff(maVh2, maVh0);

            bool nonManifoldEdge = false;
            //test non-manifold edge
            {
                auto he0 = m_mesh.find_halfedge(vh0, vh1);
                auto he1 = m_mesh.find_halfedge(vh1, vh2);
                auto he2 = m_mesh.find_halfedge(vh2, vh0);

                if (he0.is_valid() && !he0.is_boundary())
                {
                    m_mesh.status(vh0).set_tagged(true);
                    m_mesh.status(vh1).set_tagged(true);

                    nonManifoldEdge = true;
                    std::cout << "non manifold edge he0" << std::endl;
                }

                if (he1.is_valid() && !he1.is_boundary())
                {
                    m_mesh.status(vh1).set_tagged(true);
                    m_mesh.status(vh2).set_tagged(true);
                    nonManifoldEdge = true;
                    std::cout << "non manifold edge he1" << std::endl;

                }

                if (he2.is_valid() && !he2.is_boundary())
                {
                    m_mesh.status(vh2).set_tagged(true);
                    m_mesh.status(vh0).set_tagged(true);
                    nonManifoldEdge = true;
                    std::cout << "non manifold edge he2" << std::endl;

                }
            }

            bool nonManifoldVertex = false;
            //test non-manifold vertex
            for (auto vhToTest : { vh0, vh1, vh2 })
            {
                if (!m_mesh.is_boundary(vhToTest))
                {
                    std::cout << "ma indices: " <<maVh0 <<" " << maVh1 << " " << maVh2 << std::endl;
                    std::cout << "non manifold vertex encountered: " << vhToTest.idx() <<" " << operandToMeshArrangementIndex[vhToTest] << std::endl;
                    nonManifoldVertex = true;
                    m_mesh.status(vhToTest).set_tagged(true);
                }
            }

#if 0
            if (m_mesh.status(vh0).tagged()  )
            {
                if (!nonManifoldVertexMap[vh0].is_valid())
                {
                    auto newVh = m_mesh.add_vertex(m_mesh.point(vh0));
                    nonManifoldVertexMap[vh0] = newVh;
                    vh0 = newVh;
                    maVh0 = addMeshArrangementVertex(ma, maVh0, currentTriangleId * 3 + 0);
                    m_meshArrangementToVhs.push_back(newVh);
                }
                else
                {
                    vh0 = nonManifoldVertexMap[vh0];
                    maVh0 = operandToMeshArrangementIndex[vh0];
                    ma.m_triangles[currentTriangleId * 3 + 0] = maVh0;
                }
            }

            if (m_mesh.status(vh1).tagged())
            {
                if (!nonManifoldVertexMap[vh1].is_valid())
                {
                    auto newVh = m_mesh.add_vertex(m_mesh.point(vh1));
                    nonManifoldVertexMap[vh1] = newVh;
                    vh1 = newVh;
                    maVh1 = addMeshArrangementVertex(ma, maVh1, currentTriangleId * 3 + 1);
                    m_meshArrangementToVhs.push_back(newVh);
                }
                else
                {
                    vh1 = nonManifoldVertexMap[vh1];
                    maVh1 = operandToMeshArrangementIndex[vh1];
                    ma.m_triangles[currentTriangleId * 3 + 1] = maVh1;
                }
            }

            if (m_mesh.status(vh2).tagged())
            {
                if (!nonManifoldVertexMap[vh2].is_valid())
                {
                    auto newVh = m_mesh.add_vertex(m_mesh.point(vh2));
                    nonManifoldVertexMap[vh2] = newVh;
                    vh2 = newVh;
                    maVh2 = addMeshArrangementVertex(ma, maVh2, currentTriangleId * 3 + 2);
                    m_meshArrangementToVhs.push_back(newVh);
                }
                else
                {
                    vh2 = nonManifoldVertexMap[vh2];
                    maVh2 = operandToMeshArrangementIndex[vh2];
                    ma.m_triangles[currentTriangleId * 3 + 2] = maVh2;

                }
            }
#endif
            operandToMeshArrangementIndex[vh0] = maVh0;
            operandToMeshArrangementIndex[vh1] = maVh1;
            operandToMeshArrangementIndex[vh2] = maVh2;





            //if (nonManifoldVertex || nonManifoldEdge)
            //{
            //    std::cout << "Adding face will fail" << std::endl;
            //}

            tFaceHandle fh = m_mesh.add_face(vh0, vh1, vh2);

            if (fh.is_valid())
            {
                if (isCoplanar)
                {
                    coplanarFace[fh] = true;
                }
                //std::cout << "Face added" << std::endl;

            }
            else
            {
                success = false;
                std::cout << "Could not add face" << std::endl;

                //MeshType debugMesh;

                //auto dbVh0 = debugMesh.add_vertex(m_mesh.point(vh0));
                //auto dbVh1 = debugMesh.add_vertex(m_mesh.point(vh1));
                //auto dbVh2 = debugMesh.add_vertex(m_mesh.point(vh2));
                //auto dbFh = debugMesh.add_face(dbVh0, dbVh1, dbVh2);
                //debugMesh.set_color(dbFh, { 255,0,0 });
                //debugMesh.update_normals();
                //m_mesh.update_normals();

                //for (auto cfh : m_mesh.faces())
                //{
                //    m_mesh.set_color(cfh, { 192,192,192 });
                //}
                //m_mesh.update_normals();

                //DMB::saveMesh(m_mesh, "C:/T3D/Samples/Booleans/m.obj");
                //DMB::saveMesh(debugMesh, "C:/T3D/Samples/Booleans/dbgM.obj");

                //T3D_VDT_STORE_NAMED_MESH_AND_WAIT("", m_mesh);
                //T3D_VDT_STORE_NAMED_MESH_AND_WAIT("", debugMesh);


                //return false;
            }

            //this triangle was processed
            faceVisited[currentTriangleId] = true;
            pFhToMaFh[fh] = currentTriangleId;



            if (edgeFaceAdjacency[e0].size() <= 2)
            {
                for (auto adjacentTriangle : edgeFaceAdjacency[e0])
                {
                    trianglesToProcess.push(adjacentTriangle);
                }
            }

            if (edgeFaceAdjacency[e1].size() <= 2)
            {
                for (auto adjacentTriangle : edgeFaceAdjacency[e1])
                {
                    trianglesToProcess.push(adjacentTriangle);
                }
            }

            if (edgeFaceAdjacency[e2].size() <= 2)
            {
                for (auto adjacentTriangle : edgeFaceAdjacency[e2])
                {
                    trianglesToProcess.push(adjacentTriangle);
                }
            }
        }
               
        for (auto vh : m_mesh.vertices())
        {
            m_mesh.status(vh).set_tagged(true);
            nonManifoldVertexMap[vh].invalidate();
        }

    }

    //reset tagged flag
    for (auto vh : m_mesh.vertices())
    {
        m_mesh.status(vh).set_tagged(false);
    }

    m_mesh.update_normals();
    DMB::saveMesh(m_mesh, "C:/T3D/Samples/Booleans/operandMesh" + std::to_string(m_meshOperandLabel) + ".ply");
    T3D_VDT_STORE_NAMED_MESH_AND_WAIT("after creating operand", m_mesh);

    return success;
}

template<typename MeshType>
inline void DMB::MeshOperand<MeshType>::handleCoplanarFaces()
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
            m_mesh.delete_face(fh);
        }
    }

}

template<typename MeshType>
inline void DMB::MeshOperand<MeshType>::detectBoundaries()
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
inline void DMB::MeshOperand<MeshType>::detectIntersectionEdges(const MeshArrangement& ma)
{
    OpenMesh::EProp<bool> pIntersectionEdge(false, m_mesh, "bo_pIntersectionEdge");
    m_pIntersectionEdge = pIntersectionEdge.getRawProperty();

    //TODO: debug
    pIntersectionEdge.set_persistent(true);

    std::bitset<NBIT> meshArrangementLabel;
    meshArrangementLabel[m_meshOperandLabel] = 1;
    std::bitset<NBIT> meshArrangementLabelComplement = meshArrangementLabel;
    meshArrangementLabelComplement.flip();

    //for (auto vh : m_mesh.vertices())
    //{
    //    m_mesh.set_color(vh, { 0,0,0 });
    //}

    //for (int k = 0; k < m_edgeValanceMatrix.outerSize(); ++k)
    //{
    //    for (tSparseIntMatrix::InnerIterator it(m_edgeValanceMatrix, k); it; ++it)
    //    {
    //        //TODO: check if edge exists
    //        if (it.value() > 2 )
    //        {
    //            size_t r = it.row();
    //            size_t c = it.col();

    //            //if (r < 0 || r >= m_meshArrangementToVhs.size() ||
    //            //    c < 0 || c >= m_meshArrangementToVhs.size())
    //            //{
    //            //    continue;
    //            //}

    //            auto he = m_mesh.find_halfedge(m_meshArrangementToVhs[r], m_meshArrangementToVhs[c]);

    //            if (!he.is_valid())
    //            {
    //                continue;
    //            }

    //            m_mesh.set_color(he.from(), { 255,0,0 });
    //            m_mesh.set_color(he.to(), { 255,0,0 });

    //            pIntersectionEdge[he.edge()] = true;
    //        }
    //    }
    //}

    OpenMesh::VProp<uint> pOperandToMeshArrangementIndex(m_mesh, m_pOperandToMeshArrangementIndex);

    for (auto eh : m_mesh.edges())
    {
        auto v0 = eh.v0();
        auto v1 = eh.v1();

        //m_mesh.set_color(v0, { 0,0,0 });
        //m_mesh.set_color(v1, { 0,0,0 });

        auto edgeValence = m_edgeValanceMatrix.coeff(pOperandToMeshArrangementIndex[v0], pOperandToMeshArrangementIndex[v1]);

        if (edgeValence > 2)
        {
            //m_mesh.set_color(v0, { 255,0,0 });
            //m_mesh.set_color(v1, { 255,0,0 });
            pIntersectionEdge[eh] = true;

        }
    }


    //{

    //    OpenMesh::VProp<bool> pIntersectionVertex(false, m_mesh);

    //    for (auto eh : m_mesh.edges())
    //    {
    //        if (pIntersectionEdge[eh])
    //        {
    //            auto v0 = eh.v0();
    //            auto v1 = eh.v1();
    //            pIntersectionVertex[v0] = true;
    //            pIntersectionVertex[v1] = true;
    //        }
    //    }

    //    std::map<int, int> meshToObjVertexMap;
    //    int objVertexId = 1;
    //    std::ofstream lineObj("C:/T3D/Samples/Booleans/intEdges" + std::to_string(m_meshOperandLabel) + ".obj");

    //    for (auto vh : m_mesh.vertices())
    //    {
    //        if (pIntersectionVertex[vh])
    //        {
    //            meshToObjVertexMap[vh.idx()] = objVertexId;

    //            auto p = m_mesh.point(vh);
    //            lineObj << "v " << std::to_string(p[0]) << " " << std::to_string(p[1]) << " " << std::to_string(p[2]) << std::endl;

    //            ++objVertexId;
    //        }
    //    }

    //    for (auto eh : m_mesh.edges())
    //    {
    //        if (pIntersectionEdge[eh])
    //        {
    //            auto v0 = eh.v0();
    //            auto v1 = eh.v1();
    //            lineObj << "l " << std::to_string(meshToObjVertexMap[v0.idx()]) << " " << std::to_string(meshToObjVertexMap[v1.idx()]) << std::endl;
    //        }
    //    }



    //}

    //{
    //    m_mesh.update_normals();
    //    OpenMesh::IO::Options opt = OpenMesh::IO::Options::Default;
    //    opt += OpenMesh::IO::Options::FaceColor;
    //    opt += OpenMesh::IO::Options::EdgeColor;
    //    opt += OpenMesh::IO::Options::VertexColor;
    //    //DMB::saveMesh(m_mesh, "C:/T3D/Samples/Booleans/operandMesh" + std::to_string(m_meshOperandLabel) + "intEdges.om", true, opt);
    //    DMB::saveMesh(m_mesh, "C:/T3D/Samples/Booleans/operandMesh" + std::to_string(m_meshOperandLabel) + "intEdges.ply", true, opt);
    //}
    //T3D_VDT_STORE_NAMED_MESH_AND_WAIT("", m_mesh);


}

template<typename MeshType>
inline void DMB::MeshOperand<MeshType>::discardUnusableComponents()
{
    //delete components with unclosed intersection boundaries
    for (auto intersectionCurve : m_intersectionCurves)
    {
        std::vector<OpenMesh::SmartFaceHandle> intersectionRegion;

        bool loopClosed = intersectionCurve.closed;

        if (loopClosed)
        {
            continue;
        }

        //get seed face
        auto seedFh = intersectionCurve.loopHalfedges.front().face();

        std::queue<tFaceHandle> facesToVisit;
        facesToVisit.push(seedFh);
        m_mesh.status(seedFh).set_selected(true);

        while (!facesToVisit.empty())
        {
            auto topFace = OpenMesh::make_smart(facesToVisit.front(), m_mesh);
            facesToVisit.pop();

            for (auto ffh : topFace.faces())
            {
                if (!ffh.selected())
                {
                    m_mesh.status(ffh).set_selected(true);
                    facesToVisit.push(ffh);
                }
            }
        }
    }

    for (auto fh : m_mesh.faces())
    {
        if (fh.selected())
        {
            m_mesh.delete_face(fh);
        }
    }
}

template<typename MeshType>
inline uint DMB::MeshOperand<MeshType>::addMeshArrangementVertex(MeshArrangement& ma, uint vertexHandle, uint triangleVertexHandle)
{
    uint newMaVh = ma.m_coordinatesImplicit.size();

    auto* gp = ma.m_coordinatesImplicit[vertexHandle];
    //auto* gpCopy = copyGenericPoint(gp);
    ma.m_coordinatesImplicit.push_back(gp);

    ma.m_triangles[triangleVertexHandle] = newMaVh;

    ma.m_coordinates.push_back(ma.m_coordinates[vertexHandle * 3 + 0]);
    ma.m_coordinates.push_back(ma.m_coordinates[vertexHandle * 3 + 1]);
    ma.m_coordinates.push_back(ma.m_coordinates[vertexHandle * 3 + 2]);

    return newMaVh;
}

template<typename MeshType>
inline int DMB::MeshOperand<MeshType>::findHighestSetBit(const std::bitset<NBIT> bitset)
{
    for (std::size_t pos = (NBIT - 1); pos >= 0; --pos)
    {
        if (bitset[pos] == 1)
        {
            return pos;
        }
    }

    return NBIT;
}

template<typename MeshType>
inline int DMB::MeshOperand<MeshType>::findLowestSetBit(const std::bitset<NBIT> bitset)
{
    for (std::size_t pos = 0; pos < (NBIT - 1); ++pos)
    {
        if (bitset[pos] == 1)
        {
            return pos;
        }
    }

    return NBIT;
}

