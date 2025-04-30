//TODO: License

#pragma once

//stl
#include <algorithm>

//Eigen
#include <Eigen/Dense>

//local
#include "AABBTree.h"
#include "Utils.h"

namespace DMB
{
    template<typename MeshType>
    class FastWindingNumber
    {
    public: //types
        using tTree = typename AABBTree<MeshType>;
        using tNode = tTree::AABBTreeNode;
        using tPoint = typename MeshType::Point;
        using tNormal = typename MeshType::Normal;
        using tScalar = typename MeshType::Scalar;
        using tFace = typename MeshType::FaceHandle;
        using tEigenMat3 = Eigen::Matrix<tScalar, 3, 3>;

        struct Dipole
        {
            Dipole() :
                m_weightedCentroid(0, 0, 0),
                m_radiusSquared(0),
                m_order1Vec(0, 0, 0),
                m_order2Mat(tEigenMat3::Zero())
            {
            }

            Dipole(const MeshType& mesh, const std::span<tFace>& faces, const OpenMesh::FPropHandleT<tScalar>& faceAreas, const OpenMesh::FPropHandleT<tPoint>& faceCentroids, tScalar beta) :
                m_weightedCentroid(0,0,0),
                m_radiusSquared(0),
                m_order1Vec(0,0,0),
                m_order2Mat(tEigenMat3::Zero())
            {
                // compute area-weighted centroid of triangles, we use this as the expansion point
                tScalar sumArea = 0;

                for (const auto& fh : faces)
                {
                    tScalar faceArea = mesh.property(faceAreas, fh);

                    sumArea += faceArea;
                    m_weightedCentroid += faceArea * mesh.property(faceCentroids, fh);
                }

                m_weightedCentroid /= sumArea;

                tScalar radiusSquared = 0.0f;
                // compute first and second-order coefficients of FWN Taylor expansion, as well as
                // 'radius' value r, which is max dist from any tri vertex to p
                //for (It it = begin; it != end; ++it)
                for(const auto& fh : faces)
                {
                    //auto fh = *it;

                    auto fv_it = mesh.cfv_iter(fh);
                    const auto& P0 = mesh.point(*fv_it); ++fv_it;
                    const auto& P1 = mesh.point(*fv_it); ++fv_it;
                    const auto& P2 = mesh.point(*fv_it);

                    tNormal faceNormal = mesh.normal(fh);
                    tScalar faceArea = mesh.property(faceAreas, fh);

                    m_order1Vec += faceArea * faceNormal;

                    tPoint dcp = mesh.property(faceCentroids, fh) - m_weightedCentroid;

                    /// Construct outer-product of u*transpose(v) of u and v
                    /// result is that Mij = u_i * v_j
                    // u = dcp, v = faceNormal
                    tEigenMat3 order2Mat;
                    order2Mat << dcp[0] * faceNormal[0], dcp[0] * faceNormal[1], dcp[0] * faceNormal[2],
                        dcp[1] * faceNormal[0], dcp[1] * faceNormal[1], dcp[1] * faceNormal[2],
                        dcp[2] * faceNormal[0], dcp[2] * faceNormal[1], dcp[2] * faceNormal[2];

                    m_order2Mat += faceArea * order2Mat;

                    radiusSquared = std::max({ radiusSquared, DMB::length2(m_weightedCentroid - P0), DMB::length2(m_weightedCentroid - P1), DMB::length2(m_weightedCentroid - P2) });
                }

                m_radiusSquared = beta * beta * radiusSquared;
            }

            bool canUseFWNApproximation(const tPoint& p) const
            {
                return DMB::length2(m_weightedCentroid - p) > m_radiusSquared;
            }

            double order2Approx(const tPoint& p) const
            {
                static constexpr const tScalar FOUR_PI = static_cast<tScalar>(4.f * M_PI);

                tPoint dpq = (m_weightedCentroid - p);
                tScalar len = dpq.norm();
                tScalar len3 = len * len * len;
                tScalar fourPi_len3 = tScalar(1.0) / (FOUR_PI * len3);

                double order1 = fourPi_len3 * DMB::dot(m_order1Vec, dpq);

                // second-order hessian \grad^2(G)
                tScalar c = -tScalar(3.0) / (FOUR_PI * len3 * len * len);

                // expanded-out version below avoids extra constructors
                tEigenMat3 hessian;
                hessian << fourPi_len3 + c * dpq[0] * dpq[0], c* dpq[0] * dpq[1], c* dpq[0] * dpq[2],
                    c* dpq[1] * dpq[0], fourPi_len3 + c * dpq[1] * dpq[1], c* dpq[1] * dpq[2],
                    c* dpq[2] * dpq[0], c* dpq[2] * dpq[1], fourPi_len3 + c * dpq[2] * dpq[2];

                //inner product of matrices
                double order2Mat = 0;

                for (int r = 0; r < 3; ++r)
                {
                    order2Mat += hessian.row(r).dot(m_order2Mat.row(r));
                }

                return order1 + order2Mat;

            }

            tPoint     m_weightedCentroid;
            tScalar    m_radiusSquared;
            tPoint     m_order1Vec;
            tEigenMat3 m_order2Mat;
        };

    public: //methods

        void build(const MeshType& mesh)
        {
            m_mesh = &mesh;

            m_tree.build(mesh);

            buildWindingNumberCache(mesh);
        }

        bool inside(const tPoint& p) const
        {
            return windingNumber(p) > 0.5;
        }

        double windingNumber(const tPoint& p) const
        {
            return evaluateWindingNumber(p, 0);
        }

    private: //methods

        void buildWindingNumberCache(const MeshType& mesh)
        {
            m_dipoles.clear();
            m_dipoles.resize(m_tree.nodes().size());

            auto centroidProp = OpenMesh::makeTemporaryProperty<tFace, tPoint>(const_cast<MeshType&>(mesh));
            auto faceAreaProp = OpenMesh::makeTemporaryProperty<tFace, tScalar>(const_cast<MeshType&>(mesh));

            // compute face centroids
            for (auto fh : mesh.faces())
            {
                centroidProp[fh] = mesh.calc_face_centroid(fh);
                faceAreaProp[fh] = mesh.calc_face_area(fh);
            }

            //TODO: do not compute for leaf nodes?
            for (int nodeIndex = 0; nodeIndex < m_tree.nodes().size(); ++nodeIndex)
            {
                if (m_tree.node(nodeIndex).isLeaf())
                {
                    continue;
                }

                m_dipoles[nodeIndex] = Dipole(mesh, m_tree.nodeFaces(nodeIndex), faceAreaProp.getRawProperty(), centroidProp.getRawProperty(), 2.0f);
            }
        }


        double evaluateWindingNumber(const tPoint& p, std::size_t nodeIdx) const
        {
            double wn = 0.0;

            const auto& node = m_tree.node(nodeIdx);

            //can we use approximation?
            if (m_dipoles[nodeIdx] && m_dipoles[nodeIdx]->canUseFWNApproximation(p))
            {
                wn += m_dipoles[nodeIdx]->order2Approx(p);
            }
            else
            {

                if (node.isLeaf())
                {
                    for (auto fh : node.faces)
                    {
                        wn += DMB::solidAngle(*m_mesh, fh, p);
                    }
                }
                else
                {
                    wn += evaluateWindingNumber(p, node.leftNode);
                    wn += evaluateWindingNumber(p, node.rightNode);
                }
            }

            return wn;
        }

    private: //properties
        AABBTree<MeshType> m_tree;

        std::vector<std::optional<Dipole>> m_dipoles;

        //std::weak_ptr<MeshType> m_mesh;
        const MeshType* m_mesh; //TODO: rework to weak_ptr?
    };
}