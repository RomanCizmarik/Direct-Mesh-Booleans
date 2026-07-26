//TODO: License

#pragma once

#include "BoundingBox.h"

#include <random>
#include <algorithm>

#include <OpenMesh/Core/Mesh/PolyConnectivity.hh>

#include <numerics.h> //indirect predicates - for bigfloat

namespace DMB
{
    //! Dot product.
    template<typename VectorType>
    inline auto dot(const VectorType& v0, const VectorType& v1) -> decltype(v0[0] * v1[0])
    {
        using ScalarType = decltype(v0[0] * v1[0]);
        return std::inner_product(v0.begin(), v0.end(), v1.begin(), ScalarType(0));
    }

    //! Squared Euclidian norm.
    template<typename VectorType>
    inline auto length2(const VectorType& v) -> decltype(v[0] * v[0])
    {
        using ScalarType = decltype(v[0] * v[0]);

        return std::accumulate(v.begin(), v.end(),
            ScalarType(0),
            [](const ScalarType& l, const ScalarType& r) { return l + r * r; });
    }

    //! Euclidian norm.
    template<typename VectorType>
    inline auto length(const VectorType& v) -> decltype(v[0] * v[0])
    {
        return std::sqrt(DMB::length2<VectorType>(v));
    }

    template<typename T>
    void makeUniqueVector(T& vec)
    {
        std::sort(vec.begin(), vec.end());

        auto it = std::unique(vec.begin(), vec.end());

        vec.resize(std::distance(vec.begin(), it));
    }

	template<typename MeshType>
	BoundingBoxT<typename MeshType::Point> calcMeshBoundingBox(const MeshType& mesh)
	{
		BoundingBoxT<typename MeshType::Point> bb;

		for (auto vh : mesh.vertices())
		{
			bb.expandBy(mesh.point(vh));
		}

		return bb;
	}
    template<typename MeshType>
    BoundingBoxT<typename MeshType::Point> calcComponentBoundingBox(const MeshType& mesh, const std::vector<OpenMesh::SmartFaceHandle>& component)
    {
        BoundingBoxT<typename MeshType::Point> bb;

        for (auto fh : component)
        {
            for (auto vh : fh.vertices())
            {
                bb.expandBy(mesh.point(vh));
            }

        }

        return bb;
    }

	template<typename MeshType, typename FaceIterator>
	BoundingBoxT<typename MeshType::Point> calcFacesBoundingBox(const MeshType& mesh, FaceIterator facesBegin, FaceIterator facesEnd)
	{
		BoundingBoxT<typename MeshType::Point> bb;

		for (auto it = facesBegin; it != facesEnd; ++it)
		{
			auto fh = OpenMesh::make_smart(*it, &mesh);

			for (auto vh : fh.vertices())
			{
				bb.expandBy(mesh.point(vh));
			}

		}

		return bb;
	}

    template<typename MeshType, typename BoundingBoxType>
    MeshType boundingBoxToMesh(const BoundingBoxType& bb)
    {
        MeshType m;
        /*
           6 --- 7
          /|    /|
         2 --- 3 |
         | 4 --| 5
         |/    |/
         0 --- 1

         */

         /*
              +y   +z
               | /
         -x____|/___ +x
              /|
             / |
           -z  -y
         */

        using tPoint = typename MeshType::Point;

        auto min = bb.m_min;
        auto max = bb.m_max;

        auto v0 = m.add_vertex(tPoint(min[0], min[1], min[2]));
        auto v1 = m.add_vertex(tPoint(max[0], min[1], min[2]));
        auto v2 = m.add_vertex(tPoint(min[0], max[1], min[2]));
        auto v3 = m.add_vertex(tPoint(max[0], max[1], min[2]));
        auto v4 = m.add_vertex(tPoint(min[0], min[1], max[2]));
        auto v5 = m.add_vertex(tPoint(max[0], min[1], max[2]));
        auto v6 = m.add_vertex(tPoint(min[0], max[1], max[2]));
        auto v7 = m.add_vertex(tPoint(max[0], max[1], max[2]));

        //front
        m.add_face(v0, v3, v1);
        m.add_face(v0, v2, v3);

        //right
        m.add_face(v1, v7, v5);
        m.add_face(v1, v3, v7);

        //back
        m.add_face(v5, v6, v4);
        m.add_face(v5, v7, v6);

        //left
        m.add_face(v4, v2, v0);
        m.add_face(v4, v6, v2);

        //bottom          
        m.add_face(v0, v1, v5);
        m.add_face(v0, v5, v4);

        //top            
        m.add_face(v2, v7, v3);
        m.add_face(v2, v6, v7);

        return m;
    }

    template<typename Mesh>
    void OpenMesh2Matrix(const Mesh& mesh, Eigen::MatrixXd& V, Eigen::MatrixXi& F)
    {
        assert(mesh.is_trimesh());

        V.resize(mesh.n_vertices(), 3);
        F.resize(mesh.n_faces(), 3);

        for (auto vh : mesh.vertices())
        {
            const auto& p = mesh.point(vh);
            V(vh.idx(), 0) = p[0];
            V(vh.idx(), 1) = p[1];
            V(vh.idx(), 2) = p[2];
        }

        for (auto fh : mesh.faces())
        {
            int vi = 0;
            for (auto vh : fh.vertices())
            {
                F(fh.idx(), vi) = vh.idx();
                vi++;
            }
        }

    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    //!\brief   Calculates solid angle of a triangle projected onto a sphere.
    //! This is used for computing winding numbers generalized to three dimensions.
    //!
    //!\param   A        First vertex of a triangle
    //!\param   B        Second vertex of a triangle
    //!\param   C        Third vertex of a triangle
    //!\param   p        Center of the projection sphere
    //!
    //!\return  Solid angle scalar normalized by angle (so sphere sums to 1 or -1)
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    template<typename VectorType>
    double solidAngle(const VectorType& A, const VectorType& B, const VectorType& C,
        const VectorType& p)
    {
        VectorType v0 = A - p;
        VectorType v1 = B - p;
        VectorType v2 = C - p;

        auto vl0 = length(v0);
        auto vl1 = length(v1);
        auto vl2 = length(v2);

        // Compute determinant
        auto det = v0[0] * v1[1] * v2[2] +
            v1[0] * v2[1] * v0[2] +
            v2[0] * v0[1] * v1[2] -
            v2[0] * v1[1] * v0[2] -
            v1[0] * v0[1] * v2[2] -
            v0[0] * v2[1] * v1[2];

        //    Eigen::Matrix<SType, 1, 3> dp;
        auto dp0 = v1[0] * v2[0]
            + v1[1] * v2[1]
            + v1[2] * v2[2];

        auto dp1 = v2[0] * v0[0]
            + v2[1] * v0[1]
            + v2[2] * v0[2];

        auto dp2 = v0[0] * v1[0]
            + v0[1] * v1[1]
            + v0[2] * v1[2];

        // Compute winding number
        // Only divide by TWO_PI instead of 4*pi because there was a 2 out front
        return std::atan2(det, vl0 * vl1 * vl2 +
            dp0 * vl0 +
            dp1 * vl1 +
            dp2 * vl2) / (2. * M_PI);
    }

    template<typename MeshType>
    double solidAngle(const MeshType& mesh, typename MeshType::FaceHandle fh, const typename MeshType::Point& p)
    {
        assert(fh.is_valid() && !mesh.status(fh).deleted());

        auto fv_it = mesh.cfv_iter(fh);

        auto A = mesh.point(fv_it); ++fv_it;
        auto B = mesh.point(fv_it); ++fv_it;
        auto C = mesh.point(fv_it);

        return DMB::solidAngle(A, B, C, p);
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    //!\brief   Calculates generalized winding number for a whole mesh from one point.
    //! This projects every triangle onto a sphere and sums the winding numbers.
    //!
    //!\param   mesh     Tested mesh
    //!\param   p        Center of the projection sphere
    //!
    //!\return  If the mesh is watertight, it returns 1 when the point is inside and 0 when the point is outside the mesh.
    //!         For non-watertight but normal meshes, the number is inside the <0, 1> interval indicating how much inside it is.
    //!         For broken meshes, god help you.
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    template<typename MeshType>
    double windingNumber(const MeshType& mesh, const typename MeshType::Point& p)
    {
        double w = 0;

        for (auto fh : mesh.faces())
        {
            w += solidAngle(mesh, fh, p);
        }

        return w;
    }


    template<typename MeshType>
    void copyMeshPart(const MeshType& srcMesh, MeshType& destMesh, std::vector<typename MeshType::FaceHandle> facesToCopy, bool doGarbageCollection = false)
    {
        destMesh = srcMesh;

        auto leaveThisPrimitive = OpenMesh::makeTemporaryProperty<typename MeshType::FaceHandle, bool>(destMesh, "leaveThisPrimitive");

        //handles and indices should be the same in copied mesh
        for (auto fh : facesToCopy)
        {
            leaveThisPrimitive[fh] = true;
        }

        for (auto fh : destMesh.faces())
        {
            if (!leaveThisPrimitive[fh])
            {
                destMesh.delete_face(fh, true);
            }
        }

        if (doGarbageCollection)
        {
            destMesh.garbage_collection();
        }
    }

    
    template<typename MeshType>
    void copyMeshPart(const MeshType& srcMesh, MeshType& destMesh, std::vector<typename MeshType::EdgeHandle> edgesToCopy, bool doGarbageCollection = false)
    {
        destMesh = srcMesh;

        auto leaveThisPrimitive = OpenMesh::makeTemporaryProperty<typename MeshType::EdgeHandle, bool>(destMesh, "leaveThisPrimitive");

        //handles and indices should be the same in copied mesh
        for (auto eh : edgesToCopy)
        {
            leaveThisPrimitive[eh] = true;
        }

        for (auto eh : destMesh.edges())
        {
            if (!leaveThisPrimitive[eh])
            {
                destMesh.delete_edge(eh, true);
            }
        }

        if (doGarbageCollection)
        {
            destMesh.garbage_collection();
        }
    }

    
    template<typename MeshType>
    void copyMeshPart(const MeshType& srcMesh, MeshType& destMesh, std::vector<typename MeshType::VertexHandle> verticesToCopy, bool doGarbageCollection = false)
    {
        destMesh = srcMesh;

        auto leaveThisPrimitive = OpenMesh::makeTemporaryProperty<typename MeshType::VertexHandle, bool>(destMesh, "leaveThisPrimitive");

        //handles and indices should be the same in copied mesh
        for (auto vh : verticesToCopy)
        {
            leaveThisPrimitive[vh] = true;
        }

        for (auto vh : destMesh.vertices())
        {
            if (!leaveThisPrimitive[vh])
            {
                destMesh.delete_vertex(vh, true);
            }
        }

        if (doGarbageCollection)
        {
            destMesh.garbage_collection();
        }
    }

    template<typename MeshType>
    typename MeshType::Scalar calcAverageEdgeLength(const MeshType& mesh)
    {
        using tScalar = typename MeshType::Scalar;

        tScalar avgLength(0);
        std::size_t edgeCount = 0;

        if (mesh.n_edges() == 0)
        {
            return avgLength;
        }

        for (auto eh : mesh.edges())
        {
            avgLength += mesh.calc_edge_length(eh);
            edgeCount++;
        }

        if (edgeCount)
        {
            avgLength /= edgeCount;
        }

        return avgLength;
    }

    template<typename MeshType>
    double cotan(const MeshType& mesh, typename MeshType::HalfedgeHandle heI)
    {

        if (mesh.is_boundary(mesh.edge_handle(heI)))
        {
            return typename MeshType::Scalar(0);
        }

        auto sHe = OpenMesh::make_smart(heI, mesh);

        auto he0 = sHe.next();
        auto he1 = he0.next();

        auto v0 = he0.from();
        auto v1 = he1.from(); //compute cotan at this vertex
        auto v2 = sHe.from();

        const auto& p0 = mesh.point(v0);
        const auto& p1 = mesh.point(v1);
        const auto& p2 = mesh.point(v2);

        auto edgeVector0 = p0 - p1;
        auto edgeVector1 = p2 - p1;

        return  dot(edgeVector0, edgeVector1) / DMB::length(cross(edgeVector0, edgeVector1));
    }

    template<typename MeshType>
    double  edgeCotanWeight(const MeshType& m, typename MeshType::EdgeHandle eh)
    {
        double sum = 0;
        for (auto he : OpenMesh::make_smart(eh, m).halfedges())
        {
            sum += DMB::cotan(m, he);
        }
        return sum;
    };

    template<typename MeshType>
    bool edgeIsDelaunay(const MeshType& m, typename MeshType::EdgeHandle eh, double delaunayCotanLimit = 1e-6)
    {
        double cWeight = DMB::edgeCotanWeight(m, eh);
        return (cWeight > -delaunayCotanLimit);
    };

    // Uses mesh's axis-aligned bounding box, scaled around its center.
    // percent = 100.0 -> exact bbox, 110.0 -> 10% larger (can produce outside points).
    template <typename MeshType>
    std::vector<typename MeshType::Point> samplePointsInBoundingBox(const MeshType& mesh, std::size_t sampleCount, double percent)
    {
        using Point = typename MeshType::Point;
        using Scalar = typename MeshType::Scalar;

        if (mesh.n_vertices() == 0) 
        {
            return {};
        }
        if (percent <= 0.0) 
        {
            return {};
        }

        std::mt19937_64 rng(std::random_device{}());

        auto bb = DMB::calcMeshBoundingBox<MeshType>(mesh);

        Point bbMin = bb.m_min;
        Point bbMax = bb.m_max;

        const Scalar scale = static_cast<Scalar>(percent / 100.0);
        const Point center = (bbMin + bbMax) * Scalar(0.5);
        const Point half = (bbMax - bbMin) * (Scalar(0.5) * scale);

        std::uniform_real_distribution<Scalar> dx(center[0] - half[0], center[0] + half[0]);
        std::uniform_real_distribution<Scalar> dy(center[1] - half[1], center[1] + half[1]);
        std::uniform_real_distribution<Scalar> dz(center[2] - half[2], center[2] + half[2]);

        std::vector<Point> samples;
        samples.reserve(sampleCount);

        for (std::size_t i = 0; i < sampleCount; ++i) {
            Point p;
            p[0] = dx(rng);
            p[1] = dy(rng);
            p[2] = dz(rng);
            samples.push_back(p);
        }

        return samples;
    }

    template<typename MeshType>
    double calcComponentSignedVolume(const MeshType& mesh, const std::vector<OpenMesh::SmartFaceHandle>& component)
    {

        //http://chenlab.ece.cornell.edu/Publication/Cha/icip01_Cha.pdf
        //https://www.ams.org/journals/mcom/1986-46-173/S0025-5718-1986-0815838-7/S0025-5718-1986-0815838-7.pdf
        //https://dsp.stackexchange.com/questions/7856/calculating-the-volume-of-a-triangular-mesh

        double volume = 0;

        using tPoint = typename MeshType::Point;
        tPoint centroid(0, 0, 0);

        {
            std::unordered_set<typename MeshType::VertexHandle> vertices;

            for (auto fh : component)
            {
                for (auto vh : fh.vertices())
                {
                    vertices.insert(vh);
                }
            }

            if (vertices.size() == 0)
            {
                return volume;
            }

            for (auto vh : vertices)
            {
                centroid += mesh.point(vh);
            }

            centroid /= vertices.size();
        }

        for (auto fh : component)
        {
            std::vector<tPoint> pts{};

            for (auto vh : fh.vertices())
            {
                pts.push_back(mesh.point(vh) - centroid);
            }

            typename MeshType::Scalar v = (
                -pts[2][0] * pts[1][1] * pts[0][2]
                + pts[1][0] * pts[2][1] * pts[0][2]
                + pts[2][0] * pts[0][1] * pts[1][2]
                - pts[0][0] * pts[2][1] * pts[1][2]
                - pts[1][0] * pts[0][1] * pts[2][2]
                + pts[0][0] * pts[1][1] * pts[2][2]
                );

            volume += v;
        }

        return volume * (1.0 / 6.0); // should be multiplied by (1.0 / 6.0), but it does not make a difference for my use case
    }
 
    template<typename MeshType>
    bigfloat calcComponentSignedVolumeExact(const MeshType& mesh, const std::vector<OpenMesh::SmartFaceHandle>& component)
    {
        using tExactPoint = std::array<bigfloat, 3>;

            bigfloat volume = 0;

            tExactPoint centroid;
            centroid[0] = bigfloat(0);
            centroid[1] = bigfloat(0);
            centroid[2] = bigfloat(0);


            {
                std::unordered_set<typename MeshType::VertexHandle> vertices;

                for (auto fh : component)
                {
                    for (auto vh : fh.vertices())
                    {
                        vertices.insert(vh);
                    }
                }


                if (vertices.size() == 0)
                {
                    return volume;
                }

                for (auto vh : vertices)
                {
                    centroid[0] = centroid[0] + bigfloat(mesh.point(vh)[0]);
                    centroid[1] = centroid[1] + bigfloat(mesh.point(vh)[1]);
                    centroid[2] = centroid[2] + bigfloat(mesh.point(vh)[2]);
                }

                centroid[0] = centroid[0] * bigfloat((1.0 / (double)vertices.size()));
                centroid[1] = centroid[1] * bigfloat((1.0 / (double)vertices.size()));
                centroid[2] = centroid[2] * bigfloat((1.0 / (double)vertices.size()));
            }


            for (auto fh : component)
            {
                std::vector<tExactPoint> pts;

                for (auto vh : fh.vertices())
                {
                    tExactPoint p;
                    p[0] = bigfloat(mesh.point(vh)[0]) - centroid[0];
                    p[1] = bigfloat(mesh.point(vh)[1]) - centroid[1];
                    p[2] = bigfloat(mesh.point(vh)[2]) - centroid[2];

                    pts.push_back(p);
                }

                bigfloat v = (
                    -pts[2][0] * pts[1][1] * pts[0][2]
                    + pts[1][0] * pts[2][1] * pts[0][2]
                    + pts[2][0] * pts[0][1] * pts[1][2]
                    - pts[0][0] * pts[2][1] * pts[1][2]
                    - pts[1][0] * pts[0][1] * pts[2][2]
                    + pts[0][0] * pts[1][1] * pts[2][2]
                    );

                volume = volume + v;
            }

            return volume * bigfloat(1.0 / 6.0); // should be multiplied by (1.0 / 6.0), but it does not make a difference for my use case
        
    }
    
    template<typename MeshType>
    int calcComponentVolumeSignExact(const MeshType& mesh, const std::vector<OpenMesh::SmartFaceHandle>& component)
    {
        return calcComponentSignedVolumeExact(mesh, component).sgn();
    }

    template<typename MeshType>
    double calcComponentArea(const MeshType& mesh, const std::vector<OpenMesh::SmartFaceHandle>& component)
    {
        double area = 0.0;

        for (auto fh : component)
        {
            area += mesh.calc_face_area(fh);
        }

        return area;
    }

    template<typename MeshType>
    double calcMeshSignedVolume(const MeshType& mesh)
    {
        double volume = 0;

        if (mesh.n_vertices() == 0)
        {
            return volume;
        }

        using tPoint = typename MeshType::Point;
        tPoint centroid(0, 0, 0);

        {
            for (auto vh : mesh.vertices())
            {
                centroid += mesh.point(vh);
            }

            centroid /= mesh.n_vertices();
        }

        for (auto fh : mesh.faces())
        {
            std::vector<tPoint> pts{};

            for (auto vh : fh.vertices())
            {
                pts.push_back(mesh.point(vh) - centroid);
            }

            typename MeshType::Scalar v = (
                -pts[2][0] * pts[1][1] * pts[0][2]
                + pts[1][0] * pts[2][1] * pts[0][2]
                + pts[2][0] * pts[0][1] * pts[1][2]
                - pts[0][0] * pts[2][1] * pts[1][2]
                - pts[1][0] * pts[0][1] * pts[2][2]
                + pts[0][0] * pts[1][1] * pts[2][2]
                );

            volume += v;
        }

        return volume * (1.0 / 6.0); 
    }

    template<typename MeshType>
    bigfloat calcMeshSignedVolumeExact(const MeshType& mesh)
    {
        using tExactPoint = std::array<bigfloat, 3>;

        bigfloat volume = 0;

        if (mesh.n_vertices() == 0)
        {
            return volume;
        }

        tExactPoint centroid;
        centroid[0] = bigfloat(0);
        centroid[1] = bigfloat(0);
        centroid[2] = bigfloat(0);


        {
            for (auto vh : mesh.vertices())
            {
                centroid[0] = centroid[0] + bigfloat(mesh.point(vh)[0]);
                centroid[1] = centroid[1] + bigfloat(mesh.point(vh)[1]);
                centroid[2] = centroid[2] + bigfloat(mesh.point(vh)[2]);
            }

            centroid[0] = centroid[0] * bigfloat((1.0 / (double)mesh.n_vertices()));
            centroid[1] = centroid[1] * bigfloat((1.0 / (double)mesh.n_vertices()));
            centroid[2] = centroid[2] * bigfloat((1.0 / (double)mesh.n_vertices()));
        }


        for (auto fh : mesh.faces())
        {
            std::vector<tExactPoint> pts;

            for (auto vh : fh.vertices())
            {
                tExactPoint p;
                p[0] = bigfloat(mesh.point(vh)[0]) - centroid[0];
                p[1] = bigfloat(mesh.point(vh)[1]) - centroid[1];
                p[2] = bigfloat(mesh.point(vh)[2]) - centroid[2];

                pts.push_back(p);
            }

            bigfloat v = (
                -pts[2][0] * pts[1][1] * pts[0][2]
                + pts[1][0] * pts[2][1] * pts[0][2]
                + pts[2][0] * pts[0][1] * pts[1][2]
                - pts[0][0] * pts[2][1] * pts[1][2]
                - pts[1][0] * pts[0][1] * pts[2][2]
                + pts[0][0] * pts[1][1] * pts[2][2]
                );

            volume = volume + v;
        }

        return volume * bigfloat(1.0 / 6.0); 
    }

    template<typename MeshType>
    int calcMeshVolumeSignExact(const MeshType& mesh)
    {
        return calcMeshSignedVolumeExact(mesh).sgn();
    }

    template<typename MeshType>
    double calcMeshArea(const MeshType& mesh)
    {
        double area = 0.0;

        for (auto fh : mesh.faces())
        {
            area += mesh.calc_face_area(fh);
        }

        return area;
    }

} //namespace