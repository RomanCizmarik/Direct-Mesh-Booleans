//TODO: License

#pragma once

#include "BoundingBox.h"

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
} //namespace