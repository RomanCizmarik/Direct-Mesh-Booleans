//TODO: License

#pragma once

#include "BoundingBox.h"

namespace DMB
{
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
}