//TODO: License

#pragma once

#include "BoundingBox.h"

namespace DMB
{
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

}