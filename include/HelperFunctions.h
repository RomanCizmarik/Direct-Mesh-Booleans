//TODO: License

#pragma once

#include "vector"
#include <cstdint>
#include "HelperStrucutres.h"

namespace DMB
{


    ////////////////////////////////////////////////////////////////////////////////////////////////////
    //! \brief   Adds whole mesh to the triangle soup.
    //! Could be a methos inside the TriangleSoup but I like to use simple structs to indicate there is no additional magic happening.
    //! 
    //! \param   soup        The triangle soup.
    //! \param   ingredient  Mesh that shall be added.
    //! \param   transform   Some transformation of the mesh. The triangle soup puts them all into one space.
    //! \param   label       Label of the mesh. This is used to remember which trangle belongs to which mesh.
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    //template <typename MeshType>
    //inline void addMesh(TriangleSoup& soup, const MeshType& ingredient, const TransformMatrix& transform, uint label)
    //{
    //    uint vertex_start = uint(std::size(soup.coordinates) / 3);

    //    for (auto vh : ingredient.vertices())
    //    {
    //        auto worldSpacePoint = transform * ingredient.point(vh);

    //        soup.coordinates.push_back(worldSpacePoint[0]);
    //        soup.coordinates.push_back(worldSpacePoint[1]);
    //        soup.coordinates.push_back(worldSpacePoint[2]);
    //    }

    //    for (auto fh : ingredient.faces())
    //    {
    //        for (auto vh : fh.vertices_ccw())
    //        {
    //            soup.triangles.push_back(vertex_start + vh.idx());
    //        }

    //        soup.labels.push_back(label);
    //    }
    //}

    inline void addMesh(TriangleSoup& soup, const InputTriangleMesh& ingredient)
    {
        uint vertex_start = uint(std::size(soup.coordinates) / 3);

        soup.coordinates.insert(soup.coordinates.end(), ingredient.coordinates.begin(), ingredient.coordinates.end());

        //soup.triangles.insert(soup.triangles.end(), ingredient.triangles.begin(), ingredient.triangles.end());

        for (int i = 0; i < ingredient.triangles.size(); i += 3)
        {
            soup.triangles.push_back(vertex_start + ingredient.triangles[i + 0]);
            soup.triangles.push_back(vertex_start + ingredient.triangles[i + 1]);
            soup.triangles.push_back(vertex_start + ingredient.triangles[i + 2]);
            soup.labels.push_back(ingredient.label);
        }

    }
} //namespace
