//OpenMesh
//MeshIO header needs to be included before mesh type
#pragma warning(push)
#pragma warning(disable : 4723)
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#pragma warning(pop)

template <typename ScalarType>
struct OMTraitsImpl : public OpenMesh::DefaultTraits
{
    typedef OpenMesh::VectorT<ScalarType, 3> Point; // use float-values points
    typedef OpenMesh::VectorT<ScalarType, 3> Normal; // use float-values normals

    VertexAttributes(OpenMesh::Attributes::Status | OpenMesh::Attributes::Normal | OpenMesh::Attributes::Color | OpenMesh::Attributes::TexCoord2D);
    FaceAttributes(OpenMesh::Attributes::Status | OpenMesh::Attributes::Normal | OpenMesh::Attributes::Color);
    EdgeAttributes(OpenMesh::Attributes::Status | OpenMesh::Attributes::Normal);
    HalfedgeAttributes(OpenMesh::Attributes::Status);
};

using OMFloatTraits = OMTraitsImpl<float>;

using SimpleMesh = OpenMesh::TriMesh_ArrayKernelT<OMFloatTraits>;