//TODO: License

#pragma once

//STL
#include <vector>
#include <span>
#include <algorithm>
#include <stack>

//OpenMesh
#include <OpenMesh/Core/Utils/PropertyManager.hh>


//Local
#include <BoundingBox.h>
#include <Utils.h>


namespace DMB
{
#if 0
	template<typename MeshType>
	class AABBTree
	{
	public: //types

		using tScalar = typename MeshType::Scalar;
		using tPoint = typename MeshType::Point;
		using tBoundingBox = BoundingBoxT <tPoint>;
		using tFace = typename MeshType::FaceHandle;

		struct AABBNode
		{
			tBoundingBox bb;
			std::vector<tFace> faceHandles;
			AABBNode* left = nullptr;
			AABBNode* right = nullptr;
            int depth = -1;

			~AABBNode()
			{
				delete left;
				delete right;
			}

			bool isLeaf() const
			{
				return left == nullptr && right == nullptr;
			}
		};

	public: //methods
        AABBTree()
        {
        }

        ~AABBTree()
        {
            delete m_root;
        }

        void build(const MeshType& mesh)
        {
            m_root = buildRecursive(mesh, mesh.faces_begin(), mesh.faces_end(), 0);
        }

        template<typename FunctionType>
        void recurse(const FunctionType& f) const
        {
            recurseInternal(m_root, f);
        }

        void debug_outputAsMesh() const
        {
            int cnt = 0;

            auto f = [&cnt](const AABBNode* node)
                {
                    OpenMesh::IO::write_mesh(boundingBoxToMesh<MeshType, tBoundingBox>(node->bb), "C:/skola/PhD/Samples/booleans/tree/bb_" + std::to_string(cnt) + "_" +   std::to_string(node->depth) + ".obj");
                    ++cnt;
                };

            recurse(f);
        }

    public: //properties 
        AABBNode* m_root; //TODO: make accessor or something

    private: //methods

        template<typename FunctionType>
        void recurseInternal(const AABBNode* node, const FunctionType& f) const
        {
            
            if (node->left != nullptr)
            {
                recurseInternal(node->left, f);
            }

            if (node->right != nullptr)
            {
                recurseInternal(node->right, f);
            }

            f(node);
        }


        template<typename FaceIterator>
        AABBNode* buildRecursive(const MeshType& mesh, FaceIterator facesBegin, FaceIterator facesEnd, int depth)
        {
            std::cout << "depth: " << depth << std::endl;

            constexpr int facesPerNodeLimit = 4;

            auto node = new AABBNode();
            auto bb = calcFacesBoundingBox(mesh, facesBegin, facesEnd);
            node->bb = bb;
            node->depth = depth;

            node->faceHandles = std::vector<tFace>(facesBegin, facesEnd);

            if (std::distance(facesBegin, facesEnd) <= facesPerNodeLimit)
            {
                return node;
            }

            //// Split along longest axis
            //int axis;
            //auto maxCoeff = [&axis](const tPoint& diagonal)
            //    {
            //        tScalar max = std::numeric_limits<tScalar>::lowest();

            //        for (int i = 0; i < 3; ++i)
            //        {
            //            if (std::abs(diagonal[i]) > max)
            //            {
            //                axis = i;
            //                max = std::abs(diagonal[i]);
            //            }
            //        }
            //    };
            //
            //auto extents = node->bb.diagonal();
            //maxCoeff(extents);

            //const axis
            int axis = depth % 3;


            constexpr double eps = 1e-6;
            tScalar mid = node->bb.center()[axis];

            //if the extents of the bounding box in the given axis are too small, just take half the faces
            auto partitionIt = node->faceHandles.begin() + size_t(std::distance(node->faceHandles.begin(), node->faceHandles.end()) / 2);
                

            if (std::abs(bb.m_max[axis] - bb.m_min[axis]) > eps)
            {
                //otherwise, partition the faces based on mid point in given axis
                partitionIt = std::partition(node->faceHandles.begin(), node->faceHandles.end(),
                    [mid, axis, &mesh](tFace fh)
                    {
                        return mesh.calc_face_centroid(fh)[axis] <= mid; //TODO: precompute center
                    });
            }

            // Handle bad splits
            //if (leftSet.empty() || rightSet.empty())
            if(std::distance(partitionIt, node->faceHandles.begin()) == 0 || partitionIt == node->faceHandles.end())
            {
                return node;
            }

            node->left = buildRecursive(mesh, node->faceHandles.begin(), partitionIt, depth + 1);
            node->right = buildRecursive(mesh, partitionIt, node->faceHandles.end(), depth + 1);


            return node;
        }
	};
#endif

    template<typename MeshType>
    class AABBTree
    {
    public: //types

        using tScalar = typename MeshType::Scalar;
        using tPoint = typename MeshType::Point;
        using tBoundingBox = BoundingBoxT <tPoint>;
        using tFace = typename MeshType::FaceHandle;

        struct AABBTreeNode
        {
            std::span<tFace> faces;
            tBoundingBox bb;
            int depth = 0;
            int leftNode = -1;
            int rightNode = -1;

            bool isLeaf() const
            {
                return leftNode == -1 && rightNode == -1;
            }
        };


    public: //methods

        AABBTree(int maxNodeCount = 4) : 
            m_maxNodeCount(maxNodeCount)
        {

        }

        ~AABBTree()
        { }

        void build(const MeshType& mesh)
        {
            m_nodes.clear();

            m_faces.clear();
            m_faces.resize(mesh.n_faces());
            std::copy(mesh.faces_begin(), mesh.faces_end(), m_faces.begin());

            buildInternal(mesh);
        }

        const AABBTreeNode& node(size_t idx) const
        {
            assert(idx < m_nodes.size());

            return m_nodes[idx];
        }

        const std::vector<AABBTreeNode>& nodes() const
        {
            return m_nodes;
        }

        const std::span<tFace>& nodeFaces(size_t idx) const
        {
            assert(idx < m_nodes.size());

            return m_nodes[idx].faces;
        }

    private: // methods

        void buildInternal(const MeshType& mesh)
        {

            auto createNewNode = []()
                {

                };

            auto centroidProp = OpenMesh::makeTemporaryProperty<tFace, tPoint>(const_cast<MeshType&>(mesh));

            // compute face centroids
            for (auto fh : mesh.faces())
            {
                centroidProp[fh] = mesh.calc_face_centroid(fh);
            }

            AABBTreeNode rootNode;
            rootNode.faces = std::span<tFace>(m_faces.begin(), m_faces.end());
            rootNode.depth = 0;
            rootNode.bb = calcFacesBoundingBox(mesh, m_faces.begin(), m_faces.end());

            m_nodes.push_back(std::move(rootNode));

            std::stack< size_t > nodesToProcess;
            nodesToProcess.push(m_nodes.size() - 1);

            while (!nodesToProcess.empty())
            {
                auto nodeIdx = nodesToProcess.top();

                auto& node = m_nodes[nodeIdx];

                //node count is lower than requested max node count -> this is a leaf node
                if (node.faces.size() < m_maxNodeCount)
                {
                    //node processed
                    nodesToProcess.pop();
                    continue;
                }

                //parition faces based on given axis
                int axis = node.depth % 3;
                tScalar mid = 0;
                tScalar min =  std::numeric_limits<tScalar>::max();
                tScalar max = std::numeric_limits<tScalar>::lowest();

                for (auto fh : node.faces)
                {
                    const auto& centroid = centroidProp[fh];
                    mid += centroid[axis];
                    min = std::min(min, centroid[axis]);
                    max = std::max(max, centroid[axis]);
                }

                mid /= node.faces.size();

                constexpr auto eps = 1e-6;

                //if the extents of the bounding box in the given axis are too small, just take half the faces
                auto partitionIt = node.faces.begin() + size_t(std::distance(node.faces.begin(), node.faces.end()) / 2);

                if (std::abs(max - min) > eps)
                {
                    //otherwise, partition the faces based on mid point in given axis
                    partitionIt = std::partition(node.faces.begin(), node.faces.end(),
                        [axis, mid, &centroidProp](tFace fh)
                        {
                            return centroidProp[fh][axis] <= mid;
                        });
                }

                AABBTreeNode nodeLeft;
                nodeLeft.faces = std::span<tFace>(node.faces.begin(), partitionIt);
                nodeLeft.depth = node.depth + 1;
                nodeLeft.bb = calcFacesBoundingBox(mesh, nodeLeft.faces.begin(), nodeLeft.faces.end());

                AABBTreeNode nodeRight;
                nodeRight.faces = std::span<tFace>(partitionIt, node.faces.end());
                nodeRight.depth = node.depth + 1;
                nodeRight.bb = calcFacesBoundingBox(mesh, nodeRight.faces.begin(), nodeRight.faces.end());

                //this node processed
                nodesToProcess.pop();

                //add new nodes
                m_nodes.push_back(std::move(nodeLeft));
                nodesToProcess.push(m_nodes.size() - 1);
                m_nodes[nodeIdx].leftNode = m_nodes.size() - 1;

                m_nodes.push_back(std::move(nodeRight));
                nodesToProcess.push(m_nodes.size() - 1);
                m_nodes[nodeIdx].rightNode = m_nodes.size() - 1;
            }

        }

    private: //properties
        std::vector<AABBTreeNode> m_nodes;
        std::vector<tFace> m_faces;

        int m_maxNodeCount;
    };
}

#include "AABBTree.hxx"