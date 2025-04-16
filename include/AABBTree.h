//TODO: License

namespace DMB
{
	template<typename MeshType>
	class AABBTree
	{
	public:
		AABBTree();
		~AABBTree();

		void build(const MeshType& mesh);

	private:

	};
}

#include "AABBTree.hxx"