//TODO: License

#pragma once

//stl
#include <memory> //unique_ptr
#include <mutex> 

//local
#include "FastWindingNumber.h"

namespace DMB
{
    template<typename MeshType>
    class AcceleratorAccessor
    {
    public:

        using tFWN = typename DMB::FastWindingNumber<MeshType>;

    public:
        AcceleratorAccessor() = default;
        AcceleratorAccessor(AcceleratorAccessor&&) = default;

        void setMesh(MeshType* mesh /*non-const, because this class updates normals*/)
        {
            m_mesh = mesh;

            setDirty();
        }

        tFWN* getFastWindingNumber()
        {
            updateFWNTree();

            return m_FWNTree.get();
        }

        int getAcceleratorVersion() const
        {
            return m_acceleratorVersion;
        }

        void setDirty()
        {
            m_FWNTreeMutex.lock();
            m_FWNTree.reset();
            m_FWNTreeMutex.unlock();

            m_acceleratorVersion++;
        }

    private:

        void updateFWNTree()
        {
            m_FWNTreeMutex.lock();

            if (!m_FWNTree)
            {
                //make sure normals are updated 
                //this could be optimized: we could rely on the fact that the normals are already updated, but then it would need to be gauranteed everywhere in the code and
                //we would need to pay attention to always have consistently updated normals
                m_mesh->update_normals();
                m_FWNTree = std::make_unique<tFWN>();
                m_FWNTree->build(*m_mesh);
            }

            m_FWNTreeMutex.unlock();
        }

    private:
        std::unique_ptr<tFWN> m_FWNTree;
        mutable std::recursive_mutex m_FWNTreeMutex;

        std::atomic<long> m_acceleratorVersion = -1;

        MeshType* m_mesh = nullptr; //non-const, because this class updates normals
    };
}//namespace