//TODO: License

#pragma once

#include <array>
#include <cassert>

//inspired by OpenSceneGraph implementation

namespace DMB
{
        template<typename VectorType>
        class BoundingBoxT
        {
        public:

            /** Minimum extent. (Smallest X, Y, and Z values of all coordinates.) */
            VectorType m_min;
            /** Maximum extent. (Greatest X, Y, and Z values of all coordinates.) */
            VectorType m_max;

        public:
            using ValueType = decltype(m_min[0] * m_min[0]);
            using tSelf = BoundingBoxT<VectorType>;

            //This is not redundancy it's used to get internal vector type.
            using tVectorType = VectorType;

        public: //methods

            inline BoundingBoxT() :
                m_min(std::numeric_limits<ValueType>::max(),
                    std::numeric_limits<ValueType>::max(),
                    std::numeric_limits<ValueType>::max()),
                m_max(std::numeric_limits<ValueType>::lowest(),
                    std::numeric_limits<ValueType>::lowest(),
                    std::numeric_limits<ValueType>::lowest())
            {
            }

            inline BoundingBoxT(const BoundingBoxT& bb) :
                m_min(bb.m_min),
                m_max(bb.m_max)
            {
            }

            /** Creates a bounding box initialized to the given extents. */
            inline BoundingBoxT(ValueType xmin, ValueType ymin, ValueType zmin,
                ValueType xmax, ValueType ymax, ValueType zmax) :
                m_min(xmin, ymin, zmin),
                m_max(xmax, ymax, zmax)
            {
            }

            /** Creates a bounding box initialized to the given extents. */
            inline BoundingBoxT(const VectorType& min, const VectorType& max) :
                m_min(min),
                m_max(max)
            {
            }

            ~BoundingBoxT()
            {
            }

            //! Assignment
            BoundingBoxT& operator=(const BoundingBoxT& other)
            {
                m_min = other.m_min;
                m_max = other.m_max;

                return *this;
            }
            
            inline void init()
            {
                m_min = { std::numeric_limits<ValueType>::max(),
                         std::numeric_limits<ValueType>::max(),
                         std::numeric_limits<ValueType>::max() };
                m_max = { std::numeric_limits<ValueType>::lowest(),
                         std::numeric_limits<ValueType>::lowest(),
                         std::numeric_limits<ValueType>::lowest() };
            }

            inline bool operator == (const BoundingBoxT& rhs) const { return m_min == rhs.m_min && m_max == rhs.m_max; }
            inline bool operator != (const BoundingBoxT& rhs) const { return m_min != rhs.m_min || m_max != rhs.m_max; }

            
            inline bool valid() const
            {
                return m_max[0] >= m_min[0] && m_max[1] >= m_min[1] && m_max[2] >= m_min[2];
            }

            inline ValueType& xMin() { return m_min[0]; }
            inline ValueType xMin() const { return m_min[0]; }

            inline ValueType& yMin() { return m_min[1]; }
            inline ValueType yMin() const { return m_min[1]; }

            inline ValueType& zMin() { return m_min[2]; }
            inline ValueType zMin() const { return m_min[2]; }

            inline ValueType& xMax() { return m_max[0]; }
            inline ValueType xMax() const { return m_max[0]; }

            inline ValueType& yMax() { return m_max[1]; }
            inline ValueType yMax() const { return m_max[1]; }

            inline ValueType& zMax() { return m_max[2]; }
            inline ValueType zMax() const { return m_max[2]; }

            /** Calculates and returns the bounding box center. */
            inline const VectorType center() const
            {
                assert(valid());
                return  VectorType((m_min[0] + m_max[0]) * ValueType(0.5),
                    (m_min[1] + m_max[1]) * ValueType(0.5),
                    (m_min[2] + m_max[2]) * ValueType(0.5));
            }

            /** Calculates and returns the bounding box radius. */
            inline ValueType radius() const
            {
                assert(valid());
                return sqrt(radius2());
            }

            /** Calculates and returns the squared length of the bounding box radius.
            * Note, radius2() is faster to calculate than radius(). */
            inline ValueType radius2() const
            {
                assert(valid());
                return ValueType(0.25) * length2(m_max - m_min);
            }

            /** Returns a specific corner of the bounding box.
            * pos specifies the corner as a number between 0 and 7.
            * Each bit selects an axis, X, Y, or Z from least- to
            * most-significant. Unset bits select the minimum value
            * for that axis, and set bits select the maximum. */
            inline const VectorType corner(unsigned int pos) const
            {
                assert(valid());
                return VectorType(pos & 1 ? m_max[0] : m_min[0], pos & 2 ? m_max[1] : m_min[1], pos & 4 ? m_max[2] : m_min[2]);
            }

            //   y
            //
            //   6 --- 7
            //  /|  z /|
            // 2 --- 3 |
            // | 4 --| 5
            // |/    |/
            // 0 --- 1     x
            std::array<VectorType, 8> corners() const
            {
                assert(valid());
                std::array<VectorType, 8> vertices;

                for (int i = 0; i < 8; ++i)
                {
                    vertices[i] = corner(i);
                }

                return vertices;
            }

            /** Expands the bounding box to include the given coordinate.
            * If the box is uninitialized, set its min and max extents to v. */
            inline void expandBy(const VectorType& v)
            {
                if (v[0] < m_min[0]) m_min[0] = v[0];
                if (v[0] > m_max[0]) m_max[0] = v[0];

                if (v[1] < m_min[1]) m_min[1] = v[1];
                if (v[1] > m_max[1]) m_max[1] = v[1];

                if (v[2] < m_min[2]) m_min[2] = v[2];
                if (v[2] > m_max[2]) m_max[2] = v[2];
            }

            /** Expands the bounding box to include the given coordinate.
            * If the box is uninitialized, set its min and max extents to
            * Vec3(x,y,z). */
            inline void expandBy(ValueType x, ValueType y, ValueType z)
            {
                if (x < m_min[0]) m_min[0] = x;
                if (x > m_max[0]) m_max[0] = x;

                if (y < m_min[1]) m_min[1] = y;
                if (y > m_max[1]) m_max[1] = y;

                if (z < m_min[2]) m_min[2] = z;
                if (z > m_max[2]) m_max[2] = z;
            }

            /** Expands this bounding box to include the given bounding box.
            * If this box is uninitialized, set it equal to bb. */
            void expandBy(const tSelf& bb)
            {
                assert(bb.valid());

                if (bb.m_min[0] < m_min[0]) m_min[0] = bb.m_min[0];
                if (bb.m_max[0] > m_max[0]) m_max[0] = bb.m_max[0];

                if (bb.m_min[1] < m_min[1]) m_min[1] = bb.m_min[1];
                if (bb.m_max[1] > m_max[1]) m_max[1] = bb.m_max[1];

                if (bb.m_min[2] < m_min[2]) m_min[2] = bb.m_min[2];
                if (bb.m_max[2] > m_max[2]) m_max[2] = bb.m_max[2];
            }

            tSelf expandBy(const tSelf& bb) const
            {
                tSelf out{ *this };

                out.expandBy(bb);

                return out;
            }


            /** Returns the intersection of this bounding box and the specified bounding box. */
            tSelf intersect(const tSelf& bb) const
            {
                assert(valid());
                assert(bb.valid());

                return tSelf(std::max(xMin(), bb.xMin()), std::max(yMin(), bb.yMin()), std::max(zMin(), bb.zMin()),
                    std::min(xMax(), bb.xMax()), std::min(yMax(), bb.yMax()), std::min(zMax(), bb.zMax()));
            }

            /** Return true if this bounding box intersects the specified bounding box. */
            bool intersects(const tSelf& bb) const
            {
                assert(valid());
                assert(bb.valid());

                return  (std::max(xMin(), bb.xMin()) <= std::min(xMax(), bb.xMax()) &&
                    std::max(yMin(), bb.yMin()) <= std::min(yMax(), bb.yMax()) &&
                    std::max(zMin(), bb.zMin()) <= std::min(zMax(), bb.zMax()));

            }

            /** Returns true if this bounding box contains the specified coordinate. */
            inline bool contains(const VectorType& v) const
            {
                assert(valid());

                return (v[0] >= m_min[0] && v[0] <= m_max[0]) &&
                    (v[1] >= m_min[1] && v[1] <= m_max[1]) &&
                    (v[2] >= m_min[2] && v[2] <= m_max[2]);
            }

            /** Returns true if this bounding box contains the specified bounding box. */
            inline bool contains(const tSelf& bb) const
            {
                assert(valid());

                return (bb.m_min[0] >= m_min[0] && bb.m_max[0] <= m_max[0]) &&
                    (bb.m_min[1] >= m_min[1] && bb.m_max[1] <= m_max[1]) &&
                    (bb.m_min[2] >= m_min[2] && bb.m_max[2] <= m_max[2]);
            }

            /** Returns true if this bounding box contains the specified coordinate allowing for specific epsilon. */
            inline bool contains(const VectorType& v, ValueType epsilon) const
            {
                assert(valid());

                return ((v[0] + epsilon) >= m_min[0] && (v[0] - epsilon) <= m_max[0]) &&
                    ((v[1] + epsilon) >= m_min[1] && (v[1] - epsilon) <= m_max[1]) &&
                    ((v[2] + epsilon) >= m_min[2] && (v[2] - epsilon) <= m_max[2]);
            }

            inline VectorType diagonal() const
            {
                assert(valid());

                return m_max - m_min;
            }
        };
} //namepsace