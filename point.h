#ifndef POINT_H
#define POINT_H

#include <stddef.h>
#include <math.h>

#include <algorithm>
#include <vector>

template <size_t N>
class PointV
{
    /***/
    float nPoints[N];
    //std::vector<float> nPoints;

public:
    //PointV(float)

    /***/
    float& operator[](const size_t index)
    {   return nPoints[index];  }

    float operator [](const size_t index) const
    {   return nPoints[index];  }

    /***/
    size_t size() const
    {   return N;   }

    typedef float* iterator;
    typedef const float* const_iterator;

    /***/
    iterator begin()
    {   return nPoints; }

    /***/
    iterator end()
    {   return begin() + size();    }

    const_iterator begin() const;

    const_iterator end() const;
};

template <size_t N>
typename PointV<N>::const_iterator PointV<N>::begin() const
{
    return nPoints;
}

template <size_t N>
typename PointV<N>::const_iterator PointV<N>::end() const
{
    return begin() + size();
}

template <size_t N>
float Distance(const PointV<N>&p1, const PointV<N>&p2)
{
    float dist = 0.f;
    for(size_t i = 0; i<N; ++i)
        dist += (p2[i] - p1[i]) * (p2[i] - p1[i]);

    return sqrtf(dist);
}

template <size_t N>
bool operator==(const PointV<N>&p1, const PointV<N>&p2)
{
    return std::equal(p1.begin(), p1.end(), p2.begin());
}

template <size_t N>
bool operator!=(const PointV<N>&p1, const PointV<N>&p2)
{
    return !(p1 == p2);
}


#endif // POINT_H

