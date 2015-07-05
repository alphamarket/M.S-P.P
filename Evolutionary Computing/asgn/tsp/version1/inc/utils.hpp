#ifndef UTILS_HPP
#define UTILS_HPP

#include <cmath>
#include "population.hpp"
#include <boost/serialization/base_object.hpp>

namespace TSPGA {
    struct coordinate
    {
      double longitude, latitude;
      /**
       * @brief Init coordinate
       */
      coordinate() : longitude(-1), latitude(-1) { /* ctor */ }
      coordinate(double longitude, double latitude) : longitude(longitude), latitude(latitude) { /* ctor */ }
      /**
       * Serialization point
       */
#pragma GCC diagnostic ignored "-Wunused-parameter"
      friend class boost::serialization::access;
      template<class Archive> void serialize(Archive & ar, const unsigned int version) { ar  & longitude & latitude; }
#pragma GCC diagnostic pop
    };
    typedef vector<coordinate> coordinates;
}

#endif // UTILS_HPP
