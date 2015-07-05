#ifndef DATA_HPP
#define DATA_HPP
#include "utils.hpp"
#include <string>
#include <fstream>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
using namespace std;
namespace TSPGA {
    typedef unsigned long long distance;
    class data
    {
        /**
         * @brief The data coordinates
         */
        coordinates* _data;
        /**
         * @brief The data type
         */
        std::string _type;
    public:
        /**
         * @brief construct an empty data
         */
        data() { this->_data = new coordinates; }
        /**
         * @brief construct a data of coordinates
         */
        data(coordinates* _data) __nonnull() : _data(_data) { }
        /**
         * @brief ~data
         */
        virtual ~data() { delete this->_data; }
        /**
         * @brief Set data type
         */
        void type(const char* const _type) __nonnull()  { this->_type = _type; boost::to_upper(this->_type); }
        /**
         * @brief Get data type
         */
        std::string type()                      const   { return this->_type;  }
        /**
         * @brief Get a data at an index
         * @param index The index of datas
         */
        inline TSPGA::coordinate operator[] (size_t index)  const { return this->_data->at(index); }
        /**
         * @brief Set data (Note that this function WILL discard any previously binded data)
         * @param _data The target data to bind
         */
        inline void set     (coordinates*const _data) __nonnull() { delete this->_data; this->_data = _data; }
        /**
         * @brief Get a data at an index
         * @param index The index of datas
         */
        inline TSPGA::coordinate        get (size_t index)  const { return (*this)[index]; }
        /**
         * @brief Get size of data
         */
        inline size_t size() const { return this->_data->size(); }
        /**
         * @brief Get distance according to the data type between two coordinatess
         * @return The distance between coordinate i & j
         */
        distance getDistance(const size_t i, const size_t j) const {
#           define FOR_TYPE(x) if(this->type() == x)
            FOR_TYPE("GEO")
                return this->getDistance_GEO(i, j);
            FOR_TYPE("EUC_2D")
                return this->getDistance_EUC_2D(i, j);
            throw std::invalid_argument("There is no equivalent distance measure for type `" + this->_type + "`");
#           undef FOR_TYPE
        }
        /**
         * @brief Get geometric distance between two coordinatess
         * @return The distance between coordinate i & j
         */
        distance getDistance_GEO(const size_t i, const size_t j) const {
             double lati, latj, longi, longj;
             double q1, q2, q3, q4, q5;

             lati = M_PI * this->get(i).latitude / 180.0;
             latj = M_PI * this->get(j).latitude / 180.0;

             longi = M_PI * this->get(i).longitude / 180.0;
             longj = M_PI * this->get(j).longitude / 180.0;

             q1 = cos (latj) * sin(longi - longj);
             q3 = sin((longi - longj)/2.0);
             q4 = cos((longi - longj)/2.0);
             q2 = sin(lati + latj) * q3 * q3 - sin(lati - latj) * q4 * q4;
             q5 = cos(lati - latj) * q4 * q4 - cos(lati + latj) * q3 * q3;
             return (distance) (6378.3880 * atan2(sqrt(q1*q1 + q2*q2), q5) + 1.0);
        }
        /**
         * @brief Get euclidean distance between two coordinatess
         * @return The distance between coordinate i & j
         */
        distance getDistance_EUC_2D(const size_t i, const size_t j) const {
            auto xd = this->get(i).longitude - this->get(j).longitude;
            auto yd = this->get(i).latitude  - this->get(j).latitude;
            return (distance) (sqrt(xd * xd + yd * yd));
        }
        /**
         * Serialization point
         */
#pragma GCC diagnostic ignored "-Wunused-parameter"
        friend class boost::serialization::access;
        template<class Archive> void serialize(Archive & ar, const unsigned int version) { ar  & _data & _type; }
#pragma GCC diagnostic pop
        /**
         * @brief Deserializes a data instance from a file
         * @param file The data file
         * @return The data instance
         */
        static data* deserialize(const char* file) __nonnull() {
            // create data instance
            data* _data = new data();
            // open the input file
            std::ifstream ifs( file );
            // bind the input file to an archiver
            boost::archive::text_iarchive ar( ifs );
            // deserialize the data
            ar & _data;
            // close input file
            ifs.close();
            // return the deserialized data
            return _data;
        }
        /**
         * @brief Serializes passed data into a file
         * @param _data The data to serialize
         * @param file The data file
         */
        static void serialize(const data* const _data, const char* file) __nonnull((2)) {
            // open the output file
            std::ofstream ofs(file);
            // bind the output file to an archiver
            boost::archive::text_oarchive ar(ofs);
            // deserialize the data
            ar & _data;
            // close output file
            ofs.close();
        }
    };
}
#endif // DATA_HPP
