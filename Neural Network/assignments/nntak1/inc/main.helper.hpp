#ifndef MAIN_HELPER_HPP
#define MAIN_HELPER_HPP

#include "stdafx.hpp"
#include <iostream>
#include "Eigen/Core"

template<typename T, typename = typename std::enable_if<!std::is_same<T, Eigen::VectorXf>::value>::type>
ostream& operator<<(ostream& os, vector<T> v) {
    for(auto& i : v) cout<<i<<" ";
    return os;
}
ostream& operator<<(ostream& os, vector<Eigen::VectorXf> v) {
    int c = 0;
    for(auto& i : v) cout<<"["<<c++<<"] "<<i.transpose()<<endl;
    return os;
}
ostream& operator<<(ostream& os, vector<Eigen::VectorXf*> v) {
    int c = 0;
    for(auto& i : v) cout<<"["<<c++<<"] "<<i->transpose()<<endl;
    return os;
}

bool operator<(const Eigen::VectorXf& v, scalar s) {
    for(int i = 0; i < v.rows(); i++) if(v[i] >= s) return false;
    return true;
}

inline string get_date_time() {
  time_t rawtime;
  struct tm * timeinfo;
  char buffer[80];

  time (&rawtime);
  timeinfo = localtime(&rawtime);

  strftime(buffer,80,"%d-%m-%Y %H-%M-%S",timeinfo);
  return std::string(buffer);
}

#endif // MAIN_HELPER_HPP

