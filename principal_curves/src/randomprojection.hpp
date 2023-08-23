#ifndef RANDOMPROJECTION_HPP
#define RANDOMPROJECTION_HPP

#include <boost/random/taus88.hpp>
#include <boost/random/normal_distribution.hpp>

#include "projection2d.hpp"
#include "dataset.hpp"

class ProjectionGenerator {
    boost::random::taus88 generator;
    boost::random::normal_distribution<float> distribution;

    public:
    ProjectionMatrix2D generate(const DataSet& dataSet);
};

#endif
