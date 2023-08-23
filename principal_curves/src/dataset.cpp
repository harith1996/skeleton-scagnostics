
#include "dataset.hpp"
#include <sstream>

DataSet::DataSet(size_t d) : dimensions { d } {
    stats.resize(d);
}

float DataSet::operator()(size_t dimension, size_t index) const {
    return data[index*dimensions + dimension];
}

void DataSet::append(const std::vector<float>& row) {
    if(row.size() != dimensions) { fprintf(stderr, "Row of wrong size"); return; }

    data.insert(data.end(), std::begin(row), std::end(row));

    for(int i = 0; i < dimensions; i++) {
        stats[i].put(row[i]);
    }
}

void DataSet::modifyRangeInStat(int dimension, float minRange, float maxRange){
    stats[dimension].setNewRange(minRange, maxRange);
}


void DataSet::DoubleRangeInStat(int dimension, bool debug){

    double maxV = stats[dimension].maxValue;
    double minV = stats[dimension].minValue;

    double range = maxV - minV;

    if ( debug){
       // std::cout << "Range " << range << std::endl;
    }
    double maxRange = maxV + 0.5*range;
    double minRange = minV - 0.5*range;
    stats[dimension].setNewRange(minRange, maxRange);
}

void DataSet::CheckRange(int dimension){
    double maxV = stats[dimension].maxValue;
    double minV = stats[dimension].minValue;
    double range = maxV - minV;

    if (range < 0.000000001){

        float val = data[0*dimensions + dimension];
        stats[dimension].setNewRange(val-0.5, val + 0.5);

    }

}

Stat<float> DataSet::getStat(size_t index) const {
    return stats[index];
}

size_t DataSet::size() const {
    if(dimensions == 0) return 0;
    return data.size() / dimensions;
}

size_t DataSet::dimension() const {
    return dimensions;
}


void DataSet::addRandomPoints(float percent){

    int totalElements = size();
    int toAdd = percent*totalElements;

    std::random_device rd;

    std::mt19937 e2(rd());
    std::uniform_real_distribution<> dist(0, 1);


    for (int i = 0; i < toAdd; i++) {
         std::vector<float> row;
         for(int k = 0; k < dimensions; k++){
              float r = getStat(k).maxValue - getStat(k).minValue;

              if ( r < 0.0001){
                  row.push_back(getStat(k).maxValue);// doesn't matter which
              }
              else {
                  float v = r*dist(e2) + getStat(k).minValue;
                  row.push_back(v);
              }
          }
         append(row);
    }
}


void DataSet::jitterRandomPoints(float percent){
    //We just need to remove random elements
    int totalElements = size();
    int toRemove = percent*totalElements;

    //first get all the indices of elements to remove
    //then create a temporary data vector with those removed...
    std::set<int> indicesToJitter;

    std::random_device rd;

    std::mt19937 e2(rd());
    std::uniform_real_distribution<> dist(0, 1);
    std::uniform_real_distribution<> dist2(-1, 1);

    while(indicesToJitter.size() < toRemove){
          int newIdx = dist(e2)*totalElements;
          indicesToJitter.insert(newIdx);
    }

    std::vector< std::vector<float> > withJitteredPoints;

    for(int i = 0; i < totalElements; i++){
        if ( count(indicesToJitter.begin(), indicesToJitter.end(), i) == 0){

            std::vector<float> newElement;
            for(int k = 0; k < dimensions; k++)
            {
                float val = data[i*dimensions + k];
                newElement.push_back(val);
            }
            withJitteredPoints.push_back(newElement);
        }
        else {
            //the new element is the previous plus jitter

            std::vector<float> newElement;
            for(int k = 0; k < dimensions; k++)
            {
                float originalval = data[i*dimensions + k];
                float r = getStat(k).maxValue - getStat(k).minValue;
                if ( r < 0.0001) r = 0.0001;
                float maxShift = 0.01*r;
                newElement.push_back(originalval + maxShift*dist2(e2));
            }
            withJitteredPoints.push_back(newElement);
        }
    }
    stats.clear();
    stats.resize(dimensions);

    for(auto newPoint:withJitteredPoints){
        append(newPoint);
    }







}

void DataSet::removeRandomPoints(float percent){
    //We just need to remove random elements
    int totalElements = size();
    int toRemove = percent*totalElements;

    //first get all the indices of elements to remove
    //then create a temporary data vector with those removed...
    std::set<int> indicesToRemove;

    std::random_device rd;

    std::mt19937 e2(rd());
    std::uniform_real_distribution<> dist(0, 1);

    while(indicesToRemove.size() < toRemove){
          int newIdx = dist(e2)*totalElements;
          indicesToRemove.insert(newIdx);
    }

    std::vector< std::vector<float> > notRemovedPoints;

    for(int i = 0; i < totalElements; i++){
        if ( count(indicesToRemove.begin(), indicesToRemove.end(), i) == 0){

            std::vector<float> newElement;
            for(int k = 0; k < dimensions; k++)
            {
                float val = data[i*dimensions + k];
                newElement.push_back(val);
            }
            notRemovedPoints.push_back(newElement);
        }
    }
    stats.clear();
    stats.resize(dimensions);

    for(auto notRemoved: notRemovedPoints){
        append(notRemoved);
    }

}


#include <iostream>

void DataSet::SaveDataset(const char* filename){
   // each line is the number of dimensions....comma separated values
   int numElements = size();

   std::stringstream ss;
   for(int i = 0; i < numElements; i++){
       for(int j = 0; j < dimensions; j++){

           float val = data[i*dimensions + j];

           ss << val;
           if ( j < dimensions -1)
               ss << ",";
       }

       if ( i < numElements -1)
           ss << "\n";
   }

   //std::cout << "Saved? " << filename << "," << numElements << std::endl;
   std::ofstream out(filename);
   out << ss.str();
   out.close();
}

void DataSet::SaveNormalizedDataset(const char* filename){
   // each line is the number of dimensions....comma separated values
   int numElements = size();

   std::stringstream ss;
   for(int i = 0; i < numElements; i++){
       for(int j = 0; j < dimensions; j++){

           auto currentStat = stats[j];
           float val = data[i*dimensions + j];
           val = (val - currentStat.minValue)/(currentStat.maxValue - currentStat.minValue);

           ss << val;
           if ( j < dimensions -1)
               ss << ",";
       }

       if ( i < numElements -1)
           ss << "\n";
   }

   //std::cout << "Saved? " << filename << "," << numElements << std::endl;
   std::ofstream out(filename);
   out << ss.str();
   out.close();
}

void DataSet::SaveNormalizedRotatedDataset(const char* filename, int angle, int firstAxis, int secondAxis){
   // each line is the number of dimensions....comma separated values
   int numElements = size();

   std::stringstream ss;
   float min1 = 10;
   float min2 = 10;


   float values[numElements][2];


   for(int i = 0; i < numElements; i++){

      auto stat1 = stats[firstAxis];
      auto stat2 = stats[secondAxis];

      float val1 = data[i*dimensions + firstAxis];
      float val2 = data[i*dimensions + secondAxis];

       val1 = (val1 - stat1.minValue)/(stat1.maxValue - stat1.minValue);
       val2 = (val2 - stat1.minValue)/(stat2.maxValue - stat2.minValue);

       float newV1 = val1*cos(angle) - val2*sin(angle);
       float newV2 = val1*sin(angle) + val2*cos(angle);

       val1 = newV1;
       val2 = newV2;


       values[i][0] = val1;
       values[i][1] = val2;

       if(min1 > val1)
           min1 = val1;

       if(min2 > val2)
           min2 = val2;


       /*
       ss << val1 <<"," << val2;

       if ( i < numElements -1)
           ss << "\n";*/
   }



   for(int i = 0; i < numElements; i++){
       float val1 = values[i][0] -min1;
       float val2 = values[i][1] -min2;
       ss << val1 <<"," << val2;

       if ( i < numElements -1)
           ss << "\n";
   }

  //std::cout << "Saved? " << filename << "," << numElements << std::endl;
   std::ofstream out(filename);
   out << ss.str();
   out.close();
}
