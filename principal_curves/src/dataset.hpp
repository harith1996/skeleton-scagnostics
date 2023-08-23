#ifndef DATASET_HPP
#define DATASET_HPP

#include <fstream>
#include <vector>
#include <random>

#include "stat.hpp"
#include "multival.hpp"
#include "curve2d.hpp"
#include "plotanalysis.hpp"
#include <iostream>
#include <set>

struct SkeletonNode {
   int label;
   Int2 location;
   int idxInGraph = -1;
   //std::vector<int> neighborsNodes; // labels of neighbor nodes with the counting
   // > technically the same info so there's no need to repeat ...
   std::vector<int> indexNeighborNode; //index in graph with the counting
   //std::vector<bool> antecesor; //whether the node is in the parents direction or not
   std::vector<SkeletonNode*> children;
   std::vector<SkeletonNode*> parents;
   RegionStat stat;
   RegionStat stat2; //using full float resolution

   int countInPath = 0;
   bool endPoint = false;
   bool bifurcation = false;
   Float2 direction = 0;
   float amount = 0;
   float initialDt = 0;
   float amountPerLength = 0;
   float orthoVariance = 0;

   void Clean(){
       children.clear();
       parents.clear();
   }

   void PrintInfo(){
        std::cout << "Label: "<< label << "\t Idx in Graph  " << idxInGraph << std::endl;

        std::cout <<"Neighbor labels and indices >" <<std::endl;
        /*for(int k = 0; k < neighborsNodes.size(); k++){
            std::cout << "(" << neighborsNodes.at(k) << "," << indexNeighborNode.at(k) << ") " ;
        }*/
        std::cout << std::endl;
   }
};

class DataSet {
    size_t dimensions;
    std::vector<float> data;
    std::vector<Stat<float>> stats;

    public:
    DataSet(size_t d);

    float operator()(size_t dimension, size_t index) const;

    void append(const std::vector<float>& row);

    void SaveDataset(const char* filename);
    void SaveNormalizedDataset(const char* filename);
    void SaveNormalizedRotatedDataset(const char* filename, int angle, int firstAxis, int secondAxis);

    void addRandomPoints(float percent);
    void removeRandomPoints(float percent);
    void jitterRandomPoints(float percent);

    //void extend();

    Stat<float> getStat(size_t index) const;

    void modifyRangeInStat(int dimension, float minRnage, float maxRange);
    void DoubleRangeInStat(int dimension, bool debug=false);
    void CheckRange(int dimension);
    size_t size() const;

    size_t dimension() const;
};
#endif
