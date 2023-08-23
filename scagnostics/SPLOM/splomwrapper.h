#ifndef SPLOMWRAPPER_H
#define SPLOMWRAPPER_H

#include "dataset.h"
#include "Reader.h"
#include "statistics.h"
#include "splomthread.h"

#include "transferfunction.h"

#include <QObject>
#include <QThread>

#include <QTime>
#include <vector>
#include <stdio.h>
#include <iostream>
#include <utility>


class SPLOMWrapper : public QObject
{
    Q_OBJECT
public:
    explicit SPLOMWrapper(QObject *parent = 0);

    ~SPLOMWrapper();
    void Initialize(Reader* reader, Dataset* data, Statistics* stats, int numThreads, bool useAsDrawer);

    void ProcessData(int attributeToCluster);
    void ProcessData(int attributeToCluster, vector<pair<int, int> > filteredAttributes);

    void SetTransferFunction(TransferFunction* tf_){ tf = tf_; }

    void GetImage(int attr1, int attr2, int clusteringAttribute, int valuetoFilter, short texture[]);
    void GetDensityImage(int attr1, int attr2, int clusteringAttribute, int valuetoFilter, short texture[]);

    void CleanMemory();
    void ChangePointSize(float scale);

    float GetPointSize(){ return (3.0/128.0)*scale;  }

    bool IsAttributeNumerical(int clusteringAttribute){
        return _reader->IsNumerical(clusteringAttribute);
    }
    bool HasData(){ return imageData != NULL; }
    bool IsRunning(){ return isRunning;}
    bool IsFinished(){ return !isRunning;}
    int GetTotalNumericalAttributes();
    void GetListOfNumericalAttributes(vector<int>* list);
    Dataset* GetData(){ return _data;}
    Statistics* GetStats(){ return _stats;}

    int GetAttributeToCluster(){ return _attributeToCluster; }
    void CreateNumericalAttributes();
    void ChangeTypeBlend(int type){
        _typeBlend = type;
    }
    int GetTypeBlend(){ return _typeBlend;}

    void SetTextureSize(int size){ sizeT = size;}
    int GetSizeTexture(){ return sizeT;}
signals:
   void FinishedImageGeneration();

public slots:
    void FinishedGeneration();

private:
    Reader* _reader;
    Dataset* _data;
    Statistics* _stats;
    int _numThreads;
    short* imageData; // YUUUGE Image data array...
    // LESS YUUUGE Image data array, this one contains the number of elements per location
    // TODO- modify so it is always calculated (?)
    short* densityData;


    int _attributeToCluster;

    vector<QThread*> threads;
    vector<SPLOMThread*> splomThreads;
    int threadsFinished;

    bool isRunning;
    bool isInitialized;

    int _typeBlend;

    int sizeT;

    vector<int> numTypeAttributes; //
    vector< pair<int, int> > attributesToDraw;

    float scale;

    TransferFunction* tf;

    bool useAsDrawer;
};

#endif // SPLOMWRAPPER_H
