#ifndef COOCCURRENCECALCULATOR_H
#define COOCCURRENCECALCULATOR_H

/*
   Instead of making a wrapper
   We use the an original cooccurence calculator
   To create threads for other co-occurence data .
*/

#include <QObject>
#include <QThread>
#include "dataset.h"
#include "Reader.h"
#include "statistics.h"

class CooccurrenceCalculator : public QObject {

    Q_OBJECT
public:
    explicit CooccurrenceCalculator(QObject *parent = 0);
    ~CooccurrenceCalculator();

    // We just need the reader/mapper for figuring out the
    // type of data & the dataset itself
    void InitializeAsParent(Reader* reader, Dataset* data, int numThreads);
    void CalculateCramerV(Reader *reader, Statistics *stats);

    float GetCramerV(int attr1, int attr2);

    // This assumes a linear search
    // where each categorical dimension goes in order
    // and only there index are set
    int GetOccurenceNumberOfAttrAndCategory(int pos){ return totalOccurences.at(pos);}
    string GetNameOfAttrAndCategory(int pos){ return names.at(pos); }
    string GetNameOfAttr(int pos, bool direct = false);
    void GetCategoricalIndices(vector<int>* indices);
    void GetAssociationsInfo(vector< pair<int,int> >* nominalPairs, vector<float> * associations);

    float GetMaxNumberOccurrences(){ return max_numOccurences; }
    int GetNumberOfCategoricalDims();
    int GetTotalSumCategoricalSize(int sizes[]);

    void Initialize(Reader* reader, Dataset* data, int numThreads);
    void SetInfo(vector<pair<int, int> > attributes);
    void GetInfo(int data[], int from, int to);
    void Clean();
    void CleanMemory();
public slots:
    void process();
    void FinishedGeneration();
signals:
    void allThreadsDone();
    void finished();
    void error(QString err);

private:
    Reader* _reader;
    Dataset* _data;

    vector<int> totalOccurences;
    vector<string> names;

    bool isRunning;
    int* co_occurences;

    bool isParentThread;
    bool isAssociationCalculated;

    vector< pair<int, int> > attributesToProcess;

    int threadsFinished;
    int numThreads;
    int max_numOccurences;
    vector<QThread*> threads;
    vector<CooccurrenceCalculator*> cooccurenceThreads;

    vector< pair<int,int> > _nominalAttrs;
    vector<float> associations;
};

#endif // COOCCURRENCECALCULATOR_H
