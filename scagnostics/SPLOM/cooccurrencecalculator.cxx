#include "cooccurrencecalculator.h"

CooccurrenceCalculator::CooccurrenceCalculator(QObject *parent) : QObject(parent){
  isRunning = false;
  co_occurences = NULL;
  isAssociationCalculated = false;

  max_numOccurences = 0;
}

CooccurrenceCalculator::~CooccurrenceCalculator(){

    //Clean();
    if(co_occurences != NULL)
        free(co_occurences);
}


int CooccurrenceCalculator::GetNumberOfCategoricalDims(){

    int n = _reader->GetTotalAttributes(); //Total attributes
    int total = 0;
    for(int i = 0; i < n; i++){
        bool attr1 = _reader->IsCategorical(i) || _reader->IsOrdinal(i);
        if ( attr1){
            total++;
        }
    }

    return total;
}

void CooccurrenceCalculator::GetAssociationsInfo(vector< pair<int,int> >* nominalPairs, vector<float> * associations){

     for(unsigned int i = 0; i < _nominalAttrs.size(); i++)
         nominalPairs->push_back( _nominalAttrs.at(i));

     for(unsigned int i = 0; i < this->associations.size(); i++){
         associations->push_back( this->associations.at(i) );
     }
}


int CooccurrenceCalculator::GetTotalSumCategoricalSize(int sizes[]){
    int n = _reader->GetTotalAttributes(); //Total attributes
    int total = 0;
    int pos = 0;
    for(int i = 0; i < n; i++){
        bool attr1 = _reader->IsCategorical(i) || _reader->IsOrdinal(i);
        if ( attr1){

            total +=  _reader->GetCategoricalSize(i);
            sizes[pos] = _reader->GetCategoricalSize(i);
            pos++;
        }
    }
    return total;
}

string CooccurrenceCalculator::GetNameOfAttr(int pos, bool direct){

    if ( direct){

        return _reader->GetAttributeName(pos);
    }
    int n = _reader->GetTotalAttributes(); //Total attributes

    int index = 0;
    for(int i = 0; i < n; i++){
        bool attr1 = _reader->IsCategorical(i) || _reader->IsOrdinal(i);
        if ( attr1){

            if ( index == pos){
                return _reader->GetAttributeName(i);
            }
            index++;
        }
    }

    return "N/A";
}

void CooccurrenceCalculator::GetCategoricalIndices(vector<int>* indices){
  int n = _reader->GetTotalAttributes(); //Total attributes
  for(int i = 0; i < n; i++){
      bool attr1 = _reader->IsCategorical(i) || _reader->IsOrdinal(i);
      if ( attr1){
          indices->push_back(i);
      }
  }
}


void CooccurrenceCalculator::CalculateCramerV(Reader* reader,Statistics* stats){

    int n = reader->GetTotalAttributes(); //Total attributes
    //
    _nominalAttrs.clear();

    vector<int> actualCatAttrs;

    for(int i = 0; i < n; i++){
        bool attr1 = reader->IsCategorical(i) || reader->IsOrdinal(i);
        if ( attr1){

            actualCatAttrs.push_back(i);
            for(int j = i+1; j < n;j++ ){
                if ( reader->IsCategorical(j) || reader->IsOrdinal(j)){
                    _nominalAttrs.push_back( make_pair(i,j));
                }
            }
        }
    }
    std::cout << "Calculating Cramer V " << std::endl;
    if( _nominalAttrs.empty()){ return; }

    totalOccurences.clear();
    names.clear();
    // Let a sample of size n of the simultaneously distributed variables
    // A and B  for i = 1...r, j = 1...k
    // n_ij = number of times the values (A_i, B_j) were observed

    max_numOccurences =0 ;
    int stPos = 0;
    double totalElements = stats->GetNumberElements();

    for(unsigned int i = 0; i< actualCatAttrs.size(); i++){
        int n = reader->GetCategoricalSize(actualCatAttrs.at(i));
        for(int id1 = 0; id1 < n ; id1++){

            //float freq1 = stats->GetFrequency(currentPair.first, id1);
            double n_i = stats->GetTotalOfCategoricalValue(actualCatAttrs.at(i), id1); //freq1*totalElements;
            totalOccurences.push_back(n_i);
            string nameAttr = reader->GetAttributeName(actualCatAttrs.at(i));
            string nameCat = reader->GetCategoryAttributeName(actualCatAttrs.at(i),id1);
            string completeName = nameAttr + " : " + nameCat;
            names.push_back(completeName);
        }
    }

    for(unsigned int i = 0; i < _nominalAttrs.size(); i++){

         pair<int, int> currentPair = _nominalAttrs.at(i);
         // We get the size of the local, co-occurence matrix
         int n =  reader->GetCategoricalSize(currentPair.first);
         int m =  reader->GetCategoricalSize(currentPair.second);
         double X = 0;
         for(int id1 = 0; id1 < n ; id1++){
             double n_i = stats->GetTotalOfCategoricalValue(currentPair.first, id1); //freq1*totalElements;
             if ( n_i > max_numOccurences)
                max_numOccurences = n_i;

             for(int id2 = 0; id2 < m; id2++){
                 double n_ij =  co_occurences[stPos + id1*m + id2];

                  double n_j = stats->GetTotalOfCategoricalValue(currentPair.second, id2);
                 double den = (n_i*n_j)/totalElements;
                 double num = (n_ij - den)*(n_ij -den);
                 X += (num/den);
             }
         }
         double V = sqrt( (X/totalElements)/min(n-1, m-1) );
         associations.push_back(V);
         stPos += n*m;
    }
    isAssociationCalculated = true;
}

float CooccurrenceCalculator::GetCramerV(int attr1, int attr2){

    for(unsigned int i = 0; i < _nominalAttrs.size(); i++){
        if (attr1 == _nominalAttrs.at(i).first && attr2 == _nominalAttrs.at(i).second){
            return associations.at(i);
        }
    }

    return 0.0f;
}


void CooccurrenceCalculator::InitializeAsParent(Reader* reader, Dataset* data, int numThreads)
{
    _reader = reader;
    _data =  data;
    // First we get all the possible pairs of categorical + ordinal data
    // We separete them for each thread
    // and send them to process...
    isRunning = true;
    isParentThread = true;
    this->numThreads = numThreads;

    //
    int n = reader->GetTotalAttributes(); //Total attributes
    //
    vector< pair<int,int> > nominalAttrs;
    for(int i = 0; i < n; i++){
        bool attr1 = reader->IsCategorical(i) || reader->IsOrdinal(i);
        if ( attr1){
            for(int j = i+1; j < n;j++ ){
                if ( reader->IsCategorical(j) || reader->IsOrdinal(j)){
                    nominalAttrs.push_back( make_pair(i,j));
                }
            }
        }
    }
    //Got all the attributes that are categorical

    this->threadsFinished = 0;


    if( nominalAttrs.empty()){
        return;
    }

    // Attributes to be processed by each thread
    vector< pair<int, int> >  attributesPerThread[numThreads];
    for(unsigned int i = 0; i < nominalAttrs.size(); i++){
        int addTo = i % numThreads;
        attributesPerThread[addTo].push_back( nominalAttrs.at(i));
    }
    //
    for(int i = 0; i < numThreads;i++){
           threads.push_back(new QThread);
           cooccurenceThreads.push_back(new CooccurrenceCalculator());

           CooccurrenceCalculator* concurrencyThread = cooccurenceThreads.at(i);
           concurrencyThread->Initialize(reader, data,numThreads);
           concurrencyThread->SetInfo(attributesPerThread[i]);
           concurrencyThread->moveToThread(threads.at(i));
           connect( threads.at(i), SIGNAL(started()), concurrencyThread, SLOT(process()));
           connect( concurrencyThread, SIGNAL(finished()), this, SLOT(FinishedGeneration()));
           // Connect worker finished signal to trigger thread quit, then delete.

           // Lets check this one
           connect(concurrencyThread, SIGNAL(finished()), threads.at(i), SLOT(quit()));

           //TODO- for some reason the thread is being deleted
           //-- connect(latestSplomThread, SIGNAL(finished()), latestSplomThread, SLOT(deleteLater()));
           connect(threads.at(i), SIGNAL(finished()), threads.at(i),SLOT(deleteLater()));
    }

    for(int i = 0; i < numThreads; i++){
        threads.at(i)->start();
    }
}


void CooccurrenceCalculator::FinishedGeneration(){
    threadsFinished++;

    //std::cout << "New thread finished " << std::endl;
    if (threadsFinished == numThreads){


        // std::cout << "Loading everything now in memory " << std::endl;
        // Create a new array that can contain all the other memory
        int n = _reader->GetTotalAttributes(); //Total attributes
        //
        int totalSize = 0;
        vector< pair<int,int> > nominalAttrs;
        for(int i = 0; i < n; i++){
            bool attr1 = _reader->IsCategorical(i) || _reader->IsOrdinal(i);
            if ( attr1){
                for(int j = i+1; j < n;j++ ){
                    if ( _reader->IsCategorical(j) || _reader->IsOrdinal(j)){
                        nominalAttrs.push_back( make_pair(i,j));
                        int n = _reader->GetCategoricalSize(i);
                        int m = _reader->GetCategoricalSize(j);
                        totalSize += n*m;
                    }
                }
            }
        }

        //std::cout << "Total size A" << std::endl;

        if (co_occurences != NULL )
            free(co_occurences);


        co_occurences = (int*) malloc(sizeof(int)*totalSize);

        int stPos = 0;
        int stPosInThread[numThreads];
        for(int i = 0; i < numThreads; i++) stPosInThread[i] = 0;
        for(unsigned int i = 0; i < nominalAttrs.size(); i++){
            int extractFrom = i % numThreads;
            CooccurrenceCalculator* cthread = cooccurenceThreads.at(extractFrom);

            pair<int, int> currentPair = nominalAttrs.at(i);
            // We get the size of the local, co-occurence matrix
            int n =  _reader->GetCategoricalSize(currentPair.first);
            int m =  _reader->GetCategoricalSize(currentPair.second);

            int tmpData[n*m];

            cthread->GetInfo(tmpData, stPosInThread[extractFrom], stPosInThread[extractFrom] + n*m);

            for(int k = 0; k < n*m;k++){
                co_occurences[stPos + k] = tmpData[k];
            }

            stPos += n*m;
            stPosInThread[extractFrom] += n*m;
        }

        emit allThreadsDone();
        CleanMemory();

        // and then clean memory
        /*
        std::cout << "Call to clean memory from " << isParentThread << std::endl;
        std::cout <<"Finish cleaning " << std::endl;*/
    }


}

void CooccurrenceCalculator::Initialize(Reader* reader, Dataset* data, int numThreads){
    _reader = reader;
    _data =  data;
    isParentThread = false;
    this->numThreads = 0;
}

void CooccurrenceCalculator::SetInfo(vector<pair<int, int> > attributes){
    // Here we initialize the memory of which attributes im generating
    // the co-occurence ...

    attributesToProcess.clear();
    int totalSize = 0;
    for(unsigned int i = 0; i < attributes.size(); i++){
        int n =  _reader->GetCategoricalSize(attributes.at(i).first);
        int m =  _reader->GetCategoricalSize(attributes.at(i).second);

        attributesToProcess.push_back(attributes.at(i));
       totalSize  += n*m;
    }
    // The size of this is nominal size of attribute 1 x nominal size of attribute 2
    // summing all over the range

    co_occurences = (int*) malloc(sizeof(int)*totalSize*totalSize);
    // set memory to zeros
    memset(co_occurences, 0, sizeof(int)*totalSize*totalSize);
}

void CooccurrenceCalculator::GetInfo(int data[], int from, int to){
    int pos = 0;
    for(int i = from; i < to; i++){
         data[pos] = static_cast<int>(co_occurences[i]);
         pos++;
     }
}

void CooccurrenceCalculator::Clean(){

    if(co_occurences != NULL)
        free(co_occurences);
    co_occurences = NULL;
}

void CooccurrenceCalculator::CleanMemory(){

    if (cooccurenceThreads.empty()) return;

    int n = cooccurenceThreads.size();
    for(int i = 0; i < n; i++){
        CooccurrenceCalculator* sthread = cooccurenceThreads.at(i);
        if (sthread != NULL)
            delete sthread;
    }
    cooccurenceThreads.clear();
    threads.clear();

   // Clean();
}

void CooccurrenceCalculator::process(){
    //
    int stPos = 0;
    for(unsigned int i = 0; i < attributesToProcess.size(); i++){
        //For each pair of attributes ....

        pair<int, int> currentPair = attributesToProcess.at(i);
        // We get the size of the local, co-occurence matrix
        int n =  _reader->GetCategoricalSize(currentPair.first);
        int m =  _reader->GetCategoricalSize(currentPair.second);

        Attribute attr1 = _reader->GetAttribute(currentPair.first);
        Attribute attr2 = _reader->GetAttribute(currentPair.second);

        for(int j = 0; j < n*m ; j++){
            co_occurences[stPos + j] = 0;
        }

        int maxCor = 0;
        // for each element in the dataset
        for(int j = 0; j < this->_data->GetTotalNumberOfElements(); j++){
              int v1 = _data->GetElementValue(j, currentPair.first);
              int v2 = _data->GetElementValue(j, currentPair.second);

              int id1 = attr1.GetIndexValue(v1);
              int id2 = attr2.GetIndexValue(v2);

              if ( id1 == Attribute::MISSING_INDEX || id2 == Attribute::MISSING_INDEX){
                  continue;
              }
              co_occurences[stPos + id1*m + id2] += 1;
              //co_occurences[stPos + id2*m + id1] += 1;

              if (co_occurences[stPos + id1*n + id2] > maxCor )
                  maxCor = co_occurences[stPos + id1*n + id2];
        }
        stPos += n*m;
    }


    emit finished();
}
