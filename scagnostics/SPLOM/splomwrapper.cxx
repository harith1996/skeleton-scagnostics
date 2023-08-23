#include "splomwrapper.h"

SPLOMWrapper::SPLOMWrapper(QObject *parent) : QObject(parent)
{
  isRunning = false;
  sizeT = 58;
  // 64 is the minimum size that we should use
  // it crashes otherwise
  _typeBlend = 0;
  imageData = NULL;
  densityData = NULL;
  scale = 1.0;
  _reader = NULL;
  useAsDrawer = false;

}

void SPLOMWrapper::Initialize(Reader* reader, Dataset* data, Statistics* stats, int numThreads, bool useAsDrawer){
  _reader = reader;
  _data =  data;
  _stats = stats;
  _numThreads = numThreads;
  isRunning = false;
  this->useAsDrawer = useAsDrawer;
  CreateNumericalAttributes();
}

void SPLOMWrapper::ChangePointSize(float scale){
   this->scale = scale;
}

SPLOMWrapper::~SPLOMWrapper()
{
  if (imageData != NULL)
      free(imageData);

  if ( densityData != NULL)
      free(densityData);
}

void SPLOMWrapper::ProcessData(int attributeToCluster){

    if (isRunning) {
        return;
    }
    _attributeToCluster = attributeToCluster;
    isRunning = true;

    vector< pair<int, int> >  attributesPerThread[_numThreads];

    for(unsigned int i = 0; i < attributesToDraw.size(); i++){
        int addTo = i % _numThreads;
        attributesPerThread[addTo].push_back(attributesToDraw.at(i));
    }
    // separated, almost equally between threads...
    // need to take them into account...

    threadsFinished = 0;

    //std::cout << "Calling to process data with banded? " << tf->IsBandingColorMap() << std::endl;
    //std::cout << "Total Nodes " << tf->GetTotalNodes() << std::endl


    for(int j = 0; j < tf->GetTotalNodes(); j++){
        float tc[3];
        tf->GetColor(j, tc);
        //std::cout << "Color of value " << j << " is " << tc[0] << ", " << tc[1] << " ," << tc[2] << std::endl;
    }



    // std::cout << "num threads? " << _numThreads << std::endl;
    //......
    for(int i = 0; i < _numThreads; i++){
        threads.push_back(new QThread);
        splomThreads.push_back(new SPLOMThread());

        SPLOMThread* latestSplomThread = splomThreads.at(i);
        latestSplomThread->SetTransferFunction(tf);

        latestSplomThread->Initialize(_reader,_data, _stats);
        latestSplomThread->SetUseAsDrawer(this->useAsDrawer);
        latestSplomThread->ChangeTypeBlend(_typeBlend);
        latestSplomThread->SetSize(sizeT);
        latestSplomThread->ChangePointSize(scale*0.95);
        latestSplomThread->SetInfo(attributesPerThread[i], attributeToCluster);
        latestSplomThread->moveToThread(threads.at(i));
        connect( threads.at(i), SIGNAL(started()), latestSplomThread, SLOT(process()));
        connect(latestSplomThread,SIGNAL(finished()), this, SLOT(FinishedGeneration()));
        // Connect worker finished signal to trigger thread quit, then delete.
        connect(latestSplomThread, SIGNAL(finished()), threads.at(i), SLOT(quit()));

        //TODO- for some reason the thread is being deleted
        //-- connect(latestSplomThread, SIGNAL(finished()), latestSplomThread, SLOT(deleteLater()));
        connect(threads.at(i), SIGNAL(finished()), threads.at(i),SLOT(deleteLater()));
    }

    for(int i = 0; i < _numThreads;i++){
        threads.at(i)->start();
    }
}

void SPLOMWrapper::ProcessData(int attributeToCluster, vector<pair<int, int> > filteredAttributes){


    std::cout << "Filtered attribute process data " << std::endl;
    _attributeToCluster = attributeToCluster;
     isRunning = true;

    int n = attributesToDraw.size();
    vector< pair<int, int> >  attributesPerThread[_numThreads];

    for(int i = 0; i < n; i++){
        int addTo = i % _numThreads;
        // if the attribute pair to draw is in filtered attributes then we add it
        // otherwise we add  -1,-1 ....
        bool exists = false;
        for(unsigned int k = 0; k < filteredAttributes.size(); k++){
            int v1 = attributesToDraw.at(i).first;
            int v2 = filteredAttributes.at(k).first;

             if (v1 == v2){

                   if ( attributesToDraw.at(i).second == filteredAttributes.at(k).second ){
                       exists = true;
                       break;
                   }
             }
        }

        if ( exists)
            attributesPerThread[addTo].push_back(attributesToDraw.at(i));
         else
            attributesPerThread[addTo].push_back( make_pair<int, int>(-1,-1));

    }
    // separated, almost equally between threads...
    // need to take them into account...

    threadsFinished = 0;

    for(int i = 0; i < _numThreads; i++){
        threads.push_back(new QThread);
        splomThreads.push_back(new SPLOMThread());

        SPLOMThread* latestSplomThread = splomThreads.at(i);
        latestSplomThread->SetTransferFunction(tf);
        latestSplomThread->Initialize(_reader,_data, _stats);
        latestSplomThread->SetUseAsDrawer(this->useAsDrawer);
        latestSplomThread->ChangeTypeBlend(_typeBlend);
        latestSplomThread->SetSize(sizeT);
        latestSplomThread->ChangePointSize(scale*0.95);
        latestSplomThread->SetInfo(attributesPerThread[i], attributeToCluster);
        latestSplomThread->moveToThread(threads.at(i));
        connect( threads.at(i), SIGNAL(started()), latestSplomThread, SLOT(process()));
        connect(latestSplomThread,SIGNAL(finished()), this, SLOT(FinishedGeneration()));
        // Connect worker finished signal to trigger thread quit, then delete.
        connect(latestSplomThread, SIGNAL(finished()), threads.at(i), SLOT(quit()));

        //TODO- for some reason the thread is being deleted
        //-- connect(latestSplomThread, SIGNAL(finished()), latestSplomThread, SLOT(deleteLater()));
        connect(threads.at(i), SIGNAL(finished()), threads.at(i),SLOT(deleteLater()));
        threads.at(i)->start();
    }
}

void SPLOMWrapper::GetImage(int attr1, int attr2, int clusteringAttribute, int valuetoFilter, short texture[]){
    //Lets find the start of the image
    int catSize = 0;

    if (_attributeToCluster != -1)
        catSize = _reader->GetCategoricalSize(clusteringAttribute);

    if (_attributeToCluster != -1 && _reader->IsNumerical(clusteringAttribute))
        catSize = tf->GetNumberOfClases();

    int sizeOfSingle = sizeT*sizeT*(3 + catSize);

    int stPos = 0;
    int checkIndex = 0;
    while (true){

       pair<int, int> currentPair = attributesToDraw.at(checkIndex);

       if( currentPair.first == attr1 && currentPair.second == attr2) break;
       if( currentPair.first == attr2 && currentPair.second == attr1) break;
       stPos += sizeOfSingle;
       checkIndex++;
    }

    // found the starting position at stPos..
    // texture as a rgb array...
    if( clusteringAttribute == -1 || valuetoFilter == -1){ // it is the original image that we are searching for
        for(int i = 0; i < sizeT*sizeT*3;i++)
            texture[i] = imageData[stPos + i];
    }
    else if (clusteringAttribute != -1 && valuetoFilter != -1 && _typeBlend == 4 ){
       // We are actually talking about density here ...

       short tmpDensity[sizeT*sizeT];
       GetDensityImage(attr1, attr2, clusteringAttribute, valuetoFilter, tmpDensity);
       int maxPoints = 0;
       for(int i = 0; i < sizeT*sizeT;i++){ if (tmpDensity[i] > maxPoints) maxPoints = tmpDensity[i]; }

       for(int i = 0; i < sizeT*sizeT;i++){
           float weight = static_cast<float>( tmpDensity[i]) / static_cast<float>(maxPoints);
           float tmpColor[3];
           tf->GetColor(weight, tmpColor);
           texture[i*3 + 0] = tmpColor[0]*255;
           texture[i*3 + 1] = tmpColor[1]*255;
           texture[i*3 + 2] = tmpColor[2]*255;
       }

    }
    else {// it is one of the filtered....
         stPos += sizeT*sizeT*3;
         stPos += sizeT*sizeT*valuetoFilter; //
         float color[3] = {0,0,0};

         if (_reader->IsNumerical(clusteringAttribute)){

             tf->GetClassColor(valuetoFilter, color);
         }
         else
             tf->GetColor(valuetoFilter, color);


         for(int i = 0; i < sizeT*sizeT;i++){
             texture[i*3 + 0] = imageData[stPos +i]*color[0];
             texture[i*3 + 1] = imageData[stPos +i]*color[1];
             texture[i*3 + 2] = imageData[stPos +i]*color[2];
         }
    }

}

void SPLOMWrapper::GetDensityImage(int attr1, int attr2, int clusteringAttribute, int valuetoFilter, short texture[]){
    //Lets find the start of the image
    int catSize = 0;

    if (_attributeToCluster != -1)
        catSize = _reader->GetCategoricalSize(clusteringAttribute);
    if (_attributeToCluster != -1 && _reader->IsNumerical(clusteringAttribute))
        catSize = tf->GetNumberOfClases();


    int sizeOfSingle = sizeT*sizeT*(1 + catSize);

    int stPos = 0;
    int checkIndex = 0;
    while (true){

       pair<int, int> currentPair = attributesToDraw.at(checkIndex);

       if( currentPair.first == attr1 && currentPair.second == attr2) break;
       if( currentPair.first == attr2 && currentPair.second == attr1) break;
       stPos += sizeOfSingle;
       checkIndex++;
    }

    // found the starting position at stPos..
    // texture as a rgb array...
    if( clusteringAttribute == -1 || valuetoFilter == -1){ // it is the original image that we are searching for
        for(int i = 0; i < sizeT*sizeT;i++)
            texture[i] = densityData[stPos + i];
    }
    else { // it is one of the filtered....
         stPos += sizeT*sizeT*1;
         stPos += sizeT*sizeT*valuetoFilter; //

         for(int i = 0; i < sizeT*sizeT;i++){
             texture[i + 0] = densityData[stPos +i];
         }
    }

}

void SPLOMWrapper::CreateNumericalAttributes(){
     numTypeAttributes.clear();
     GetListOfNumericalAttributes(&numTypeAttributes);
     attributesToDraw.clear();
     for(unsigned int i = 0; i < numTypeAttributes.size(); i++)
         for(unsigned int j = i+1; j < numTypeAttributes.size(); j++){
             unsigned int attr1 = numTypeAttributes.at(i);
             unsigned int attr2 = numTypeAttributes.at(j);

             attributesToDraw.push_back( std::make_pair<int, int>(attr1, attr2));
         }

}

int SPLOMWrapper::GetTotalNumericalAttributes(){


    int qty = 0;
    //
    for(int i = 0; i < _reader->GetTotalAttributes(); i++)
        if ( _reader->IsNumerical(i)) qty++;

    return qty;
}

void SPLOMWrapper::GetListOfNumericalAttributes(vector<int>* list){
    list->clear();

    for(int i = 0; i < _reader->GetTotalAttributes(); i++)
        if ( _reader->IsNumerical(i))
            list->push_back(i);

}

void SPLOMWrapper::FinishedGeneration(){

   //clear them  vector<SPLOMThread*> splomThreads;
    threadsFinished++;
    if (threadsFinished == _numThreads){

        // all the threads finish, we need to extract the info
        //and place it so it can be easily queried...
        int catSize = 0;
        if (_attributeToCluster != -1)
            catSize = _reader->GetCategoricalSize(_attributeToCluster);

        if (_attributeToCluster != -1 && _reader->IsNumerical(_attributeToCluster))
            catSize = tf->GetNumberOfClases();


        int sizeOfSingle = sizeT*sizeT*(3 + catSize);
        int sizeOfSingleDensity = sizeT*sizeT*(1+ catSize);

        long int totalSize = sizeT*sizeT*(3 + catSize)*attributesToDraw.size();

        if (imageData != NULL)
            free(imageData);

        if ( densityData != NULL)
            free(densityData);

        imageData = (short*)malloc(totalSize*sizeof(short));
        densityData = (short*)malloc(sizeOfSingleDensity*attributesToDraw.size()*sizeof(short));

        int attributesPerThread[_numThreads];
        for(int i = 0; i < _numThreads; i++){
            attributesPerThread[i] = 0;
        }

        int stPos = 0;
        int stDPos = 0;

        short tmpData[sizeOfSingle];
        short tmpDData[sizeOfSingleDensity];

        for(unsigned int i = 0; i < attributesToDraw.size(); i++){
            int addTo = i % _numThreads;
            //addTo which, and how many have I seen from it ...
            SPLOMThread* cthread = splomThreads.at(addTo);
            cthread->GetInfo(tmpData, attributesPerThread[addTo]*sizeOfSingle,   (attributesPerThread[addTo] + 1)*sizeOfSingle);
            cthread->GetDensityInfo(tmpDData,attributesPerThread[addTo]*sizeOfSingleDensity, (attributesPerThread[addTo]+1)*sizeOfSingleDensity   );
            attributesPerThread[addTo] += 1;

            for(int k = 0; k < sizeOfSingle;k++){
                imageData[stPos + k] = tmpData[k];
            }
            for(int k = 0; k < sizeOfSingleDensity;k++){
                densityData[stDPos + k] = tmpDData[k];
            }

            stPos += sizeOfSingle;
            stDPos += sizeOfSingleDensity;

        }

        isRunning = false;
        emit FinishedImageGeneration();
    }
}



void SPLOMWrapper::CleanMemory(){

    // all the threads finished, we clear the splom vector,
    // and then we clear the threads vector

    if (splomThreads.empty()) return;

    for(int i = 0; i < _numThreads; i++){
        SPLOMThread* sthread = splomThreads.at(i);
        delete sthread;
    }
    splomThreads.clear();
    threads.clear();

   // if (imageData != NULL)
   //     free(imageData);
}
