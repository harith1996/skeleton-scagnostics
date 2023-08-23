#include "statistics.h"

Statistics::Statistics()
{
   currentSimilarityMeasure = Statistics::Lin;
   useBestOrdering = false;
   useCustomOrdering = false;
   customOrdering = NULL;
   bestOrdering = NULL;
   SimilarityMatrix = NULL;
   maxNumericalDistance = -1;
   covarianceMatrix = NULL;
   equalSpaceRepresentation = false; //TODO default false
   nonBestOrdering = true;
   dataNotSet = true;
}

void Statistics::SetToEqualSpace(){
    equalSpaceRepresentation = true;
    for(int j = 0; j < _mapper->GetTotalAttributes(); j++){
        int currentSize = _mapper->GetCategoricalSize(j);

        int order = GetBestOrdering(j);

        for(int i = 0; i < currentSize; i++){

            float mid = GetMidPointPDF(j, i, order, true);

        }

    }
}
void Statistics::ToggleRepresentationToEqualSpaced(){


    equalSpaceRepresentation = !equalSpaceRepresentation;


    for(int j = 0; j < _mapper->GetTotalAttributes(); j++){
        int currentSize = _mapper->GetCategoricalSize(j);

        int order = GetBestOrdering(j);

        for(int i = 0; i < currentSize; i++){

            float mid = GetMidPointPDF(j, i, order, true);

        }

    }
    /*int low = 0, up = 0;

    low = _mapper->GetMinCategoricalSize();   up = _mapper->GetMaxCategoricalSize();

    equalSpaceRepresentation = !equalSpaceRepresentation;
    possibleOrderings.clear();
    ranges.clear();

   for(int i = low; i <= up ; i++ ){
       Ordering newOrder;
       newOrder.SetNumberOfElements(i);
       possibleOrderings.push_back(newOrder);
   }
    // For each non ordinal variable, calculate which is the best ordering...

    CalculateBestOrdering(true);
    int attrs = _mapper->GetTotalAttributes();

    for(int i = 0; i < attrs; i++){

        AttributeRanges range;

       if (_mapper->IsNumerical(i)){
           //
            range = GetNumericalVariableRange(i);
       }
       else if (_mapper->IsOrdinal(i) || _mapper->IsCategorical(i)){
            range = GetCategoricalVariableRange(i);
       }

       ranges.push_back(range);
    }
    dataset->SetMaxsAndMins(&ranges);

    CreateCovarianceMatrix();*/
}

Statistics::~Statistics(){


    for (std::map<const int, float*>::iterator it=attributeSimilarity.begin(); it!=attributeSimilarity.end(); ++it){
           free(it->second);
    }
    attributeSimilarity.clear();
}

float Statistics::GetFrequency(int attributeIndex, int valueIndex, bool force){

    /*Return the frequency according to a single attribute, and how
      many times it appears*/
    if ( valueIndex == Attribute::MISSING_INDEX) return 0;
    if (_mapper->IsNumerical(attributeIndex) ) return 0;


   //pair<int,int> key =  make_pair<int, int>(attributeIndex, valueIndex);
   pair<int,int> key =  make_pair(attributeIndex, valueIndex);

    if ( frequencies.count(key) != 0 && !force){

        return frequencies[key];
    }
    else {

        int totalNonMissing = 0;
        int totalOfValue = 0;
        //

        Attribute att = _mapper->GetAttribute(attributeIndex);

        for(int i = 0 ; i < dataset->GetTotalNumberOfElements(); i++){

            if ( dataset->IsElementValueMissing(i, attributeIndex))
                continue;

            int value = static_cast<int>( dataset->GetElementValue(i, attributeIndex));

            int vIndex =  att.GetIndexValue(value); //_mapper->GetIndexOfValue(attributeIndex, value);

            if ( vIndex == valueIndex)
                totalOfValue++;

            // We only take into account those elements that are not missing
            // to calculate the frequency of the element
            totalNonMissing++;
        }

        float frequency = static_cast<float>(totalOfValue)/static_cast<float>(totalNonMissing);

        if ( equalSpaceRepresentation)
            frequency = 1.0/att.GetNumCatDimensions();
        frequencies[key] = frequency;
        categoricalValueTotals[key] = totalOfValue;
        return frequencies[key];
    }
}

int Statistics::GetTotalOfCategoricalValueSkippable(int attributeIndex, int valueIndex,bool skip[]){
     Attribute att = _mapper->GetAttribute(attributeIndex);

     int totalOfValue = 0;
     for(int i = 0 ; i < dataset->GetTotalNumberOfElements(); i++){

         if ( skip[i]) continue;

         int value = static_cast<int>( dataset->GetElementValue(i, attributeIndex));

         int vIndex =  att.GetIndexValue(value); //_mapper->GetIndexOfValue(attributeIndex, value);

         if ( vIndex == valueIndex)
             totalOfValue++;
     }

     return totalOfValue;
}

int Statistics::GetTotalOfTwoCategoricalValueSkippable(int attributeIndex1, int valueIndex1,int attributeIndex2, int valueIndex2, bool skip[]){
     Attribute att1 = _mapper->GetAttribute(attributeIndex1);
     Attribute att2 = _mapper->GetAttribute(attributeIndex2);

     int totalOfValue = 0;
     for(int i = 0 ; i < dataset->GetTotalNumberOfElements(); i++){

         if ( skip[i]) continue;

         int value1 = static_cast<int>( dataset->GetElementValue(i, attributeIndex1));
         int value2 = static_cast<int>( dataset->GetElementValue(i, attributeIndex2));



         int vIndex1 =  att1.GetIndexValue(value1); //_mapper->GetIndexOfValue(attributeIndex, value);
         int vIndex2 =  att2.GetIndexValue(value2); //_mapper->GetIndexOfValue(attributeIndex, value);

         if ( vIndex1 == valueIndex1 && vIndex2 == valueIndex2)
             totalOfValue++;
     }
     return totalOfValue;
}

int Statistics::GetTotalOfCategoricalValue(int attributeIndex, int valueIndex){
    pair<int,int> key =  make_pair(attributeIndex, valueIndex);
    return  categoricalValueTotals[key];
}


void Statistics::SetCustomOrdering(int attributeIndex, int newOrdering){
    customOrdering[attributeIndex] = newOrdering;
}


int Statistics::GetCustomOrdering(int attributeIndex){
 return customOrdering[attributeIndex];
}

void Statistics::GetTextInfo(vector<QString> *info, int elementIndex){
  // First for each attribute we generate an array of strings
  int n = _mapper->GetTotalAttributes();
  //QString values[n];

  for(int i = 0; i < n; i++){
      if ( _mapper->IsNumerical(i)){
          // if its numerical, we get just the value.
          info->push_back( QString::number(  dataset->GetElementValue(elementIndex, i)  ) );
      }
      else if (_mapper->IsCategorical(i) || _mapper->IsOrdinal(i)){
          // if its ordinal or categorical we get the name
          int val = static_cast<int>(GetElementValue(elementIndex,i));
          // value is not the same as the index...
          Attribute attr = _mapper->GetAttribute(i);
          int idx = attr.GetIndexValue(val);

          if (idx == Attribute::MISSING_INDEX ){
              info->push_back(QString::fromStdString("Missing"));
          }
          else {
              std::string name = _mapper->GetCategoryAttributeName(i, idx);
              info->push_back(QString::fromStdString(name));
          }
      }
      else if (_mapper->IsDescriptive(i) || _mapper->IsDate(i)){
          info->push_back( "_tmp_" );
      }
  }

}

float Statistics::GetMaximumInAttribute(int attributeIndex){
    return ranges.at(attributeIndex).maximum;//  maximumInAttribute.at(attributeIndex);
}
float Statistics::GetMinimumInAttribute(int attributeIndex){
    return ranges.at(attributeIndex).minimum; //minimumInAttribute.at(attributeIndex);
}

float Statistics::GetMeanInAttribute(int attributeIndex){

   return ranges.at(attributeIndex).mean;
}

void Statistics::GetHeader(vector<QString>* header){
   int n = _mapper->GetTotalAttributes();
   for(int i = 0; i < n; i++)
       header->push_back(QString::fromStdString(_mapper->GetAttributeName(i)) );
}

void Statistics::GetOrder(int orderIndex, int attributeIndex, vector<int> *order){

   int currentSize = _mapper->GetCategoricalSize(attributeIndex);
   if (possibleOrderings.empty())
       CreateOrderings();
   Ordering orders = possibleOrderings.at(currentSize -_mapper->GetMinCategoricalSize());
   int desiredOrder = 0;
   if ( _mapper->IsCategorical(attributeIndex)){
          desiredOrder = orderIndex;
     }

    orders.GetOrder(desiredOrder, order);
}

float  Statistics::GetPDF(int attributeIndex, int valueIndex, int order, bool force, float* freq){
    if ( valueIndex < 0) return 0;

    // This is actually the important part...
    // on reordering the frequencies...
    // There is a total of getcategoricalsize(valueIndex)!
    // different ordering possibilities...


    QString key = QString::number(attributeIndex) + "," + QString::number(valueIndex) + "," + QString::number(order) + ".";

    std::string ckey = key.toStdString();

    if ( pdfWithOrdering.count(ckey) != 0 && !force){

        if (freq != nullptr){
            float f = GetFrequency(attributeIndex, valueIndex, force);
            *freq = f;
        }

        return pdfWithOrdering[ckey];
    }
    else {

        vector<int> currentOrder;
        GetOrder(order, attributeIndex,&currentOrder);


        float sum = 0;
        int i =0;
        while (true){

            float f = GetFrequency(attributeIndex, currentOrder.at(i), force);
            sum += f;
            if ( currentOrder.at(i) == valueIndex){
                if (freq != nullptr) *freq = f;
                break;
            }
            i++;
        }
        pdfWithOrdering[ckey] = sum;
        return sum;
    }
}


int Statistics::GetPrevValueInOrder(int attributeIndex, int valueIndex, int order){
    vector<int> currentOrder;
    GetOrder(order, attributeIndex, &currentOrder);

    int prev = -1;
    for(int i = 0; i < static_cast<int>(currentOrder.size())-1;i++){
        if ( currentOrder.at(i+1) == valueIndex){
            prev = currentOrder.at(i);
        }
    }
    return prev;
}

int Statistics::GetValueInOrder(int attributeIndex, int valueIndex, int order){
    vector<int> currentOrder;
    GetOrder(order, attributeIndex, &currentOrder);

    int next = -1;
    for(int i = 0; i < static_cast<int>(currentOrder.size());i++){
        if ( currentOrder.at(i) == valueIndex){
            next = currentOrder.at(i);
        }
    }
    return next;
}

int Statistics::GetNextValueInOrder(int attributeIndex, int valueIndex, int order){
    vector<int> currentOrder;
    GetOrder(order, attributeIndex, &currentOrder);

    int next = -1;
    for(int i = 1; i < static_cast<int>(currentOrder.size());i++){
        if ( currentOrder.at(i-1) == valueIndex){
            next = currentOrder.at(i);
        }
    }
    return next;
}


void Statistics::ClearMaps(){

    for (std::map<const int, float*>::iterator it=attributeSimilarity.begin(); it!=attributeSimilarity.end(); ++it){
           free(it->second);
    }
    attributeSimilarity.clear();

    categoricalValueTotals.clear();
    frequencies.clear();
    pdfWithOrdering.clear();
    midpointpdfWithOrdering.clear();

}
float Statistics::GetMidPointPDF(int attributeIndex, int valueIndex, int order, bool force){

    if (_mapper->IsNumerical(attributeIndex) ) return 0;

    QString key = QString::number(attributeIndex) + "," + QString::number(valueIndex) + "," + QString::number(order) + ".";
    std::string ckey = key.toStdString();


    if (midpointpdfWithOrdering.count(ckey) != 0 && !force){

        return midpointpdfWithOrdering[ckey];
    }
    else {

        int prev = GetPrevValueInOrder(attributeIndex, valueIndex, order);

        float v1 = GetPDF(attributeIndex, prev, order, force);
        float v2 = GetPDF(attributeIndex, valueIndex, order, force);

        midpointpdfWithOrdering[ckey] = (v1 + v2)/2.0;
        return midpointpdfWithOrdering[ckey];
    }

}

int Statistics::GetDesiredOrder(int desiredOrder[], int attributeIndex){
    int currentSize = _mapper->GetCategoricalSize(attributeIndex);
    Ordering orders = possibleOrderings.at(currentSize -_mapper->GetMinCategoricalSize());

    return orders.FindOrder(desiredOrder);
}

int Statistics::GetValueIndex(int attributeIndex, float frequency, int order  ){

    int currentSize = _mapper->GetCategoricalSize(attributeIndex);


    for(int i = 0; i < currentSize; i++){
        float mid = GetMidPointPDF(attributeIndex, i, order);

        float freq = GetFrequency(attributeIndex, i) / 2.0;
        float dif = fabs(mid - frequency);

        if ( dif < freq)
            return i;
    }
    return -1;
}

void Statistics::ForceRecalculationOrdering(int attributeIndex){

  if (_mapper->IsOrdinal(attributeIndex)){
      // if the attribute is ordinal then we still need to recalculate
      // the pdf

      int currentSize = _mapper->GetCategoricalSize(attributeIndex);
      for(int i = 0; i < currentSize; i++){
          GetMidPointPDF(attributeIndex, i, 0, true);
      }
     // the ordinal is done...
  }
  if (_mapper->IsCategorical(attributeIndex))
  {
      // if is categorical, we need to recalculate the order and the pdf
      int currentSize = _mapper->GetCategoricalSize(attributeIndex);
      Ordering orders = possibleOrderings.at(currentSize - _mapper->GetMinCategoricalSize());
      int total = orders.GetTotalNumberOfOrderings();

      int indexOfMax = 0;
      float maxDist = 0;
      for(int j = 0; j < total; j++){
         vector<int> tmp;
         // We need the distance
         orders.GetOrder(j, &tmp );
         float sum = 0;
         for(int k = 0; k < currentSize -1; k++){
             float v1 = GetMidPointPDF(attributeIndex, tmp.at(k), j, true);
             float v2 = GetMidPointPDF(attributeIndex, tmp.at(k+1), j, true);
             sum += fabs(v1-v2);
         }
         if ( sum > maxDist){
             maxDist = sum;
             indexOfMax = j;
         }
      }
      bestOrdering[attributeIndex] = indexOfMax;
  }
}


void Statistics::CalculateBestOrdering(bool force = false){

    if ( customOrdering != NULL)
        free(customOrdering);

    if ( bestOrdering != NULL)
        free(bestOrdering);

    customOrdering = (int*) malloc(sizeof(int)*_mapper->GetTotalAttributes());
    bestOrdering = (int*) malloc(sizeof(int)*_mapper->GetTotalAttributes());
    for(int category = 0; category < _mapper->GetTotalAttributes(); category++){
        // For each category.

        customOrdering[category] = 0; // initialize the custom ordering to 0

        if (!_mapper->IsCategorical(category)){
            bestOrdering[category] = 0;
        }
        else {
            // if its not ordinal, we have to check all possible ordering...
            int currentSize = _mapper->GetCategoricalSize(category);
            Ordering orders = possibleOrderings.at(currentSize - _mapper->GetMinCategoricalSize());
            int total = orders.GetTotalNumberOfOrderings();

            int indexOfMax = 0;
            float maxDist = 0;

            if (!nonBestOrdering){
                for(int j = 0; j < total; j++){
                   vector<int> tmp;
                   // We need the distance
                   orders.GetOrder(j, &tmp );
                   float sum = 0;
                   for(int k = 0; k < currentSize -1; k++){
                       float v1 = GetMidPointPDF(category, tmp.at(k), j, force);
                       float v2 = GetMidPointPDF(category, tmp.at(k+1), j, force);
                       sum += fabs(v1-v2);
                   }
                   if ( sum > maxDist){
                       maxDist = sum;
                       indexOfMax = j;
                   }
                }
            }

            bestOrdering[category] = indexOfMax;
        }
    }
}

int Statistics::GetBestOrdering(int attributeIndex){

    if (!useBestOrdering) return 0;
    else {
        return bestOrdering[attributeIndex];
    }
}


AttributeRanges Statistics::GetNumericalVariableRange(int index){
    int i =  index;
    AttributeRanges range;
    int n = dataset->GetTotalNumberOfElements();

    float val = dataset->GetElementValue(0, i);
    float mean = 0;
    int nonEmpty = 0;
    float variance = 0;
    float max = val, min = val;

    for(int j = 1; j < n; j++){
         val = dataset->GetElementValue(j, i);


         if ( dataset->IsElementValueMissing(j,i)) continue; // Missing data element ...
         if ( max == -99999 ) max = val;
         if ( min == -99999 ) min = val;

         if ( val > max) max = val;
         if ( val < min) min = val;
         mean += val;
         nonEmpty += 1;
    }

    if (nonEmpty > 1){
        mean /= nonEmpty;

        variance = 0;
        for(int j = 0; j < n; j++){
             val = dataset->GetElementValue(j, i);
             if ( val == -99999) continue; // Missing data element ...
             variance += (val- mean)*(val - mean);
        }

        variance = variance/(nonEmpty-1);
    }
    else {
        variance = -1;// It cannot be negative, so basically and error-check
    }

    range.maximum = max;
    range.minimum = min;
    range.mean = mean;
    range.variance = variance;
    float stdDev = sqrt(variance);
    range.stdDev = stdDev;

    float minStandard = (val - mean)/stdDev;
    float maxStandard = minStandard;

    for(int j = 0; j < n; j++){
         val = dataset->GetElementValue(j, i);
         if ( val == -99999) continue; // Missing data element ...

         float newVal = (val - mean)/stdDev;


         if ( newVal > maxStandard) maxStandard = newVal;
         if ( newVal < minStandard) minStandard = newVal;
    }

    range.minStandarized = minStandard;
    range.maxStandarized = maxStandard;
    return range;
}


AttributeRanges Statistics::GetNumericalVariableRange(int index, vector<int>* elementIndices, bool skip[]){
    int i =  index;
    AttributeRanges range;
    int n = elementIndices->size();

    float val = dataset->GetElementValue(0, elementIndices->at(0));
    float mean = 0;
    int nonEmpty = 0;
    float variance = 0;
    float max = val, min = val;

    for(unsigned int j = 1; j < elementIndices->size(); j++){
        if ( skip[ elementIndices->at(j)]) continue;
         val = dataset->GetElementValue( elementIndices->at(j), i);

         if ( dataset->IsElementValueMissing(elementIndices->at(j),i)) continue; // Missing data element ...

         if ( max == -99999 ) max = val;
         if ( min == -99999 ) min = val;

         if ( val > max) max = val;
         if ( val < min) min = val;
         mean += val;
         nonEmpty += 1;
    }

    if (nonEmpty > 1){
        mean /= nonEmpty;

        variance = 0;
        for(int j = 0; j < n; j++){
            if ( skip[ elementIndices->at(j)]) continue;

             val = dataset->GetElementValue(elementIndices->at(j), i);
             if ( val == -99999) continue; // Missing data element ...
             variance += (val- mean)*(val - mean);
        }

        variance = variance/(nonEmpty-1);
    }
    else {
        variance = -1;// It cannot be negative, so basically and error-check
    }

    range.maximum = max;
    range.minimum = min;
    range.mean = mean;
    range.variance = variance;
    float stdDev = sqrt(variance);
    range.stdDev = stdDev;

    float minStandard = (val - mean)/stdDev;
    float maxStandard = minStandard;

    for(int j = 0; j < n; j++){
        if ( skip[ elementIndices->at(j)]) continue;

         val = dataset->GetElementValue(elementIndices->at(j), i);
         if ( val == -99999) continue; // Missing data element ...

         float newVal = (val - mean)/stdDev;


         if ( newVal > maxStandard) maxStandard = newVal;
         if ( newVal < minStandard) minStandard = newVal;
    }

    range.minStandarized = minStandard;
    range.maxStandarized = maxStandard;
    return range;
}

AttributeRanges Statistics::GetCategoricalVariableRange(int index, vector<int>* elementIndices, bool skip[]){
    AttributeRanges range;
    int i = index;
    int order = GetBestOrdering(i);

    float maxFreq = 0;
    float maxMidPoint = 0, minMidPoint = 1.0;

    int sz = _mapper->GetCategoricalSize(i);
    int idxOfMax = 0;
    for(int j = 0; j < sz; j++){
        float frequency = 1.0;
         GetPDF(i, j, order, false, &frequency );
         float midPoint = GetMidPointPDF(i, j, order);
         if (midPoint > maxMidPoint)
             maxMidPoint = midPoint;
         if (midPoint < minMidPoint)
             minMidPoint = midPoint;

         if ( frequency > maxFreq){
             maxFreq = frequency;
             idxOfMax = j;
         }
    }




    range.mean = GetMidPointPDF(i, idxOfMax, order);
    range.maximum = maxMidPoint;
    range.minimum = minMidPoint;
    // variance
    // std-dev
    // minStandard
    // maxStandard
    int n = elementIndices->size();

    Attribute attr = GetAttr(index);

    float variance = 0;
    int totPoints = 0;
    for(int j = 0; j < n ; j++){
        int val = GetElementValue(elementIndices->at(j), index);

        if ( val == -99999) continue; // Missing data element ...

        if (skip[elementIndices->at(j)]) continue;


        totPoints += 1;
        int idx = attr.GetIndexValue(val);
        float midPoint = GetMidPointPDF(i, idx, order);
        variance += (midPoint - range.mean)*(midPoint - range.mean);
    }

    if ( totPoints < 1)
        variance = 0;
    else
       variance /= (totPoints-1);

    range.variance = variance;
    float stdDev = sqrt(variance);
    range.stdDev = stdDev;

    range.maxStandarized = (maxMidPoint - range.mean)/stdDev;
    range.minStandarized = (minMidPoint - range.mean)/stdDev;


   return range;
}


AttributeRanges Statistics::GetCategoricalVariableRange(int index){
    AttributeRanges range;
    int i = index;
    int order = GetBestOrdering(i);

    float maxFreq = 0;
    float maxMidPoint = 0, minMidPoint = 1.0;

    int sz = _mapper->GetCategoricalSize(i);
    int idxOfMax = 0;
    for(int j = 0; j < sz; j++){
        float frequency = 1.0;
         GetPDF(i, j, order, false, &frequency );
         float midPoint = GetMidPointPDF(i, j, order);
         if (midPoint > maxMidPoint)
             maxMidPoint = midPoint;
         if (midPoint < minMidPoint)
             minMidPoint = midPoint;

         if ( frequency > maxFreq){
             maxFreq = frequency;
             idxOfMax = j;
         }
    }


    range.mean = GetMidPointPDF(i, idxOfMax, order);
    range.maximum = maxMidPoint;
    range.minimum = minMidPoint;
    // variance
    // std-dev
    // minStandard
    // maxStandard
    int n = dataset->GetTotalNumberOfElements();

    Attribute attr = GetAttr(index);

    float variance = 0;
    int totPoints = 0;
    for(int j = 0; j < n ; j++){
        int val = GetElementValue(j, index);
        if ( val == -99999) continue; // Missing data element ...

        totPoints += 1;
        int idx = attr.GetIndexValue(val);
        float midPoint = GetMidPointPDF(i, idx, order);
        variance += (midPoint - range.mean)*(midPoint - range.mean);
    }

    variance /= (totPoints-1);

    range.variance = variance;
    float stdDev = sqrt(variance);
    range.stdDev = stdDev;

    range.maxStandarized = (maxMidPoint - range.mean)/stdDev;
    range.minStandarized = (minMidPoint - range.mean)/stdDev;


   return range;

}


void Statistics::CreateOrderings(){


    int low = 0, up = 0;
    low = _mapper->GetMinCategoricalSize();   up = _mapper->GetMaxCategoricalSize();


    std::cout << "Creating ordering? " << std::endl;
     possibleOrderings.clear();
    // ranges.clear();
     //maximumInAttribute.clear();
     //minimumInAttribute.clear();

    for(int i = low; i <= up ; i++ ){
        std::cout << "i " <<i << " -- limit " << up << std::endl;
        Ordering newOrder;
        newOrder.SetNumberOfElements(i);
        possibleOrderings.push_back(newOrder);
    }
     // For each non ordinal variable, calculate which is the best ordering...
     CalculateBestOrdering(false);
}

void Statistics::CreateRanges(){

    int attrs = _mapper->GetTotalAttributes();

    for(int i = 0; i < attrs; i++){

        //std::cout << "Creating range for " << i  << std::endl;
        AttributeRanges range;

       if (_mapper->IsNumerical(i)){
           //
              range = GetNumericalVariableRange(i);
              //std::cout << "Numerical ? " << std::endl;
       }
       else if (_mapper->IsOrdinal(i) || _mapper->IsCategorical(i)){
              range = GetCategoricalVariableRange(i);
              //std::cout << "Cate" << std::endl;
       }
       else {
           //std::cout << "Variable  " << i << " has no range? " << std::endl;
       }

       ranges.push_back(range);
       //std::cout << "Now " << ranges.size() << std::endl;
    }


    dataset->SetMaxsAndMins(&ranges);

    //Not really needed now... the covariance...
   // CreateCovarianceMatrix();
}

void Statistics::SetDataset(Dataset *data,  Reader *mapper)
{
    dataset = data;
    _mapper = mapper;
     dataNotSet = false;
}


void Statistics::CreateWeightedCovarianceMatrix(int attributeToJump, vector<int>* indices, float weights[], gsl_matrix* covarianceMatrix, bool skip[]){


    int attrs = _mapper->GetTotalAttributes();
    vector<int> localAttrs;

    for(int i = 0; i < attrs; i++){
        if (i == attributeToJump) continue;
        if ( weights[i] < 0.0001 ) continue;

        localAttrs.push_back(i);
    }

    // Steps:
    // 1. Compute the weighted  mean of the data
    // 2. Center the data by the weighted mean
    // 3. Compute the weighted covariance

    int d = localAttrs.size();
    float mean[d];
    bool numerical[d];
    Attribute locals[d];
    for(int i = 0; i < d; i++) {
        mean[i] = 0;
        numerical[i] = _mapper->IsNumerical(localAttrs.at(i));
        locals[i] = _mapper->GetAttribute(localAttrs.at(i));
    }
    int total = 0;

    for(unsigned int k = 0; k < indices->size(); k++){

        int element =  indices->at(k);
        if ( skip[element]) continue;
        // given that we have the element... now we add to the mean
        for(int h = 0; h < d; h++){

            if ( numerical[h]){
                mean[h] += dataset->GetElementValue(element, localAttrs.at(h) );
            }
            else {
                int val = dataset->GetElementValue(element, localAttrs.at(h));
                if ( val == -99999) continue; // Missing data element ...
                int idx =  locals[h].GetIndexValue(val);
                mean[h] += GetMidPointPDF( localAttrs.at(h), idx, 0);
            }


        }
        total++;
    }

    for( int i = 0; i < d; i++) mean[i] /= total;


    //2. Compute the mean-mean matrices

   gsl_matrix* W    = gsl_matrix_alloc(d, d);
   gsl_matrix* C    = gsl_matrix_alloc(d, d);
   gsl_matrix* Tmp    = gsl_matrix_alloc(1, d);


   float sum = 0;

   for(int k = 0; k < d; k++){
       for(int h = 0; h < d; h++){
            gsl_matrix_set(W, k, h, 0 );
            gsl_matrix_set(C, k, h, 0 );
       }
       gsl_matrix_set(W, k, k, weights[localAttrs.at(k)] );
       sum += weights[localAttrs.at(k)];
   }

   gsl_matrix* tmp1 = gsl_matrix_alloc(1, d);
   gsl_matrix* tmp2 = gsl_matrix_alloc(d, 1);


   for(int i= 0; i <  indices->size(); i++){
       int element = indices->at(i);
       if ( skip[element]) continue;

       for(int h = 0; h < d; h++){
          float val = 0;

          if ( numerical[h]){
              val = dataset->GetElementValue(element, localAttrs.at(h) );
          }
          else {
              int val = dataset->GetElementValue(element, localAttrs.at(h));
              if ( val == -99999) continue; // Missing data element ...
              int idx =  locals[h].GetIndexValue(val);
              val = GetMidPointPDF( localAttrs.at(h), idx, 0);
          }

          gsl_matrix_set(tmp1, 0, h, (val - mean[h]));
          gsl_matrix_set(tmp2, h, 0, (val  - mean[h]));
       }

       gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, tmp1, W, 0.0, Tmp);
       gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, tmp2, Tmp, 1.0, C);
   }

   PrintMatrix("Tmp Covariance Matrix ", C);

   sum = 1.0 / sum;

   for(int k = 0; k < covarianceMatrix->size1; k++){
       for(int h = 0; h < covarianceMatrix->size2; h++){
             gsl_matrix_set(covarianceMatrix, h, k, 0);
       }
   }

   for(int k = 0; k < d; k++){
       for(int h = 0; h < d; h++){
           gsl_matrix_set(covarianceMatrix, localAttrs.at(h), localAttrs.at(k), sum* gsl_matrix_get(C, h,k));
       }
   }

   PrintMatrix("Weighted Covariance Matrix ", covarianceMatrix);

   gsl_matrix_free(W);
   gsl_matrix_free(C);
   gsl_matrix_free(Tmp);
   gsl_matrix_free(tmp1);
   gsl_matrix_free(tmp2);

}


void  Statistics::CreateCovarianceMatrix(int attributeToJump, vector<int>* indices, gsl_matrix* covarianceMatrix, bool skip[]){
    int attrs = _mapper->GetTotalAttributes();

    vector<AttributeRanges> ranges;
    vector<int> localAttrs;
    for(int i = 0; i < attrs; i++){

        if (i == attributeToJump) continue;
        localAttrs.push_back(i);
        AttributeRanges range;

        if (_mapper->IsNumerical(i)){
           range = GetNumericalVariableRange(i,indices, skip);
        }
        else if (_mapper->IsOrdinal(i) || _mapper->IsCategorical(i)){
           range = GetCategoricalVariableRange(i, indices, skip);
        }
        ranges.push_back(range);
    }

    int n = localAttrs.size();
    int numElements = indices->size();

    for(int attr1 = 0; attr1 < n; attr1++){

        gsl_matrix_set(covarianceMatrix, attr1, attr1, ranges.at(attr1).variance);

        bool isFirstNumerical = _mapper->IsNumerical(localAttrs.at(attr1));

        Attribute attrib1 = GetAttr( localAttrs.at(attr1) );
        int order1 = GetBestOrdering(localAttrs.at(attr1));

        for(int attr2 = attr1+1; attr2 < n ; attr2++){

              Attribute attrib2 = GetAttr(localAttrs.at(attr2));
              int order2 = GetBestOrdering(localAttrs.at(attr2));
              bool isSecondNumerical = _mapper->IsNumerical(localAttrs.at(attr2));

              float covariance = 0;
              int totPoints = 0;

              for(int j = 0; j < numElements; j++){

                      if (skip[j] ) continue;
                      float v1,v2;

                      if ( isFirstNumerical){
                          v1 = dataset->GetElementValue(indices->at(j), localAttrs.at(attr1));
                          if ( v1 == -99999) continue; // Missing data element ...
                      }
                      else {
                          int val = GetElementValue(indices->at(j), localAttrs.at(attr1));
                          if ( val == -99999) continue; // Missing data element ...
                          int idx =  attrib1.GetIndexValue(val);
                           v1 = GetMidPointPDF( localAttrs.at(attr1), idx, order1);
                      }

                      if ( isSecondNumerical){
                          v2 = dataset->GetElementValue(indices->at(j), localAttrs.at(attr2));
                          if ( v2 == -99999) continue; // Missing data element ...
                      }
                      else {
                          int val = dataset->GetElementValue(indices->at(j), localAttrs.at(attr2));
                          if ( val == -99999) continue; // Missing data element ...
                          int idx =  attrib2.GetIndexValue(val);
                          v2 = GetMidPointPDF(localAttrs.at(attr2), idx, order2);
                      }

                      totPoints += 1;
                      covariance += (v1- ranges.at(attr1).mean )*(v2 - ranges.at(attr2).mean);
              }

              if (totPoints == 0)
                  covariance = 0;
              else
                 covariance /= (totPoints -1);

              gsl_matrix_set(covarianceMatrix,attr1, attr2, covariance);
              gsl_matrix_set(covarianceMatrix,attr2, attr1, covariance);
        }
    }
    PrintMatrix("Covariance Matrix ", covarianceMatrix);
}


float Statistics::GetSimilarityFromMem(int i, int j){

    if (SimilarityMatrix == NULL)
        return 0;

    int n = dataset->GetTotalNumberOfElements();
    return SimilarityMatrix[i*n + j];
}

float Statistics::GetSimilarity(int i, int j, int attributeToIgnore ){
    float sum = 0;
    int cAttrs = 0;

    int numAttrs = dataset->GetTotalNumberOfAttributes();
    for(int k = 0; k < numAttrs; k++){
        if ( attributeToIgnore == k) continue;

        if (_mapper->IsOrdinal(k) || _mapper->IsCategorical(k)){
          sum += GetCategoricalSimilarityFromMemory(i, j, k);
          cAttrs += 1;
        }
    }


    sum += GetNumericalSimilarity(i, j);
    sum /= (cAttrs +1);
    return sum;
}


float Statistics::GetDistance(vector<float> dataPoint, int j){
    float sum = 0;
    int cAttrs = 0;

    int numAttrs = dataset->GetTotalNumberOfAttributes();
    for(int k = 0; k < numAttrs; k++){
        if (_mapper->IsOrdinal(k) || _mapper->IsCategorical(k)){
          sum += (1.0 - GetCategoricalSimilarity(dataPoint, j, k));
          cAttrs += 1;
        }
    }

    //TODO- remove assumption that there's numerical attributes
    sum += GetNumericalDistance(dataPoint, j);
    sum /= (cAttrs +1);
    return sum;

}

float Statistics::GetNumericalDistance(vector<float>& dataPoint, int j){

    int numAttrs = dataset->GetTotalNumberOfAttributes();
    float distanceV = 0;
    int totalNumerical = 0;
    for(int k = 0; k < numAttrs; k++){
        if (_mapper->IsNumerical(k)){
            float maxV  =  GetMaximumInAttribute(k);
            float minV  =  GetMinimumInAttribute(k);
            float rangeV = maxV - minV;

            float v1 = (dataPoint.at(k) - minV) /(rangeV) ;
            float v2 = (dataset->GetElementValue(j, k) - minV)/(rangeV);

            float d = (v1 -v2)*(v1-v2);
            if ( rangeV < 0.00000001)
                d = 0;
            distanceV += d;
            totalNumerical++;
        }
    }

    if ( totalNumerical == 0)
        return 0;

    return sqrt(distanceV);
}

float Statistics::GetNumericalSimilarity(int i, int j){

    float val = 0;
    if ( attributeSimilarity.count(-1) != 0){
        // if its one then the key exists
        float* numericalSimilarityMatrix = attributeSimilarity[-1];
        val = numericalSimilarityMatrix[i*dataset->GetTotalNumberOfElements() + j];
    }
    return val;
}

float Statistics::GetCategoricalSimilarityFromMemory(int i, int j, int attributeIndex){
    float val = 0;
    if ( attributeSimilarity.count(-1) != 0){
        // if its one then the key exists
        float* SimilarityMatrix = attributeSimilarity[attributeIndex];
        val = SimilarityMatrix[i*dataset->GetTotalNumberOfElements() + j];
    }
    return val;
}

void Statistics::CreateCategoricalSimilarity(){

    int numAttrs = dataset->GetTotalNumberOfAttributes();
    for(int i = 0; i < numAttrs; i++){

        if (_mapper->IsOrdinal(i) || _mapper->IsCategorical(i)){
            CreateCategoricalSimilarityForAttribute(i);
        }
    }
}


void Statistics::CreateCPUSimilarityMatrix(int attributeToIgnore ){

    if ( SimilarityMatrix != NULL)
        free(SimilarityMatrix);
    int n = dataset->GetTotalNumberOfElements();

    SimilarityMatrix = (float*)malloc(sizeof(float)*n*n);
    for(int i = 0; i < n; i++){
        for(int j = i; j < n; j++){
            float similarity = GetSimilarity(i,j, attributeToIgnore);

            SimilarityMatrix[i*n +j] = similarity;
            SimilarityMatrix[j*n +i] = similarity;
        }
    }
}

void Statistics::CreateSimilarityMatrix(){
    if ( SimilarityMatrix != NULL)
        free(SimilarityMatrix);

    int n = dataset->GetTotalNumberOfElements();


    int cAttrs = 0; // This assumes numerical attributes....

    int numAttrs = dataset->GetTotalNumberOfAttributes();
    for(int k = 0; k < numAttrs; k++){
        if (_mapper->IsOrdinal(k) || _mapper->IsCategorical(k)){
          cAttrs += 1;
        }
    }

    SimilarityMatrix = (float*)malloc(sizeof(float)*n*n);

     if ( attributeSimilarity.count(-1) > 0){
         cAttrs += 1;
     }



    //Allocate the memory for all the similarity matrices together...
    float* tmpData = (float*) malloc(sizeof(float)*cAttrs*n*n);

    int pos = 0;
    for(int k = 0; k < numAttrs; k++){
        if (_mapper->IsOrdinal(k) || _mapper->IsCategorical(k)){
              float* mem = attributeSimilarity[k];
              memcpy(  &(tmpData[pos]), mem, sizeof(float)*n*n);
              pos += n*n;
        }
    }

    if ( attributeSimilarity.count(-1) > 0){
        float* mem = attributeSimilarity[-1];
        memcpy(  &(tmpData[pos]), mem, sizeof(float)*n*n);
    }

    CombineSimilarity(tmpData, cAttrs, n, &SimilarityMatrix);

    std::cout << "Similarity matrix created " << std::endl;
    free(tmpData);
}

void Statistics::CreateNumericalSimilarity(){


    if ( attributeSimilarity.count(-1) > 0){
        // if the key already exists
        float* inf = attributeSimilarity[-1];
        free(inf);
        attributeSimilarity.erase(-1);
    }

   int numAttrs = dataset->GetTotalNumberOfAttributes();
   bool isNumerical[numAttrs];
   int numNumericalAttrs = 0;

   for(int i = 0; i < numAttrs; i++){
        isNumerical[i] = _mapper->IsNumerical(i);
        if (isNumerical[i]) numNumericalAttrs++;
   }

   if (numNumericalAttrs == 0) return; // No need to return if there are no numerical values

   int n = dataset->GetTotalNumberOfElements();

   float* numericalSimilarityMatrix = (float*)malloc(sizeof(float)*n*n);
   CreateCudaNumericalSimilarity( dataset->GetTotalNumberOfAttributes(), n, isNumerical,
                                  dataset->GetNormalizedDataPointer(), &numericalSimilarityMatrix);
   //find the maximum and re-scale the matrix
   float maxDist = 0;
   for(int i = 0; i < n; i++){
       for(int j = i; j < n; j++){
           if ( numericalSimilarityMatrix[i*n + j] > maxDist)
               maxDist = numericalSimilarityMatrix[i*n + j];
       }
   }

   maxNumericalDistance = maxDist;
   for(int i = 0; i < n; i++){
       for(int j = i; j < n; j++){

           float v =  numericalSimilarityMatrix[i*n + j] / maxDist;
           numericalSimilarityMatrix[i*n + j] = 1.0 - v;
       }
   }

   // Currently it is the dissimilarity matrix... so this is the distance
   attributeSimilarity[-1] = numericalSimilarityMatrix;
}

void Statistics::CreateCategoricalSimilarityForAttribute(const int attributeIndex){
    if ( attributeSimilarity.count(attributeIndex) > 0){
        // if the key already exists
        float* inf = attributeSimilarity[attributeIndex];
        free(inf);
        attributeSimilarity.erase(attributeIndex);
    }

    bool nominal =  _mapper->IsCategorical(attributeIndex) || _mapper->IsOrdinal(attributeIndex);

    if ( !nominal ) return;


    int numAttrs = dataset->GetTotalNumberOfAttributes();
    int n = dataset->GetTotalNumberOfElements();
    int catSize = _mapper->GetCategoricalSize(attributeIndex);

    float freq[catSize];
    for(int i = 0; i < catSize; i++){
        freq[i] = GetFrequency(attributeIndex, i);
    }
    float* catSimilarityMatrix = (float*)malloc(sizeof(float)*n*n);
    CreateCudaCategoricalSimilarity(numAttrs, n, attributeIndex, catSize, freq, dataset->GetDataPointer(),
                                    &catSimilarityMatrix, currentSimilarityMeasure);

    attributeSimilarity[attributeIndex] = catSimilarityMatrix;
}


void Statistics::CreateLongProjectionMatrix(float* prevProjectionMatrix, float *newProjectionMatrix, int clusteringAttribute,
                                            vector<pair<float, float> > numericalClassesRanges){
    //1. Compute the d-dimensional mean vectors for the different classes from the dataset.
    // The prev-projection matrix is used to basically remove the attributes that are not needed...

    // First let's get the number of classes...
    int numClasses = 0;
    bool numerical = false;

    if ( _mapper->IsNumerical(clusteringAttribute)){
        numClasses = numericalClassesRanges.size();
        numerical = true;
    }
    else
        numClasses = _mapper->GetCategoricalSize(clusteringAttribute);

    // by default, whatever clustering attribute we are using will have a weight of 0 ...
    int non_zero_features = 0;
    int n = _mapper->GetTotalAttributes();
    vector<int> nonZeroAttrs;
    for(int i = 0; i < n ; i++){
     if ( prevProjectionMatrix[i*3+ 2]  > 0.00000001) {
         non_zero_features++;
         nonZeroAttrs.push_back(i);
     }
   }

    // non_zero_features is the dimensionality that we have now

    int d = non_zero_features;
    vector<int> groupings[numClasses]; // group the indices of elements that belong to a class...

    float means[numClasses][d];
    float mean[d];
    int totalSize = 0;
    for(int i = 0; i < d; i++) mean[i] = 0;


    Attribute attr = _mapper->GetAttribute(clusteringAttribute);

    for(int i = 0; i < numClasses; i++){

        vector<int> indices;
        if (numerical){
            dataset->FilterDatasetByNumericalValue( clusteringAttribute, numericalClassesRanges.at(i).first, numericalClassesRanges.at(i).second*1.01, &indices);
        }
        else {
            std::string val = attr.GetDims().at(i).GetValue();
            QString tmp = QString::fromStdString(val);
            int value = tmp.toInt();
            dataset->FilterDatasetByCategoricalValue(clusteringAttribute, value, &indices);
        }
        // Create the mean of it...
        for(int k = 0; k < d; k++)
            means[i][k] = 0; // initialize to 0..0

        for(unsigned int k = 0; k < indices.size(); k++){

            int element =  indices.at(k);
            // given that we have the element... now we add to the mean
            for(int h = 0; h < d; h++){
                means[i][h] += dataset->GetElementValue(element, nonZeroAttrs.at(h) );
                mean[h] += dataset->GetElementValue(element, nonZeroAttrs.at(h) );
            }
            groupings[i].push_back(indices.at(k));
        }

        for(int k = 0; k < d; k++){
            means[i][k] /= groupings[i].size();
        }
        totalSize += groupings[i].size();
    }

    for(int i = 0; i < d; i++) mean[i] /= totalSize;

     //2. Compute the mean-mean matrices
    gsl_matrix* tmp1 = gsl_matrix_alloc(1, d);
    gsl_matrix* tmp2 = gsl_matrix_alloc(d, 1);
    gsl_matrix* withinClassScatter = gsl_matrix_alloc(d, d);

    gsl_matrix* C    = gsl_matrix_alloc(d, d);

    for(int i= 0; i < numClasses; i++){

          for(int k = 0; k < d; k++)
              for(int h = 0; h < d; h++)
                   gsl_matrix_set(C, k, h, 0 );

          for(int j = 0; j < numClasses; j++){

              for(int h = 0; h < d; h++){
                  gsl_matrix_set(tmp1, 0, h, means[j][h]  - means[i][h]);
                  gsl_matrix_set(tmp2, h, 0, means[j][h]  - means[i][h]);
              }

              gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, tmp2, tmp1, 1.0, C);
          }

          gsl_matrix_add(withinClassScatter, C);
         //
    }

    // Compute Eigen-Vectors of it
    gsl_vector *eval = gsl_vector_alloc (d);
    gsl_matrix *evec = gsl_matrix_alloc (d, d);

    gsl_eigen_symmv_workspace * w =
      gsl_eigen_symmv_alloc (d);

    gsl_eigen_symmv (C, eval, evec, w);

    gsl_eigen_symmv_free (w);

    gsl_eigen_symmv_sort (eval, evec,
                          GSL_EIGEN_SORT_VAL_DESC);



    gsl_vector_view evec_first = gsl_matrix_column (evec, 0);
    gsl_vector_view evec_second = gsl_matrix_column (evec, 1);

    int nonZero = 0;

    double first_eigenVector[non_zero_features];
    double second_eigenVector[non_zero_features];

    for(int k = 0; k < n; k++){
       // for each feature
       if ( prevProjectionMatrix[k*3+ 2]  > 0.00000001){

           first_eigenVector[nonZero] = gsl_vector_get(&evec_first.vector,nonZero );
           second_eigenVector[nonZero] = gsl_vector_get(&evec_second.vector,nonZero);
           nonZero++;
       }
    }

    nonZero = 0;

    for(int k = 0; k < n; k++){
      // for each feature
         if ( prevProjectionMatrix[k*3+ 2]  > 0.00000001){
             // non-zero
             double v1[3] = {first_eigenVector[nonZero], second_eigenVector[nonZero],0};
             double w = MathHelper::Norm(v1,3);
             MathHelper::Normalize(v1);

             newProjectionMatrix[k*3+ 0] = v1[0];
             newProjectionMatrix[k*3+ 1] = v1[1];
             newProjectionMatrix[k*3+ 2] = w;

             nonZero++;
         }
         else {
             newProjectionMatrix[k*3+ 0] = 0;
             newProjectionMatrix[k*3+ 1] = 0;
             newProjectionMatrix[k*3+ 2] = 0;
         }

    }

    gsl_matrix_free(tmp1);
    gsl_matrix_free(tmp2);
    gsl_matrix_free(withinClassScatter);
    gsl_matrix_free(C);
}

void Statistics::CreateLDAProjectionMatrix(float* prevProjectionMatrix, float *newProjectionMatrix, int clusteringAttribute, float regularizationParameter,
                                           vector<pair<float, float> > numericalClassesRanges){
   //1. Compute the d-dimensional mean vectors for the different classes from the dataset.
   // The prev-projection matrix is used to basically remove the attributes that are not needed...

   // First let's get the number of classes...
   int numClasses = 0;
   bool numerical = false;

   if ( _mapper->IsNumerical(clusteringAttribute)){
       numClasses = numericalClassesRanges.size();
       numerical = true;
   }
   else
       numClasses = _mapper->GetCategoricalSize(clusteringAttribute);

   // by default, whatever clustering attribute we are using will have a weight of 0 ...
   int non_zero_features = 0;
   int n = _mapper->GetTotalAttributes();
   vector<int> nonZeroAttrs;
   for(int i = 0; i < n ; i++){
    if ( prevProjectionMatrix[i*3+ 2]  > 0.00000001) {
        non_zero_features++;
        nonZeroAttrs.push_back(i);
    }
  }

   // non_zero_features is the dimensionality that we have now

   int d = non_zero_features;
   vector<int> groupings[numClasses]; // group the indices of elements that belong to a class...

   float means[numClasses][d];
   float mean[d];
   int totalSize = 0;
   for(int i = 0; i < d; i++) mean[i] = 0;


   Attribute attr = _mapper->GetAttribute(clusteringAttribute);

   for(int i = 0; i < numClasses; i++){

       vector<int> indices;
       if (numerical){
           dataset->FilterDatasetByNumericalValue( clusteringAttribute, numericalClassesRanges.at(i).first, numericalClassesRanges.at(i).second*1.01, &indices);
       }
       else {
           std::string val = attr.GetDims().at(i).GetValue();
           QString tmp = QString::fromStdString(val);
           int value = tmp.toInt();
           dataset->FilterDatasetByCategoricalValue(clusteringAttribute, value, &indices);
       }
       // Create the mean of it...
       for(int k = 0; k < d; k++)
           means[i][k] = 0; // initialize to 0..0

       for(unsigned int k = 0; k < indices.size(); k++){

           int element =  indices.at(k);
           // given that we have the element... now we add to the mean
           for(int h = 0; h < d; h++){
               means[i][h] += dataset->GetElementValue(element, nonZeroAttrs.at(h) );
               mean[h] += dataset->GetElementValue(element, nonZeroAttrs.at(h) );
           }
           groupings[i].push_back(indices.at(k));
       }

       for(int k = 0; k < d; k++){
           means[i][k] /= groupings[i].size();
       }
       totalSize += groupings[i].size();
   }

   for(int i = 0; i < d; i++) mean[i] /= totalSize;


   //2. Compute the scatter matrices (in-between-class and within-class scatter matrix).

   gsl_matrix* withinClassScatter = gsl_matrix_alloc(d, d);
   gsl_matrix* regularizationScatter = gsl_matrix_alloc(d, d);

   gsl_matrix* invWithinClassScatter = gsl_matrix_alloc(d, d);

   gsl_matrix* betweenClassScatter = gsl_matrix_alloc(d, d);

   for(int j = 0; j < d; j++){
       for(int k =0; k < d; k++){
              gsl_matrix_set( withinClassScatter, j, k, 0 );
              gsl_matrix_set( regularizationScatter, j, k, 0 );
              gsl_matrix_set( invWithinClassScatter, j, k, 0 );
              gsl_matrix_set( betweenClassScatter, j, k, 0 );
       }

       gsl_matrix_set( regularizationScatter, j, j, regularizationParameter*1.0);

   }

   gsl_matrix* tmp1 = gsl_matrix_alloc(1, d);
   gsl_matrix* tmp2 = gsl_matrix_alloc(d, 1);

   gsl_matrix* C    = gsl_matrix_alloc(d, d);

   for(int i= 0; i < numClasses; i++){

         vector<int>* elements = &(groupings[i]);

         for(int k = 0; k < d; k++)
             for(int h = 0; h < d; h++)
                  gsl_matrix_set(C, k, h, 0 );

         for(unsigned int elementIdx = 0; elementIdx < elements->size(); elementIdx++){
             int element = elements->at(elementIdx);

             for(int h = 0; h < d; h++){
                 gsl_matrix_set(tmp1, 0, h, dataset->GetElementValue(element, nonZeroAttrs.at(h) )  - means[i][h]);
                 gsl_matrix_set(tmp2, h, 0, dataset->GetElementValue(element, nonZeroAttrs.at(h) )  - means[i][h]);
             }

             gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, tmp2, tmp1, 1.0, C);
         }

         gsl_matrix_add(withinClassScatter, C);
        //
   }

   gsl_matrix_add(withinClassScatter, regularizationScatter);

   PrintMatrix("within-class Scatter-matrix", withinClassScatter);

   for(int i= 0; i < numClasses; i++){

       for(int h = 0; h < d; h++){
           gsl_matrix_set(tmp1, 0, h, means[i][h] - mean[h]);
           gsl_matrix_set(tmp2, h, 0, means[i][h] - mean[h]);
       }

       gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, groupings[i].size(), tmp2, tmp1, 1.0, betweenClassScatter);
   }

   PrintMatrix("between-class Scatter-matrix", betweenClassScatter);
   //3.Compute the eigenvectors (ee1,ee2,...,eed ) and corresponding eigenvalues (λλ1,λλ2,...,λλd ) for the scatter matrices.
   MathHelper::InverseMatrixGSL2(withinClassScatter, invWithinClassScatter, d);

   gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, invWithinClassScatter, betweenClassScatter, 0.0, C);


   gsl_vector *eval = gsl_vector_alloc (d);
   gsl_matrix *evec = gsl_matrix_alloc (d, d);

   gsl_eigen_symmv_workspace * w =
     gsl_eigen_symmv_alloc (d);

   gsl_eigen_symmv (C, eval, evec, w);

   gsl_eigen_symmv_free (w);

   gsl_eigen_symmv_sort (eval, evec,
                         GSL_EIGEN_SORT_VAL_DESC);

   //4. Sort the eigenvectors by decreasing eigenvalues and choose k eigenvectors with the largest eigenvalues to form a d×k dimensional matrix WW (where every column represents an eigenvector).


   gsl_vector_view evec_first = gsl_matrix_column (evec, 0);
   gsl_vector_view evec_second = gsl_matrix_column (evec, 1);

   int nonZero = 0;

   double first_eigenVector[non_zero_features];
   double second_eigenVector[non_zero_features];

   for(int k = 0; k < n; k++){
      // for each feature
      if ( prevProjectionMatrix[k*3+ 2]  > 0.00000001){

          first_eigenVector[nonZero] = gsl_vector_get(&evec_first.vector,nonZero );
          second_eigenVector[nonZero] = gsl_vector_get(&evec_second.vector,nonZero);
          nonZero++;
      }
   }


   //5. Use this d×k eigenvector matrix to transform the samples onto the new subspace. This can be summarized by the matrix multiplication: YY=XX×WW
   // (where XX is a n×d-dimensional matrix representing the n samples, and yy are the transformed n×k-dimensional samples in the new subspace).


   nonZero = 0;

   for(int k = 0; k < n; k++){
     // for each feature

        if ( prevProjectionMatrix[k*3+ 2]  > 0.00000001){
            // non-zero

            double v1[3] = {first_eigenVector[nonZero], second_eigenVector[nonZero],0};
            double w = MathHelper::Norm(v1,3);
            MathHelper::Normalize(v1);


            newProjectionMatrix[k*3+ 0] = v1[0];
            newProjectionMatrix[k*3+ 1] = v1[1];
            newProjectionMatrix[k*3+ 2] = w;

            nonZero++;
        }
        else {
            newProjectionMatrix[k*3+ 0] = 0;
            newProjectionMatrix[k*3+ 1] = 0;
            newProjectionMatrix[k*3+ 2] = 0;
        }

   }



   gsl_matrix_free(tmp1);
   gsl_matrix_free(tmp2);
   gsl_matrix_free(withinClassScatter);
   gsl_matrix_free(betweenClassScatter);
   gsl_matrix_free(regularizationScatter);
   gsl_matrix_free(invWithinClassScatter);
   gsl_matrix_free(evec);
   gsl_vector_free(eval);
   gsl_matrix_free(C);
}


void Statistics::CreateSkippableLDAProjectionMatrix(float* prevProjectionMatrix, float *newProjectionMatrix, int clusteringAttribute, float regularizationParameter,
                                           vector<pair<float, float> > numericalClassesRanges, bool skip[]){
   //1. Compute the d-dimensional mean vectors for the different classes from the dataset.
   // The prev-projection matrix is used to basically remove the attributes that are not needed...

   // First let's get the number of classes...
   int numClasses = 0;
   bool numerical = false;

   if ( _mapper->IsNumerical(clusteringAttribute)){
       numClasses = numericalClassesRanges.size();
       numerical = true;
   }
   else
       numClasses = _mapper->GetCategoricalSize(clusteringAttribute);

   // by default, whatever clustering attribute we are using will have a weight of 0 ...
   int non_zero_features = 0;
   int n = _mapper->GetTotalAttributes();
   vector<int> nonZeroAttrs;
   for(int i = 0; i < n ; i++){
    if ( prevProjectionMatrix[i*3+ 2]  > 0.00000001) {
        non_zero_features++;
        nonZeroAttrs.push_back(i);
    }
  }

   // non_zero_features is the dimensionality that we have now

   int d = non_zero_features;
   vector<int> groupings[numClasses]; // group the indices of elements that belong to a class...

   float means[numClasses][d];
   float mean[d];
   int totalSize = 0;
   for(int i = 0; i < d; i++) mean[i] = 0;


   Attribute attr = _mapper->GetAttribute(clusteringAttribute);

   for(int i = 0; i < numClasses; i++){

       vector<int> indices;
       if (numerical){
           dataset->FilterDatasetByNumericalValue( clusteringAttribute, numericalClassesRanges.at(i).first, numericalClassesRanges.at(i).second*1.01, &indices);
       }
       else {
           std::string val = attr.GetDims().at(i).GetValue();
           QString tmp = QString::fromStdString(val);
           int value = tmp.toInt();
           dataset->FilterDatasetByCategoricalValue(clusteringAttribute, value, &indices);
       }
       // Create the mean of it...
       for(int k = 0; k < d; k++)
           means[i][k] = 0; // initialize to 0..0

       for(unsigned int k = 0; k < indices.size(); k++){

           int element =  indices.at(k);
           if ( skip[element]) continue;

           // given that we have the element... now we add to the mean
           for(int h = 0; h < d; h++){
               means[i][h] += dataset->GetElementValue(element, nonZeroAttrs.at(h) );
               mean[h] += dataset->GetElementValue(element, nonZeroAttrs.at(h) );
           }
           groupings[i].push_back(indices.at(k));
       }

       for(int k = 0; k < d; k++){
           means[i][k] /= groupings[i].size();
       }
       totalSize += groupings[i].size();
   }

   for(int i = 0; i < d; i++) mean[i] /= totalSize;


   //2. Compute the scatter matrices (in-between-class and within-class scatter matrix).

   gsl_matrix* withinClassScatter = gsl_matrix_alloc(d, d);
   gsl_matrix* regularizationScatter = gsl_matrix_alloc(d, d);

   gsl_matrix* invWithinClassScatter = gsl_matrix_alloc(d, d);

   gsl_matrix* betweenClassScatter = gsl_matrix_alloc(d, d);

   for(int j = 0; j < d; j++){
       for(int k =0; k < d; k++){
              gsl_matrix_set( withinClassScatter, j, k, 0 );
              gsl_matrix_set( regularizationScatter, j, k, 0 );
              gsl_matrix_set( invWithinClassScatter, j, k, 0 );
              gsl_matrix_set( betweenClassScatter, j, k, 0 );
       }

       gsl_matrix_set( regularizationScatter, j, j, regularizationParameter*1.0);

   }

   gsl_matrix* tmp1 = gsl_matrix_alloc(1, d);
   gsl_matrix* tmp2 = gsl_matrix_alloc(d, 1);

   gsl_matrix* C    = gsl_matrix_alloc(d, d);

   for(int i= 0; i < numClasses; i++){

         vector<int>* elements = &(groupings[i]);

         for(int k = 0; k < d; k++)
             for(int h = 0; h < d; h++)
                  gsl_matrix_set(C, k, h, 0 );

         for(unsigned int elementIdx = 0; elementIdx < elements->size(); elementIdx++){
             int element = elements->at(elementIdx);

             for(int h = 0; h < d; h++){
                 gsl_matrix_set(tmp1, 0, h, dataset->GetElementValue(element, nonZeroAttrs.at(h) )  - means[i][h]);
                 gsl_matrix_set(tmp2, h, 0, dataset->GetElementValue(element, nonZeroAttrs.at(h) )  - means[i][h]);
             }

             gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, tmp2, tmp1, 1.0, C);
         }

         gsl_matrix_add(withinClassScatter, C);
        //
   }

   gsl_matrix_add(withinClassScatter, regularizationScatter);

   PrintMatrix("within-class Scatter-matrix", withinClassScatter);

   for(int i= 0; i < numClasses; i++){

       for(int h = 0; h < d; h++){
           gsl_matrix_set(tmp1, 0, h, means[i][h] - mean[h]);
           gsl_matrix_set(tmp2, h, 0, means[i][h] - mean[h]);
       }

       gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, groupings[i].size(), tmp2, tmp1, 1.0, betweenClassScatter);
   }

   PrintMatrix("between-class Scatter-matrix", betweenClassScatter);
   //3.Compute the eigenvectors (ee1,ee2,...,eed ) and corresponding eigenvalues (λλ1,λλ2,...,λλd ) for the scatter matrices.
   MathHelper::InverseMatrixGSL2(withinClassScatter, invWithinClassScatter, d);

   gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, invWithinClassScatter, betweenClassScatter, 0.0, C);


   gsl_vector *eval = gsl_vector_alloc (d);
   gsl_matrix *evec = gsl_matrix_alloc (d, d);

   gsl_eigen_symmv_workspace * w =
     gsl_eigen_symmv_alloc (d);

   gsl_eigen_symmv (C, eval, evec, w);

   gsl_eigen_symmv_free (w);

   gsl_eigen_symmv_sort (eval, evec,
                         GSL_EIGEN_SORT_VAL_DESC);

   //4. Sort the eigenvectors by decreasing eigenvalues and choose k eigenvectors with the largest eigenvalues to form a d×k dimensional matrix WW (where every column represents an eigenvector).


   gsl_vector_view evec_first = gsl_matrix_column (evec, 0);
   gsl_vector_view evec_second = gsl_matrix_column (evec, 1);

   int nonZero = 0;

   double first_eigenVector[non_zero_features];
   double second_eigenVector[non_zero_features];

   for(int k = 0; k < n; k++){
      // for each feature
      if ( prevProjectionMatrix[k*3+ 2]  > 0.00000001){

          first_eigenVector[nonZero] = gsl_vector_get(&evec_first.vector,nonZero );
          second_eigenVector[nonZero] = gsl_vector_get(&evec_second.vector,nonZero);
          nonZero++;
      }
   }


   //5. Use this d×k eigenvector matrix to transform the samples onto the new subspace. This can be summarized by the matrix multiplication: YY=XX×WW
   // (where XX is a n×d-dimensional matrix representing the n samples, and yy are the transformed n×k-dimensional samples in the new subspace).


   nonZero = 0;

   for(int k = 0; k < n; k++){
     // for each feature

        if ( prevProjectionMatrix[k*3+ 2]  > 0.00000001){
            // non-zero

            double v1[3] = {first_eigenVector[nonZero], second_eigenVector[nonZero],0};
            double w = MathHelper::Norm(v1,3);
            MathHelper::Normalize(v1);


            newProjectionMatrix[k*3+ 0] = v1[0];
            newProjectionMatrix[k*3+ 1] = v1[1];
            newProjectionMatrix[k*3+ 2] = w;

            nonZero++;
        }
        else {
            newProjectionMatrix[k*3+ 0] = 0;
            newProjectionMatrix[k*3+ 1] = 0;
            newProjectionMatrix[k*3+ 2] = 0;
        }

   }



   gsl_matrix_free(tmp1);
   gsl_matrix_free(tmp2);
   gsl_matrix_free(withinClassScatter);
   gsl_matrix_free(betweenClassScatter);
   gsl_matrix_free(regularizationScatter);
   gsl_matrix_free(invWithinClassScatter);
   gsl_matrix_free(evec);
   gsl_vector_free(eval);
   gsl_matrix_free(C);
}



float Statistics::GetCategoricalSimilarity(vector<float> dataPointA, vector<float> dataPointB, int attributeIndex){

    float similarity = 0;
    int v1, v2;
    int t1, t2;
    int n = dataset->GetTotalNumberOfElements();

    switch( currentSimilarityMeasure){
        case Gower:
            v1 = static_cast<int>(dataPointA.at(attributeIndex)); //dataset->GetElementValue(i, attributeIndex);
            v2 = static_cast<int>(dataPointB.at(attributeIndex));
            similarity = (v1 == v2)? 1:0;
            break;
        case Lin:
            v1 = static_cast<int>(dataPointA.at(attributeIndex));
            v2 = static_cast<int>(dataPointB.at(attributeIndex));
            if ( v1 == v2) {
                similarity = 2* log( GetFrequency(attributeIndex, v1));
            }
            else{
                similarity = 2* log( GetFrequency(attributeIndex, v1) +GetFrequency(attributeIndex, v2)  );
            }
            break;
        case Eskin:
            v1 = static_cast<int>(dataPointA.at(attributeIndex));
            v2 = static_cast<int>(dataPointB.at(attributeIndex));
            if ( v1 == v2) similarity = 1.0;
            else similarity = pow(_mapper->GetCategoricalSize(attributeIndex), 2.0) /( pow(_mapper->GetCategoricalSize(attributeIndex), 2.0) + 2.0);
            break;

        case Own:
            v1 = static_cast<int>(dataPointA.at(attributeIndex));
            v2 = static_cast<int>(dataPointB.at(attributeIndex));
            if ( v1 == v2) similarity = 0;
            else {
                t1 = GetFrequency(attributeIndex, v1)*n;
                t1 = (t1 - 1 >= 0)? t1 -1 : 0;

                t2 = GetFrequency(attributeIndex, v2)*n;
                t2 = (t2 - 1 >= 0)? t2 -1 : 0;

                similarity = sqrt( pow(static_cast<float>(t1)/(n-2),2.0)  +  pow(static_cast<float>(t2)/(n-2),2.0)  );
            }
            similarity = 1.0 - similarity;
            break;

    }
    return similarity;
}

float Statistics::GetCategoricalSimilarity(  vector<float> dataPoint, int j, int attributeIndex){
    // i and j are the elements in the dataset

    float similarity = 0;
    int v1, v2;
    int t1, t2;
    int n = dataset->GetTotalNumberOfElements();

    switch( currentSimilarityMeasure){
        case Gower:
            v1 = static_cast<int>(dataPoint.at(attributeIndex)); //dataset->GetElementValue(i, attributeIndex);
            v2 = dataset->GetElementValue(j, attributeIndex);
            similarity = (v1 == v2)? 1:0;
            break;
        case Lin:
            v1 = static_cast<int>(dataPoint.at(attributeIndex));
            v2 = dataset->GetElementValue(j, attributeIndex);
            if ( v1 == v2) {
                similarity = 2* log( GetFrequency(attributeIndex, v1));
            }
            else{
                similarity = 2* log( GetFrequency(attributeIndex, v1) +GetFrequency(attributeIndex, v2)  );
            }
            break;
        case Eskin:
            v1 = static_cast<int>(dataPoint.at(attributeIndex));
            v2 = dataset->GetElementValue(j, attributeIndex);
            if ( v1 == v2) similarity = 1.0;
            else similarity = pow(_mapper->GetCategoricalSize(attributeIndex), 2.0) /( pow(_mapper->GetCategoricalSize(attributeIndex), 2.0) + 2.0);
            break;

        case Own:
            v1 = static_cast<int>(dataPoint.at(attributeIndex));
            v2 = dataset->GetElementValue(j, attributeIndex);
            if ( v1 == v2) similarity = 0;
            else {
                t1 = GetFrequency(attributeIndex, v1)*n;
                t1 = (t1 - 1 >= 0)? t1 -1 : 0;

                t2 = GetFrequency(attributeIndex, v2)*n;
                t2 = (t2 - 1 >= 0)? t2 -1 : 0;

                similarity = sqrt( pow(static_cast<float>(t1)/(n-2),2.0)  +  pow(static_cast<float>(t2)/(n-2),2.0)  );
            }
            similarity = 1.0 - similarity;
            break;

    }
    return similarity;
}



float Statistics::GetCategoricalSimilarity(int i, int j, int attributeIndex){
    // i and j are the elements in the dataset
    if ( i == j) return 1;

    float similarity = 0;
    int v1, v2;
    int t1, t2;
    int n = dataset->GetTotalNumberOfElements();

    switch( currentSimilarityMeasure){
        case Gower:
            v1 = dataset->GetElementValue(i, attributeIndex);
            v2 = dataset->GetElementValue(j, attributeIndex);
            similarity = (v1 == v2)? 1:0;
            break;
        case Lin:
            v1 = dataset->GetElementValue(i, attributeIndex);
            v2 = dataset->GetElementValue(j, attributeIndex);
            if ( v1 == v2) {
                similarity = 2* log( GetFrequency(attributeIndex, v1));
            }
            else{
                similarity = 2* log( GetFrequency(attributeIndex, v1) +GetFrequency(attributeIndex, v2)  );
            }
            break;
        case Eskin:
            v1 = dataset->GetElementValue(i, attributeIndex);
            v2 = dataset->GetElementValue(j, attributeIndex);
            if ( v1 == v2) similarity = 1.0;
            else similarity = pow(_mapper->GetCategoricalSize(attributeIndex), 2.0) /( pow(_mapper->GetCategoricalSize(attributeIndex), 2.0) + 2.0);
            break;

        case Own:
            v1 = dataset->GetElementValue(i, attributeIndex);
            v2 = dataset->GetElementValue(j, attributeIndex);
            if ( v1 == v2) similarity = 0;
            else {
                t1 = GetFrequency(attributeIndex, v1)*n;
                t1 = (t1 - 1 >= 0)? t1 -1 : 0;

                t2 = GetFrequency(attributeIndex, v2)*n;
                t2 = (t2 - 1 >= 0)? t2 -1 : 0;

                similarity = sqrt( pow(static_cast<float>(t1)/(n-2),2.0)  +  pow(static_cast<float>(t2)/(n-2),2.0)  );
            }
            similarity = 1.0 - similarity;
            break;

    }
    return similarity;
}


void Statistics::CreatePCAProjectionMatrix(float* prevProjectionMatrix, float* newProjectionMatrix, pair<int,int> eigenVectorsToUse, gsl_matrix *covarianceMatrix){
  // Create a PCA Projection Matrix (?) possibly...
  // Define the plane by the two most important eigenvectors
  // Project each of the axis unto the plane (?)
  // Use that as the new projection, but don't use the ones whose weight is 0

  int non_zero_features = 0;
  int n = _mapper->GetTotalAttributes();
  vector<int> nonZeroAttrs;
  for(int i = 0; i < n ; i++){
    std::cout << "Attribute i " << _mapper->GetAttributeName(i) << ":" << i <<" :: " << prevProjectionMatrix[i*3+ 2]  << std::endl;

   if ( prevProjectionMatrix[i*3+ 2]  > 0.00000001) {
       non_zero_features++;
       nonZeroAttrs.push_back(i);
   }
 }

  std::cout << "Non-zero attribs " << non_zero_features << std::endl;
  for(int k = 0; k < non_zero_features; k++)
      std::cout << nonZeroAttrs.at(k) << ",";
  std::cout << std::endl;

    if ( covarianceMatrix == nullptr)
      covarianceMatrix = this->covarianceMatrix;



     // Use the covariance matrix now to generate...
     gsl_matrix* tmpMatrix = gsl_matrix_alloc(non_zero_features, non_zero_features);

     for(unsigned int j = 0; j < nonZeroAttrs.size(); j++){
         for(unsigned int k =0; k < nonZeroAttrs.size(); k++){
                gsl_matrix_set(tmpMatrix, j, k,  gsl_matrix_get(covarianceMatrix,nonZeroAttrs.at(j), nonZeroAttrs.at(k) ) );
         }
     }

   std::cout << "Total non_zero features " << non_zero_features << std::endl;

   PrintMatrix("Covariance Matrix ", tmpMatrix);

   int M = non_zero_features;
   gsl_vector *eval = gsl_vector_alloc (M);
   gsl_matrix *evec = gsl_matrix_alloc (M, M);

   gsl_eigen_symmv_workspace * w =
     gsl_eigen_symmv_alloc (M);

   gsl_eigen_symmv (tmpMatrix, eval, evec, w);

   gsl_eigen_symmv_free (w);

   gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_VAL_DESC);


   // We can have the two biggest eigenvectors based on the points selected...
   // Now for each vector we assume the point e_i which is 1 when i==j, otherwise 0s
   // and perform the projection using those axis vectors... what would happen?

   gsl_vector_view evec_first = gsl_matrix_column (evec, eigenVectorsToUse.first);
   gsl_vector_view evec_second = gsl_matrix_column (evec, eigenVectorsToUse.second);



   double first_eigenVector[non_zero_features];
   double second_eigenVector[non_zero_features];

   int nonZero = 0;
   for(int k = 0; k < n; k++){
      // for each feature
      if ( prevProjectionMatrix[k*3+ 2]  > 0.00000001){
          first_eigenVector[nonZero] = gsl_vector_get(&evec_first.vector,nonZero );
          second_eigenVector[nonZero] = gsl_vector_get(&evec_second.vector,nonZero);
          nonZero++;
      }
   }
   nonZero = 0;

   for(int k = 0; k < n; k++){
     // for each feature

        if ( prevProjectionMatrix[k*3+ 2]  > 0.00000001){
            // non-zero

            double v1[3] = {first_eigenVector[nonZero], second_eigenVector[nonZero],0};
            double w = MathHelper::Norm(v1,3);
            MathHelper::Normalize(v1);
            newProjectionMatrix[k*3+ 0] = v1[0];
            newProjectionMatrix[k*3+ 1] = v1[1];
            newProjectionMatrix[k*3+ 2] = w;
            nonZero++;
        }
        else {
            newProjectionMatrix[k*3+ 0] = 0;
            newProjectionMatrix[k*3+ 1] = 0;
            newProjectionMatrix[k*3+ 2] = 0;
        }
   }

   gsl_vector_free (eval);
   gsl_matrix_free (evec);
   gsl_matrix_free ( tmpMatrix);
}


void Statistics::CreateCovarianceMatrix(){
    if ( covarianceMatrix == NULL)
        gsl_matrix_free(covarianceMatrix);


    int n = _mapper->GetTotalAttributes();
    int numElements = dataset->GetTotalNumberOfElements();
    covarianceMatrix = gsl_matrix_alloc(n,n);

    for(int attr1 = 0; attr1 < n; attr1++){
        gsl_matrix_set(covarianceMatrix,attr1, attr1, ranges.at(attr1).variance);

        bool isFirstNumerical = _mapper->IsNumerical(attr1);

        Attribute attrib1 = GetAttr(attr1);
        int order1 = GetBestOrdering(attr1);

        for(int attr2 = attr1+1; attr2 < n ; attr2++){
               Attribute attrib2 = GetAttr(attr2);
               int order2 = GetBestOrdering(attr2);

              bool isSecondNumerical = _mapper->IsNumerical(attr2);

              float covariance = 0;
              int totPoints = 0;

              for(int j = 0; j < numElements; j++){

                      float v1,v2;


                      if ( isFirstNumerical){
                          v1 = dataset->GetElementValue(j, attr1);
                          if ( v1 == -99999) continue; // Missing data element ...
                      }
                      else {
                          int val = GetElementValue(j, attr1);
                          if ( val == -99999) continue; // Missing data element ...
                          int idx =  attrib1.GetIndexValue(val);
                           v1 = GetMidPointPDF(attr1, idx, order1);

                      }

                      if ( isSecondNumerical){
                          v2 = dataset->GetElementValue(j, attr2);
                          if ( v2 == -99999) continue; // Missing data element ...
                      }
                      else {
                          int val = dataset->GetElementValue(j, attr2);
                          if ( val == -99999) continue; // Missing data element ...
                          int idx =  attrib2.GetIndexValue(val);
                          v2 = GetMidPointPDF(attr2, idx, order2);
                      }


                      totPoints += 1;
                      //covariance += (v1- ranges.at(attr1).mean )*(v2 - ranges.at(attr2).mean);
                      // This next is how Molchanov calculated... by dividing by mean
                      covariance += (v1- ranges.at(attr1).mean )*(v2 - ranges.at(attr2).mean) / ranges.at(attr1).mean /ranges.at(attr2).mean ;



              }

              if (totPoints == 0)
                  covariance = 0;
              else
                  covariance /= (totPoints -1);

              gsl_matrix_set(covarianceMatrix,attr1, attr2, covariance);
              gsl_matrix_set(covarianceMatrix,attr2, attr1, covariance);
        }


    }
}

void  Statistics::CalculateMDS(vector<pair<float, float> > *normalizedMDSpoints){

    // use the similarity matrix that was defined previously to the class
    //float* SimilarityMatrix;

    int n = dataset->GetTotalNumberOfElements();

    double* SquaredSimilarityMatrix = (double*)malloc(sizeof(double)*n*n);
    for(int i = 0; i < n; i++){
        for(int j = i; j < n; j++){

            SquaredSimilarityMatrix[i*n +j] = SimilarityMatrix[i*n +j]*SimilarityMatrix[i*n +j];
            SquaredSimilarityMatrix[j*n +i] = SimilarityMatrix[j*n +i]*SimilarityMatrix[j*n +i];
        }
    }

    CalculateMDS(SquaredSimilarityMatrix, n, normalizedMDSpoints);

    free(SquaredSimilarityMatrix);
}

void Statistics::CalculateMDSWithPointIndices(vector<int>* indices, vector<pair<float, float> > *normalizedMDSpoints){
    int n = indices->size();
    double* SquaredSimilarityMatrix = (double*)malloc(sizeof(double)*n*n);
    for(int h = 0; h < n; h++){
        int i =  indices->at(h);
        for(int k = i; k < n; k++){
            int j = indices->at(k);
            SquaredSimilarityMatrix[h*n +k] = SimilarityMatrix[i*n +j]*SimilarityMatrix[i*n +j];
            SquaredSimilarityMatrix[k*n +h] = SimilarityMatrix[j*n +i]*SimilarityMatrix[j*n +i];
        }
    }

    CalculateMDS(SquaredSimilarityMatrix, n, normalizedMDSpoints);
    free(SquaredSimilarityMatrix);
}


void Statistics::CalculateMDS(double* similarityMatrixSquared, int n, vector<pair<float,float>>* normalizedMDSpoints){/*
    1. Set up the matrix of squared proximities P(2) = [p2].
    2. Apply the double centering:     ,
    where n is the number of objects.
    3. Extract the m largest positive eigenvalues λ1 . . . λm of B and the corresponding
    m eigenvectors e1 . . . em.
   */
   // squared proximities is in the weightedSimilarity measure
    gsl_matrix_view psquared = gsl_matrix_view_array (similarityMatrixSquared, n, n);

    gsl_matrix* J = gsl_matrix_alloc(n, n);

    for(int  i = 0; i < n; i++){
       for(int j = 0; j < n; j++){
           gsl_matrix_set(J, i, j, (-1.0/static_cast<float>(n))*1.0);
       }
       gsl_matrix_set(J, i, i, 1.0 - (1.0/static_cast<float>(n))*1.0);
    }


    gsl_matrix* t1 =  gsl_matrix_alloc(n,n);
    gsl_matrix* B = gsl_matrix_alloc(n,n);
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, &psquared.matrix, J, 0.0, t1);
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, -0.5, J, t1, 0.0, B);

     gsl_vector *eval = gsl_vector_alloc (n);
     gsl_matrix *evec = gsl_matrix_alloc (n, n);

     gsl_eigen_symmv_workspace * w =  gsl_eigen_symmv_alloc (n);

     gsl_eigen_symmv (B, eval, evec, w);
     gsl_eigen_symmv_free (w);
     // positive eigenvalues
     gsl_eigen_symmv_sort (eval, evec,
                            GSL_EIGEN_SORT_VAL_DESC);

     /*4. A m-dimensional spatial configuration of the n objects is derived from the
     coordinate matrix X = EmΛ1/2
     m , where Em is the matrix of m eigenvectors
     and Λm is the diagonal matrix of m eigenvalues of B, respectively.*/

     // eigen values has the diagonal matrix
     gsl_matrix* eigenvalues = gsl_matrix_alloc(2,2);
     gsl_matrix_set(eigenvalues, 0,1, 0);
     gsl_matrix_set(eigenvalues, 1,0, 0);
     for(int i =0; i < 2; i++)
     {
         // std::cout << "eigen val " << gsl_vector_get(eval,i) << std::endl;
         gsl_matrix_set(eigenvalues, i,i, sqrt(gsl_vector_get(eval,i)) );
     }
     //

     gsl_matrix* Em = gsl_matrix_alloc(n,2);
     for(int i = 0; i < 2; i++){
         // Get the first two eigen vectors
         gsl_vector_view evec_i = gsl_matrix_column (evec, i);

         // set each of them as a column
         for(int j = 0; j < n; j++){
             gsl_matrix_set(Em, j, i,  gsl_vector_get(&evec_i.vector,j) );
         }
     }

     gsl_matrix* finalPoints = gsl_matrix_alloc(n,2);
     gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, Em, eigenvalues, 0.0, finalPoints);

     vector<pair<float,float>> MDSpoints;

     float MDSmin_dim1 = 1, MDSmax_dim1 = 1, MDSmin_dim2 = 1,MDSmax_dim2 = 1;

     for(int i = 0; i < n ; i++){

         float x = gsl_matrix_get(finalPoints, i, 0);
         float y = gsl_matrix_get(finalPoints, i, 1);

         if ( i == 0){
            MDSmin_dim1 = x;
            MDSmax_dim1 = x;
            MDSmin_dim2 = y;
            MDSmax_dim2 = y;
         }
         MDSmin_dim1 = min(MDSmin_dim1, x);
         MDSmax_dim1 = max(MDSmax_dim1, x);
         MDSmin_dim2 = min(MDSmin_dim2, y);
         MDSmax_dim2 = max(MDSmax_dim2, y);


         MDSpoints.push_back( make_pair(x,y));
     }
     normalizedMDSpoints->clear();

     for(int index = 0; index < n; index++){
         pair<float, float> nonNormalizedPoint = MDSpoints.at(index);
         // for range 0  to 1
         float x = ((nonNormalizedPoint.first - MDSmin_dim1 )/ ( MDSmax_dim1 - MDSmin_dim1));
         float y = ((nonNormalizedPoint.second - MDSmin_dim2 )/ ( MDSmax_dim2 - MDSmin_dim2));
         normalizedMDSpoints->push_back( make_pair(x,y));
     }

     gsl_matrix_free(Em);
     gsl_matrix_free(eigenvalues);
     gsl_matrix_free(t1);
     gsl_matrix_free(B);
     gsl_matrix_free(evec);
     gsl_vector_free(eval);
}

double Statistics::GetCorrelation(int attrib1, int attrib2, bool skip[]){
    double value = 0;

    if (_mapper->IsNumerical(attrib1) && _mapper->IsNumerical(attrib2) ){
        // both numerical... we use Pearson correlation (?)
        value = fabs(CalculatePearsonCoefficient(attrib1, attrib2, skip));
    }

    if (!_mapper->IsNumerical(attrib1) && !_mapper->IsNumerical(attrib2) ){

        if ( _mapper->IsCategorical(attrib1) || _mapper->IsOrdinal(attrib1))
        {
            if ( _mapper->IsCategorical(attrib2) || _mapper->IsOrdinal(attrib2))
                value = CalculateCramersV(attrib1, attrib2, skip);
        }
    }

    if ( _mapper->IsNumerical(attrib1) && _mapper->IsOrdinal(attrib2)){
        value = CalculatePearsonCoefficientWithOrdinal(attrib1, attrib2, skip);
    }

    if ( _mapper->IsOrdinal(attrib1) && _mapper->IsNumerical(attrib2)){
        value = CalculatePearsonCoefficientWithOrdinal(attrib2, attrib1, skip);

    }
    return value;
}

double Statistics::CalculatePearsonCoefficientWithOrdinal(int attrib1, int ordinalAttrib2, bool skip[]){
    double meanX = 0, meanY = 0;
    float total = 0;

    Attribute attr = _mapper->GetAttribute(ordinalAttrib2);
    for(int j = 0; j < dataset->GetTotalNumberOfElements(); j++){

        if ( skip[j]) continue;

        double xi = dataset->GetElementValue(j, attrib1);
        int val =  static_cast<int>(dataset->GetElementValue(j, ordinalAttrib2));
        int idx = attr.GetIndexValue(val);

        meanX += xi;
        meanY += idx;
        total++;

    }

    meanX /= total;
    meanY /= total;

    double up = 0;
    double den1 = 0;
    double den2 = 0;

    for(int j = 0; j < dataset->GetTotalNumberOfElements(); j++){
        if ( skip[j]) continue;

        double xi = dataset->GetElementValue(j, attrib1);
        int val =  static_cast<int>(dataset->GetElementValue(j, ordinalAttrib2));
        int idx = attr.GetIndexValue(val);

         up += (xi - meanX)*(idx - meanY);
         den1 += (xi - meanX)*(xi- meanX);
         den2 += (idx - meanY)*(idx- meanY);
    }

    double den = sqrt(den1)*sqrt(den2);

    double r = up/ den;

    return r;

}


double Statistics::CalculatePearsonCoefficient(int attrib1, int attrib2, bool skip[]){
    double meanX = 0, meanY = 0;
    float total = 0;
    for(int j = 0; j < dataset->GetTotalNumberOfElements(); j++){

        if ( skip[j]) continue;

        double xi = dataset->GetElementValue(j, attrib1);
        double yi = dataset->GetElementValue(j, attrib2);
        meanX += xi;
        meanY += yi;
        total++;

    }

    meanX /= total;
    meanY /= total;

    double up = 0;
    double den1 = 0;
    double den2 = 0;

    for(int j = 0; j < dataset->GetTotalNumberOfElements(); j++){
        if ( skip[j]) continue;

        double xi = dataset->GetElementValue(j, attrib1);
        double yi = dataset->GetElementValue(j, attrib2);

         up += (xi - meanX)*(yi - meanY);
         den1 += (xi - meanX)*(xi- meanX);
         den2 += (yi - meanY)*(yi- meanY);
    }

    double den = sqrt(den1)*sqrt(den2);

    double r = up/ den;

    return r;

}
double Statistics::CalculateCramersV(int attrib1, int attrib2, bool skip[]){

    int totalElements = 0;

    double max_numOccurences = 0;
    int n =  _mapper->GetCategoricalSize(attrib1);
    int m =  _mapper->GetCategoricalSize(attrib2);

    int co_occurences[n*m];
    for(int j = 0; j < n*m ; j++){
        co_occurences[j] = 0;
    }

    Attribute attr1 = _mapper->GetAttribute(attrib1);
    Attribute attr2 = _mapper->GetAttribute(attrib2);

    for(int j = 0; j < dataset->GetTotalNumberOfElements(); j++){

        if ( skip[j]) continue;
        totalElements++;

        int v1 = dataset->GetElementValue(j, attrib1);
        int v2 = dataset->GetElementValue(j, attrib2);

        int id1 = attr1.GetIndexValue(v1);
        int id2 = attr2.GetIndexValue(v2);
        if ( id1 == Attribute::MISSING_INDEX || id2 == Attribute::MISSING_INDEX){
            continue;
        }

        co_occurences[id1*m + id2] += 1;
    }



    double X = 0;
    for(int id1 = 0; id1 < n ; id1++){
        double n_i = GetTotalOfCategoricalValueSkippable(attrib1, id1, skip); //freq1*totalElements;
        if ( n_i > max_numOccurences)
           max_numOccurences = n_i;

        for(int id2 = 0; id2 < m; id2++){
            double n_ij =  co_occurences[id1*m + id2];

            double n_j = GetTotalOfCategoricalValueSkippable(attrib2, id2, skip);
            double den = (n_i*n_j)/totalElements;
            double num = (n_ij - den)*(n_ij -den);
            X += (num/den);
        }
    }
    double V = sqrt( (X/totalElements)/min(n-1, m-1) );

    return V;

}


double Statistics::SpearmanCorrelation(vector<LocalPoint> *points){
    // first we need the mean of x and y
    double meanX = 0, meanY = 0;
    double total = 0;
    for(unsigned int i = 0; i < points->size(); i++){
       meanX += points->at(i).p[0];
       meanY += points->at(i).p[1];
       total +=1;
    }
    meanX /= total;
    meanY /= total;


    double up = 0;
    double den1 = 0;
    double den2 = 0;

    for(unsigned int i = 0; i < points->size(); i++){
         double xi = points->at(i).p[0];
         double yi = points->at(i).p[1];

         up += (xi - meanX)*(yi - meanY);
         den1 += (xi - meanX)*(xi- meanX);
         den2 += (yi - meanY)*(yi- meanY);
    }

    double den = sqrt(den1)*sqrt(den2);

    double r = up/ den;

    return r;
}

double Statistics::NumericalMutualInformation(int numericalAttribute1, int numericalAttribute2, bool skip[], int k ){
   // "Estimating Mutual Information"  Kraskov et al. 2004

   // The estimate of I(X,Y) = digamma(k) - 1/N sum_1^N(  digamma(n_x+1) + digamma(ny+1)  ) + digamma(N)
   // For each point
   //   1. Find the closest k_neighbors
   //   2. Project the k-th neighbor to each dimension X, Y ( i.e. just get its value at X or Y)
   //   3. Define distance e ( twice the maximum distance to either X or Y)
   //   4. Get n_x, n_y ... number of points lower or equal than distance e
   //   5. apply the formula

   // Note 1. in order for e to work, we need to normalize the values... after our skip phase....
    //     2. Ranking doesn't necessarily have to be done with euclidean distance
   double RANGE_X[2] = {-DBL_MAX ,DBL_MAX};// max and min
   double RANGE_Y[2] = {-DBL_MAX ,DBL_MAX};// max and min

   int n = dataset->GetTotalNumberOfElements();
   // int k = 2;


   //find the range for normalization
   for(int i = 0; i < n; i++){
       if ( skip[i]) continue;

       double xvalue = dataset->GetElementValue(i, numericalAttribute1);
       double yvalue = dataset->GetElementValue(i, numericalAttribute2);

       if ( xvalue > RANGE_X[0]) RANGE_X[0] = xvalue;
       if ( yvalue > RANGE_Y[0]) RANGE_Y[0] = yvalue;

       if ( xvalue < RANGE_X[1]) RANGE_X[1] = xvalue;
       if ( yvalue < RANGE_Y[1]) RANGE_Y[1] = yvalue;
   }

   // For each sample find the k-nearest neighbors
   // The distance matrix can be calculated to speed this up
   int nonSkipped = 0;
   double sum = 0;
   for(int i = 0; i < n; i++){
       if ( skip[i]) continue;

       vector<pair<int, double>> distances;

       // Normalized location
       double cx = (dataset->GetElementValue(i, numericalAttribute1) - RANGE_X[1]) /(RANGE_X[0] - RANGE_X[1]);
       double cy = (dataset->GetElementValue(i, numericalAttribute2) - RANGE_Y[1]) /(RANGE_Y[0] - RANGE_Y[1]);

       for(int j =0 ; j < n; j++){
           if ( skip[j] ) continue;
           if (i == j) continue;

           // normalized other
           double ox = (dataset->GetElementValue(j, numericalAttribute1) - RANGE_X[1]) /(RANGE_X[0] - RANGE_X[1]);
           double oy = (dataset->GetElementValue(j, numericalAttribute2) - RANGE_Y[1]) /(RANGE_Y[0] - RANGE_Y[1]);
           double d =sqrt( (cx-ox)*(cx-ox) + (cy-oy)*(cy-oy));
           distances.push_back(make_pair(j, d));
       }


       std::sort( distances.begin(), distances.end(),
              boost::bind(&std::pair<int, double>::second, _1) <
              boost::bind(&std::pair<int, double>::second, _2 ));

       // Get the kth value

       if (k >= distances.size())
           k = distances.size() -1;
       int kthIndex = distances.at(k).first;

       double ox = (dataset->GetElementValue(kthIndex, numericalAttribute1) - RANGE_X[1]) /(RANGE_X[0] - RANGE_X[1]);
       double oy = (dataset->GetElementValue(kthIndex, numericalAttribute2) - RANGE_Y[1]) /(RANGE_Y[0] - RANGE_Y[1]);

       ox = fabs(ox - cx);
       oy = fabs(oy - cy);

       double e = max(ox, oy); // distance to nearest projection, times two

       double nx = 0;
       double ny = 0;
       for(int j =0 ; j < n; j++){
           if ( skip[j] ) continue;
           if (i == j) continue;

           double ox = (dataset->GetElementValue(j, numericalAttribute1) - RANGE_X[1]) /(RANGE_X[0] - RANGE_X[1]);
           double oy = (dataset->GetElementValue(j, numericalAttribute2) - RANGE_Y[1]) /(RANGE_Y[0] - RANGE_Y[1]);
           ox = fabs(ox - cx);
           oy = fabs(oy - cy);

           if ( ox <= e) nx++;
           if ( oy <= e) ny++;
       }


       double tmp = MathHelper::DigammaFunction(nx+1) + MathHelper::DigammaFunction(ny+1);

       sum += tmp;
       nonSkipped++;
   }

   double avg = sum / nonSkipped;
   //Could this be negative....
   double I = MathHelper::DigammaFunction(k) + MathHelper::DigammaFunction(nonSkipped) - avg;
   return I;
}

double Statistics::MixedMutualInformation(int categoricalAttribute, int numericalAttribute, bool skip[], int k){

     // There are four things that we need to calculate
     // N = total points
     // N_xi = totalPoints that discrete value is the same as x_i
     // k = given
     // mi = number of neighbors that lie within the distance d

     int m =  _mapper->GetCategoricalSize(categoricalAttribute);
     Attribute attr = _mapper->GetAttribute(categoricalAttribute);

     vector<vector<int>*> samplesInClass;

     for(int i = 0; i < m; i++){

         std::string val = attr.GetDims().at(i).GetValue();
         QString tmp = QString::fromStdString(val);
         int value = tmp.toInt();

         vector<int>* indices = new vector<int>();
         dataset->FilterDatasetByCategoricalValue(categoricalAttribute, value, indices, skip);

         samplesInClass.push_back(indices);
     }


     int n = dataset->GetTotalNumberOfElements();

     double averageI = 0;
     int totalElements = 0;
     for(int i = 0; i < n; i++){
         if ( skip[i]) continue;
         if ( dataset->IsElementValueMissing(i, categoricalAttribute)) continue;
         if ( dataset->IsElementValueMissing(i, numericalAttribute)) continue;

         double I = 0;

         int value = static_cast<int>( dataset->GetElementValue(i, categoricalAttribute));
         int vIndex =  attr.GetIndexValue(value);

         // Now we need to use calculate the m

         // First we calculate the absolute difference between the current point
         double nvalue = dataset->GetElementValue(i, numericalAttribute);

         vector<double> differences;

         // and all other from the same class
         int sz = samplesInClass[vIndex]->size();
         for(int g = 0; g < sz; g++)
         {
             if ( samplesInClass[vIndex]->at(g) == i) continue;
             if ( dataset->IsElementValueMissing(samplesInClass[vIndex]->at(g), numericalAttribute)) continue;

             double oVal = dataset->GetElementValue(samplesInClass[vIndex]->at(g), numericalAttribute);

             double diff = fabs(nvalue- oVal);
             differences.push_back(diff);
         }

         totalElements++;

         sort(differences.begin(), differences.end());
         if ( differences.empty() ) continue;

         I += MathHelper::DigammaFunction(samplesInClass[vIndex]->size());

        double d = 0;
        if ( differences.size() >= k)
            d = differences.at(k-1); // because we start from 0-index
         else continue; //only if there's enough elements?
           // d = differences.back();

         // now how many are within that d...

         int m = 0;
         for(int g = 0; g < n; g++){
             if (skip[g]) continue;
             if ( dataset->IsElementValueMissing(g, numericalAttribute)) continue;
             double oVal = dataset->GetElementValue(g, numericalAttribute);
             double diff = fabs(nvalue- oVal);

             if (diff <= d){
                 m += 1;
             }
         }

         // This is actually the N_xi,
         I += MathHelper::DigammaFunction(m);

         averageI += I;
     }

     averageI /= totalElements;


     averageI = (MathHelper::DigammaFunction(totalElements) + MathHelper::DigammaFunction(k))  - averageI;
     return averageI;
}


double Statistics::CategoricalMutualInformation(int categoricalAttribute1, int categoricalAttribute2, bool skip[]){
    // For categorical mutual information
    // We calculate the joint probability and the marginal probability of each category

    int x_n = _mapper->GetCategoricalSize(categoricalAttribute1);
    int y_n = _mapper->GetCategoricalSize(categoricalAttribute2);


    double total = 0;

    double frequencies1[x_n];
    double frequencies2[y_n];
    int total1[x_n];
    int total2[y_n];

    double totalB = 0;
    for(int i = 0; i < x_n; i++){
        double n_j = GetTotalOfCategoricalValueSkippable(categoricalAttribute1, i, skip);
        frequencies1[i] = n_j;
        total1[i] = n_j;
        // We skip the same quantity of values here and the next one, so we only need to
        // calculate one total
        total += n_j;
    }

    for(int i = 0; i < y_n; i++){
        double n_j = GetTotalOfCategoricalValueSkippable(categoricalAttribute2, i, skip);
        total2[i] = n_j;
        frequencies2[i] = n_j;
        totalB += n_j;
    }
    //The total should be the same for both...

    double Entropy = 0;

    for(int x = 0; x < x_n; x++){

        double p_x = frequencies1[x] / total;
        if ( total1[x] == 0) continue;
        // if there's no element, the probability is 0, which will cause the
        // division in the log to be undefined

        for(int y = 0; y < y_n; y++){
            double p_y = frequencies2[y] / total;
            if ( total2[y] == 0) continue;

            int n = GetTotalOfTwoCategoricalValueSkippable(categoricalAttribute1,x, categoricalAttribute2, y, skip);

            if ( n == 0 ) continue; // undefined value otherwise

            double jointProbability = static_cast<float>(n)/total;
            double add = log2(jointProbability/(p_x*p_y));


            Entropy += jointProbability*add;
        }
    }
    return Entropy;
}

double Statistics::MixedNormalizedMutualInformation(int categoricalAttribute, int numericalAttribute, bool skip[]){
   // We need a measure of association between the categorical and numerical attribute
   // Strength of association based on covariance and correlation doesn't cut it as we
   // have no definition of variance in the categorical attributes,
   // but still we could use things like ANOVA and Kruskal-willis test to figure out
   // whether the numerical attributes means and variances differ given the categorical attribute
   // but they do not have a normalized value that can be used as a source of comparision
   // So we look at the mutual information. Which can be done with binning or nearest neighbor
   // the binning is not nice, as there is an overreliance of bin size, so we use the
   // Discrete-Continuous Mutual Information as defined by  Brian C. Ross in
   // Mutual Information between Discrete and Continuous Data Sets
   // as for there is no upper bound for I(X,Y), a normalized version could be quite useful
   // in this case, we use the normalized mutual information as defined by Strehl and Gosh, Cluster Ensembles
   // - A framework for Combining Multiple partitions.
   // The skip is in case of stratification, to filter out non-wanted samples.


    int amountNotSkip = 0;
    for(int i = 0; i < dataset->GetTotalNumberOfElements(); i++)
        if ( !skip[i]) amountNotSkip++;

    double I = 0;
    for(int j = 0; j < 10; j++){
        int k =  j +1;
        double tI =  MixedMutualInformation(categoricalAttribute, numericalAttribute, skip, k);
        I = max(I,tI);
    }


    float percent = 0.01;
    while( percent < 0.1){
        int k = percent*amountNotSkip;
        if ( k <= 0) k = 2;
        double tI =  MixedMutualInformation(categoricalAttribute, numericalAttribute, skip, k);
        I = max(I,tI);
        percent += 0.01;
    }

   // by def. Mutual Info cannot be negative.. if it does, it is because of the substraction terms in
   // the estimator
   if ( I < 0){
       //I = MixedMutualInformation(categoricalAttribute, numericalAttribute, skip, k);
       //k += 2;
       //std::cout << "Mutual info is negative " << I << ".." << std::endl;
       I = 0;// by definition Mutual Information should be positive, if its negative then
       // there's an issue with estimation... a problem is that a high K value is basically forcing a positive
       // value, and would be dependent on the attribute instead of being used through all the application
   }


   double HX = CategoricalEntropy(categoricalAttribute, skip);

   double HY = NumericalEntropy(numericalAttribute, skip);
   double HYY = NumericalKernelEntropy(numericalAttribute, skip, 1000);


   if (I > max(HX,HY))
       HY = HYY;
   //std::cout << "I " << I  << std::endl;
   //std::cout << HX << " <---------> " << HY  << std::endl;


   //double NMI = I / sqrt(HX*HY);
   double NMI = I / max(HX,HY);

   if (isnan(NMI ) || isinf(NMI))
       NMI = 0;
   if ( NMI > 1.0 || isinf(NMI)){
       std::cout << "Mixed Info " << I << ".." << HX << " .."  <<  HYY << ":: " << NMI << "-> 1,2,4 " << MixedMutualInformation(categoricalAttribute, numericalAttribute, skip, 1) << "::" <<MixedMutualInformation(categoricalAttribute, numericalAttribute, skip, 2) << "::" << MixedMutualInformation(categoricalAttribute, numericalAttribute, skip, 4) << std::endl;

   }

   return NMI;
}

double Statistics::NumericalNormalizedMutualInformation(int numericalAttribute1, int numericalAttribute2, bool skip[]){
    // There are several measures of association for numerical data, mostly we measure the relationship
    // either linear or monotonic between the variables. Another way is cosine similarity measure, Li distance dissimilarity.
    // (examples of the calculation can be seen in Laplacian Star Coordinates).
    // A more general approach is looking at the Mutual Information, ... Information gain is helpful, but it is not a true metric
    // ... its relavite to Kullback-Leibler divergence...
    // Also, given that there's we might not know the data distribution of numerical info, we have to calculate it
    // in a different fashion. 1) Bins... create bins and apply categorical, but then we have issues
    // of the sampling rate. 2)Adaptive Bins... they perform much better, but they have several issues.
    // 3) Kernel density estimators
    // 4) Nearest-neighbors methods
    // Here we apply the nearest neighbor method as defined by Kraskov et al. (2004), Estimating mutual information

    int amountNotSkip = 0;
    for(int i = 0; i < dataset->GetTotalNumberOfElements(); i++)
        if ( !skip[i]) amountNotSkip++;


    double I = NumericalMutualInformation(numericalAttribute1, numericalAttribute2, skip, 5);

    /*
    double I = 0;
    for(int j = 0; j < 10; j++){
        double tI = NumericalMutualInformation(numericalAttribute1, numericalAttribute2, skip, j*5 + 5);
        I = max(I,tI);
    }

    float percent = 0.02;
    while( percent <0.2){
        int k = percent*amountNotSkip;
        double tI = NumericalMutualInformation(numericalAttribute1, numericalAttribute2, skip,k);
        I = max(I,tI);
        percent += 0.02;
    }*/


    // by def. Mutual Info cannot be negative.. if it does, it is because of the substraction terms in
    // the estimator
    if ( I < 0) {
        //k++;
        //I = NumericalMutualInformation(numericalAttribute1, numericalAttribute2, skip, k) ;
        I = 0;
    }

    //std::cout << "k > " << k << std::endl;
    double HX = NumericalEntropy(numericalAttribute1, skip);
    double HXX = NumericalKernelEntropy(numericalAttribute1, skip, 50);

    double HY = NumericalEntropy(numericalAttribute2, skip);
    double HYY = NumericalKernelEntropy(numericalAttribute2, skip, 50);

    if ( I > max(HX, HY)){
        HX = HXX;
        HY = HYY;
    }

    //double NMI = I / sqrt(HX*HY);
    double NMI = I / max(HX,HY);
    //if(isnan(NMI)) NMI = 0;

    if (isnan(NMI))
        NMI = 0;

    if (NMI > 1.0){
        std::cout << "I" <<I << " .. " << HX << "," << HY << std::endl;
        std::cout << "Numerical Mutual Information " << NMI << std::endl;

    }

    return NMI;
}


double Statistics::CategoricalNormalizedMutualInformation(int categoricalAttribute1, int categoricalAttribute2, bool skip[]){

  double I = CategoricalMutualInformation(categoricalAttribute1, categoricalAttribute2, skip);
  double HX = CategoricalEntropy(categoricalAttribute1, skip);

  // by def. Mutual Info cannot be negative.. if it does, it is because of the substraction terms in
  // the estimator



  if ( I < 0) I = 0;
  if ( categoricalAttribute1 == 2 && HX < 0.001){
      std::cout << "Categorical " << categoricalAttribute1 <<" __ " << categoricalAttribute2 << " Var names " <<_mapper->GetAttributeName(categoricalAttribute1) << " , " << _mapper->GetAttributeName(categoricalAttribute2 ) << std::endl;
      I = 0;
  }

  if (isnan(I)){
      std::cout << "MI" <<_mapper->GetAttributeName(categoricalAttribute1) << " , " << _mapper->GetAttributeName(categoricalAttribute2 ) << std::endl;
      I = 0;
  }
  if (isnan(HX)){
      std::cout << "Categorical Entropy ? " << categoricalAttribute1 << std::endl;
  }
  double HY = CategoricalEntropy(categoricalAttribute2, skip);

  if ( categoricalAttribute2 == 2 && HY < 0.001){
      std::cout << "Categorical " << categoricalAttribute1 <<" __ " << categoricalAttribute2 << " Var names " <<_mapper->GetAttributeName(categoricalAttribute1) << " , " << _mapper->GetAttributeName(categoricalAttribute2 ) << std::endl;
  }
  // std::cout << "I? " << I  << "..HX.." << HX << " , " << HY << " " <<  std::endl;
 // double NMI = I / sqrt(HX*HY);
  double NMI = I / max(HX,HY);

  if (isnan(NMI))
      NMI = 0;

  if ( NMI > 1){
      std::cout << "Categorical Mutual Inf " << NMI << std::endl;
  }

  return NMI;
}

double  Statistics::CategoricalEntropy(int attrib, bool skip[]){
    // Categorical entropy is defined as by Claude Shannon > Information Entropy
    // as  sum_i=1^n p(x_i) log p(x_i)
    // so first, for each category
    int m =  _mapper->GetCategoricalSize(attrib);
    double total = 0;
    double frequencies[m];
    int amount[m];
    Attribute attr = _mapper->GetAttribute(attrib);
    // calculate the number of times it appears
    int amountSkip = 0;
    for(int i = 0; i < dataset->GetTotalNumberOfElements(); i++)
        if ( skip[i]) amountSkip++;


    for(int i = 0; i < m; i++){

        std::string v = attr.GetDims().at(i).GetValue();
        int value = QString::fromStdString(v).toInt();
        vector<int> tindices;
        dataset->FilterDatasetByCategoricalValue(attrib, value, &tindices, skip );
        double n_j = tindices.size(); //GetTotalOfCategoricalValueSkippable(attrib, , skip);
        frequencies[i] = n_j;
        amount[i] = static_cast<int>(n_j);
        total += n_j;
    }
    // Now that we have the total amount it appears
    // and the total number we calculate the sum

    double Entropy = 0;

    // Find the probabilibty p(x_i)
    for(int i = 0; i < m; i++){

        double p = frequencies[i] / total;
        // it may be possible that the category doesnt appear
        // then it generates a log(0) which is undefined
        // therefore we remove that
        if ( amount[i] > 0 )
            Entropy += p*log2(p);
    }

    Entropy *= (-1);

    //Error checking... it should not print this ....
    if ( isnan(Entropy) ){


        std::cout << "Amount skipped were "  << amountSkip << "  // " << dataset->GetTotalNumberOfElements() << std::endl;

        std::cout <<"Attribute " << _mapper->GetAttributeName(attrib) << " ..has nan entropy " << total << "::" << attrib << std::endl;
        for(int i = 0; i < m; i++){
            std::cout << "F " <<i << " :: " << frequencies[i] << std::endl;
        }

        Entropy = 0;

        std::cout << std::endl << std::endl;
    }

    if ( Entropy < 0){
        std::cout << "Negative entropy should not exist?" << std::endl;
    }

    return Entropy;
}

double  Statistics::NumericalEntropy(int attrib, bool skip[]){

    // Estimates of entropy based on nearest neigbor distance,
    // Given the formula in non-parametric entropy: an overview -> Beirlant
    // Equation 18. Nearest neighbor estimate

    // Let ro_n,i be the nearest neighbor distance of X_i and the other X_j
    // r_o,i = min j != i, j <= n || X_i - X_j || ... i.e. the distance to the nearest neighbor
    // Then the nearest neighbor estimate is
    //  [ (1/n)*sum_i=1^n  ln( n*ro_n_i) ] + ln2 + Ce
    // where Ce is EulerMascheroni constant

    vector<double> values;

    // So first we get all the values from the dataset, that are not skipped
    int n = dataset->GetTotalNumberOfElements();
    for(int i = 0; i < n; i++){
        if ( skip[i]) continue;
        if (dataset->IsElementValueMissing(i, attrib)) continue;

        double val = dataset->GetElementValue(i, attrib);

        values.push_back(val);
    }

    // Given that we are looking at a single variable, this is the same as
    // they being a sorted, and the closest one located at either side

    sort(values.begin(), values.end());

    int totalElements = values.size();

    double sum = 0;
    for(int i = 0; i < values.size(); i++){
        // distance to the nearest neighbor
        double valueLeft = DBL_MAX;
        double valueRight = DBL_MAX;

        if ( i - 1 >= 0) valueLeft = fabs(values.at(i-1) - values.at(i));
        if ( i + 1 < values.size() -1) valueRight = fabs(values.at(i+1)- values.at(i));

        // Get the mininmum distance, either to the left or right
        double ro = std::min(valueLeft, valueRight);

        // std::cout << "i:: " << i << " (" << valueLeft << ", " << valueRight << ")" << ro << std::endl;
        //log returns the natural logarithm

        double nro_n_i = totalElements*ro;

        // Ro should always be positive, there cannot be a negative distance
        if ( ro < 0)
            std::cout << "Error " << std::endl;


        //It may be that the distance to the
        // The distance may be 0 as well
        if ( fabs(ro) > 0.00000001 )
             sum += log(nro_n_i);

    }


    // Now we take the average...
    sum /= totalElements;
    // basically now that it is sorted

    // std::cout << "Average " << sum  << " total elements  " << totalElements << std::endl;

    // and sum the other two values...
    double Entropy = sum + log(2.0 ) + MathHelper::EulerMascheroniConstant;
    return Entropy;
}


double Statistics::GetDimensionSimilarity(int attrib1, int attrib2, bool skip[]){

    // If both attributes are numerical, we calculate the spearman correlation, re-ranged between 0..1

    double similarity = 0;
    if (_mapper->IsNumerical(attrib1) && _mapper->IsNumerical(attrib2) ){
        //std::cout << "Dims " << attrib1 << ","<< attrib2 << " NN " <<std::endl;
        similarity = NumericalNormalizedMutualInformation(attrib1, attrib2, skip);

        // std::cout << "Numerical Normalized ?"  << similarity << " btw" << attrib1 << " / "<< attrib2 << std::endl;
    }

    // if both attributes are categorical or ordinal we use Cramer's V
    if (!_mapper->IsNumerical(attrib1) && !_mapper->IsNumerical(attrib2) ){

        if ( _mapper->IsCategorical(attrib1) || _mapper->IsOrdinal(attrib1))
        {
            if ( _mapper->IsCategorical(attrib2) || _mapper->IsOrdinal(attrib2))
                //std::cout << "Dims " << attrib1 << ","<< attrib2 << " CC " <<std::endl;
                similarity = CategoricalNormalizedMutualInformation(attrib1, attrib2, skip);

        }
        //similarity = CalculateCramersV(attrib1, attrib2, skip);
    }


    // if one attribute is numerical and the other one is categorical/ordinal then we
    // use normalized mutual information
    if (!_mapper->IsNumerical(attrib1) && _mapper->IsNumerical(attrib2) ){
        // attrib2 is numerical, attrib1 is categorical or ordinal

        if ( _mapper->IsCategorical(attrib1) || _mapper->IsOrdinal(attrib1)){
            //std::cout << "Dims " << attrib1 << ","<< attrib2 << " CN " <<std::endl;

            similarity = MixedNormalizedMutualInformation(attrib1, attrib2, skip);

        }
       // std::cout << "Mixed similarity " << similarity  << " btw" << attrib1 << " / "<< attrib2 << std::endl;

    }

    if (_mapper->IsNumerical(attrib1) && !_mapper->IsNumerical(attrib2) ){
        // attrib2 is numerical, attrib1 is categorical or ordinal
        if ( _mapper->IsCategorical(attrib2) || _mapper->IsOrdinal(attrib2)){
            // std::cout << "Dims " << attrib1 << ","<< attrib2 << " NC " <<std::endl;

            similarity = MixedNormalizedMutualInformation(attrib2, attrib1, skip);
        }

       // std::cout << "Mixed similarity " << similarity  << " btw" << attrib1 << " / "<< attrib2 << std::endl;

    }

    if ( similarity > 1.0 || isinf(similarity)){

        std::cout << "Error " << similarity << ".." << _mapper->GetAttributeName(attrib1) << ", " << _mapper->GetAttributeName(attrib2) << std::endl;
    }


    return similarity;
}

void Statistics::CreateLaplacianProjectionMatrix(float* prevProjectionMatrix, float *newProjectionMatrix, bool normalize){
   // Based on Tran Van Long Laplacian Star Coordinates

    std::cout << "Creating Laplacian Projection Matrix " << std::endl;
    int non_zero_features = 0;
    int n = _mapper->GetTotalAttributes();
    vector<int> nonZeroAttrs;
    for(int i = 0; i < n ; i++){
     if ( prevProjectionMatrix[i*3+ 2]  > 0.00000001) {
         non_zero_features++;
         nonZeroAttrs.push_back(i);
     }
    }

    // We define the Similarity Matrix and the Degree Matrix
    gsl_matrix* S = gsl_matrix_alloc(non_zero_features, non_zero_features);

    bool skip[dataset->GetTotalNumberOfElements()];
    for(int i = 0; i < dataset->GetTotalNumberOfElements(); i++)
        skip[i] = false;


    // Symmetric, could basically avoid half the computations here ....
    for(unsigned int j = 0; j < nonZeroAttrs.size(); j++){
        for(unsigned int k =0; k < nonZeroAttrs.size(); k++){

               if ( j == k)
                   gsl_matrix_set(S, j, k, 1.0 );

               else {
                   double Sim = GetDimensionSimilarity(nonZeroAttrs.at(j), nonZeroAttrs.at(k) , skip );
                   std::cout << "Similarity? "<< j << "," << k << " :: " << Sim << std::endl;
                   gsl_matrix_set(S, j, k,  Sim);
               }
        }
    }

    gsl_matrix* D = gsl_matrix_alloc(non_zero_features, non_zero_features);


    for(unsigned int j = 0; j < nonZeroAttrs.size(); j++){

        double sum = 0;

        for(unsigned int k =0; k < nonZeroAttrs.size(); k++){
            gsl_matrix_set(D, j, k,  0 );
            sum += gsl_matrix_get(S, j, k);
        }
        gsl_matrix_set(D, j, j,  sum );
    }

    gsl_matrix* C  = gsl_matrix_alloc(nonZeroAttrs.size(), nonZeroAttrs.size());
    gsl_matrix* invD  = gsl_matrix_alloc(nonZeroAttrs.size(), nonZeroAttrs.size());

    int d = nonZeroAttrs.size();

    MathHelper::InverseMatrixGSL2(D,invD, d);

    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, invD, S, 0.0, C);

    gsl_vector *eval = gsl_vector_alloc (d);
    gsl_matrix *evec = gsl_matrix_alloc (d, d);

    gsl_eigen_symmv_workspace * w =
      gsl_eigen_symmv_alloc (d);

    gsl_eigen_symmv (C, eval, evec, w);

    gsl_eigen_symmv_free (w);

    gsl_eigen_symmv_sort (eval, evec,
                          GSL_EIGEN_SORT_VAL_DESC);



    gsl_vector_view evec_first = gsl_matrix_column (evec, 0);
    gsl_vector_view evec_second = gsl_matrix_column (evec, 1);

    int nonZero = 0;

    double first_eigenVector[non_zero_features];
    double second_eigenVector[non_zero_features];

    std::cout << "Eigen Vectors " << std::endl;
    for(int k = 0; k < n; k++){
       // for each feature
       if ( prevProjectionMatrix[k*3+ 2]  > 0.00000001){

           first_eigenVector[nonZero] = gsl_vector_get(&evec_first.vector,nonZero );
           second_eigenVector[nonZero] = gsl_vector_get(&evec_second.vector,nonZero);
           std::cout <<  first_eigenVector[nonZero]  << " , " <<  second_eigenVector[nonZero] << std::endl;

           nonZero++;
       }
    }


    nonZero = 0;

    for(int k = 0; k < n; k++){
      // for each feature

         if ( prevProjectionMatrix[k*3+ 2]  > 0.00000001){
             // non-zero

              //for (i = 0; i < 2; i++)
              //  {
             //double first_eigenvalue = gsl_vector_get (eval, 0);
             //double second_eigenvalue = gsl_vector_get (eval, 1);


             // According to Molchanovś  code, we use the normalized eigenvectors to calculate it
             //double v1[3] = {first_eigenvalue*gsl_vector_get(&evec_first.vector,nonZero ), second_eigenvalue*gsl_vector_get(&evec_second.vector,nonZero ), 0};
             //double w = MathHelper::Norm(v1,3);
             //MathHelper::Normalize(v1);

             double v1[3] = {first_eigenVector[nonZero], second_eigenVector[nonZero],0};
             double w = MathHelper::Norm(v1,3);
             MathHelper::Normalize(v1);


             newProjectionMatrix[k*3+ 0] = v1[0];
             newProjectionMatrix[k*3+ 1] = v1[1];

             newProjectionMatrix[k*3+ 2] = w;
             if ( normalize)
                 newProjectionMatrix[k*3+ 2] = 1.0;

                 /* printf ("eigenvalue = %g\n", eval_i);
                  printf ("eigenvector = \n");


                  gsl_vector_fprintf (stdout,
                                      &evec_i.vector, "%g");*/
             //   }



             nonZero++;
         }
         else {
             newProjectionMatrix[k*3+ 0] = 0;
             newProjectionMatrix[k*3+ 1] = 0;
             newProjectionMatrix[k*3+ 2] = 0;
         }

    }



    gsl_matrix_free(S);
    gsl_matrix_free(D);
    gsl_matrix_free(evec);
    gsl_vector_free(eval);
    gsl_matrix_free(C);
}


void Statistics::GetDensityDistributionViaBinning(int attributeIndex, int numberOfBins, vector<double>* bins, bool skip[]){

    // For now this is only taking into account the full distribution
    if ( numberOfBins == 0) return;

    double step = 1.0 / numberOfBins;
    double start = step / 2.0;

    double accumulate[numberOfBins];
    for(int i = 0; i < numberOfBins; i++) accumulate[i] = 0;

    float max = GetMaximumInAttribute(attributeIndex);
    float min = GetMinimumInAttribute(attributeIndex);

    int n = dataset->GetTotalNumberOfElements();
    for(int i = 0; i < n; i++){
        if ( skip[i]) continue;
        double val = dataset->GetElementValue(i, attributeIndex);

        double nval = (val - min)/(max-min);

        // normalized value

        start = step /2.0;
        for(int k =0; k < numberOfBins; k++){
            accumulate[k] += MathHelper::NormalDistribution(start - nval, 0.2 );
            start += step;
        }
    }

    /*
    float maxVal = 0;
    for(int k = 0; k < numberOfBins; k++)
        if ( accumulate[k] > maxVal) maxVal = accumulate[k];
    */
    // The sum over all the area should be 1
    float maxVal = 0;
    for(int k = 0; k < numberOfBins; k++)
         maxVal += accumulate[k];


    for(int k =0; k < numberOfBins; k++){
        double normalized = accumulate[k] / maxVal;
        bins->push_back(normalized);
    }
}



void Statistics::CreateSkippableLaplacianProjectionMatrix(float* prevProjectionMatrix, float *newProjectionMatrix, bool normalize, bool skip[]){
   // Based on Tran Van Long Laplacian Star Coordinates

    std::cout << "Creating Laplacian Projection Matrix " << std::endl;
    int non_zero_features = 0;
    int n = _mapper->GetTotalAttributes();
    vector<int> nonZeroAttrs;
    for(int i = 0; i < n ; i++){
     if ( prevProjectionMatrix[i*3+ 2]  > 0.00000001) {
         non_zero_features++;
         nonZeroAttrs.push_back(i);
     }
    }

    // We define the Similarity Matrix and the Degree Matrix
    gsl_matrix* S = gsl_matrix_alloc(non_zero_features, non_zero_features);

    // Symmetric, could basically avoid half the computations here ....
    for(unsigned int j = 0; j < nonZeroAttrs.size(); j++){
        for(unsigned int k =0; k < nonZeroAttrs.size(); k++){

               if ( j == k)
                   gsl_matrix_set(S, j, k, 1.0 );

               else {
                   double Sim = GetDimensionSimilarity(nonZeroAttrs.at(j), nonZeroAttrs.at(k) , skip );
                   std::cout << "Similarity? "<< j << "," << k << " :: " << Sim << std::endl;
                   gsl_matrix_set(S, j, k,  Sim);
               }
        }
    }

    gsl_matrix* D = gsl_matrix_alloc(non_zero_features, non_zero_features);


    for(unsigned int j = 0; j < nonZeroAttrs.size(); j++){

        double sum = 0;

        for(unsigned int k =0; k < nonZeroAttrs.size(); k++){
            gsl_matrix_set(D, j, k,  0 );
            sum += gsl_matrix_get(S, j, k);
        }
        gsl_matrix_set(D, j, j,  sum );
    }

    gsl_matrix* C  = gsl_matrix_alloc(nonZeroAttrs.size(), nonZeroAttrs.size());
    gsl_matrix* invD  = gsl_matrix_alloc(nonZeroAttrs.size(), nonZeroAttrs.size());

    int d = nonZeroAttrs.size();

    MathHelper::InverseMatrixGSL2(D,invD, d);

    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, invD, S, 0.0, C);

    gsl_vector *eval = gsl_vector_alloc (d);
    gsl_matrix *evec = gsl_matrix_alloc (d, d);

    gsl_eigen_symmv_workspace * w =
      gsl_eigen_symmv_alloc (d);

    gsl_eigen_symmv (C, eval, evec, w);

    gsl_eigen_symmv_free (w);

    gsl_eigen_symmv_sort (eval, evec,
                          GSL_EIGEN_SORT_VAL_DESC);



    gsl_vector_view evec_first = gsl_matrix_column (evec, 0);
    gsl_vector_view evec_second = gsl_matrix_column (evec, 1);

    int nonZero = 0;

    double first_eigenVector[non_zero_features];
    double second_eigenVector[non_zero_features];

    std::cout << "Eigen Vectors " << std::endl;
    for(int k = 0; k < n; k++){
       // for each feature
       if ( prevProjectionMatrix[k*3+ 2]  > 0.00000001){

           first_eigenVector[nonZero] = gsl_vector_get(&evec_first.vector,nonZero );
           second_eigenVector[nonZero] = gsl_vector_get(&evec_second.vector,nonZero);
           std::cout <<  first_eigenVector[nonZero]  << " , " <<  second_eigenVector[nonZero] << std::endl;

           nonZero++;
       }
    }


    nonZero = 0;

    for(int k = 0; k < n; k++){
      // for each feature

         if ( prevProjectionMatrix[k*3+ 2]  > 0.00000001){
             // non-zero

              //for (i = 0; i < 2; i++)
              //  {
             //double first_eigenvalue = gsl_vector_get (eval, 0);
             //double second_eigenvalue = gsl_vector_get (eval, 1);


             // According to Molchanovś  code, we use the normalized eigenvectors to calculate it
             //double v1[3] = {first_eigenvalue*gsl_vector_get(&evec_first.vector,nonZero ), second_eigenvalue*gsl_vector_get(&evec_second.vector,nonZero ), 0};
             //double w = MathHelper::Norm(v1,3);
             //MathHelper::Normalize(v1);

             double v1[3] = {first_eigenVector[nonZero], second_eigenVector[nonZero],0};
             double w = MathHelper::Norm(v1,3);
             MathHelper::Normalize(v1);


             newProjectionMatrix[k*3+ 0] = v1[0];
             newProjectionMatrix[k*3+ 1] = v1[1];

             newProjectionMatrix[k*3+ 2] = w;
             if ( normalize)
                 newProjectionMatrix[k*3+ 2] = 1.0;

                 /* printf ("eigenvalue = %g\n", eval_i);
                  printf ("eigenvector = \n");


                  gsl_vector_fprintf (stdout,
                                      &evec_i.vector, "%g");*/
             //   }



             nonZero++;
         }
         else {
             newProjectionMatrix[k*3+ 0] = 0;
             newProjectionMatrix[k*3+ 1] = 0;
             newProjectionMatrix[k*3+ 2] = 0;
         }

    }



    gsl_matrix_free(S);
    gsl_matrix_free(D);
    gsl_matrix_free(evec);
    gsl_vector_free(eval);
    gsl_matrix_free(C);
}

double Statistics::CategoricalKullBackDivergence(int attributeIndex, vector<int>& P, vector<int>& Q){
   //
   // - sum_i P(i) log(Q(i)/P(i))
    if ( P.empty() && Q.empty()) // weird case
        return -99999;

    if ( P.empty() && !Q.empty()){
        // P is the combined, i.e. those that have missing
        // in the other attribute, if P is empty
        // it means that everything is complete
        return -99; // green
    }
    else if ( !P.empty() && Q.empty()){
        //Q is the complete, this means then
        // that there are no variables where the mix are complete
        return -999;
    }
    int currentSize = _mapper->GetCategoricalSize(attributeIndex);
    Attribute att = _mapper->GetAttribute(attributeIndex);
    int frequenciesP[currentSize];
    int frequenciesQ[currentSize];
    for(int i =0; i < currentSize; i++){
        frequenciesP[i] = 0;
        frequenciesQ[i] = 0;
    }

    for(unsigned int i = 0; i < P.size(); i++){
        int value = static_cast<int>( dataset->GetElementValue( P.at(i), attributeIndex));
        int vIndex =  att.GetIndexValue(value);
        if ( vIndex == Attribute::MISSING_INDEX){
           //std::cout << "This should not happen  " << att.GetName() << "?? " << value << " element " << i << " . " <<  attributeIndex << ",P::" << P.at(i) << std::endl;
        }
        else
           frequenciesP[vIndex] += 1;
    }

    for(unsigned int i = 0; i < Q.size(); i++){
        int value = static_cast<int>( dataset->GetElementValue( Q.at(i), attributeIndex));
        int vIndex =  att.GetIndexValue(value);
        if ( vIndex == Attribute::MISSING_INDEX){
           //std::cout << "This should not happen  " << att.GetName() << "?? " << value << " element " << i << " . " <<  attributeIndex << ",Q::" << Q.at(i) << std::endl;
        }
        else
            frequenciesQ[vIndex] += 1;
    }


    double divergence = 0;

    for(int i = 0; i < currentSize; i++){
        double P_i = static_cast<double>(frequenciesP[i]) / static_cast<double>(P.size());
        double Q_i = static_cast<double>(frequenciesQ[i]) / static_cast<double>(Q.size());

        if ( P_i > 0.000001 && Q_i > 0.00000001){
           divergence += P_i*log2( Q_i/P_i);
           // divergence += P_i*log( Q_i/P_i);

        }

    }
    divergence *= -1;




    return divergence;
}


double  Statistics::CategoricalKullBackDivergenceGivenValues(int attributeIndex, vector<int>& P, vector<int>& Q){
    if ( P.empty() && Q.empty()) // weird case
        return -99999;

    if ( P.empty() && !Q.empty()){
        // P is the combined, i.e. those that have missing
        // in the other attribute, if P is empty
        // it means that everything is complete
        return -99; // green
    }
    else if ( !P.empty() && Q.empty()){
        //Q is the complete, this means then
        // that there are no variables where the mix are complete
        return -999;
    }

    int currentSize = _mapper->GetCategoricalSize(attributeIndex);
    int frequenciesP[currentSize];
    int frequenciesQ[currentSize];
    for(int i =0; i < currentSize; i++){
        frequenciesP[i] = 0;
        frequenciesQ[i] = 0;
    }

    for(unsigned int i = 0; i < P.size(); i++){
        frequenciesP[P.at(i)] += 1;
    }
    for(unsigned int i = 0; i < Q.size(); i++){
        frequenciesQ[Q.at(i)] += 1;
    }

    double divergence = 0;

    for(int i = 0; i < currentSize; i++){
        double P_i = static_cast<double>(frequenciesP[i]) / static_cast<double>(P.size());
        double Q_i = static_cast<double>(frequenciesQ[i]) / static_cast<double>(Q.size());

        if ( P_i > 0.000001 && Q_i > 0.00000001){
            divergence += P_i*log2( Q_i/P_i);
            //divergence += P_i*log( Q_i/P_i);

        }

    }
    divergence *= -1;


    if ( divergence < 0){
          double tmp = 0;
          double tmp2 = 0;
           for(int i = 0; i < currentSize; i++){

               double P_i = static_cast<double>(frequenciesP[i]) / static_cast<double>(P.size());
               double Q_i = static_cast<double>(frequenciesQ[i]) / static_cast<double>(Q.size());
               double v = 0;
               if ( P_i > 0.000001 && Q_i > 0.00000001){
                   v += P_i*log( Q_i/P_i);
                  tmp2 += P_i*log( P_i/Q_i);
               }
               tmp += v;
               // std::cout << i  << ":: " << P_i << ", " << Q_i <<  " v " << v << " tmp " << tmp << "/" << tmp2 << std::endl;

           }

    }

    return divergence;
}



double Statistics::NumericalKernelEntropy(int attrib, bool skip[], int numberOfBins){
    vector<double> values;


    int n = dataset->GetTotalNumberOfElements();
    for(int i = 0; i < n; i++){
        if ( skip[i]) continue;
        if (dataset->IsElementValueMissing(i, attrib)) continue;

        double val = dataset->GetElementValue(i, attrib);
        values.push_back(val);
    }


    double meanP = 0;
    double minV = DBL_MAX;
    double maxV = -DBL_MAX;

    for(int i = 0; i < values.size(); i++){
        double v = values.at(i);
        if (v < minV) minV = v;
        if (v > maxV) maxV = v;
        meanP += v;
    }

    meanP /= values.size();
    meanP = (meanP - minV)/(maxV - minV);
    double stdDevP = 0;
    for(int i = 0; i < values.size(); i++){
        double v = values.at(i);
        v =  (v-minV)/(maxV - minV);
        stdDevP +=  pow(meanP - v, 2.0);
    }
    stdDevP = sqrt(stdDevP/(values.size()-1));//-1 as is the sample standard deviation
    float optP = 1.06*stdDevP*pow(values.size(), -0.2);

    //int numberOfBins = 30;
    double step = 1.0 / numberOfBins;
    double start = step / 2.0;

    double accumulateP[numberOfBins];
    for(int i = 0; i < numberOfBins; i++){
        accumulateP[i] = 0;
    }
    for(int i = 0; i < values.size(); i++){
        double v = values.at(i);
        v =  (v-minV)/(maxV - minV);
        start = step /2.0;
        for(int k =0; k < numberOfBins; k++){
            accumulateP[k] += MathHelper::NormalDistribution(start - v, optP );
            start += step;
        }
    }

    double sumP = 0;

    for(int i = 0; i < numberOfBins; i++){
        sumP += accumulateP[i];
    }
    double entropy = 0;

    float sumP_i = 0;

    for(int i = 0; i < numberOfBins; i++){
        double P_i = accumulateP[i]/ sumP;
        sumP_i += P_i;
        if ( P_i > 0.000001 )
             entropy += P_i*log2(P_i);

    }

    entropy = entropy*(-1.0);
    return entropy;
}


double Statistics::NumericalKullBackDivergence(int attributeIndex, vector<int>& P, vector<int>& Q, int numBins){
   // The numerical kullback leibler divergence is performed doing kernel density estimation
   // so first...
   // calculate the std dev of the elements in P and Q, in order to get the optimum bandwith
    if ( P.empty() && Q.empty()) // weird case
        return -99999;

    if ( P.empty() && !Q.empty()){
        // P is the combined, i.e. those that have missing
        // in the other attribute, if P is empty
        // it means that everything is complete
        return -99; // green
    }
    else if ( !P.empty() && Q.empty()){
        //Q is the complete, this means then
        // that there are no variables where the mix are complete
        return -999;
    }

   double meanP = 0;
   double meanQ = 0;

   double minV = DBL_MAX;
   double maxV = -DBL_MAX;

   for(int i = 0; i < P.size(); i++){
       double v = dataset->GetElementValue( P.at(i), attributeIndex);
       if (v < minV) minV = v;
       if (v > maxV) maxV = v;
       meanP += v;
   }

   for(int i = 0; i < Q.size(); i++){
       double v = dataset->GetElementValue( Q.at(i), attributeIndex);
       if (v < minV) minV = v;
       if (v > maxV) maxV = v;
       meanQ += v;
   }

   meanP /= P.size();
   meanQ /= Q.size();

   meanP = (meanP - minV)/(maxV - minV);
   meanQ = (meanQ - minV)/(maxV - minV);

   double stdDevP = 0;
   double stdDevQ = 0;

   for(int i = 0; i < P.size(); i++){
       double v = dataset->GetElementValue( P.at(i), attributeIndex);
       v =  (v-minV)/(maxV - minV);
       stdDevP +=  pow(meanP - v, 2.0);
   }
   stdDevP = sqrt(stdDevP/P.size());

   for(int i = 0; i < Q.size(); i++){
       double v = dataset->GetElementValue( Q.at(i), attributeIndex);
       v =  (v-minV)/(maxV - minV);
       stdDevQ +=  pow(meanQ - v, 2.0);
   }
   stdDevQ = sqrt(stdDevQ/Q.size());

   float optP = 1.06*stdDevP*pow(P.size(), -0.2);
   float optQ = 1.06*stdDevQ*pow(Q.size(), -0.2);

   if (attributeIndex == 1)
          std::cout << "P" << optP << " , Q" << optQ << ".." << numBins<< std::endl;

   int numberOfBins = numBins;
   double step = 1.0 / numberOfBins;
   double start = step / 2.0;

   double accumulateP[numberOfBins];
   double accumulateQ[numberOfBins];

   for(int i = 0; i < numberOfBins; i++){
       accumulateP[i] = 0;
       accumulateQ[i] = 0;
   }

   for(int i = 0; i < P.size(); i++){
       double v = dataset->GetElementValue( P.at(i), attributeIndex);
       v =  (v-minV)/(maxV - minV);
       start = step /2.0;
       for(int k =0; k < numberOfBins; k++){
           accumulateP[k] += MathHelper::NormalDistribution(start - v, optP );
           start += step;
       }
   }

   for(int i = 0; i < Q.size(); i++){
       double v = dataset->GetElementValue( Q.at(i), attributeIndex);
       v =  (v-minV)/(maxV - minV);
       start = step /2.0;
       for(int k =0; k < numberOfBins; k++){
           accumulateQ[k] += MathHelper::NormalDistribution(start - v, optQ );
           start += step;
       }
   }

   double sumP = 0;
   double sumQ = 0;

   for(int i = 0; i < numberOfBins; i++){
       sumP += accumulateP[i];
       sumQ += accumulateQ[i];
   }

   double divergence = 0;

   double sumPi = 0;
   double sumQi = 0;
   for(int i = 0; i < numberOfBins; i++){
       double P_i = accumulateP[i]/ sumP;
       double Q_i = accumulateQ[i]/ sumQ;
       sumPi += P_i;
       sumQi += Q_i;
       if (attributeIndex == 1)
             std::cout << i<< "P,Q "<< P_i << "," << Q_i << divergence  << ".." << std::endl;
      if ( P_i > 0.000001 && Q_i > 0.00000001){
         // divergence += P_i*log( Q_i/P_i);

         divergence += P_i*log2( Q_i/P_i);
      }
   }
   if (attributeIndex == 1)
   {
       std::cout << "Sum " << sumPi << ", " << sumQi <<std::endl;
   }
   divergence *= -1;


   return divergence;

}



double Statistics::NumericalKullBackDivergenceGivenValues(int attributeIndex, vector<float>& P, vector<float>& Q){

    // The numerical kullback leibler divergence is performed doing kernel density estimation
    // so first...
    // calculate the std dev of the elements in P and Q, in order to get the optimum bandwith
     if ( P.empty() && Q.empty()) // weird case
         return -99999;

     if ( P.empty() && !Q.empty()){
         // P is the combined, i.e. those that have missing
         // in the other attribute, if P is empty
         // it means that everything is complete
         return -99; // green
     }
     else if ( !P.empty() && Q.empty()){
         //Q is the complete, this means then
         // that there are no variables where the mix are complete
         return -999;
     }
     double meanP = 0;
     double meanQ = 0;

     double minV = DBL_MAX;
     double maxV = -DBL_MAX;

     for(int i = 0; i < P.size(); i++){
         double v = P.at(i);
         if (v < minV) minV = v;
         if (v > maxV) maxV = v;
         meanP += v;
     }

     for(int i = 0; i < Q.size(); i++){
         double v = Q.at(i);
         if (v < minV) minV = v;
         if (v > maxV) maxV = v;
         meanQ += v;
     }

     meanP /= P.size();
     meanQ /= Q.size();

     meanP = (meanP - minV)/(maxV - minV);
     meanQ = (meanQ - minV)/(maxV - minV);

     double stdDevP = 0;
     double stdDevQ = 0;

     for(int i = 0; i < P.size(); i++){
         double v = P.at(i);
         v =  (v-minV)/(maxV - minV);
         stdDevP +=  pow(meanP - v, 2.0);
     }
     stdDevP = sqrt(stdDevP/P.size());

     for(int i = 0; i < Q.size(); i++){
         double v = Q.at(i);
         v =  (v-minV)/(maxV - minV);
         stdDevQ +=  pow(meanQ - v, 2.0);
     }
     stdDevQ = sqrt(stdDevQ/Q.size());

     float optP = 1.06*stdDevP*pow(P.size(), -0.2);
     float optQ = 1.06*stdDevQ*pow(Q.size(), -0.2);
     int numberOfBins = 100;
     double step = 1.0 / numberOfBins;
     double start = step / 2.0;

     double accumulateP[numberOfBins];
     double accumulateQ[numberOfBins];

     for(int i = 0; i < numberOfBins; i++){
         accumulateP[i] = 0;
         accumulateQ[i] = 0;
     }

     for(int i = 0; i < P.size(); i++){
         double v = P.at(i);
         v =  (v-minV)/(maxV - minV);
         start = step /2.0;
         for(int k =0; k < numberOfBins; k++){
             accumulateP[k] += MathHelper::NormalDistribution(start - v, optP );
             start += step;
         }
     }

     for(int i = 0; i < Q.size(); i++){
         double v =Q.at(i);
         v =  (v-minV)/(maxV - minV);
         start = step /2.0;
         for(int k =0; k < numberOfBins; k++){
             accumulateQ[k] += MathHelper::NormalDistribution(start - v, optQ );
             start += step;
         }
     }

     double sumP = 0;
     double sumQ = 0;

     for(int i = 0; i < numberOfBins; i++){
         sumP += accumulateP[i];
         sumQ += accumulateQ[i];
     }

     double divergence = 0;

     for(int i = 0; i < numberOfBins; i++){
         double P_i = accumulateP[i]/ sumP;
         double Q_i = accumulateQ[i]/ sumQ;
        if ( P_i > 0.000001 && Q_i > 0.00000001){
             divergence += P_i*log2( Q_i/P_i);
            //divergence += P_i*log( Q_i/P_i);
        }
     }

     divergence *= -1;
     return divergence;
}


