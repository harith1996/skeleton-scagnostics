#include "transferfunction.h"
#include <fstream>
#include <unordered_map>
TransferFunction::TransferFunction()
{
   isBanding = true;
   deuteranopiaColor = false;
}



void TransferFunction::GetColor(float value, float* r, float* g, float *b){
    if (nodeValues.empty())
    {
        *r = 0; *g = 0; *b = 0;
        return;
    }
    float color[3];
    GetColor(value, color);
    *r = color[0];
    *g = color[1];
    *b = color[2];
}

void TransferFunction::GetNodeColor(int i, float color[]){
    if (nodeValues.empty())
    {
        for(int k = 0; k < 3; k++) color[k] = 0;
        return;
    }
    color[0] = colors.at(i).r;
    color[1]= colors.at(i).g;
    color[2] = colors.at(i).b;

}

void TransferFunction::ReMap(vector<pair<float, float> > ranges, vector<QColor> ccolors){

    nodeValues.clear();
    colors.clear();
    // isBanding = true;
    for(int i = 0; i < ranges.size(); i++){

        nodeValues.push_back(ranges.at(i).first);
        nodeValues.push_back(ranges.at(i).second);

        auto current = ccolors.at(i);
         ColorRGB c;
         c.r = current.red()/255.0f;
         c.g = current.green()/255.0f;
         c.b = current.blue()/255.0f;

         colors.push_back(c);
         colors.push_back(c);
    }

}

void TransferFunction::LoadFromData(vector<std::string> data){

    int bands = QString::fromStdString(data.at(0)).toInt();
    isBanding = (bands == 1);
    // now the node values

    nodeValues.clear();
    vector<string> nodes;
    boost::split(nodes,data.at(1), boost::is_any_of(","));

    std::cout << "Size? " << nodes.size() << std::endl;
    for(unsigned int k = 0; k < nodes.size(); k++)
        nodeValues.push_back(  QString::fromStdString(nodes.at(k)).toFloat() );

    if (data.at(1).compare("") == 0 )
    {
        ClearTF();
    }
    else {
    // now the rgb

    colors.clear();
    nodes.clear();
    boost::split(nodes,data.at(2), boost::is_any_of(","));

    for(unsigned int k = 0; k < nodes.size(); k+= 3){
        ColorRGB c;
        c.r = QString::fromStdString(nodes.at(k)).toFloat();
        c.g = QString::fromStdString(nodes.at(k+1)).toFloat();
        c.b = QString::fromStdString(nodes.at(k+2)).toFloat();
        std::cout << "New color " << c.r << " , " << c.g << " , " << c.b <<std::endl;
        colors.push_back(c);
        }
    }
}
std::string TransferFunction::GetAsDataStream(){

    // ; separated.. is banding... node values and colors rgb

    std::string res = "";
    if ( isBanding) res = "1;";
    else res = "0;";

    //now the node values
    for(int i = 0; i < nodeValues.size(); i++){

        res += QString::number(nodeValues.at(i)).toStdString();

        if ( i == nodeValues.size() -1)
            res += ";";
        else
            res += ",";
    }

    // now the rgb

    for(int i = 0; i < colors.size(); i++){

        ColorRGB current = colors.at(i);
        res += QString::number(current.r).toStdString()  + ",";
        res += QString::number(current.g).toStdString()  + ",";
        res += QString::number(current.b).toStdString();

        if ( i == colors.size() -1)
            res += ";";
        else
            res += ",";
    }

   return res;
}

void TransferFunction::ChangeNodeLocation(int i, float newLoc){

    if (isBanding){

    }
    else {

       int size = nodeValues.size();

       if (i == 0) return;
       if (i == size -1 ) return;


       vector<float> tmpValues;

       for(int j = 0; j < size; j++){

           if (i ==j){ //i + 1 at most will be j

               if (newLoc <= nodeValues.at(i-1) ) {  newLoc = nodeValues.at(i-1) + nodeValues.at(i-1)*0.005; }
               if (newLoc >= nodeValues.at(i+1) ) {  newLoc = nodeValues.at(i+1) - nodeValues.at(i+1)*0.005;}
               tmpValues.push_back(newLoc);
           }
           else
               tmpValues.push_back(nodeValues.at(j));
       }

       nodeValues.clear();
       for(unsigned int j = 0; j < tmpValues.size(); j++){
           nodeValues.push_back(tmpValues.at(j));
       }
    }
}

void TransferFunction::ChangeNodeColor(int i, float color[]){
    if (isBanding){
         // *******************
        vector<ColorRGB> tmpColors;

        for(int j = 0; j < GetTotalNodes(); j++){

            if ( i == j){
                ColorRGB c;
                c.r = color[0];
                c.g = color[1];
                c.b = color[2];
                tmpColors.push_back(c);
                tmpColors.push_back(c);

            }
            else {
                tmpColors.push_back(colors.at(2*j));
                tmpColors.push_back(colors.at(2*j+1));

            }
        }

        colors.clear();
        for(unsigned int j = 0; j < tmpColors.size(); j++){
            colors.push_back(tmpColors.at(j));
        }

    }
    else {
        vector<ColorRGB> tmpColors;
        for(int j = 0; j < (int)nodeValues.size(); j++){

            if (i ==j){
                ColorRGB c;
                c.r = color[0];
                c.g = color[1];
                c.b = color[2];

                tmpColors.push_back(c);
            }
            else
                tmpColors.push_back(colors.at(j));
        }
        colors.clear();
        for(unsigned int j = 0; j < tmpColors.size(); j++){
            colors.push_back(tmpColors.at(j));
        }
    }
}

float TransferFunction::GetNodeInnerLoc(int i){

    float min = GetMinNodeVal();
    float max = GetMaxNodeVal();

    float v = nodeValues.at(i);

    float r =  (v - min)/(max - min);
    return r;
}


int TransferFunction::GetClass(float value){
   int cls = -1;
   if (nodeValues.empty())
   {
       return cls;
   }

   if ( GetNumberOfClases() > 1){

       int n = GetTotalNodes();
       for(int i = 0; i < n;  i+=1){

           if (nodeValues.at(2*i) <= value && value <= nodeValues.at(2*i+1) ){
               return i;
           }
       }
   }

   return cls;

}

void TransferFunction::GetColor(float value, float color[]){
    if (nodeValues.empty())
    {
        for(int k = 0; k < 3; k++) color[k] = 0;
        return;
    }

    const bool currentlyBanding = isBanding;



    if ( currentlyBanding){

        int n = GetTotalNodes();
        for(int i = 0; i < n;  i+=1){

            if (nodeValues.at(2*i) <= value && value <= nodeValues.at(2*i+1) ){

                ColorRGB c = colors.at(2*i);
                color[0] = c.r;
                color[1] = c.g;
                color[2] = c.b;
                return;
            }
        }
    }
    else {
        for(unsigned int i = 0; i< nodeValues.size()-1; i++){
            if ( nodeValues.at(i) <= value && value <= nodeValues.at(i+1) ){
               // Linearly interpolate
               float r = nodeValues.at(i+1) -nodeValues.at(i);
               if ( r < 0.00001 ) continue;
               //---------------
               float s = value - nodeValues.at(i);
               float alpha = s / r;

               color[0] = colors.at(i).r*(1.0 -alpha)  + colors.at(i+1).r*alpha;
               color[1]= colors.at(i).g*(1.0 -alpha)  + colors.at(i+1).g*alpha;
               color[2] = colors.at(i).b*(1.0 -alpha)  + colors.at(i+1).b*alpha;
               return;
            }
        }
    }


}

void TransferFunction::ClearTF(){
   colors.clear();
   nodeValues.clear();

}

void TransferFunction::SetAsBanded(){
  isBanding = true;
  ClearTF();
}

void TransferFunction::SetAsPiecewise(){
  isBanding = false;
  ClearTF();
}


void TransferFunction::RemoveNodeLoc(int loc){

    int size = nodeValues.size();
    if ( size < 3) return;
    if (isBanding) return;

    if (loc == 0) return;
    if (loc ==  size - 1) return;
    // If its banding, we cannot remove nodes
    vector<ColorRGB> tmpColors;
    vector<float> tmpNodeValues;

    for(int j = 0; j < size; j++){

        if (loc ==j){

        }
        else{
            tmpColors.push_back(colors.at(j));
            tmpNodeValues.push_back(nodeValues.at(j));
        }
    }
    nodeValues.clear();
    colors.clear();
    for(unsigned int j = 0; j < tmpColors.size(); j++){
        colors.push_back(tmpColors.at(j));
        nodeValues.push_back(tmpNodeValues.at(j));
    }

}

int TransferFunction::GetNumberOfClases(){

    int classes = 0;
    // Depends
    if (isBanding){
         classes = GetTotalNodes();
    }
    else { // on the distance b/w nodes
         float max = GetMaxNodeVal();
         float min = GetMinNodeVal();

         // epsilon
         double epsilon = 0.03;

         if ( !nodeValues.empty()){
             for(unsigned int i = 0; i < nodeValues.size()-1; i++){
                  float current = nodeValues.at(i);
                  float next = nodeValues.at(i+1);
                  float cP = (current - min)/(max - min);
                  float nP = (next - min)/(max - min);


                  if ( fabs(nP-cP) < epsilon )
                      continue;
                  classes++;
             }
         }
    }
    return classes;
}


void TransferFunction::GetClassColor(int idx, float color[]){
    color[0] = 0;
    color[1] = 0;
    color[2] = 0;

    float max = GetMaxNodeVal();
    float min = GetMinNodeVal();

    // epsilon
    double epsilon = 0.03;
    int classes = 0;

    for(unsigned int i = 0; i < nodeValues.size()-1; i++){
         float current = nodeValues.at(i);
         float next = nodeValues.at(i+1);
         float cP = (current - min)/(max - min);
         float nP = (next - min)/(max - min);


         if ( fabs(nP-cP) < epsilon )
             continue;

         if ( classes == idx){
              auto c1 = colors.at(i);
              auto c2 = colors.at(i+1);

              color[0] = (c1.r +c2.r)/2.0;
              color[1] = (c1.g +c2.g)/2.0;
              color[2] = (c1.b +c2.b)/2.0;
              break;
         }
         classes++;
    }

}

pair<float, float> TransferFunction::GetClassRange(int idx){

    float max = GetMaxNodeVal();
    float min = GetMinNodeVal();

    auto currentRange = make_pair(min, max);

    // epsilon
    double epsilon = 0.03;
    int classes = 0;

    for(unsigned int i = 0; i < nodeValues.size()-1; i++){
         float current = nodeValues.at(i);
         float next = nodeValues.at(i+1);
         float cP = (current - min)/(max - min);
         float nP = (next - min)/(max - min);


         if ( fabs(nP-cP) < epsilon )
             continue;

         if ( classes == idx){
              currentRange = make_pair(current, next);
              break;
         }
         classes++;
    }

    return currentRange;
}

void TransferFunction::AddNode(float value, float r, float g, float b){
   // value is between 0 and 1

    float max = GetMaxNodeVal();
    float min = GetMinNodeVal();

    float desired = (max-min)*value + min;

    vector<float> tmpValues;
    vector<ColorRGB> tmpColors;


    for(unsigned int i = 0; i < nodeValues.size(); i++){
        tmpValues.push_back(nodeValues.at(i));
        tmpColors.push_back(colors.at(i));

        if ( i < nodeValues.size() -1)
            if ( nodeValues.at(i) < desired && desired < nodeValues.at(i+1)  ){
                tmpValues.push_back(desired);
                ColorRGB c;
                c.r = r;             c.g = g;
                c.b = b;
                tmpColors.push_back(c);
            }
    }

    nodeValues.clear();
    colors.clear();

    for(unsigned int i = 0; i < tmpValues.size(); i++){
        nodeValues.push_back(tmpValues.at(i));
        colors.push_back(tmpColors.at(i));
    }
}

void TransferFunction::CreateDefaultBandedMap(int qty, bool ordinal){
     // Base on the index

    isBanding = true;
    ClearTF();

    if(!deuteranopiaColor){ // Assumes that not red-green color friendly
         // is the normal, not handling the other two types of color

        float color[3];
        for(int i = 0; i < qty; i++){

            // if ordinal
            float v1 = i - 0.4;
            float v2 = i + 0.4;

            nodeValues.push_back(v1);
            nodeValues.push_back(v2);

            ColorRGB c;
            if ( !ordinal) // Categorical... lower the saturation
                GetHSVColor(i, qty, 0.6, color);
            else{
                GetHSVColor(0, 1, (static_cast<float>(i)/static_cast<float>(qty))*0.5 + 0.5, color);
            }
            c.r = color[0];
            c.g = color[1];
            c.b = color[2];
            colors.push_back(c);
            colors.push_back(c);
        }
    }
    else {
        // need to create qty different colors

        std::ifstream file;
        file.open("redgreen.txt");

        if (!file.is_open()){
            std::cerr << "Error opening the file in the specified directory" << std::endl;
            exit(1);
        }


        std::cout << "Going to generate the colors " << std::endl;
        int colorsGenerated = 0;

        vector<QColor> original;
        vector<QColor> visible;

        std::string line;
        while(std::getline(file, line)){
             std::vector<std::string> strings;
             boost::split(strings, line, boost::is_any_of(" "));
             //0 is deuteranpoia label
             //1-2-3 original rgb
             //4-5-6 what is visible as
             //so first we need to
             QColor o;
             o.setRed( QString::fromStdString(strings.at(1)).toInt());
             o.setGreen(QString::fromStdString(strings.at(2)).toInt());
             o.setBlue( QString::fromStdString(strings.at(3)).toInt());

             QColor v;
             v.setRed( QString::fromStdString(strings.at(4)).toInt());
             v.setGreen(QString::fromStdString(strings.at(5)).toInt());
             v.setBlue( QString::fromStdString(strings.at(6)).toInt());


             if ( o.red() == o.green() && o.green() == o.blue()) continue;

             float maxC = std::max(std::max( o.red(), o.green()), o.blue());
             float minC = std::min(std::min( o.red(), o.green()), o.blue());

             float dif = maxC - minC;
             float saturation = dif / maxC;

             if (saturation < 0.5) continue;
             if (saturation > 0.8) continue;

             float value = maxC / 255.0f;

             if ( value < 0.45 ) continue;

             original.push_back(o); visible.push_back(v);

        }
        // Now
        float color[3];
        GetHSVColor(2,3, 0.6, color);

        vector<QColor> insertedColors;

        ColorRGB c; QColor qc;


        float v1 = colorsGenerated - 0.4;
        float v2 = colorsGenerated + 0.4;

        nodeValues.push_back(v1);
        nodeValues.push_back(v2);



        c.r = color[0];
        c.g = color[1];
        c.b = color[2];
        qc.setRgb(color[0]*255, color[1]*255, color[2]*255);
        insertedColors.push_back(qc);
        colors.push_back(c);
        colors.push_back(c);
        colorsGenerated++;

        vector<QColor> visibilityOfInsertedColors;
        int idx = -1;
        float minDistance = 999999999999999;//? there must be a better way for this...
        for(int i = 0; i < original.size(); i++){
            // find the closest color index....

           float distance = fabs(original.at(i).red() - qc.red()) +  fabs(original.at(i).green() - qc.green())  + fabs(original.at(i).blue() - qc.blue()) ;

           if ( distance < minDistance ){
               minDistance = distance;
               idx = i;
           }
        }

        std::cout << "Added? " << qc.red() << ", " << qc.green() << " ," << qc.black() << std::endl;


        visibilityOfInsertedColors.push_back(visible.at(idx));

        while( colorsGenerated < qty){
           // We have one color generated, so we need to generate new colors
           // such as the distance to the ones inserted works...
            float v1 = colorsGenerated - 0.4;
            float v2 = colorsGenerated + 0.4;

            nodeValues.push_back(v1);
            nodeValues.push_back(v2);
            QColor qc = GetNewColorblindColor(&original, &visible, &visibilityOfInsertedColors);
            ColorRGB c;

            std::cout << "Adding " << qc.red() << ", " << qc.green() << " ," << qc.blue() << std::endl;

            c.r = qc.red()/255.0f;  c.g = qc.green()/255.0f; c.b = qc.blue()/255.0f;
            colors.push_back(c);
            colors.push_back(c);

            std::cout << "Added? " << c.r << ", " << c.g << " ," << c.b << std::endl;
            colorsGenerated++;
        }
    }
}

QColor TransferFunction::GetNewColorblindColor(vector<QColor>* originalColors, vector<QColor>* visibleColors, vector<QColor>* sofarInserted){


    // find the color that has the max distance to the last two, inserted
    // if its  the last one, its just going to toggle in between
    int idx = -1;
    float maxDistance = 0;

    int to = 2;
    if ( to > sofarInserted->size())
        to = 1;

    std::cout << "Comparing agains? " << to << std::endl;

    for(int i = 0; i < visibleColors->size(); i++){
        QColor c = visibleColors->at(i);

        float d = 0;
        for(int j = 0; j < sofarInserted->size(); j++){
             QColor o = sofarInserted->at(j);//sofarInserted->at( sofarInserted->size() - j -1);
             float cd = fabs(c.red() - o.red()) + fabs(c.green() - o.green()) + fabs(c.blue() - o.blue());

             if ( cd < 110 ) d -= 99999999999;
             d+= cd;
        }

        if ( d > maxDistance ){
            maxDistance  = d;
            idx = i;
        }
    }


    QColor c = originalColors->at(idx);
    sofarInserted->push_back(visibleColors->at(idx));
    std::cout << "num colors " << visibleColors->size() << std::endl;
    std::cout << idx << " , " << maxDistance << " :: " << c.red() << "," << c.green() << "," << c.blue() << std::endl;

    return c;
}

int TransferFunction::GetTotalNodes(){

    if ( isBanding)
        return nodeValues.size()/2;
    else
        return nodeValues.size();
}

float TransferFunction::GetMaxNodeVal(){
   if (nodeValues.empty()) return -9999999;

   return nodeValues.at(nodeValues.size()-1);
}
float TransferFunction::GetMinNodeVal(){
    if (nodeValues.empty()) return -9999999;
    return nodeValues.at(0);
}

void  TransferFunction::CreateDefaultColorMap(float max, float min){
    ClearTF();
    isBanding = false;

    nodeValues.push_back(min);
    nodeValues.push_back(max);

    std::cout << "Max? "<<  max << " , min?"  << min << std::endl;
    ColorRGB c, c2;
    c.r = 0;   c.g = 0;  c.b = 0;
    colors.push_back(c);
    c2.r = 1.0;   c2.g = 0;  c2.b = 0;
    colors.push_back(c2);

}

void TransferFunction::CreateDefaultDensityMap(){

    ClearTF();

    isBanding = false;
    float colorsDensity[6][3] = {{0.95, 0.95, 0.95}, {0.95, 0.95, 0.95},{0.55, 0.55, 0},
                                 {1.0, 1.0, 0}, {1.0,0.62, 0 },{1.0, 0,0}};

    nodeValues.push_back(0);
    nodeValues.push_back(0.001);
    nodeValues.push_back(0.001);
    nodeValues.push_back(1.0f/3.0f);
    nodeValues.push_back(2.0f/3.0f);
    nodeValues.push_back(1.0f);

    for(int i = 0; i < 6; i++){
        ColorRGB c;
        c.r = colorsDensity[i][0];
        c.g = colorsDensity[i][1];
        c.b = colorsDensity[i][2];
        colors.push_back(c);
    }
}

void TransferFunction::GetHSVColor(int index, int total, float s, float color[]){

    // the idea is to change from hsv to rgb...
    float angleValue =  index * ( (360.0 )/ static_cast<float>(total ));
    //float s = 1.0;
    float v = 0.8;

    int Hi =  static_cast<int>( (angleValue / 60) ) % 6;
    // f is how far it is in the location
    float f =  (angleValue - Hi*60) / 60.0;

    float p = v*(1.0 - s);
    float q = v*(1.0 - f*s);
    float t = v*(1.0 -  (1.0-f)*s);


    switch(Hi){
       case  0:
           color[0] = v;  color[1] = t; color[2] = p;
           break;
        case  1:
            color[0] = q;  color[1] = v; color[2] = p;
            break;
        case  2:
            color[0] = p;  color[1] = v; color[2] = t;
            break;
        case  3:
            color[0] = p;  color[1] = q; color[2] = v;
            break;
        case  4:
            color[0] = t;  color[1] = p; color[2] = v;
            break;
        case  5:
            color[0] = v;  color[1] = p; color[2] = q;
            break;
    }
}
