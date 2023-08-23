#ifndef TRANSFERFUNCTION_H
#define TRANSFERFUNCTION_H

// Transfer function for
// Bandmap or Piecewise Transfer functions
// Default Piecewise Transfer functions

#include <vector>
#include <iostream>
#include <math.h>
#include <QString>
#include <QColor>

#include <boost/algorithm/string.hpp>
#include <boost/regex.hpp>

using namespace std;
class ColorRGB {
    public:
       float r;
       float g;
       float b;
};

class TransferFunction
{
public:
    TransferFunction();

    void GetColor(float value, float* r, float* g, float *b);
    void GetColor(float value, float color[]);
    int GetClass(float value);


    void GetNodeColor(int i, float color[]);
    void GetClassColor(int idx, float color[]);

    void ChangeNodeLocation(int i, float newLoc);
    void ChangeNodeColor(int i, float color[]);

    void RemoveNodeLoc(int loc);

    float GetNodeInnerLoc(int i);


    void EnableRedGreenColorBlindness(){ deuteranopiaColor = true;}
    void DisableRedGreenColorBlindness(){ deuteranopiaColor = false;}
    bool IsColorBlindFriendly(){ return deuteranopiaColor; }
    int GetTotalNodes();
    float GetMaxNodeVal();
    float GetMinNodeVal();

    void ClearTF();
    void AddNode(float value, float r, float g, float b);
    bool IsBandingColorMap(){ return isBanding;}
    bool IsPieceWieseLinear(){ return !isBanding; }

    int GetNumberOfClases();
    pair<float, float> GetClassRange(int idx);


    void SetAsBanded();
    void SetAsPiecewise();
    void CreateDefaultBandedMap(int qty, bool ordinal); // Base on the index
    void CreateDefaultDensityMap();
    void CreateDefaultColorMap(float max, float min);

    void LoadFromData(vector<std::string> data);
    std::string GetAsDataStream();
    static void GetHSVColor(int index, int total, float s, float color[]);

    void ReMap(vector<pair<float, float> > ranges, vector<QColor> colors);


    vector<float> nodeValues; // In the range from 0...1 where are they located
    // For the bandmap... repeated values are possible, also for the regular

    vector<ColorRGB> colors;
private:

    QColor GetNewColorblindColor(vector<QColor>* originalColors, vector<QColor>* visibleColors, vector<QColor>* sofarInserted);
    bool isBanding;
    bool deuteranopiaColor; // apply color blindness, for default colormapping
};

#endif // TRANSFERFUNCTION_H
