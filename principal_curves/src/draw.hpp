#ifndef DRAW_HPP
#define DRAW_HPP

#include <vector>
#include <functional>
#include <cairo.h>

#include "curve2d.hpp"
#include "image2d.hpp"
#include "dataset.hpp"
#include "pipeline.hpp"

class DrawContext {
    cairo_surface_t* _surface;
    cairo_t* cr;
    int width, height;
    int baseWidth, baseHeight;
    float colorR = 1, colorG = 0, colorB = 0;

    const float colorMin = 0.1;

    public:

    DrawContext(int bw, int bh, int w, int h);
    DrawContext(int bw, int bh, int w, int h, int hplus);
    ~DrawContext();

    void setColor(float r, float g, float b);

    void drawPlotBackground(const Image2D<float>& background);
    void drawSkeleton(const Image2D<char>& background);
    void drawPlot(const Image2D<unsigned int>& background);
    void drawVoronoiRegions(const Image2D<int>& labels, const Image2D<unsigned> &plot, int numRegions);
    void SaveBitmapTemplate(const Image2D<unsigned int>& plot, const char* filename);
    void SaveBitmapTemplate(const Image2D<char>& plot, const char* filename);


    void drawBand(const std::vector<AugmentedPoint>& curve);
    void drawGraphBand( std::vector<SkeletonNode*>* nodes);

    void drawMainLine(const std::vector<AugmentedPoint>& curve);
    void drawMainLine(const std::vector<Int2>& curve);

    void drawDashLines(const std::vector<AugmentedPoint>& curve);
    void drawControlPoints(const std::vector<AugmentedPoint>& curve, float rad = 1.0);
    void drawControlPoints(const std::vector<Int2>& curve, float rad =1.0);

    void drawLine(const Int2 p1, const Int2 p2);
    void drawVoronoiGraph(PipelineResult* result);
    void drawVoronoiGraph( std::vector<SkeletonNode*>* nodes);

    void drawInfo(const std::vector<AugmentedPoint>& curve);

    void writeToFile(const std::string& fileName) const;
    cairo_surface_t* surface() const;
};

cairo_surface_t* drawDots(const std::vector<Float2> dots, int width, int height, int highlight, std::function<bool(int)> mark);

#endif
