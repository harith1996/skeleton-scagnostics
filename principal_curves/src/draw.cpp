#include "draw.hpp"
#include <iostream>

DrawContext::DrawContext(int bw, int bh, int w, int h)
: DrawContext(bw, bh, w, h, 0) {
}

DrawContext::DrawContext(int bw, int bh, int w, int h, int hplus)
: width {w}, height {h}, baseWidth{bw}, baseHeight{bh} {
    _surface = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, width, height + hplus);
    cr = cairo_create(_surface);

    cairo_set_source_rgb(cr, colorMin, colorMin, colorMin);
    cairo_paint(cr);

    cairo_scale(cr, width * 1.0 / baseWidth, height * 1.0 / baseHeight);
}


DrawContext::~DrawContext() {
    cairo_surface_destroy(_surface);
    cairo_destroy(cr);
}

void DrawContext::setColor(float r, float g, float b) {
    colorR = r;
    colorG = g;
    colorB = b;
}

void DrawContext::drawPlotBackground(const Image2D<float>& background) {
    for(int x = 0; x < background.width(); x++) {
        for(int y = 0; y < background.height(); y++) {
            auto val = std::max(0.0f, std::min(background(x, y), 1.0f));
            cairo_set_source_rgba(cr, colorR, colorG, colorB, val);
            cairo_rectangle(cr, x, y, 1, 1);
            cairo_fill(cr);
        }
    }
}

void DrawContext::drawVoronoiRegions(const Image2D<int>& labels, const Image2D<unsigned>& plot, int numRegions){



    for(int x = 0; x < labels.width(); x++) {
        for(int y = 0; y < labels.height(); y++) {
            auto val = labels(x,y);
            auto dense = plot(x,y);

            float opacity = 1.0;

            if (dense==0)
                opacity = 0;

            //if ( dense != 0)
            //    opacity *=  static_cast<float>(val)/static_cast<float>(numRegions);
            if ( val % 3 == 0)
               cairo_set_source_rgba(cr,  1.0, 0, 0, opacity);
            else if (val % 3 == 1){
                cairo_set_source_rgba(cr, 0, 1.0, 0, opacity);
            }
            else if (val % 3 == 2){
                cairo_set_source_rgba(cr, 0, 0, 1.0, opacity);

            }

            cairo_rectangle(cr, x, y, 1, 1);
            cairo_fill(cr);
        }
    }


}


void  DrawContext::drawPlot(const Image2D<unsigned int>& background){
    for(unsigned int x = 0; x < background.width(); x++) {
        for(unsigned int y = 0; y < background.height(); y++) {
            float val = 0.0;
            if (background(x, y))
                 val = 1.0;

            cairo_set_source_rgba(cr, colorR, colorG, colorB, val);
            cairo_rectangle(cr, x, y, 1, 1);
            cairo_fill(cr);
        }
    }
}

void DrawContext::SaveBitmapTemplate(const Image2D<unsigned int>& plot, const char* filename){
    // First line is the width and height of template
    stringstream ss;
    ss << plot.width() << " " << plot.height() << "\n";
    for(unsigned int y = 0; y < plot.height(); y++) {
        for(unsigned int x = 0; x < plot.width(); x++) {
            int val = 0;
            if (plot(x, y))
                val = 1;
            ss << val;
            if ( x < plot.width() -1 )
                ss << " ";

        }
        if ( y < plot.height() -1);
            ss << "\n";
    }

    std::ofstream out(filename);
    out << ss.str();
    out.close();
}


void  DrawContext::SaveBitmapTemplate(const Image2D<char>& plot, const char* filename){
    // First line is the width and height of template
    stringstream ss;
    ss << plot.width() << " " << plot.height() << "\n";
    for(unsigned int y = 0; y < plot.height(); y++) {
        for(unsigned int x = 0; x < plot.width(); x++) {
            int val = 0;
            //if (plot(x, y) )
            //    val = 1;*/
            if ( plot(x,y) == 1 )
                val = 1;
            ss << val;
            if ( x < plot.width() -1 )
                ss << " ";

        }
        if ( y < plot.height() -1);
            ss << "\n";
    }

    std::ofstream out(filename);
    out << ss.str();
    out.close();
}
void DrawContext::drawSkeleton(const Image2D<char>& background) {
    int tot = 0;
    for(unsigned int x = 0; x < background.width(); x++) {
        for(unsigned int y = 0; y < background.height(); y++) {
            float val = 0.0;
            if (background(x, y))
                 val = 1.0;

            if ( val > 0.5) tot++;

            cairo_set_source_rgba(cr, colorR, colorG, colorB, val);
            cairo_rectangle(cr, x, y, 1, 1);
            cairo_fill(cr);
        }
    }
  //  std::cout <<"Total Skeleton points (Drawing) " << tot << "/" << background.width()*background.height() <<  std::endl;
}



void DrawContext::drawBand(const std::vector<AugmentedPoint>& curve) {
    cairo_save(cr);

    //offset to move curve points into center of "pixel"
    float subpixelX = 0.5;
    float subpixelY = 0.5;
    cairo_translate(cr, subpixelX, subpixelY);

    const float sigmaFactor = 1.2;
    cairo_pattern_t* pattern =  cairo_pattern_create_mesh();
    float curveLength = 0;

    for(int i = 0; i < curve.size(); i++) {
        curveLength += curve[i].amount / curve[i].amountPerLength;
    }

    for(int i = 1; i < curve.size(); i++) {
        Float2 vectorAB {  curve[i-1].direction[1], -curve[i-1].direction[0] };
        Float2 vectorCD {  curve[i].direction[1], -curve[i].direction[0] };

        float variance1 = curve[i-1].orthoVariance;
        float variance2 = curve[i].orthoVariance;
        float sigma1 = sqrt(std::max(0.0f, variance1));
        float sigma2 = sqrt(std::max(0.0f, variance2));

        Float2 pointA = curve[i-1].point - sigmaFactor * sigma1 * vectorAB;
        Float2 pointB = curve[i-1].point + sigmaFactor * sigma1 * vectorAB;
        Float2 pointC = curve[i].point + sigmaFactor * sigma2 * vectorCD;
        Float2 pointD = curve[i].point - sigmaFactor * sigma2 * vectorCD;

        cairo_mesh_pattern_begin_patch(pattern);
        cairo_mesh_pattern_move_to(pattern, pointA[0], pointA[1]);
        cairo_mesh_pattern_line_to(pattern, pointB[0], pointB[1]);
        cairo_mesh_pattern_line_to(pattern, pointC[0], pointC[1]);
        cairo_mesh_pattern_line_to(pattern, pointD[0], pointD[1]);

        const float lr = 0.92, lf = 50;
        float l1 = std::atan(lf * (curve[i-1].amountPerLength - 1 / curveLength)) / (0.5 * M_PI) + 0.5;
        float l2 = std::atan(lf * (curve[i  ].amountPerLength - 1 / curveLength)) / (0.5 * M_PI) + 0.5;

        /*
        l1 /= sigma1;
        l2 /= sigma2;

        l1 = 1 - (1 - l1)/1.5;
        l2 = 1 - (1 - l2)/1.5;
        */

        l1 = 1 - lr + lr * l1;
        l2 = 1 - lr + lr * l2;

        cairo_mesh_pattern_set_corner_color_rgba(pattern, 0, colorR, colorG, colorB, l1);
        cairo_mesh_pattern_set_corner_color_rgba(pattern, 1, colorR, colorG, colorB, l1);
        cairo_mesh_pattern_set_corner_color_rgba(pattern, 2, colorR, colorG, colorB, l2);
        cairo_mesh_pattern_set_corner_color_rgba(pattern, 3, colorR, colorG, colorB, l2);
        cairo_mesh_pattern_end_patch(pattern);
    }
    for(int i = 0; i < curve.size(); i++) {
        Float2 vectorCD {  curve[i].direction[1], -curve[i].direction[0] };
        float variance2 = curve[i].orthoVariance;
        Float2 pointC = curve[i].point + sigmaFactor * sqrt(std::max(0.0f, variance2)) * vectorCD;
        cairo_line_to(cr, pointC[0], pointC[1]);

    }
    for(int i = curve.size() - 1; i >= 0; i--) {
        Float2 vectorCD {  curve[i].direction[1], -curve[i].direction[0] };
        float variance2 = curve[i].orthoVariance;
        Float2 pointD = curve[i].point - sigmaFactor * sqrt(std::max(0.0f, variance2)) * vectorCD;
        cairo_line_to(cr, pointD[0], pointD[1]);
    }
    cairo_set_source(cr, pattern);
    cairo_close_path(cr);
    cairo_fill(cr);

    cairo_restore(cr);
}

void DrawContext::drawMainLine(const std::vector<AugmentedPoint>& curve) {
    cairo_save(cr);

    //offset to move curve points into center of "pixel"
    float subpixelX = 0.5;
    float subpixelY = 0.5;
    cairo_translate(cr, subpixelX, subpixelY);

    cairo_set_source_rgb(cr, colorR, colorG, colorB);

    double line_width = 2.0, dummy = 0;
    cairo_device_to_user_distance(cr, &dummy, &line_width);
    cairo_set_line_width(cr, line_width);

    if(curve.size() > 0) {
        cairo_move_to(cr, curve[0].point[0], curve[0].point[1]);
    }
    for(int i = 1; i < curve.size(); i++)
    {
        cairo_line_to(cr, curve[i].point[0], curve[i].point[1]);
    }
    cairo_stroke(cr);

    cairo_restore(cr);
}

void DrawContext::drawLine(const Int2 p1, const Int2 p2) {
    cairo_save(cr);

    //offset to move curve points into center of "pixel"
    float subpixelX = 0.5;
    float subpixelY = 0.5;
    cairo_translate(cr, subpixelX, subpixelY);

    cairo_set_source_rgb(cr, colorR, colorG, colorB);

    double line_width = 2.0, dummy = 0;
    cairo_device_to_user_distance(cr, &dummy, &line_width);
    cairo_set_line_width(cr, line_width);

    cairo_move_to(cr, p1[0],p1[1]);

    cairo_line_to(cr, p2[0],p2[1]);

    cairo_stroke(cr);

    cairo_restore(cr);
}

void DrawContext::drawVoronoiGraph(PipelineResult* result){

    std::vector<SkeletonNode*>* nodes = &(result->graphNodes);
    // basically draw bunch of lines ...
    for(int i = 0; i < nodes->size(); i++){
        SkeletonNode* currentNode = nodes->at(i);
        Int2 p1 = currentNode->location;

        for(int j = 0; j < currentNode->indexNeighborNode.size(); j++ ){
            int o = currentNode->indexNeighborNode.at(j);
            if ( o <  nodes->size()){
                Int2 p2 = nodes->at(o)->location;
                drawLine(p1,p2);
            }
        }
    }
}

void DrawContext::drawVoronoiGraph( std::vector<SkeletonNode*>* nodes){
    // basically draw bunch of lines ...


    for(int i = 0; i < nodes->size(); i++){
        SkeletonNode* currentNode = nodes->at(i);
        Int2 p1 = currentNode->location;

        for(int j = 0; j < currentNode->indexNeighborNode.size(); j++ ){
            int o = currentNode->indexNeighborNode.at(j);
            if ( o <  nodes->size()){
                Int2 p2 = nodes->at(o)->location;
                drawLine(p1,p2);
            }
        }
    }
}

void DrawContext::drawMainLine(const std::vector<Int2>& curve) {
    cairo_save(cr);

    //offset to move curve points into center of "pixel"
    float subpixelX = 0.5;
    float subpixelY = 0.5;
    cairo_translate(cr, subpixelX, subpixelY);

    cairo_set_source_rgb(cr, colorR, colorG, colorB);

    double line_width = 2.0, dummy = 0;
    cairo_device_to_user_distance(cr, &dummy, &line_width);
    cairo_set_line_width(cr, line_width);

    if(curve.size() > 0) {
        cairo_move_to(cr, curve[0][0], curve[0][1]);
    }
    for(int i = 1; i < curve.size(); i++)
    {
        cairo_line_to(cr, curve[i][0], curve[i][1]);
    }
    cairo_stroke(cr);

    cairo_restore(cr);
}


void DrawContext::drawDashLines(const std::vector<AugmentedPoint>& curve) {
    cairo_save(cr);

    //offset to move curve points into center of "pixel"
    float subpixelX = 0.5;
    float subpixelY = 0.5;
    cairo_translate(cr, subpixelX, subpixelY);

    double line_width = 1.2, dummy = 0;
    cairo_device_to_user_distance(cr, &dummy, &line_width);
    cairo_set_line_width(cr, line_width);

    cairo_set_source_rgb(cr, colorR, colorG, colorB);

    double dashes[] = {1.0, .5};
    //cairo_device_to_user_distance(cr, &dashes[0], &dashes[1]);
    const size_t numDashes = 2;
    cairo_set_dash(cr, dashes, numDashes, 0.0);

    const float sigmaFactor = 1.2;
    // For each point in the curve we create a line connecting all the points
    // in the + sign of the orthogonal direction vector
    for(int i = 0; i < curve.size(); i++) {
        Float2 vectorCD {  curve[i].direction[1], -curve[i].direction[0] };
        float variance2 = curve[i].orthoVariance;
        Float2 pointC = curve[i].point + sigmaFactor * sqrt(std::max(0.0f, variance2)) * vectorCD;

        //division by zero in empty sections
        //check for nan, isnan does not want to work
        if(pointC[0] == pointC[0] && pointC[1] == pointC[1]) {
            cairo_line_to(cr, pointC[0], pointC[1]);
        }
    }
    cairo_stroke(cr);

    //equally we now move in the negative direction ...
    for(int i = 0; i < curve.size(); i++) {
        Float2 vectorCD {  curve[i].direction[1], -curve[i].direction[0] };
        float variance2 = curve[i].orthoVariance;
        Float2 pointD = curve[i].point - sigmaFactor * sqrt(std::max(0.0f, variance2)) * vectorCD;

        //division by zero in empty sections
        //check for nan, isnan does not want to work
        if(pointD[0] == pointD[0] && pointD[1] == pointD[1]) {
            cairo_line_to(cr, pointD[0], pointD[1]);
        }
    }
    cairo_stroke(cr);

    double noDashes[] = {1.0};
    cairo_set_dash(cr, noDashes, 1, 0.0);

    cairo_restore(cr);
}

void DrawContext::drawControlPoints(const std::vector<AugmentedPoint>& curve, float rad) {
    cairo_save(cr);

    //offset to move curve points into center of "pixel"
    float subpixelX = 0.5;
    float subpixelY = 0.5;
    cairo_translate(cr, subpixelX, subpixelY);

    double radius = rad, dummy = 0;
    cairo_device_to_user_distance(cr, &dummy, &radius);

    cairo_set_source_rgb(cr, colorR, colorG, colorB);

    for(int i = 0; i < curve.size(); i++) {
        auto p = curve[i].point;
        cairo_arc(cr, p[0], p[1], radius, 0, 2*M_PI);
        cairo_fill(cr);
    }
}

void DrawContext::drawControlPoints(const std::vector<Int2>& curve, float rad) {
    cairo_save(cr);

    //offset to move curve points into center of "pixel"
    float subpixelX = 0.5;
    float subpixelY = 0.5;
    cairo_translate(cr, subpixelX, subpixelY);

    double radius = rad, dummy = 0;
    cairo_device_to_user_distance(cr, &dummy, &radius);

    cairo_set_source_rgb(cr, colorR, colorG, colorB);

    for(int i = 0; i < curve.size(); i++) {
        auto p = curve[i];
        cairo_arc(cr, p[0], p[1], radius, 0, 2*M_PI);
        cairo_fill(cr);
    }
}

void DrawContext::drawInfo(const std::vector<AugmentedPoint>& curve) {
    if(curve.size() <= 1) return;
    cairo_set_line_width(cr, baseHeight * 1.0 / 200);
    cairo_move_to(cr, 0, baseHeight);
    cairo_line_to(cr, baseWidth, baseHeight);
    cairo_set_source_rgb(cr, 1, 1, 1);
    cairo_stroke(cr);

    float curveLength = 0;

    for(int i = 0; i < curve.size(); i++) {
        curveLength += curve[i].amount / curve[i].amountPerLength;
    }

    for(int i = 0; i < 100; i++) {
        const float lr = 0.92, lf = 50;
        float a = (i - 50) / 100.0 * 1 / curveLength;
        float l1 = std::atan(lf * a) / (0.5 * M_PI) + 0.5;
        l1 = 1 - lr + lr * l1;
        cairo_rectangle(cr, i * baseWidth / 100.0, baseHeight + 1, baseWidth, baseHeight / 20.0);

        cairo_set_source_rgb(cr, colorMin, colorMin, colorMin);
        cairo_fill_preserve(cr);
        cairo_set_source_rgba(cr, colorR, colorG, colorB, l1);
        cairo_fill(cr);
    }

}


void DrawContext::writeToFile(const std::string& fileName) const {
    cairo_surface_write_to_png(_surface, fileName.c_str());
}

cairo_surface_t* DrawContext::surface() const {
    return _surface;
}

cairo_surface_t* drawDots(const std::vector<Float2> dots, int width, int height, int highlight, std::function<bool(int)> mark) {
    cairo_surface_t *surface;
    cairo_t *cr;

    surface = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, width, height);
    cr = cairo_create(surface);

    float colorMin = 0.1;
    cairo_set_source_rgb(cr, colorMin, colorMin, colorMin);
    cairo_paint(cr);

    cairo_scale(cr, width * 1.0, height * 1.0);

    double radius = 4, line_width = 1;
    cairo_device_to_user_distance(cr, &radius, &line_width);
    cairo_set_line_width(cr, line_width);

    for(int i = 0; i < dots.size(); i++) {
        cairo_arc(cr, dots[i][0], dots[i][1], radius, 0, 2*M_PI);
        auto m = mark(i);
        if(!m) {
            cairo_set_source_rgb(cr, 0.0, 0.0, 0.7);
        } else {
            cairo_set_source_rgb(cr, 0.0, 0.7, 0.0);
        }
        cairo_fill_preserve(cr);
        if(!m) {
            cairo_set_source_rgb(cr, 0.2, 0.2, 0.8);
        } else {
            cairo_set_source_rgb(cr, 0.2, 0.8, 0.2);
        }
        cairo_stroke(cr);
    }

    if(highlight >= 0 && highlight < dots.size()) {
        cairo_arc(cr, dots[highlight][0], dots[highlight][1], radius, 0, 2*M_PI);
        cairo_set_source_rgb(cr, 0.7, 0.0, 0.0);
        cairo_fill_preserve(cr);
        cairo_set_source_rgb(cr, 0.8, 0.2, 0.2);
        cairo_stroke(cr);
    }

    cairo_destroy(cr);

    return surface;
}

