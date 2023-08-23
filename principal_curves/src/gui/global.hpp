#ifndef GUI_GLOBAL_HPP
#define GUI_GLOBAL_HPP

#include <gtkmm.h>

#include "../dataset.hpp"
#include "../multival.hpp"
#include "../randomprojection.hpp"
#include "../pipeline.hpp"
#include "../projectionlist.hpp"
#include "imagecontroller.hpp"

struct GuiGlobal {
    DataSet dataSet { 0 };
    ImageController mainImageController;
    Gtk::FileChooserButton *dataFileChooserButton;
    Gtk::Switch *calculationSwitch;

    Gtk::CheckButton *checkDrawPlot;
    Gtk::CheckButton *checkDrawGauss;
    Gtk::CheckButton *checkDrawSkeleton;
    Gtk::CheckButton *checkDrawDash;
    Gtk::CheckButton *checkDrawBand;
    Gtk::CheckButton *checkDrawCurve;
    Gtk::CheckButton *checkDrawStartingCurve;
    Gtk::CheckButton *checkDrawControlPoints;

    Gtk::Grid* projectionMatrixGrid;
    Gtk::Button* projectionMatrixResetButton;
    Gtk::Button* projectionMatrixRandomButton;
    std::vector<Glib::RefPtr<Gtk::Adjustment>> projectionMatrixX;
    std::vector<Glib::RefPtr<Gtk::Adjustment>> projectionMatrixY;

    Gtk::Label* projectionListCountLabel;
    Gtk::Window* projectionListWindow;

    ProjectionGenerator projectionMatrixGenerator;
    ProjectionList projectionList;
    int projectionListCardinalCount = 0;
    ImageController projectionListImageController;

    Glib::RefPtr<Gtk::Adjustment> adjustmentSmooth;
    Glib::RefPtr<Gtk::Adjustment> adjustmentWidth;
    Glib::RefPtr<Gtk::Adjustment> adjustmentSegmentPix;
    Glib::RefPtr<Gtk::Adjustment> adjustmentVarianceSmooth;

    Glib::RefPtr<Gtk::Adjustment> adjustmentVarianceFactor;
    Glib::RefPtr<Gtk::Adjustment> adjustmentAmountFactor;
    Gtk::CheckButton *checkMirror;


    SkeletonGenerator* generator;
    PipelineParameters parameters;
    PipelineResult result;
    std::vector<Float2> dots;
    int dotsHighlight = -1;
};

extern GuiGlobal global;

#endif
