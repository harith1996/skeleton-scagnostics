
#include "parsedata.hpp"
#include "draw.hpp"
#include "transform2d.hpp"
#include "principalgraph.hpp"
#include "gui/global.hpp"
#include "gui/projectionmatrix.hpp"
#include "gui/plist.hpp"

GuiGlobal global;

void renderImageOn(ImageController *imgc) {
    bool asGraph = global.parameters.asGraph;


    if(global.result.augmentedCurve.size() == 0 && !asGraph) return;
    if(global.result.controlPoints.size() == 0 && asGraph) return;

    auto zoom = imgc->getZoom();
    int imgWidth = int(global.parameters.width * zoom);
    int imgHeight = int(global.parameters.height * zoom);

    DrawContext drawContext { global.parameters.width, global.parameters.height, imgWidth, imgHeight, imgHeight / 10 };
    if(global.checkDrawPlot->get_active()) {
        drawContext.setColor(0, 0, 1);
        if(global.checkDrawGauss->get_active()) {
            drawContext.drawPlotBackground(normalizeImageValues(global.result.plotGauss));
        } else {
            drawContext.drawPlotBackground(normalizeImageValues(global.result.plot));
        }
    }
    if(global.checkDrawSkeleton->get_active()) {
        drawContext.setColor(0.8, 0.8, 0.8);
        drawContext.drawPlotBackground(normalizeImageValues(global.result.skeleton));
    }

    if (!asGraph){

        if(global.checkDrawDash->get_active()) {
            drawContext.setColor(1, 1, 0);
            drawContext.drawDashLines(global.result.augmentedCurve);
        }
        if(global.checkDrawBand->get_active()) {
            drawContext.setColor(1, 1, 0);
            drawContext.drawBand(global.result.augmentedCurve);
        }
        if(global.checkDrawCurve->get_active()) {
            drawContext.setColor(1, 0, 0);
            drawContext.drawMainLine(global.result.augmentedCurve);
        }
        if(global.checkDrawControlPoints->get_active()) {
            drawContext.setColor(1, 0, 0);
            drawContext.drawControlPoints(global.result.augmentedCurve);
        }
        drawContext.setColor(1, 1, 0);
        drawContext.drawInfo(global.result.augmentedCurve);
    }
    else {
        if ( global.checkDrawBand->get_active()){
            auto tmp = GetBands(&(global.result.graphNodes), global.result.plot.width());
            drawContext.setColor(1, 1, 0);
            drawContext.drawPlotBackground(tmp);
        }
        if(global.checkDrawControlPoints->get_active()) {
            drawContext.setColor(1, 0, 0);
            drawContext.drawControlPoints(global.result.controlPoints);
        }
        if (global.checkDrawCurve->get_active()){
            drawContext.setColor(1, 0, 0);
            drawContext.drawVoronoiGraph(&(global.result));

        }


    }

    imgc->setImage(drawContext.surface());
}

void renderImage() {
    renderImageOn(&global.mainImageController);
}


void computeResult() {
    //switch controls calculation
    if(!global.calculationSwitch->get_active()) return;

    ProjectionMatrix2D matrix = projectionMatrixGet();
    Projection2D projection = calibrateMatrix(&matrix, global.dataSet);
    global.parameters.width = global.adjustmentWidth->get_value();
    global.parameters.height = global.parameters.width;
    global.parameters.smoothCount = global.adjustmentSmooth->get_value();
    global.parameters.smoothCountVariance = global.adjustmentVarianceSmooth->get_value();
    global.parameters.convergenceDistance = 0;
    global.parameters.segment = global.adjustmentSegmentPix->get_value();

    //global.result = processPipelineNormal(global.parameters, global.dataSet, projection);
    double time = 0;
    global.result = processPipeline(global.parameters, global.dataSet, projection, global.generator, &time);

    std::cout << "Time? " << time <<std::endl;

    renderImage();
}

void openData() {
    std::ifstream infile(global.dataFileChooserButton->get_filename());
    if(!infile.is_open()) return;
    global.dataSet = parseData(infile);

    createProjectionGUI();
    computeResult();
}

int main(int argc, char *argv[])
{
    auto app =
        Gtk::Application::create(argc, argv,
                "de.uni-muenster.hks");
    Glib::RefPtr<Gtk::Builder> builder = Gtk::Builder::create_from_file("guidefinition.glade");

    Gtk::Window *window;
    builder->get_widget("main_window", window);

    builder->get_widget("main_image", global.mainImageController.image);
    builder->get_widget("main_image_viewport", global.mainImageController.viewport);
    builder->get_widget("main_image_scrolledwindow", global.mainImageController.viewportScrolledWindow);
    builder->get_widget("data_chooser", global.dataFileChooserButton);
    builder->get_widget("calculation_switch", global.calculationSwitch);
    builder->get_widget("check_draw_plot", global.checkDrawPlot);
    builder->get_widget("check_draw_gauss", global.checkDrawGauss);
    builder->get_widget("check_draw_skeleton", global.checkDrawSkeleton);
    builder->get_widget("check_draw_dash", global.checkDrawDash);
    builder->get_widget("check_draw_band", global.checkDrawBand);
    builder->get_widget("check_draw_curve", global.checkDrawCurve);
    builder->get_widget("check_draw_starting_curve", global.checkDrawStartingCurve);
    builder->get_widget("check_draw_control_points", global.checkDrawControlPoints);
    builder->get_widget("check_mirror", global.checkMirror);


    builder->get_widget("projection_matrix_grid", global.projectionMatrixGrid);
    builder->get_widget("projection_matrix_reset", global.projectionMatrixResetButton);
    builder->get_widget("projection_matrix_random", global.projectionMatrixRandomButton);
    builder->get_widget("plist_count", global.projectionListCountLabel);
    builder->get_widget("plist_map_scrolledwindow", global.projectionListImageController.viewportScrolledWindow);
    builder->get_widget("plist_map_viewport", global.projectionListImageController.viewport);
    builder->get_widget("plist_map_image", global.projectionListImageController.image);
    builder->get_widget("plist_map_window", global.projectionListWindow);

    Gtk::Button *plist_reset, *plist_extend, *plist_map, *plist_map_tsne;
    builder->get_widget("plist_reset", plist_reset);
    builder->get_widget("plist_extend", plist_extend);
    builder->get_widget("plist_map", plist_map);
    builder->get_widget("plist_map_tsne", plist_map_tsne);
    plist_reset->signal_clicked().connect(sigc::ptr_fun(&projectionListReset));
    plist_extend->signal_clicked().connect(sigc::ptr_fun(&projectionListExtend));
    plist_map->signal_clicked().connect(sigc::ptr_fun(&projectionListMap));
    plist_map_tsne->signal_clicked().connect(sigc::ptr_fun(&projectionListMapTSNE));

    Gtk::Button *projectionMatrixSharpenButton;
    builder->get_widget("projection_matrix_sharpen", projectionMatrixSharpenButton);
    projectionMatrixSharpenButton->signal_clicked().connect(sigc::ptr_fun(&projectionMatrixSharpen));

    createProjectionGUI();

    global.mainImageController.adjustmentZoom = Glib::RefPtr<Gtk::Adjustment>::cast_dynamic(builder->get_object("zoom"));
    global.projectionListImageController.adjustmentZoom = Glib::RefPtr<Gtk::Adjustment>::cast_dynamic(builder->get_object("plist_zoom"));
    global.adjustmentSmooth = Glib::RefPtr<Gtk::Adjustment>::cast_dynamic(builder->get_object("smooth"));
    global.adjustmentSegmentPix = Glib::RefPtr<Gtk::Adjustment>::cast_dynamic(builder->get_object("segment_pix"));
    global.adjustmentWidth = Glib::RefPtr<Gtk::Adjustment>::cast_dynamic(builder->get_object("width_setting"));
    global.adjustmentVarianceSmooth = Glib::RefPtr<Gtk::Adjustment>::cast_dynamic(builder->get_object("variance_smooth"));
    global.adjustmentVarianceFactor = Glib::RefPtr<Gtk::Adjustment>::cast_dynamic(builder->get_object("variance_factor"));
    global.adjustmentAmountFactor = Glib::RefPtr<Gtk::Adjustment>::cast_dynamic(builder->get_object("amount_factor"));

    global.adjustmentSmooth->signal_value_changed().connect(sigc::ptr_fun(&computeResult));
    global.adjustmentSegmentPix->signal_value_changed().connect(sigc::ptr_fun(&computeResult));
    global.adjustmentWidth->signal_value_changed().connect(sigc::ptr_fun(&computeResult));
    global.adjustmentVarianceSmooth->signal_value_changed().connect(sigc::ptr_fun(&computeResult));

    global.dataFileChooserButton->signal_file_set().connect(sigc::ptr_fun(&openData));
    global.calculationSwitch->property_active().signal_changed().connect(sigc::ptr_fun(&computeResult));
    global.checkDrawPlot->property_active().signal_changed().connect(sigc::ptr_fun(&renderImage));
    global.checkDrawGauss->property_active().signal_changed().connect(sigc::ptr_fun(&renderImage));
    global.checkDrawSkeleton->property_active().signal_changed().connect(sigc::ptr_fun(&renderImage));
    global.checkDrawDash->property_active().signal_changed().connect(sigc::ptr_fun(&renderImage));
    global.checkDrawBand->property_active().signal_changed().connect(sigc::ptr_fun(&renderImage));
    global.checkDrawCurve->property_active().signal_changed().connect(sigc::ptr_fun(&renderImage));
    global.checkDrawStartingCurve->property_active().signal_changed().connect(sigc::ptr_fun(&renderImage));
    global.checkDrawControlPoints->property_active().signal_changed().connect(sigc::ptr_fun(&renderImage));

    global.projectionMatrixResetButton->signal_clicked().connect(sigc::ptr_fun(&projectionMatrixReset));
    global.projectionMatrixRandomButton->signal_clicked().connect(sigc::ptr_fun(&projectionMatrixRandom));

    global.parameters.asGraph = true;

    global.mainImageController.renderCallback = &renderImageOn;
    global.mainImageController.init();
    global.projectionListImageController.renderCallback = &projectionListMapDrawOn;
    global.projectionListImageController.clickCallback = &projectionListMapClick;
    global.projectionListImageController.init();
    global.generator = new SkeletonGenerator;

    return app->run(*window);

    delete window;
}
