
#include "global.hpp"
#include "projectionmatrix.hpp"

//forward reference
void computeResult();

void createProjectionGUI() {
    int newSize = global.dataSet.dimension();
    int oldSize = global.projectionMatrixX.size();

    if(newSize < oldSize) {
        for(int i = newSize; i < oldSize; i++) {
            global.projectionMatrixGrid->remove_row(newSize);
        }

        global.projectionMatrixX.erase(global.projectionMatrixX.begin() + newSize, global.projectionMatrixX.end());
        global.projectionMatrixY.erase(global.projectionMatrixY.begin() + newSize, global.projectionMatrixY.end());
    } else {

        for(int i = oldSize; i < newSize; i++) {
            Gtk::Label* label = Gtk::manage(new Gtk::Label(std::to_string(i)));
            global.projectionMatrixGrid->attach(*label, 0, i, 1, 1);

            Glib::RefPtr<Gtk::Adjustment> ax = Gtk::Adjustment::create(0, -1, 1, 0.05, 0, 0);
            Glib::RefPtr<Gtk::Adjustment> ay = Gtk::Adjustment::create(0, -1, 1, 0.05, 0, 0);
            ax->signal_value_changed().connect(sigc::ptr_fun(&computeResult));
            ay->signal_value_changed().connect(sigc::ptr_fun(&computeResult));
            global.projectionMatrixX.push_back(ax);
            global.projectionMatrixY.push_back(ay);

            Gtk::SpinButton* spinButtonX = Gtk::manage(new Gtk::SpinButton(ax));
            Gtk::SpinButton* spinButtonY = Gtk::manage(new Gtk::SpinButton(ay));
            spinButtonX->set_digits(2);
            spinButtonY->set_digits(2);
            global.projectionMatrixGrid->attach(*spinButtonX, 1, i, 1, 1);
            global.projectionMatrixGrid->attach(*spinButtonY, 2, i, 1, 1);
        }
    }

    global.projectionMatrixGrid->show_all();
}

void projectionMatrixSet(const ProjectionMatrix2D& matrix) {
    global.calculationSwitch->set_active(false);
    for(int i = 0; i < global.projectionMatrixX.size(); i++) {
        global.projectionMatrixX[i]->set_value(matrix.entry(i)[0]);
    }
    for(int i = 0; i < global.projectionMatrixY.size(); i++) {
        global.projectionMatrixY[i]->set_value(matrix.entry(i)[1]);
    }
    global.calculationSwitch->set_active(true);
}

ProjectionMatrix2D projectionMatrixGet() {
    ProjectionMatrix2D matrix { global.projectionMatrixX.size() };
    for(int i = 0; i < global.projectionMatrixX.size(); i++) {
        matrix.entry(i)[0] = global.projectionMatrixX.at(i)->get_value();
        matrix.entry(i)[1] = global.projectionMatrixY.at(i)->get_value();
    }
    return matrix;
}

void projectionMatrixReset() {
    ProjectionMatrix2D matrix { global.dataSet.dimension() };
    for(int i = 0; i < global.dataSet.dimension(); i++) {
        matrix.entry(i) = 0;
    }
    projectionMatrixSet(matrix);
}

void projectionMatrixRandom() {
    ProjectionMatrix2D matrix = global.projectionMatrixGenerator.generate(global.dataSet);
    projectionMatrixSet(matrix);
}

void projectionMatrixSharpen() {
    float alpha = 0.1;

    float score = (1 - alpha) * sharpenssOfAugmentedCurve(global.result.augmentedCurve);
    ProjectionMatrix2D matrix = projectionMatrixGet();
    ProjectionMatrix2D newMatrix = matrix;
    ProjectionMatrix2D bestMatrix = matrix;

    for(int i = 0; i < matrix.dimension(); i++) {
        for(int j = 0; j < 2; j++) {
            for(int down = 0; down <= 1; down++) {
                newMatrix = matrix;
                if(down == 0) {
                    matrix.entry(i)[j] *= 1.1;
                } else {
                    matrix.entry(i)[j] /= 1.1;
                }

                double time;
                auto result = processPipeline(global.parameters, global.dataSet, calibrateMatrix(&newMatrix, global.dataSet), global.generator, &time);
                //auto result = processPipelineNormal(global.parameters, global.dataSet, calibrateMatrix(&newMatrix, global.dataSet));

                auto newScore = alpha * couplingDistanceAugmented(global.result.augmentedCurve, result.augmentedCurve, 0, 0);
                newScore += (1 - alpha) * sharpenssOfAugmentedCurve(result.augmentedCurve);

                if(newScore < score) {
                    score = newScore;
                    bestMatrix = newMatrix;
                }
            }
        }
    }

    projectionMatrixSet(bestMatrix);
}
