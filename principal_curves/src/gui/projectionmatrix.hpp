#ifndef GUI_PROJECTIONMATRIX_HPP
#define GUI_PROJECTIONMATRIX_HPP

void createProjectionGUI();

ProjectionMatrix2D projectionMatrixGet();

void projectionMatrixSet(const ProjectionMatrix2D& matrix);

void projectionMatrixReset();

void projectionMatrixRandom();

void projectionMatrixSharpen();

#endif
