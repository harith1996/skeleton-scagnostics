#include "skeletonvis.h"

SkeletonVis::SkeletonVis(QWidget *parent) :
    QGLWidget( QGLFormat(QGL::SampleBuffers), parent)
{
    //setMouseTracking(true);
    //setFocusPolicy(Qt::StrongFocus);
    scale = 1.0;
    transX = 0;
    transY = 0;
    isLeftMouseActive = false; isRightMouseActive = false;
    oldMouseX = 0; oldMouseY = 0;

    _splomGenerator = NULL;
    _cooccurCalc = NULL;
    tex3DCombined = NULL;

    texture3DCombined = -1;
    //changeInTextureSkel = false;
    changeInTexture = false;
    imgSize = 300;

    latestDesiredSize = 64;
    useSubset = false;
    glGenTextures(1,&texture3DCombined);
    totalSeenBefore = -1;

    runningCenterline = false;
    pause = false;
    selectedGrid = make_pair(-1,-1);
    selectedCenterline = make_pair(-1,-1);
    selectedCenterlineTree = make_pair(-1,-1);

    areaSelectionInGrid = false;
    highlightedPoints = NULL;

    currentMultiVisMethod = SkeletonVis::SPLOM;

    hoverRow = -1;    hoverColumn = -1;
    hoverSphereLoc[0] = 0;     hoverSphereLoc[1] = 0;
    hoverRowLoc[0] = 0;  hoverRowLoc[1] = 0;
    currentRowSelected = -1; currentColSelected = -1;

    brushHistSize = 0.05;
    brushSkeletonSize = 0.005;
    tfDialog = new TransferFunctionDialog();
    //freeFromDialog = new FreeFormWidget();

    mapProperties = 0;
    useDensityContour = false;
    similarityTreshold = 0.10;
    contourOpacity = 0.5;

    colorContour.setRgb(0,0,0);
    colorCenterline.setRgb(255,0,0);

    fromBackground.setRgb(228,251, 153);
    toBackground.setRgb(255,202, 228);

    displayCenterlineLine = true;
    centerlineWidth = 2.0;
    useMDSCenterlines = false;

    typeOfMeasureMD = SkeletonVis::Wilkinson;
    scagnosticsValues = NULL;
    tatuValues = NULL;
    modelComparison = false;
    a_model = 0;
    b_model = 1.0;
    c_model = 0;

    multiResolution = false;
    doubleClickMode = 0;
    useFreeFormModel = false;

}

SkeletonVis::~SkeletonVis()
{
   ClearMemory();
}

void SkeletonVis::ChangeFromSimilarityColor(QColor c){
    fromBackground = c;
    std::cout << "From " << c.red() << " ," << c.green() << " , " << c.blue() << std::endl;

}

void SkeletonVis::ChangeToSimilarityColor(QColor c){
    toBackground = c;
    std::cout << "To " << c.red() << " ," << c.green() << " , " << c.blue() << std::endl;
}
void SkeletonVis::ChangeContourOpacity(float val){
    contourOpacity = val;
}

void SkeletonVis::ChangeContourColor(QColor val){
    colorContour = val;
}

void SkeletonVis::ChangeCenterlineColor(QColor val){
    colorCenterline = val;

}

void SkeletonVis::SetThreshParams(float start, float stride, float limit){
    this->thresStart = start;
    this->thresStride = stride;
    this->thresLimit = limit;
    generator.SetThreshParams(start, stride, limit);
}
void SkeletonVis::ReorderAccordingToColumn(int col){
    vector<int> indices;
    vector< pair<int,int> > nominalPairs;
    vector<float> associationValues;

    // We basically get the first categorical attribute
    // and use it as the basis for the order, the
    // correlation to it is used
    _cooccurCalc->GetAssociationsInfo(&nominalPairs, &associationValues);
    _cooccurCalc->GetCategoricalIndices(&indices);

    int actualCol = columnOrder.at(col);
    int firstIndex = indices.at(actualCol);

    vector< pair<int, float> > actualOrder;
    actualOrder.push_back( make_pair(firstIndex, 1.0));

    for(unsigned int i = 0; i < nominalPairs.size(); i++){
        if ( nominalPairs.at(i).first == firstIndex){
            actualOrder.push_back(make_pair(nominalPairs.at(i).second , associationValues.at(i)) );
        }
        else if ( nominalPairs.at(i).second == firstIndex){
            actualOrder.push_back(make_pair(nominalPairs.at(i).first , associationValues.at(i)) );
        }
    }

    sort(actualOrder.begin(), actualOrder.end(),
         boost::bind(&std::pair<int,float>::second,_1) > boost::bind(&std::pair<int,float>::second,_2) );

    columnOrder.clear();

    for(unsigned int i = 0; i< actualOrder.size(); i++){
        for(unsigned int j = 0; j < indices.size(); j++){
            if ( indices.at(j) == actualOrder.at(i).first){
                columnOrder.push_back(j);
            }
        }
    }
    updateGL();
}

void SkeletonVis::CreateInitialHistogramOrdering(){
    if ( _cooccurCalc == NULL) return;
    columnOrder.clear();
    rowOrder.clear();
    rowBelonging.clear();
    //
    int cols = _cooccurCalc->GetNumberOfCategoricalDims();
    int sizes[cols];
    int rows = _cooccurCalc->GetTotalSumCategoricalSize(sizes);
    // Now we sort the columns according to their associations

    vector<int> indices;
    vector< pair<int,int> > nominalPairs;
    vector<float> associationValues;

    // We basically get the first categorical attribute
    // and use it as the basis for the order, the
    // correlation to it is used
    _cooccurCalc->GetAssociationsInfo(&nominalPairs, &associationValues);
    _cooccurCalc->GetCategoricalIndices(&indices);
    //***************************

    vector< pair<int, int> > quantityRowOccurences;
    // How many samples are in each row ...
    int pos = 0, stRow = 0;
    for(int col = 0; col < cols ; col++){
        for(int row = 0; row < rows; row++){
            //  int GetTotalOfCategoricalValue(int attributeIndex, int valueIndex);
            if ( stRow <= row && row < stRow + sizes[col]){
                int totalOccur = _cooccurCalc->GetOccurenceNumberOfAttrAndCategory(pos);
                quantityRowOccurences.push_back( make_pair(row,totalOccur));
                pos++;
                rowBelonging.push_back( make_pair(indices.at(col), row- stRow));
            }
        }
        stRow += sizes[col];
    }

    // Sort according to the number of occurences

    sort(quantityRowOccurences.begin(), quantityRowOccurences.end(),
         boost::bind(&std::pair<int,int>::second,_1) > boost::bind(&std::pair<int,int>::second,_2) );

    for(unsigned int i = 0; i < quantityRowOccurences.size(); i++){
        rowOrder.push_back(quantityRowOccurences.at(i));
    }

    int firstIndex = indices.at(0);

    vector< pair<int, float> > actualOrder;
    actualOrder.push_back( make_pair(firstIndex, 1.0));

    for(unsigned int i = 0; i < nominalPairs.size(); i++){
        if ( nominalPairs.at(i).first == firstIndex){
            actualOrder.push_back(make_pair(nominalPairs.at(i).second , associationValues.at(i)) );
        }
        else if ( nominalPairs.at(i).second == firstIndex){
            actualOrder.push_back(make_pair(nominalPairs.at(i).first , associationValues.at(i)) );
        }
    }

    sort(actualOrder.begin(), actualOrder.end(),
         boost::bind(&std::pair<int,float>::second,_1) > boost::bind(&std::pair<int,float>::second,_2) );

    for(unsigned int i = 0; i< actualOrder.size(); i++){
        for(unsigned int j = 0; j < indices.size(); j++){
            if ( indices.at(j) == actualOrder.at(i).first){
                columnOrder.push_back(j);
            }
        }
    }
}

void SkeletonVis::ShowTF(){
    tfDialog->show();
}

void SkeletonVis::ChangeMultiVisMethod(int method){
    currentMultiVisMethod = method;
    brushingArea.clear();

    if ( currentMultiVisMethod == SkeletonVis::MDS && mdsPoints.empty()){
        CreateMDSProjection();
    }
    Redraw();
}

void SkeletonVis::ClearMemory(){
    for(unsigned int i = 0; i < centerlineMemory.size(); i++){
        Memoization* mem = centerlineMemory.at(i);
        delete mem;
    }
    centerlineMemory.clear();
    //changeInTextureSkel = false;
    changeInTexture = false;
    useSubset = false;
    totalSeenBefore = -1;

    if (scagnosticsValues != NULL){
       free(scagnosticsValues);
       scagnosticsValues = NULL;
    }

    if ( tatuValues != NULL){
        free(tatuValues);
        tatuValues = NULL;
    }
}

void SkeletonVis::CreateTatuMeasures(){

    if ( tatuValues != NULL) return;

    vector<int> numericalAttributes;
    _splomGenerator->GetListOfNumericalAttributes(&numericalAttributes);

    WilkinsonScagnostics scagnosticsMeasures;
    scagnosticsMeasures.SetSplomGenerator(_splomGenerator);
    int num = numericalAttributes.size();

    vector< pair<int,float> > tatuLoc;


    int pos =0;
    for(int i = 0; i < num; i++){
        for(int j = i+1; j < num; j++ ){
            double tatu = scagnosticsMeasures.GetRVM(numericalAttributes.at(i),numericalAttributes.at(j),-1,-1);
            tatuLoc.push_back( make_pair(pos, tatu) );
            pos++;
        }
    }


    tatuValues = (double*) malloc(sizeof(double)*pos);

    for(int i = 0; i < pos; i++){
        tatuValues[i] = tatuLoc.at(i).second;
    }

}

void SkeletonVis::CreateScagnosticsMeasures(){

    if (scagnosticsValues != NULL) return;

    std::cout << "Going to create MDS Scagnostics "<< std::endl;
    vector<int> numericalAttributes;
    _splomGenerator->GetListOfNumericalAttributes(&numericalAttributes);
    int num = numericalAttributes.size();
    int totalImages = ((num)*(num-1))/2; // n*(n-1) /2

    WilkinsonScagnostics scagnosticsMeasures;


    scagnosticsMeasures.SetSplomGenerator(_splomGenerator);

    scagnosticsValues = (double*)malloc(sizeof(double)*totalImages*9);

    // First create the nine measures

    int pos = 0;
    for(int i = 0; i < num; i++){
        for(int j = i+1; j < num; j++ ){
            //pairs.push_back(make_pair(numericalAttributes.at(i),numericalAttributes.at(j)));
            vector<float> edgeDistances;
            vector<LocalPoint> points;
            vtkSmartPointer<vtkTree> mst = scagnosticsMeasures.CreateMST(numericalAttributes.at(i),
                                                                         numericalAttributes.at(j), -1,-1, &edgeDistances ,&points,true);




            float w = scagnosticsMeasures.GetW(&edgeDistances);
            vector<int> outliers;
            scagnosticsMeasures.GetOutliers(mst, &outliers, w, &points );

            scagnosticsValues[pos*9 + 0] = scagnosticsMeasures.SparsenessMeasure(mst,&outliers,&points );
            scagnosticsValues[pos*9 + 1] = scagnosticsMeasures.SkewednessMeasure(mst, &outliers,&points);
            scagnosticsValues[pos*9 + 2] = scagnosticsMeasures.OutlyingMeasure(mst, &outliers,&points);
            scagnosticsValues[pos*9 + 3] = scagnosticsMeasures.StringyMeasure(mst, &outliers);
            scagnosticsValues[pos*9 + 4] = scagnosticsMeasures.StriateMeasure(mst, &outliers, &points);
            scagnosticsValues[pos*9 + 5] = scagnosticsMeasures.ClumpyMeasure(mst, &outliers, &points);


            vtkSmartPointer<vtkPolyData> alphaShape = scagnosticsMeasures.GetAlphaShape(scagnosticsValues[pos*9 +0],&outliers, &points);
            vtkSmartPointer<vtkPolyData> convexHull = scagnosticsMeasures.GetConvexHull(&outliers, &points);

            scagnosticsValues[pos*9 + 6] = scagnosticsMeasures.ConvexMeasure(alphaShape, convexHull);
            scagnosticsValues[pos*9 + 7] = scagnosticsMeasures.SkinnyMeasure(alphaShape);
            scagnosticsValues[pos*9 + 8] = scagnosticsMeasures.MonotonicMeasure(&outliers, &points);

            std::cout << "For pos " << pos << " :: " << std::endl;

            /*
            std::cout <<"Tatu " << scagnosticsMeasures.GetRVM(numericalAttributes.at(i),
                                                              numericalAttributes.at(j), -1,-1) << std::endl; */
            for(int k = 0; k  < 9; k++){
                std::cout << scagnosticsValues[pos*9 + k] << "...";
            }
            std::cout << std::endl;

            if ( pos == 0){ // initialize the maxs & mins
                   for(int k = 0; k  < 9; k++){
                       scagnosticsMaxsMins[2*k +0] = scagnosticsValues[pos*9 + k];
                       scagnosticsMaxsMins[2*k +1]  = scagnosticsValues[pos*9 + k];
                   }
            }
            else {
                for(int k = 0; k  < 9; k++){
                   if (scagnosticsValues[pos*9 +k] > scagnosticsMaxsMins[2*k +0] )  scagnosticsMaxsMins[2*k +0]  = scagnosticsValues[pos*9 + k];
                   if (scagnosticsValues[pos*9 +k] < scagnosticsMaxsMins[2*k +1] )  scagnosticsMaxsMins[2*k +1]  = scagnosticsValues[pos*9 + k];
                }
            }

            pos++;
        }
    }


    // Now we can sort them
    orderedMeasures.clear();
    for(int k = 0; k < 9; k++){
        // For each measure
        vector< pair<int, float> > measure;
        int pos = 0;
        for(int i = 0; i < num; i++){
            for(int j = i+1; j < num; j++ ){
                float v = scagnosticsValues[pos*9 +k];


                if ( k != 0 && k != 8){

                    float range = scagnosticsMaxsMins[2*k +0] - scagnosticsMaxsMins[2*k +1];
                    if ( range < 0.00001 ) v = 0;
                    else v = (v - scagnosticsMaxsMins[2*k +1])/range;

                }
                measure.push_back( make_pair(pos, v));
                pos++;
            }
        }
        //***************
        std::sort( measure.begin(), measure.end(),
               boost::bind(&std::pair<int, float>::second, _1) >
               boost::bind(&std::pair<int,float>::second, _2 ));

         orderedMeasures.push_back(measure);
         // max first


         /*if (k == 8){
             for(int h = 0; h < measure.size(); h++){
                 std::cout << measure.at(h).first << ", " << measure.at(h).second << "...";
             }
             std::cout << std::endl;
         } */

    }

}


void SkeletonVis::SetSplomGenerator(SPLOMWrapper* wrapper, SPLOMWrapper* drawer){
    _splomGenerator = wrapper;
    _splomDrawer = drawer;
}

void SkeletonVis::SetCoocurrenceInfo(CooccurrenceCalculator *cooccurCalc){
    _cooccurCalc = cooccurCalc;
}

void SkeletonVis::SkeletonPaint(){
    if (pause) return;

    glEnable( GL_LINE_SMOOTH );
    glEnable( GL_POLYGON_SMOOTH );
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    glViewport(0, 0, (GLsizei) viewport[2], (GLsizei) viewport[3]);
    glClearColor(0.95, 0.95, 0.95, 0.0);
    glClear(GL_COLOR_BUFFER_BIT);

    bool tex_interp = true;
    // Setup projection matrix
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(0.0, 1.0, 0.0, 1.0);

    // Setup modelview matrix
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glScalef(scale, scale, 1.0);
    glTranslatef(transX, transY, 0.0);
    glBindTexture(GL_TEXTURE_3D, 0);

    glLineWidth(1.0);

    HighlightSimilarCenterlines();

   if ( _splomGenerator != NULL){
        if (_splomGenerator->HasData() && _splomGenerator->IsFinished()){

            int totalNumerical = _splomGenerator->GetTotalNumericalAttributes();

            vector<int> list;
            _splomGenerator->GetListOfNumericalAttributes(&list);

            float blockSize = 1.0 / (totalNumerical);

            // We have now them combined ...
            int totalImages = (totalNumerical*(totalNumerical-1))/2;
            float planeSize = (1.0 / (totalImages  + 1.0)) * 0.5;

            int sizeT = _splomGenerator->GetSizeTexture();
            int pos = 1;
            float extra = -planeSize*0.5;

            //glActiveTexture(GL_TEXTURE0);

           if( tex3DCombined != NULL){

                // If there is a change in texture, re-set it ....
                if ( changeInTexture){
                    // glGenTextures(1,&texture3D);
                    glBindTexture(GL_TEXTURE_3D, texture3DCombined);
                    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, (tex_interp)?GL_LINEAR:GL_NEAREST);
                    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, (tex_interp)?GL_LINEAR:GL_NEAREST);
                    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
                    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
                    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
                    glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
                    glTexImage3D(GL_TEXTURE_3D, 0, GL_RGBA,sizeT, sizeT,(totalImages+1)*2, 0, GL_RGBA, GL_UNSIGNED_BYTE, tex3DCombined);
                    changeInTexture = false;
                }
                else {
                    glBindTexture(GL_TEXTURE_3D, texture3DCombined);
                }
                glEnable(GL_TEXTURE_3D);

                glBegin(GL_QUADS);
                for(int i = 0; i < totalNumerical; i++){
                    for(int j = i+1; j < totalNumerical; j++ ){

                        int x = j;
                        int y = totalNumerical - i- 1;//totalNumerical - j - 2;

                        int o = pos;

                        if ( !reOrderedSplom.empty()){
                            o =  reOrderedSplom.at(pos-1) + 1;
                        }
                        // (0,1) (1,1)  (1,0)  (0,0)
                        glTexCoord3f(0.0, 1.0,  o*planeSize + extra ); glVertex2f(x*blockSize, y*blockSize);
                        glTexCoord3f(1.0, 1.0,  o*planeSize + extra); glVertex2f((x+1)*blockSize, y*blockSize);
                        glTexCoord3f(1.0, 0.0,  o*planeSize + extra); glVertex2f((x+1)*blockSize, (y+1)*blockSize);
                        glTexCoord3f(0.0, 0.0,  o*planeSize + extra); glVertex2f(x*blockSize, (y+1)*blockSize);
                        pos += 1;
                    }
                    //break;
                }
                glEnd();
                pos += 1;

                // Now the centerlines
                glColor3f(0.0,1.0,0);
                glBegin(GL_QUADS);


                int locInSimilarity = 0;

                int stPos = pos;
                vector< pair<int, int> > originalPairs;

                for(int i = 0; i < totalNumerical; i++){
                    for(int j = i+1; j < totalNumerical; j++ ){
                        originalPairs.push_back( make_pair(i,j));
                    }
                }
                // perhaps here  similarCenterlines
                for(int i = 0; i < totalNumerical; i++){
                    // (1,1)  (0,1)  (0,0) ( 1,0)

                    for(int j = i+1; j < totalNumerical; j++ ){
                    //for(int j = 0; j < i; j++ ){
                        int x = i;
                        int y = totalNumerical -j- 1;

                        int o = pos;

                        if ( !reOrderedSplom.empty()){
                            o =  reOrderedSplom.at(pos-stPos) + stPos;
                        }



                        if (!reOrderedSplom.empty()){

                               glTexCoord3f(1.0, 0.0,  o*planeSize +extra);
                               glVertex2f(x*blockSize, y*blockSize);
                               glTexCoord3f(1.0, 1.0,  o*planeSize  +extra);
                               glVertex2f((x+1)*blockSize, y*blockSize );
                               glTexCoord3f(0.0, 1.0,  o*planeSize +extra);
                               glVertex2f((x+1)*blockSize, (y+1)*blockSize);
                               glTexCoord3f(0.0, 0.0,  o*planeSize  +extra);
                               glVertex2f(x*blockSize, (y+1)*blockSize);

                        }
                        else if (similarCenterlines.empty()){
                            glTexCoord3f(1.0, 0.0,  o*planeSize +extra);
                            glVertex2f(x*blockSize, y*blockSize);
                            glTexCoord3f(1.0, 1.0,  o*planeSize  +extra);
                            glVertex2f((x+1)*blockSize, y*blockSize );
                            glTexCoord3f(0.0, 1.0,  o*planeSize +extra);
                            glVertex2f((x+1)*blockSize, (y+1)*blockSize);
                            glTexCoord3f(0.0, 0.0,  o*planeSize  +extra);
                            glVertex2f(x*blockSize, (y+1)*blockSize);
                        }
                        else {


                               //----------
                               int c1 = i;
                               int c2 = j;

                               /*if (!reOrderedSplom.empty()){
                                    pair<int,int> actual = originalPairs.at(pos-stPos);
                                    c1 = actual.first;
                                    c2 = actual.second;
                               }*/

                               //******************
                              for(unsigned int k = 0; k < similarCenterlines.size(); k++){

                                  if ( c1 == similarCenterlines.at(k).first && c2 == similarCenterlines.at(k).second){
                                      glTexCoord3f(1.0, 0.0,  o*planeSize +extra);
                                      glVertex2f(x*blockSize, y*blockSize);
                                      glTexCoord3f(1.0, 1.0,  o*planeSize  +extra);
                                      glVertex2f((x+1)*blockSize, y*blockSize );
                                      glTexCoord3f(0.0, 1.0,  o*planeSize +extra);
                                      glVertex2f((x+1)*blockSize, (y+1)*blockSize);
                                      glTexCoord3f(0.0, 0.0,  o*planeSize  +extra);
                                      glVertex2f(x*blockSize, (y+1)*blockSize);
                                      break;
                                  }
                              }

                        }

                        locInSimilarity++;
                        pos += 1;
                    }
                }
                glEnd();
            glDisable(GL_TEXTURE_3D);

       }
            glBindTexture(GL_TEXTURE_3D, 0);
            //*******************************
            glBegin(GL_QUADS);
            for(int i = 0; i < totalNumerical; i++){
                     int x = i;
                     int y = totalNumerical - i - 1;
                     glColor3f(0.95,0.95,0.95);
                     glVertex2f(x*blockSize, y*blockSize);
                     glVertex2f((x+1)*blockSize, y*blockSize);
                     glVertex2f((x+1)*blockSize, (y+1)*blockSize);
                     glVertex2f(x*blockSize, (y+1)*blockSize);
            }
            glEnd();

            glBegin(GL_LINES);
            for(int i = 0; i < totalNumerical+1; i++){
                     glColor3f(0,0,0);
                     glVertex2f(i*blockSize, 0);
                     glVertex2f(i*blockSize, 1.0);

                     glVertex2f(0, i*blockSize);
                     glVertex2f(1.0, i*blockSize);
            }
            glEnd();
            glColor3f(0.0,0.0,0.0);

            if ( names.size() > 0){
                for(int i = 0; i < totalNumerical;i++){
                    int x = i;
                    int y = totalNumerical - i - 1;

                    if (reOrderedSplom.empty()){
                        renderText(x*blockSize + blockSize*0.15, y*blockSize+ blockSize*0.5, 0, names.at(i).c_str());
                    }
                    else {
                        renderText(x*blockSize + blockSize*0.15, y*blockSize+ blockSize*0.5, 0, names.at(orderAccordingToSimilarity.at(i)).c_str());

                    }
                }
            }

            // now the contours
            if (useDensityContour && !contours.empty()){
                int loc = 0;
                for(int i = 0; i < totalNumerical; i++){
                    for(int j = i +1; j < totalNumerical; j++ ){
                        int x = j;
                        int y = totalNumerical - i- 1;
                        float startX = x*blockSize;
                        float endX = (x+1)*blockSize;
                        float startY = y*blockSize;
                        float endY = (y+1)*blockSize;

                        vtkSmartPointer<vtkPolyData> currentContour = contours.at(loc);
                        vtkSmartPointer<vtkCellArray> cellArray = vtkSmartPointer<vtkCellArray>::New();
                        vtkSmartPointer<vtkPoints> points = currentContour->GetPoints();
                        vtkSmartPointer<vtkIdList> list = vtkSmartPointer<vtkIdList>::New();
                        cellArray = currentContour->GetLines();

                        double p[2];
                        glColor3f(0.0,0,0);
                        glBegin(GL_LINES);
                        for(int i = 0; i <  currentContour->GetNumberOfLines(); i++){
                            cellArray->GetNextCell(list);
                            for(int k = 0; k < list->GetNumberOfIds(); k++){
                                float tx, ty;
                                points->GetPoint(list->GetId(k),p);
                                tx = (endX-startX)*p[0] + startX;
                                ty = (endY-startY)*(1.0 - p[1]) + startY;

                                glVertex3f(tx, ty,-0.1);
                            }
                        }
                        glEnd();

                       loc++;
                    }
                }
            }
        }
    }


   glColor4f(1.0,1.0,1.0, 0.7);
   glBegin(GL_TRIANGLE_FAN);
       for(unsigned int i = 0; i < areaSelectionOutlineInGrid.size(); i++ ){
           // Convert the points from the real world to the image world...
           LocalPoint p = areaSelectionOutlineInGrid.at(i);
           float mx = p.p[0];
           float my = p.p[1];
           glVertex3f(mx,my, -0.3);
       }
   glEnd();


   for(unsigned int i = 0; i < brushingArea.size(); i++){
       glColor4f(1.0,1.0,0.0, 0.1);
       glPushMatrix();
          glTranslatef( brushingArea.at(i).p[0],brushingArea.at(i).p[1], 0);
          glutSolidSphere(brushSkeletonSize ,10,10);
       glPopMatrix();
   }

   glPointSize(2.0);
   glBegin(GL_POINTS);
   for(unsigned int i = 0; i < selectedCenterlinePoints.size(); i++){
       glColor3f(0,0,1);
       glVertex3f(selectedCenterlinePoints.at(i).first.p[0], selectedCenterlinePoints.at(i).first.p[1], 0);  
   }
   glEnd();


  DisplayCenterlineWithInfo();
  glBindTexture(GL_TEXTURE_3D, 0);

}

void SkeletonVis::HighlightSimilarCenterlines(){
    if  (selectedCenterlineTree.first != -1 ){
        vector<int> numericalAttributes;
        _splomGenerator->GetListOfNumericalAttributes(&numericalAttributes);
        int totalNumerical = numericalAttributes.size();
        float blockSize = 1.0 / (totalNumerical);
        float color[3] = {0,0,1};
        glLineWidth(2.0);

        float z = -0.1;

        int loc = 0;
        if (!similarityDistancesFromSelected.empty()){
            float maxFrechet = similarityDistancesFromSelected.at(0);

            for(unsigned int i = 0; i < similarityDistancesFromSelected.size();i++)
                if ( similarityDistancesFromSelected.at(i) > maxFrechet) maxFrechet = similarityDistancesFromSelected.at(i);

            for(int i = 0; i < totalNumerical; i++){
                for(int j = i +1 ; j < totalNumerical; j++ ){
                    int x = i;
                    int y = totalNumerical -j - 1;

                    int o = loc;
                    if ( !reOrderedSplom.empty()){
                        o =  reOrderedSplom.at(loc);
                    }
                           float frechetValue = similarityDistancesFromSelected.at(o) / maxFrechet;

                           float backgroundColor[4];
                           backgroundColor[0] = (fromBackground.red()*(1.0 - frechetValue) + frechetValue*toBackground.red())/ 255.0f;
                           backgroundColor[1] = (fromBackground.green()*(1.0 - frechetValue) + frechetValue*toBackground.green())/ 255.0f;
                           backgroundColor[2] = (fromBackground.blue()*(1.0 - frechetValue) + frechetValue*toBackground.blue())/ 255.0f;
                           backgroundColor[3] = 1.0;
                           glColor4fv(backgroundColor);
                           glBegin(GL_QUADS);
                               glVertex3f(x*blockSize, y*blockSize, z);
                               glVertex3f((x+1)*blockSize, y*blockSize, z );
                               glVertex3f((x+1)*blockSize, (y+1)*blockSize, z);
                               glVertex3f(x*blockSize, (y+1)*blockSize,z );
                           glEnd();
                     loc++;
               }
            }
        }


        if (doubleClickMode == 0){
              for(unsigned int i =0 ; i < similarCenterlines.size(); i++){
                    int x = similarCenterlines.at(i).first;
                    int y = totalNumerical -similarCenterlines.at(i).second - 1;


                     if ( selectedCenterline.first == similarCenterlines.at(i).first && selectedCenterline.second == similarCenterlines.at(i).second && !modelComparison)
                         glColor3f(1.0,0,0);
                     else
                      glColor3fv(color);

                    glBegin(GL_LINE_LOOP);
                        glVertex2f(x*blockSize, y*blockSize);
                        glVertex2f((x+1)*blockSize, y*blockSize);
                        glVertex2f((x+1)*blockSize, (y+1)*blockSize);
                        glVertex2f(x*blockSize, (y+1)*blockSize);
                    glEnd();

           }
        }
       glLineWidth(1.0);
    }
}

void SkeletonVis::DisplayMDSCenterlineWithInfo(){

    if (!mapProperties) return;

    //mapProperties = 3;
    //mapProperties = 4;
    if ( displayCenterlineLine && centerlineWidth > 0.01){
        if ( _splomGenerator != NULL){
                if (_splomGenerator->HasData() && _splomGenerator->IsFinished()){

                    vector<int> numericalAttributes;
                    _splomGenerator->GetListOfNumericalAttributes(&numericalAttributes);

                    int totalNumerical = numericalAttributes.size();

                    float blockSize = 1.0 / (totalNumerical);
                    blockSize *= 0.9;

                    double maxD, minD, maxDistance;

                    glLineWidth(centerlineWidth);


                    TransferFunction* tf  = new TransferFunction();


                    //tf->Initialize("color_schemes.txt");
                    tf->CreateDefaultDensityMap();


                    vector<LocalPoint> freeFormPoints;
                    if (useFreeFormModel){
                          std::cout << "Getting the curve? " << std::endl;
                          freeFromDialog->GetCurve(&freeFormPoints);
                    }

                    int posInMDS = 0;
                    for(int i = 0; i < totalNumerical; i++){
                        for(int j = i+1; j < totalNumerical; j++ ){
                            //int x = i;
                            //int y = totalNumerical -j- 1;

                            LocalPoint p = mdsPoints.at(posInMDS);
                            posInMDS++;
                            float startX = p.p[0] - blockSize/2.0f;
                            float startY = p.p[1] - blockSize/2.0f;


                            vector<SkeletonSegment*> currentTree;
                            GetSkeletonSegments(numericalAttributes.at(i), numericalAttributes.at(j), -1,-1, &currentTree, &maxD, &minD, &maxDistance);

                            // std::cout << "Max D " << maxD << " , " << minD << " .. " << maxDistance << std::endl;
                            for(unsigned int h = 0; h < currentTree.size(); h++){
                                SkeletonSegment* currentSegment = currentTree.at(h);
                                glBegin(GL_LINE_STRIP);

                                for(unsigned int k= 0; k < currentSegment->imageCoorPoints.size(); k++){
                                    double density = currentSegment->densities.at(k);
                                    double distance = currentSegment->distances.at(k);


                                    if (modelComparison){

                                        if ( mapProperties == 4) {// combined matric
                                            glColor3f(0,0,0);

                                            float x = currentSegment->imageCoorPoints.at(k).p[0];
                                            float y = currentSegment->imageCoorPoints.at(k).p[1];


                                            float desiredY = a_model*x*x + b_model*x + c_model;

                                            float d = fabs(desiredY - currentSegment->imageCoorPoints.at(k).p[1]);
                                            float de = (density - minD)/(maxD - minD);

                                            if ( useFreeFormModel){
                                                  d = 1.0;

                                                  for(unsigned int  k = 0; k <  freeFormPoints.size(); k++){
                                                        float ox = freeFormPoints.at(k).p[0];
                                                        float oy = freeFormPoints.at(k).p[1];

                                                        float newD = sqrt(  pow(x-ox,2.0) + pow(y-oy,2.0));
                                                        if (newD < d)
                                                            d = newD;

                                                  }
                                            }


                                            //max(currentSegment->imageCoorPoints.at(k).p[0],1.0 - currentSegment->imageCoorPoints.at(k).p[0]);
                                            float col[3];
                                            tf->GetColor(de, col);
                                            glColor4f(col[0], col[1], col[2], 1.0 - d);

                                        }
                                        else {
                                            float x = currentSegment->imageCoorPoints.at(k).p[0];
                                            float y = currentSegment->imageCoorPoints.at(k).p[1];

                                            float desiredY = a_model*x*x + b_model*x + c_model;
                                            float d = fabs(desiredY - currentSegment->imageCoorPoints.at(k).p[1]);

                                            if ( useFreeFormModel){
                                                d = 1.0;
                                                for(unsigned int  k = 0; k <  freeFormPoints.size(); k++){
                                                      float ox = freeFormPoints.at(k).p[0];
                                                      float oy = freeFormPoints.at(k).p[1];

                                                      float newD = sqrt(  pow(x-ox,2.0) + pow(y-oy,2.0));
                                                      if (newD < d)
                                                          d = newD;

                                                }
                                            }

                                            float col[3];
                                            tf->GetColor(d, col);
                                            glColor3fv(col);

                                        }
                                    }
                                    else if (mapProperties == 1){
                                         // Display density
                                         float v = (density - minD)/(maxD - minD);
                                         float col[3];
                                         tf->GetColor(v, col);
                                         glColor3fv(col);
                                    }
                                    else if (mapProperties == 2){
                                         // Display distance
                                        float v = distance / maxDistance;
                                        // blue to red

                                        float col[3];
                                        tf->GetColor(v, col);
                                        glColor3fv(col);

                                    }

                                    glVertex3f(startX + blockSize*currentSegment->imageCoorPoints.at(k).p[0],
                                               startY + blockSize*currentSegment->imageCoorPoints.at(k).p[1],
                                               0);
                                }
                                glEnd();
                            }

                        }
                    }
             }
        }


    }



}

void SkeletonVis::CompareModel(bool compare, float a, float b, float c ){

    modelComparison = compare;
    a_model = a;
    b_model = b;
    c_model = c;
    Redraw();
}

void SkeletonVis::DisplayCenterlineWithInfo(){

    //std::cout << "Calling to display centerline with info ? " << std::endl;
    if (!mapProperties) return;
    //mapProperties = 3;

    if ( displayCenterlineLine && centerlineWidth > 0.01){
        if ( _splomGenerator != NULL){
                if (_splomGenerator->HasData() && _splomGenerator->IsFinished()){

                    vector<int> numericalAttributes;
                    _splomGenerator->GetListOfNumericalAttributes(&numericalAttributes);

                    int totalNumerical = numericalAttributes.size();

                    float blockSize = 1.0 / (totalNumerical);

                    double maxD, minD, maxDistance;

                    glLineWidth(centerlineWidth);


                    TransferFunction* tf  = new TransferFunction();


                    //tf->Initialize("color_schemes.txt");
                    tf->CreateDefaultDensityMap();
                    vector<LocalPoint> freeFormPoints;
                    if (useFreeFormModel){
                          std::cout << "Getting the curve? " << std::endl;
                          freeFromDialog->GetCurve(&freeFormPoints);
                    }

                    for(int i = 0; i < totalNumerical; i++){
                        for(int j = i+1; j < totalNumerical; j++ ){
                            int x = i;
                            int y = totalNumerical -j- 1;


                            float startX = x*blockSize;
                            float startY = y*blockSize;

                            vector<SkeletonSegment*> currentTree;
                            GetSkeletonSegments(numericalAttributes.at(i), numericalAttributes.at(j), -1,-1, &currentTree, &maxD, &minD, &maxDistance);

                            // std::cout << "Max D " << maxD << " , " << minD << " .. " << maxDistance << std::endl;
                            for(unsigned int h = 0; h < currentTree.size(); h++){
                                SkeletonSegment* currentSegment = currentTree.at(h);
                                glBegin(GL_POINTS);

                                for(unsigned int k= 0; k < currentSegment->imageCoorPoints.size(); k++){
                                    double density = currentSegment->densities.at(k);
                                    double distance = currentSegment->distances.at(k);


                                    if (modelComparison){

                                        if ( mapProperties == 4) {// combined matric
                                            glColor3f(0,0,0);

                                            float x = currentSegment->imageCoorPoints.at(k).p[0];
                                            float py = currentSegment->imageCoorPoints.at(k).p[1];

                                            float desiredY = a_model*x*x + b_model*x + c_model;

                                            float d = fabs(desiredY - currentSegment->imageCoorPoints.at(k).p[1])/0.5;
                                            float de = (density - minD)/(maxD - minD);


                                            if ( useFreeFormModel){
                                                  d = 1.0;

                                                  for(unsigned int  k = 0; k <  freeFormPoints.size(); k++){
                                                        float ox = freeFormPoints.at(k).p[0];
                                                        float oy = freeFormPoints.at(k).p[1];

                                                        float newD = sqrt(  pow(x-ox,2.0) + pow(py-oy,2.0));
                                                        if (newD < d)
                                                            d = newD;

                                                  }
                                            }

                                            //max(currentSegment->imageCoorPoints.at(k).p[0],1.0 - currentSegment->imageCoorPoints.at(k).p[0]);
                                            float col[3];
                                            tf->GetColor(de, col);
                                            glColor4f(col[0], col[1], col[2], 1.0 - d);

                                        }
                                        else {
                                            float x = currentSegment->imageCoorPoints.at(k).p[0];
                                            float py = currentSegment->imageCoorPoints.at(k).p[1];

                                            float desiredY = a_model*x*x + b_model*x + c_model;
                                            float d = fabs(desiredY - currentSegment->imageCoorPoints.at(k).p[1])/0.5;


                                            if ( useFreeFormModel){
                                                  d = 1.0;

                                                  for(unsigned int  k = 0; k <  freeFormPoints.size(); k++){
                                                        float ox = freeFormPoints.at(k).p[0];
                                                        float oy = freeFormPoints.at(k).p[1];

                                                        float newD = sqrt(  pow(x-ox,2.0) + pow(py-oy,2.0));
                                                        if (newD < d)
                                                            d = newD;

                                                  }
                                            }

                                            float col[3];
                                            tf->GetColor(d, col);
                                            glColor3fv(col);

                                        }
                                    }
                                    else if (mapProperties == 1){
                                         // Display density
                                         float v = (density - minD)/(maxD - minD);
                                         float col[3];
                                         tf->GetColor(v, col);
                                         glColor3fv(col);
                                    }
                                    else if (mapProperties == 2){
                                         // Display distance
                                        float v = distance / maxDistance;
                                        // blue to red

                                        float col[3];
                                        tf->GetColor(v, col);
                                        glColor3fv(col);

                                    }
                                    /*else if ( mapProperties == 3){
                                        float v = distance / maxDistance;
                                        float d = (density - minD)/(maxD - minD);
                                        float col[3];
                                        tf->GetColor(d, col);
                                        glColor4f(col[0], col[1], col[2],(1.0 - v)/0.75 + 0.25 );
                                    }*/

                                    glVertex3f(startX + blockSize*currentSegment->imageCoorPoints.at(k).p[0],
                                               startY + blockSize*currentSegment->imageCoorPoints.at(k).p[1],
                                               0);
                                }
                                glEnd();
                            }

                        }
                    }
             }
        }


    }



}

void SkeletonVis::paintGL(){
    if ( currentMultiVisMethod == SkeletonVis::SPLOM)
         SkeletonPaint();
    if ( currentMultiVisMethod == SkeletonVis::CategoricalHist)
         CategoricalHistPaint();
    if ( currentMultiVisMethod == SkeletonVis::MDS){
        MDSPaint();
    }
    if ( currentMultiVisMethod == SkeletonVis::Scagnostics){
         ScagnosticsPaint();
    }
}

void SkeletonVis::ScagnosticsPaint(){

    // Basic setup
    glEnable( GL_LINE_SMOOTH );
    glEnable( GL_POLYGON_SMOOTH );
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    glViewport(0, 0, (GLsizei) viewport[2], (GLsizei) viewport[3]);
    glClearColor(0.95, 0.95, 0.95, 0.0);
    glClear(GL_COLOR_BUFFER_BIT);
    bool tex_interp = true;

    // Setup projection matrix
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(0.0, 1.0, 0.0, 1.0);

    // Setup modelview matrix
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glScalef(scale, scale, 1.0);
    glTranslatef(transX, transY, 0.0);
    glBindTexture(GL_TEXTURE_3D, 0);
    glLineWidth(1.0);
    glColor3f(1.0,0, 0);
    // ...........................

    // The was it depends on the measures used...
    //

    if ( _splomGenerator != NULL && _splomGenerator->HasData() && _splomGenerator->IsFinished()){
       if( tex3DCombined != NULL){
           /************************************************/
           int sizeT = _splomGenerator->GetSizeTexture();
           int totalNumerical = _splomGenerator->GetTotalNumericalAttributes();
           int totalImages = (totalNumerical*(totalNumerical-1))/2;

           if ( changeInTexture){
               // glGenTextures(1,&texture3D);
               glBindTexture(GL_TEXTURE_3D, texture3DCombined);
               glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, (tex_interp)?GL_LINEAR:GL_NEAREST);
               glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, (tex_interp)?GL_LINEAR:GL_NEAREST);
               glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
               glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
               glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
               glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
               glTexImage3D(GL_TEXTURE_3D, 0, GL_RGBA,sizeT, sizeT,(totalImages+1)*2, 0, GL_RGBA, GL_UNSIGNED_BYTE, tex3DCombined);
               changeInTexture = false;
           }
           else {
               glBindTexture(GL_TEXTURE_3D, texture3DCombined);
           }
           /************************************************/


            if ( typeOfMeasureMD == SkeletonVis::Wilkinson){
                // For each of the nine measures we show need to order them ...
                // and figure out where to place them if at all...
                float blockSize = 1.0 / 9.0; // The nine measures
                if ( scagnosticsValues == NULL)
                    CreateScagnosticsMeasures();


                // y size depends on the maximum close to the threshold.....
                int maxNumberUnderThreshold = 1;
                for(int measure = 0; measure < 9; measure++){
                    vector< pair<int, float> > currentMeasure = orderedMeasures.at(measure);
                    for(int c = 0; c < static_cast<int>(currentMeasure.size()); c++){
                        if ( currentMeasure.at(c).second > similarityTreshold) continue;
                        if ( c > maxNumberUnderThreshold) maxNumberUnderThreshold = c;
                    }
                }

                // blockSize = 1.0f/maxNumberUnderThreshold;
                float ySize = blockSize;

                float planeSize = (1.0 / (totalImages  + 1.0)) * 0.5;
                float extra = -planeSize*0.5;

                string names[9] = {"Spar","Skw","Out","Strgy","Ste","Clmp","Cnvx","Sky","Mono"};


                // for each measure create the column
                for(int measure = 0; measure < 9; measure++){
                      vector< pair<int, float> > currentMeasure = orderedMeasures.at(measure);
                      float stX = measure*blockSize;
                      renderText(stX, 1.0 , 0, names[measure].c_str());

                      for(int c = 0; c < static_cast<int>(currentMeasure.size()); c++){

                          float s = 1 - similarityTreshold; //if similarity threshold is 0.15, then everything below 0.85
                          if (measure == 0){
                              std::cout << "Getting for image  c " << c<<  "thrs " << s << " " <<  currentMeasure.at(c).second  << std::endl;
                          }
                          if ( currentMeasure.at(c).second >= s){

                              float stY = 1.0 -  c*ySize;
                              int pos = currentMeasure.at(c).first +1 ;
                              glEnable(GL_TEXTURE_3D);

                              glBegin(GL_QUADS);

                              glTexCoord3f(0.0, 1.0,  pos*planeSize + extra ); glVertex2f(stX, stY);
                              glTexCoord3f(1.0, 1.0,  pos*planeSize + extra);  glVertex2f(stX + blockSize, stY);
                              glTexCoord3f(1.0, 0.0,  pos*planeSize + extra);  glVertex2f(stX + blockSize, stY - ySize);
                              glTexCoord3f(0.0, 0.0,  pos*planeSize + extra);  glVertex2f(stX,  stY - ySize);

                              glEnd();

                              glDisable(GL_TEXTURE_3D);
                          }
                      }
                  }
            }
       }
    }

}

void SkeletonVis::MDSPaint(){

    glEnable( GL_LINE_SMOOTH );
    glEnable( GL_POLYGON_SMOOTH );
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    glViewport(0, 0, (GLsizei) viewport[2], (GLsizei) viewport[3]);
    glClearColor(0.95, 0.95, 0.95, 0.0);
    glClear(GL_COLOR_BUFFER_BIT);
    bool tex_interp = true;

    // Setup projection matrix
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(0.0, 1.0, 0.0, 1.0);

    // Setup modelview matrix
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glScalef(scale, scale, 1.0);
    glTranslatef(transX, transY, 0.0);
    glBindTexture(GL_TEXTURE_3D, 0);
    glLineWidth(1.0);
    glColor3f(1.0,0, 0);

    if ( _splomGenerator != NULL){
         if (_splomGenerator->HasData() && _splomGenerator->IsFinished()){

             int totalNumerical = _splomGenerator->GetTotalNumericalAttributes();

             vector<int> list;
             _splomGenerator->GetListOfNumericalAttributes(&list);

             float blockSize = 1.0 / (totalNumerical);

             blockSize *= 0.9;
             // We have now them combined ...
             int totalImages = (totalNumerical*(totalNumerical-1))/2;
             float planeSize = (1.0 / (totalImages  + 1.0)) * 0.5;

             int sizeT = _splomGenerator->GetSizeTexture();
             int pos = 1;

             int posInMDS = 0;
             float extra = -planeSize*0.5;

             //glActiveTexture(GL_TEXTURE0);
            if( tex3DCombined != NULL){

                if (!useMDSCenterlines){
                    for(int i = 0; i < totalNumerical; i++){
                        for(int j = i+1; j < totalNumerical; j++ ){


                            glLineWidth(2.0);
                            glColor3f(0,0,0);
                            LocalPoint p = mdsPoints.at(posInMDS);
                            float stX = p.p[0] - blockSize/2.0f;
                            float stY = p.p[1] - blockSize/2.0f;

                            glBegin(GL_LINE_LOOP);
                               glVertex2f(stX, stY);
                               glVertex2f(stX + blockSize, stY);;
                               glVertex2f(stX + blockSize, stY + blockSize);
                               glVertex2f(stX,  stY + blockSize);
                            glEnd();

                            posInMDS++;
                        }
                    }
                }
                posInMDS = 0;
                 // If there is a change in texture, re-set it ....
                 if ( changeInTexture){
                     // glGenTextures(1,&texture3D);
                     glBindTexture(GL_TEXTURE_3D, texture3DCombined);
                     glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, (tex_interp)?GL_LINEAR:GL_NEAREST);
                     glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, (tex_interp)?GL_LINEAR:GL_NEAREST);
                     glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
                     glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
                     glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
                     glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
                     glTexImage3D(GL_TEXTURE_3D, 0, GL_RGBA,sizeT, sizeT,(totalImages+1)*2, 0, GL_RGBA, GL_UNSIGNED_BYTE, tex3DCombined);
                     changeInTexture = false;
                 }
                 else {
                     glBindTexture(GL_TEXTURE_3D, texture3DCombined);
                 }

                 for(int i = 0; i < totalNumerical; i++){
                     for(int j = i; j < totalNumerical; j++ ){
                         if ( i == j){

                         }
                         else {

                            if (!useMDSCenterlines){
                            // (0,1) (1,1)  (1,0)  (0,0)

                                LocalPoint p = mdsPoints.at(posInMDS);
                                float stX = p.p[0] - blockSize/2.0f;
                                float stY = p.p[1] - blockSize/2.0f;

                                glEnable(GL_TEXTURE_3D);

                                glBegin(GL_QUADS);

                                glTexCoord3f(0.0, 1.0,  pos*planeSize + extra ); glVertex2f(stX, stY);
                                glTexCoord3f(1.0, 1.0,  pos*planeSize + extra);  glVertex2f(stX + blockSize, stY);;
                                glTexCoord3f(1.0, 0.0,  pos*planeSize + extra);  glVertex2f(stX + blockSize, stY + blockSize);
                                glTexCoord3f(0.0, 0.0,  pos*planeSize + extra);  glVertex2f(stX,  stY + blockSize);

                                glEnd();

                                glDisable(GL_TEXTURE_3D);

                                glLineWidth(2.0);
                                glColor3f(0,0,0);

                                glBegin(GL_LINE_LOOP);
                                   glVertex2f(stX, stY);
                                   glVertex2f(stX + blockSize, stY);;
                                   glVertex2f(stX + blockSize, stY + blockSize);
                                   glVertex2f(stX,  stY + blockSize);
                                glEnd();
                                posInMDS++;
                            }
                            pos += 1;
                         }
                     }
                 }

                 pos += 1;

                 // Now the centerlines
                 glColor3f(0.0,1.0,0);


                 // perhaps here  similarCenterlines
                 for(int i = 0; i < totalNumerical; i++){
                     // (1,1)  (0,1)  (0,0) ( 1,0)

                     for(int j = i; j < totalNumerical; j++ ){
                     //for(int j = 0; j < i; j++ ){

                         if ( i == j){

                         }
                         else {

                            if ( useMDSCenterlines){
                                LocalPoint p = mdsPoints.at(posInMDS);
                                float stX = p.p[0] - blockSize/2.0f;
                                float stY = p.p[1] - blockSize/2.0f;
                                glEnable(GL_TEXTURE_3D);
                                glBegin(GL_QUADS);

                                glTexCoord3f(1.0, 0.0,  pos*planeSize +extra);
                                glVertex2f(stX, stY);
                                glTexCoord3f(1.0, 1.0,  pos*planeSize  +extra);
                                glVertex2f(stX +blockSize, stY);
                                glTexCoord3f(0.0, 1.0,  pos*planeSize +extra);
                                glVertex2f(stX + blockSize, stY + blockSize);
                                glTexCoord3f(0.0, 0.0,  pos*planeSize  +extra);
                                glVertex2f(stX, stY + blockSize);

                                glEnd();

                                glDisable(GL_TEXTURE_3D);

                                /*glLineWidth(2.0);
                                glColor3f(0,0,0);

                                glBegin(GL_LINE_LOOP);
                                   glVertex2f(stX, stY);
                                   glVertex2f(stX + blockSize, stY);;
                                   glVertex2f(stX + blockSize, stY + blockSize);
                                   glVertex2f(stX,  stY + blockSize);
                                glEnd();*/

                                posInMDS++;


                            }

                            pos += 1;
                         }
                     }
                 }
             glDisable(GL_TEXTURE_3D);
        }
      }
    }

    if (useMDSCenterlines)
     DisplayMDSCenterlineWithInfo();
}

void SkeletonVis::CategoricalHistPaint(){
/*categorical attributes could be shown using a 2D histogram. The columns of the histogram represent the individ-
  ual categorical dimensions, while the rows represent the individual categories within each category.
  Columns shall be sorted by correlation between dimensions and shall be interactively adjustable.
  Rows shall be sorted by size (number of subjects falling into the categories) or by maximally match-
  ing neighboring columns (in terms of correlation) and, again, interactively adjustable. Number of
  occurences/subjects per category shall be color coded using a suitable transfer function.*/

   if ( _cooccurCalc == NULL) return;
  // Number of columns are the number of dimensions
  int cols = _cooccurCalc->GetNumberOfCategoricalDims();


  // Number of rows are the categorical values within each category
  // also store the categorical size per attribute, that way the size of columns
  // are based on the categorical size
  int sizes[cols];
  int rows = _cooccurCalc->GetTotalSumCategoricalSize(sizes);

  glEnable( GL_LINE_SMOOTH );
  glEnable( GL_POLYGON_SMOOTH );
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  glViewport(0, 0, (GLsizei) viewport[2], (GLsizei) viewport[3]);
  glClearColor(0.95, 0.95, 0.95, 0.0);
  glClear(GL_COLOR_BUFFER_BIT);

  if ( cols <= 1 ) return;
  // Setup projection matrix
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluOrtho2D(0.0, 1.0, 0.0, 1.0);

  // Setup modelview matrix
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  glScalef(scale, scale, 1.0);
  glTranslatef(transX, transY, 0.0);

  // From 0..1
  // I need a little less than 0..1 since I need to put as well something to move
  glColor3f(1.0, 1.0,1.0);

  float rowSize = 0.95 / rows;
  float startingColPos = 0;

  float maxOccurs = _cooccurCalc->GetMaxNumberOccurrences();
  int pos = 0; //, stRow = 0;

  vector<int> indices;
  _cooccurCalc->GetCategoricalIndices(&indices);


  glBegin(GL_QUADS);
  for(int col = 0; col < cols ; col++){
      // Drawing Column ...
      int actualCol = columnOrder.at(col);

      float percent = static_cast<float>(sizes[actualCol]) / static_cast<float>(rows);
      float colSize = percent*0.95;
      float st = startingColPos;
      float end = startingColPos + colSize;

      for(int row = 0; row < rows; row++){
          //  int GetTotalOfCategoricalValue(int attributeIndex, int valueIndex);

          pair<int,int> actualRow = rowOrder.at(row);
          int rowBelongsTo = rowBelonging.at(actualRow.first).first;

          if ( indices.at(actualCol) == rowBelongsTo){
          //if ( stRow <= actualRow.first && actualRow.first < stRow + sizes[actualCol]){
              int totalOccur = actualRow.second;  //_cooccurCalc->GetOccurenceNumberOfAttrAndCategory(pos);
              float v =  static_cast<float>(totalOccur)/ maxOccurs;
              float color[3];
              tfDialog->GetColorAsVector(v, color);

              glColor3fv(color);
              pos++;
          }
          else {
              glColor3f(0.95, 0.95, 0.95);
          }
          glVertex2f( st  +0.05, row*rowSize  +0.05);
          glVertex2f( end +0.05, row*rowSize  +0.05);
          glVertex2f( end +0.05, (row+1)*rowSize+0.05);
          glVertex2f( st  +0.05, (row+1)*rowSize+0.05);
       }

      //stRow += sizes[actualCol];
      startingColPos += colSize;
  }
  glEnd();
  // Now let's get the control circles for the columns

  startingColPos = 0;
  for(int col = 0; col < cols ; col++){
      int actualCol = columnOrder.at(col);

      float percent = static_cast<float>(sizes[actualCol]) / static_cast<float>(rows);
      float colSize = percent*0.95;

       glColor3f(0,0,0);

      glPushMatrix();
         glTranslatef( startingColPos + colSize*0.5 + 0.05, 0.02, 0);
         glutSolidSphere(0.015,10,10);
      glPopMatrix();

      if ( hoverColumn == col){
            glColor3f(1.0,1.0,1.0);
            glPushMatrix();
               glTranslatef( startingColPos + colSize*0.5 + 0.05, 0.02, 0);
               glutSolidSphere(0.012,10,10);
            glPopMatrix();
      }
      startingColPos += colSize;
  }
  // Now let's make the control for the rows
  for(int row = 0; row < rows; row++){

      if ( row == hoverRow)
         glColor3f(1.0, 1.0, 1.0);
      else
         glColor3f(0.6, 0.6, 0.6);

      glBegin(GL_QUADS);
          glVertex2f( 0 , row*rowSize  +0.05);
          glVertex2f( 0.05, row*rowSize  +0.05);
          glVertex2f( 0.05, (row+1)*rowSize+0.05);
          glVertex2f( 0 , (row+1)*rowSize+0.05);
      glEnd();

      glColor3f(0,0,0);
      glBegin(GL_LINE_LOOP);
          glVertex2f( 0 , row*rowSize  +0.05);
          glVertex2f( 0.05, row*rowSize  +0.05);
          glVertex2f( 0.05, (row+1)*rowSize+0.05);
          glVertex2f( 0 , (row+1)*rowSize+0.05);
      glEnd();      
  }


  glBegin(GL_LINE_LOOP);
      glLineWidth(2.0);
      glColor3f(0,0,0);
      glVertex2f( 0.05, 0.05);
      glVertex2f( 1.0 , 0.05);
      glVertex2f( 1.0 , 1.0 );
      glVertex2f( 0.05, 1.0 );
  glEnd();


  // We need to generate a background area for drawing the text

  if ( hoverRow != -1){
      // Display the name of the col/ attr
      pair<int,int> actualRow = rowOrder.at(hoverRow);
      string name = _cooccurCalc->GetNameOfAttrAndCategory(actualRow.first);

      int totalChars = name.length() + 1;

      if ( totalChars <= 2) totalChars += 1;

      float cS = 0.017; // Character spacing
      glColor3f(1.0, 1.0, 1.0);

      float over = 0.038;
      float under = 0.008;
      glBegin(GL_TRIANGLES);
         glVertex3f(0.025, hoverRow*rowSize + 0.05 - under, 0);
         glVertex3f(0.025, hoverRow*rowSize + 0.05 + over,0);
         glVertex3f(0.025 + totalChars*cS, hoverRow*rowSize + 0.05 + over, 0);
         glVertex3f(0.025, hoverRow*rowSize + 0.05 - under, 0);
         glVertex3f(0.025 + totalChars*cS, hoverRow*rowSize + 0.05 - under, 0);
         glVertex3f(0.025 + totalChars*cS, hoverRow*rowSize + 0.05 + over, 0);
      glEnd();

      glColor3f(0,0,0);

      glBegin(GL_LINE_LOOP);
         glVertex3f(0.025, hoverRow*rowSize + 0.05 - under, 0);
         glVertex3f(0.025, hoverRow*rowSize + 0.05 + over,0);
         glVertex3f(0.025 + totalChars*cS, hoverRow*rowSize + 0.05 + over, 0);
         glVertex3f(0.025 + totalChars*cS, hoverRow*rowSize + 0.05 - under, 0);
      glEnd();

      renderText(0.025,hoverRow*rowSize  + 0.05, 0, name.c_str());
  }

  // In the case of columns the starting location may differ
  if ( hoverColumn != -1){
      int actualCol = columnOrder.at(hoverColumn);
      string name = _cooccurCalc->GetNameOfAttr(actualCol);
      int totalChars = name.length() + 1;
      if ( totalChars <= 2) totalChars += 1;
      float cS = 0.017; // Character spacing
      float startX = hoverSphereLoc[0];
      float endX = hoverSphereLoc[0] + totalChars*cS;
      if ( endX > 1.0){
          // It goes over 1.0, so its not completely visible
          float over1 = endX - 1.0;
          startX -= over1;
      }

      glColor3f(1.0, 1.0, 1.0);

      float over = 0.038;
      float under = 0.008;
      glBegin(GL_TRIANGLES);
         glVertex3f(startX,  hoverSphereLoc[1] - under, 0);
         glVertex3f(startX,  hoverSphereLoc[1] + over,0);
         glVertex3f(startX + totalChars*cS,  hoverSphereLoc[1] + over, 0);
         glVertex3f(startX,  hoverSphereLoc[1] - under, 0);
         glVertex3f(startX + totalChars*cS,  hoverSphereLoc[1] - under, 0);
         glVertex3f(startX + totalChars*cS,  hoverSphereLoc[1] + over, 0);
      glEnd();

      glColor3f(0,0,0);

      glBegin(GL_LINE_LOOP);
         glVertex3f(startX,  hoverSphereLoc[1] - under, 0);
         glVertex3f(startX,  hoverSphereLoc[1] + over,0);
         glVertex3f(startX + totalChars*cS,  hoverSphereLoc[1] + over, 0);
         glVertex3f(startX + totalChars*cS,  hoverSphereLoc[1] - under, 0);
      glEnd();
      renderText( startX, hoverSphereLoc[1], 0, name.c_str());

      if ( currentColSelected != -1){
          glColor3f(0,0,0);
          glPushMatrix();
            glTranslatef( hoverSphereLoc[0] , hoverSphereLoc[1],  0);
            glutSolidSphere(0.010,10,10);
          glPopMatrix();
               glColor3f(1.0,1.0,1.0);
               glPushMatrix();
                  glTranslatef( hoverSphereLoc[0] , hoverSphereLoc[1], 0);
                  glutSolidSphere(0.008,10,10);
               glPopMatrix();
      }
  }

  for(unsigned int i = 0; i < brushingArea.size(); i++){
      glColor4f(1.0,1.0,1.0, 0.6);
      glPushMatrix();
         glTranslatef( brushingArea.at(i).p[0],brushingArea.at(i).p[1], 0);
         glutSolidSphere(brushHistSize,10,10);
      glPopMatrix();
  }
}

void SkeletonVis::SetNames(vector<string>* attributes){
    names.clear();
    for(int i = 0; i < static_cast<int>(attributes->size()); i++){
        names.push_back(attributes->at(i));
    }
}

void SkeletonVis::Generate3DCenterlineTextureGrouped(int clusteringVariable, int totalVariables){
    if ( tex3DCombined == NULL) return;

    // Here we dont use that... except the opacity part

    runningCenterline = true;
    // Here we are doing all the variables together...
    int size = _splomGenerator->GetSizeTexture();
    bool isNumerical = _splomGenerator->IsAttributeNumerical(clusteringVariable);
    vector<int> numericalAttributes;
    _splomGenerator->GetListOfNumericalAttributes(&numericalAttributes);
    int num = numericalAttributes.size();
    int totalImages = ((num)*(num-1))/2;

    changeInTexture = true;

    int stPos = size*size*(totalImages +1 )*4;
    int pos = stPos;
    int currentImage = 0;

    float color[3];
    short singlePlane[size*size*3];
    int type = 1;


    std::cout << "Creating centerline texture grouped " << std::endl;
    for(int i = 0; i < num;i++){
        for(int j = i+1; j < num; j++){
            bool runSkeletonAlg = false;
            if ( useSubset){
               // if we use a subset then we need to check whether we actually need it or not...
               // so we look at the centerlines to gen... and the pair should be there
               for(unsigned int k = 0; k < centerlineToGen.size(); k++){
                    if ( centerlineToGen.at(k).first == numericalAttributes.at(i) && centerlineToGen.at(k).second == numericalAttributes.at(j) ){
                         runSkeletonAlg = true;
                         break;
                    }
                }
            }
            else runSkeletonAlg = true;

            if (runSkeletonAlg){
                int startingPos = pos;
                // The complete background...
                GenerateSkeleton(numericalAttributes.at(i), numericalAttributes.at(j),clusteringVariable, -1, type, singlePlane);
                for(int k = 0; k < size*size;k++){
                    // Let's just do the red cover, change it to white... this encompasses all the centerlines ( it should as it has all the data)
                    if ( singlePlane[k*3 + 0] == 255 && singlePlane[k*3 + 1] == 0 &&singlePlane[k*3 + 2] == 0 ){
                        // it is red, so let's do the black border
                        //int opacity = 255;
                        //opacity *= contourOpacity;

                        //if (opacity == 0){

                            tex3DCombined[pos] = 255; pos++;
                            tex3DCombined[pos] = 255; pos++;
                            tex3DCombined[pos] = 255; pos++;
                            tex3DCombined[pos] = 255; pos++;
                        /*}
                        else {

                        }*/
                    }
                    else {
                        tex3DCombined[pos] = 0; pos++;
                        tex3DCombined[pos] = 0; pos++;
                        tex3DCombined[pos] = 0; pos++;
                        tex3DCombined[pos] = 255; pos++;
                    }
                }
                // now for each
                for(int h = 0; h < totalVariables; h++){
                    // Each variable...             SPLOMThread::GetColor(h, totalVariables,1.0, color);

                    if ( isNumerical)
                       mainTf->GetClassColor(h, color);
                    else
                       mainTf->GetColor(h, color);

                    GenerateSkeleton(numericalAttributes.at(i), numericalAttributes.at(j),clusteringVariable, h,type, singlePlane);
                    // if its white, that's skeleton and we aggregate it..
                    int tmpPos = startingPos;
                    for(int k = 0; k < size*size;k++){

                        // it has at least something on, and is a shade of white, the other (up) already took care
                        // of setting to black the background
                        if ( singlePlane[k*3 + 0] != 0 && singlePlane[k*3 + 0]  == singlePlane[k*3 + 1] && singlePlane[k*3 + 0]  == singlePlane[k*3 + 2] ){
                            // It sets the luminance ....TODO.. change the luminance accordingly of the color chosen....
                            tex3DCombined[tmpPos] += color[0]*255;  if ( tex3DCombined[tmpPos] > 255)  tex3DCombined[tmpPos] = 255; tmpPos++;
                            tex3DCombined[tmpPos] += color[1]*255;  if ( tex3DCombined[tmpPos] > 255)  tex3DCombined[tmpPos] = 255; tmpPos++;
                            tex3DCombined[tmpPos] += color[2]*255;  if ( tex3DCombined[tmpPos] > 255)  tex3DCombined[tmpPos] = 255; tmpPos++;
                            tex3DCombined[tmpPos] = 255; tmpPos++;
                        }
                        else  tmpPos += 4;
                     }
                }

                int tmpPos = startingPos;
                for(int k = 0; k < size*size;k++){
                    if ( tex3DCombined[tmpPos] == 0  &&  tex3DCombined[tmpPos+1] == 0 &&  tex3DCombined[tmpPos+2] == 0 ){
                        tex3DCombined[tmpPos] = 255*0.95; tmpPos++;
                        tex3DCombined[tmpPos] = 255*0.95; tmpPos++;
                        tex3DCombined[tmpPos] = 255*0.95; tmpPos++;
                        tex3DCombined[tmpPos] = 255; tmpPos++;
                    }
                    else if ( tex3DCombined[tmpPos] == 255  &&  tex3DCombined[tmpPos+1] == 255 &&  tex3DCombined[tmpPos+2] == 255){
                        tex3DCombined[tmpPos] = 0; tmpPos++;
                        tex3DCombined[tmpPos] = 0; tmpPos++;
                        tex3DCombined[tmpPos] = 0; tmpPos++;
                        tex3DCombined[tmpPos] = 255; tmpPos++;
                    }
                    else tmpPos += 4;
                }

                emit NewSkeletonMade();
            }
            else {
                emit NewSkeletonMade();
                for(int k = 0; k < size*size;k++){
                    tex3DCombined[pos] = 255*0.95; pos++;
                    tex3DCombined[pos] = 255*0.95; pos++;
                    tex3DCombined[pos] = 255*0.95; pos++;
                    tex3DCombined[pos] = 255; pos++;
                }
            }
            currentImage++;
        }
    }

    for(int k = 0; k < size*size; k++){
        tex3DCombined[pos] = 0; pos++;
        tex3DCombined[pos] = 0; pos++;
        tex3DCombined[pos] = 0; pos++;
        tex3DCombined[pos] = 255; pos++;
    }
    runningCenterline = false;
}

void SkeletonVis::Generate3DCenterlineTexture(int clusteringVariable, int filterVariable){
    if ( tex3DCombined == NULL) return;
    runningCenterline = true;

    int size = _splomGenerator->GetSizeTexture();

    vector<int> numericalAttributes;
    _splomGenerator->GetListOfNumericalAttributes(&numericalAttributes);
    int num = numericalAttributes.size();
    int totalImages = ((num)*(num-1))/2;

    changeInTexture = true;

    int stPos = size*size*(totalImages +1 )*4;
    int pos = stPos;
    int currentImage = 0;

    short singlePlane[size*size*3];

    // It is actually simpler... since we have create the images
    // we need to generate them from here

    //std::cout << "Calling to Generate Skeleton? " << clusteringVariable  << " , " << filterVariable << std::endl;
    int type = 1;
    centerlinesGen.clear();
    for(int i = 0; i < num;i++){
        for(int j = i+1; j < num; j++){
             bool runSkeletonAlg = false;
             if ( useSubset){
                // if we use a subset then we need to check whether we actually need it or not...
                // so we look at the centerlines to gen... and the pair should be there
                 for(unsigned int k = 0; k < centerlineToGen.size(); k++){
                     if ( centerlineToGen.at(k).first == numericalAttributes.at(i) && centerlineToGen.at(k).second == numericalAttributes.at(j) ){
                          runSkeletonAlg = true;
                          break;
                     }
                 }
             }
             else runSkeletonAlg = true; // if we are not using a subset, then run the skeleton algorithm


             if (runSkeletonAlg){
                //std::cout << "Creating skeleton ......." << runSkeletonAlg << std::endl;
                GenerateSkeleton(numericalAttributes.at(i), numericalAttributes.at(j),clusteringVariable, filterVariable,type , singlePlane);

                centerlinesGen.insert(make_pair(numericalAttributes.at(i), numericalAttributes.at(j) ));

                emit NewSkeletonMade();


                int nonZero = 0;

                    for(int k = 0; k < size*size;k++){
                        // if its all zeros then is background

                        int opacity = 255;

                        if ( singlePlane[k*3 + 0] == 0 && singlePlane[k*3 + 0] == singlePlane[k*3 +1] && singlePlane[k*3 + 0] == singlePlane[k*3 +2] )
                        {
                            singlePlane[k*3 + 0] = 255*0.95;
                            singlePlane[k*3 + 1] = 255*0.95;
                            singlePlane[k*3 + 2] = 255*0.95;
                            opacity = 0;
                            nonZero++;
                        }
                        //  if its all the same then it is the centerline...
                        else if ( singlePlane[k*3 + 0] != 0 && singlePlane[k*3 + 0] == singlePlane[k*3 +1] && singlePlane[k*3 + 0] == singlePlane[k*3 +2] )
                        {
                            // All together will be 255/255

                            float percent = static_cast<float>(singlePlane[k*3 + 0]) / 255.0f;
                            singlePlane[k*3 + 0] =  colorCenterline.red()*percent + (1.0 - percent)*255*0.95; //255 - singlePlane[k*3 + 0];
                            singlePlane[k*3 + 1] =  colorCenterline.green()*percent+ (1.0 - percent)*255*0.95; //255 - singlePlane[k*3 + 1];
                            singlePlane[k*3 + 2] =  colorCenterline.blue()*percent + (1.0 - percent)*255*0.95;  //255 - singlePlane[k*3 + 2];

                        }

                        // could be border now...

                        else if ( singlePlane[k*3 + 0] > 0 &&  singlePlane[k*3 +1] == 0 &&  singlePlane[k*3 +2] == 0){
                           opacity *= contourOpacity;

                           if (opacity == 0){
                               singlePlane[k*3 + 0] = 255*0.95;
                               singlePlane[k*3 + 1] = 255*0.95;
                               singlePlane[k*3 + 2] = 255*0.95;
                               opacity = 0;
                           }
                           else {
                               // color of contour
                               singlePlane[k*3 + 0] = colorContour.red();
                               singlePlane[k*3 + 1] = colorContour.green();
                               singlePlane[k*3 + 2] = colorContour.blue();

                           }
                        }

                        tex3DCombined[pos] = singlePlane[k*3 + 0]; pos++;
                        tex3DCombined[pos] = singlePlane[k*3 + 1]; pos++;
                        tex3DCombined[pos] = singlePlane[k*3 + 2]; pos++;
                        tex3DCombined[pos] = opacity; pos++;
                     }
             }
             else {
                 // if not, then just fill it up
                 emit NewSkeletonMade();
                 memset(&(tex3DCombined[pos]), 255*0.95, size*size*sizeof(GLubyte)*4);
                 pos += size*size*4;
             }
           currentImage++;
        }
    }

    for(int k = 0; k < size*size; k++){
        tex3DCombined[pos] = 0; pos++;
        tex3DCombined[pos] = 0; pos++;
        tex3DCombined[pos] = 0; pos++;
        tex3DCombined[pos] = 255; pos++;
    }

    runningCenterline = false;

    updateGL();
}

void SkeletonVis::initializeGL(){
    glutInitDisplayMode(GLUT_RGBA | GLUT_ALPHA);
    glViewport(0,0, width(), height());
    glGetIntegerv( GL_VIEWPORT, viewport );
   //std::cout << "width & height "<< width() << " , " << height() << std::endl;
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(0.0, 1.0, 0.0, 1.0);
    glMatrixMode(GL_MODELVIEW);
    glGetIntegerv( GL_VIEWPORT, viewport );
}

void SkeletonVis::resizeGL(int width, int height){
    int side = qMin(width, height);

     glViewport((width - side) / 2, (height - side) / 2, side, side);
     glGetIntegerv( GL_VIEWPORT, viewport );

    //float blockSize = 2.0;
     glMatrixMode(GL_PROJECTION);
     glLoadIdentity();
     gluOrtho2D(0.0, 1.0, 0.0, 1.0);
     glMatrixMode(GL_MODELVIEW);
     updateGL();
}


void SkeletonVis::mouseSkeletonPress(QMouseEvent *event){

    if ( _splomGenerator == NULL ) return;
    if (!_splomGenerator->HasData() ) return;
    oldMouseX = event->x(); oldMouseY = event->y();


    selectingOnSkeletonType = 0;

    if ( event->button() == Qt::RightButton){
         isRightMouseActive = true;
    }
    if ( event->button() == Qt::LeftButton){
         isLeftMouseActive = true;
         selectedCenterlineTree = make_pair(-1,-1);

         similarCenterlines.clear();

         // If we are on the grid, we have the group selection
         // if we are on the centerline's area, then we have a brush ...
         int totalNumerical = _splomGenerator->GetTotalNumericalAttributes();
         float blockSize = 1.0 / (totalNumerical);

         // If it's pressed in the middle we should return without doing anything...

         //******************
         double pX, pY;
         GetWorldCoordinates(event->x(), event->y(), &pX, &pY);

         if ( pX < 0 || pX > 1.0 || pY < 0 || pY > 1.0) {
             glMatrixMode(GL_PROJECTION);
             glLoadIdentity();
             gluOrtho2D(0.0, 1.0, 0.0, 1.0);
             std::cout <<"Outside locs " <<  pX << ", " << pY << std::endl;
             return;
         }
         // Diagonal, gets ignored right now
         for(int i = 0; i < totalNumerical; i++){
                  int x = i;
                  int y = totalNumerical - i - 1;

                  float startX = x*blockSize;
                  float endX = (x+1)*blockSize;
                  float startY = y*blockSize;
                  float endY = (y+1)*blockSize;
                  if (InRange(pX,pY, startX,endX, startY,endY)){
                      return;
                  }
         }

         pY = 1.0 - pY;

         if ( pX > pY){
             pY = 1.0 - pY;
             for(int i = 0; i < totalNumerical; i++){
                 // (1,1)  (0,1)  (0,0) ( 1,0)

                 for(int j = i; j < totalNumerical; j++ ){

                     int x = j;
                     int y = totalNumerical - i- 1;

                     float startX = x*blockSize;
                     float endX = (x+1)*blockSize;
                     float startY = y*blockSize;
                     float endY = (y+1)*blockSize;
                     if (InRange(pX,pY, startX,endX, startY,endY)){

                         // Selecting the grid
                         currentMovingDirectionInGrid[0] = 0;
                         currentMovingDirectionInGrid[1] = 0;
                         currentMovingDirectionInGrid[2] = 0;
                         areaSelectionOutlineInGrid.clear();
                         LocalPoint p;
                         p.p[0] = pX;
                         p.p[1] = pY;
                         areaSelectionOutlineInGrid.push_back(p);
                         areaSelectionInGrid = true;
                         selectedGrid.first = i;
                         selectedGrid.second = j;
                         selectingOnSkeletonType = 1;
                    }
                 }}
         }
         else {
             selectingOnSkeletonType = 2;
             pY = 1.0 - pY;
             for(int i = 0; i < totalNumerical; i++){
                 // (1,1)  (0,1)  (0,0) ( 1,0)

                 for(int j = i; j < totalNumerical; j++ ){
                 //for(int j = 0; j < i; j++ ){
                     int x = i;
                     int y = totalNumerical -j- 1;
                     if ( i == j){

                     }
                     else {

                         brushingArea.clear();
                         float startX = x*blockSize;
                         float endX = (x+1)*blockSize;
                         float startY = y*blockSize;
                         float endY = (y+1)*blockSize;
                         if (InRange(pX,pY, startX,endX, startY,endY)){
                                selectedCenterline.first = i;
                                selectedCenterline.second = j;
                                LocalPoint p;
                                p.p[0] = pX; p.p[1] = pY;

                                brushingArea.push_back(p);
                                updateGL();
                         }
                     }
                 }
             }
         }
    }
}

void SkeletonVis::mouseHistogramPress(QMouseEvent *event){


    if ( event->button() == Qt::RightButton){
         isRightMouseActive = true;
         oldMouseX = event->x(); oldMouseY = event->y();


    }
    if ( event->button() == Qt::LeftButton){
         isLeftMouseActive = true;
         actionOnHistogram = -1;

         double pX, pY;
         GetWorldCoordinates(event->x(), event->y(), &pX, &pY);
         int cols = _cooccurCalc->GetNumberOfCategoricalDims();
         if (cols <= 1) return;

         int sizes[cols];
         int rows = _cooccurCalc->GetTotalSumCategoricalSize(sizes);
         float rowSize = 0.95/rows;
         // Moving row    1
         // Moving column 2
         // or Brushing   3
         if ( pX >= 0.05 && pX <= 1.0 && pY >= 0.05 && pY <= 1.0 )
         {
             // Brushing
             actionOnHistogram = 3;

             brushingArea.clear();

             LocalPoint p;
             p.p[0] = pX;
             p.p[1] = pY;

             brushingArea.push_back(p);
             //vector<LocalPoint> brushingArea;
             //float brushSize;
             updateGL();

         }
         else if ( 0 <= pX && pX <= 0.05 && pY >= 0.05){
            // on row locations..
             actionOnHistogram = 1;

             for(int row = 0; row < rows; row++){
                 if ( row*rowSize +0.05 <= pY && pY < (row+1)*rowSize +0.05){
                     currentRowSelected = row;
                     hoverRow = row;
                     updateGL();
                     break;
                 }
             }
         }
         else if ( pX >= 0.05 && 0 <= pY && pY <= 0.05){
             // may not be on the cols, need to check...

             float startingColPos = 0;
             for(int col = 0; col < cols ; col++){
                 int actualCol = columnOrder.at(col);

                 float percent = static_cast<float>(sizes[actualCol]) / static_cast<float>(rows);
                 float colSize = percent*0.95;

                 float distance = pow(pX - (startingColPos + colSize*0.5 + 0.05), 2.0) + pow(pY -0.02,2.0 );
                 distance = sqrt(distance);

                 if ( distance < 0.015){
                     currentColSelected = col;
                     hoverColumn = col;
                     actionOnHistogram = 2;
                     hoverSphereLoc[0] = pX;
                     hoverSphereLoc[1] = 0.02;
                     updateGL();
                     return;
                 }
                 startingColPos += colSize;
             }
         }
    }
}

void SkeletonVis::mouseScagnosticsPress(QMouseEvent *event){
    if ( event->button() == Qt::RightButton){
         isRightMouseActive = true;
         oldMouseX = event->x(); oldMouseY = event->y();
    }
}

void SkeletonVis::mousePressEvent(QMouseEvent* event){
    switch(currentMultiVisMethod){
        case SkeletonVis::SPLOM:
             mouseSkeletonPress(event);
        break;
        case SkeletonVis::CategoricalHist:
             mouseHistogramPress(event);
        break;
        case SkeletonVis::Scagnostics:
             mouseScagnosticsPress(event);
        break;
        case SkeletonVis::MDS:
             mouseMDSPress(event);
        break;
    }


}

void SkeletonVis::mouseMDSMovement(QMouseEvent *event){

    int x = event->x();
    int y = event->y();

   if (isRightMouseActive) {

        transX += 2.0 * double(x - oldMouseX) / scale / imgSize;
        transY -= 2.0 * double(y - oldMouseY) / scale / imgSize;
        updateGL();
    }
    oldMouseX = x; oldMouseY = y;
}

void SkeletonVis::mouseMDSPress(QMouseEvent *event){

    if ( event->button() == Qt::RightButton){
         isRightMouseActive = true;
         oldMouseX = event->x(); oldMouseY = event->y();
    }
}

void SkeletonVis::mouseMDSRelease(QMouseEvent *event){
    isRightMouseActive = false;
}

void SkeletonVis::mouseDoubleClickEvent(QMouseEvent* event){

    similarCenterlines.clear();

    switch(currentMultiVisMethod){
        case SkeletonVis::SPLOM:
             mouseSkeletonDoubleClick(event);
        break;
        case SkeletonVis::CategoricalHist:
             mouseHistogramDoubleClick(event);
        break;
        case SkeletonVis::Scagnostics:
             mouseScagnosticsDoubleClick(event);
        break;
    }


}

void SkeletonVis::mouseSkeletonDoubleClick(QMouseEvent *event){
    //Let's use this to figure out where and which are available
    if ( _splomGenerator == NULL ) return;
    if (!_splomGenerator->HasData() ) return;
    double pX, pY;
    GetWorldCoordinates(event->x(), event->y(), &pX, &pY);
    // Let's look at detecting whether we are in a grid

    if ( pX < 0 || pX > 1.0 || pY < 0 || pY > 1.0) return;
    int totalNumerical = _splomGenerator->GetTotalNumericalAttributes();
    float blockSize = 1.0 / (totalNumerical);
    // Diagonal, gets ignored right now
    for(int i = 0; i < totalNumerical; i++){
             int x = i;
             int y = totalNumerical - i - 1;
             float startX = x*blockSize;
             float endX = (x+1)*blockSize;
             float startY = y*blockSize;
             float endY = (y+1)*blockSize;
             if (InRange(pX,pY, startX,endX, startY,endY)){
                 return;
             }
    }

    pY = 1.0 - pY;
    if ( pX > pY){

        //MovingOnGrid(pX, 1.0 - pY);
        std::cout << " Double-Clicked on grid " << std::endl;
    }
    else {

        pY = 1.0 - pY;
        selectedCenterlineTree = make_pair(-1,-1);


        vector< pair<int, int> > originalPairs;

        for(int i = 0; i < totalNumerical; i++){
            for(int j = i+1; j < totalNumerical; j++ ){
                originalPairs.push_back( make_pair(i,j));
            }
        }
        int p = 0;

        std::cout << "double click on tree" << std::endl;
        for(int i = 0; i < totalNumerical; i++){
            for(int j = i+1; j < totalNumerical; j++ ){
                int x = i;
                int y = totalNumerical -j- 1;
                float startX = x*blockSize;
                float endX = (x+1)*blockSize;
                float startY = y*blockSize;
                float endY = (y+1)*blockSize;
                if (InRange(pX,pY, startX,endX, startY,endY)){

                    selectedCenterlinePoints.clear();

                    if( reOrderedSplom.empty() ){
                        std::cout << "Clicked at area " << i << " , " << j << std::endl;

                        selectedCenterlineTree.first = i;
                        selectedCenterlineTree.second = j;
                    }
                    else {
                          std::cout << "Clicked at area " << i << " , " << j << std::endl;
                          pair<int,int> actual = originalPairs.at( reOrderedSplom.at(p) );
                          std::cout << "But" << actual.first << " , " << actual.second << std::endl;
                          selectedCenterlineTree.first = actual.first;
                          selectedCenterlineTree.second = actual.second;
                    }


                    if (doubleClickMode == 1 || doubleClickMode == 2){

                        if (typeOfMeasureMD == SkeletonVis::Frechet){
                              ChangeSPLOMOrderingFrechet();
                              if (doubleClickMode == 2)
                                  OrderSkeletonsAccordingToFrechetDistance();
                        }
                        else if ( typeOfMeasureMD == SkeletonVis::Tatu){
                              ChangeSPLOMOrderingTatu();
                              if (doubleClickMode == 2)
                                  OrderSkeletonsAccordingToTatuDistance();

                        }
                        else if ( typeOfMeasureMD == SkeletonVis::Hausdorff){
                              ChangeSPLOMOrderingHausdorff();

                              if (doubleClickMode == 2)
                                    OrderSkeletonsAccordingToHausDorffDistance();

                        }
                        updateGL();
                        break;
                    }
                    else {
                        if (typeOfMeasureMD == SkeletonVis::Frechet){
                              OrderSkeletonsAccordingToFrechetDistance();
                        }
                        else if ( typeOfMeasureMD == SkeletonVis::Tatu){
                                OrderSkeletonsAccordingToTatuDistance();
                        }
                        else if ( typeOfMeasureMD == SkeletonVis::Hausdorff){
                             OrderSkeletonsAccordingToHausDorffDistance();
                        }
                    }

                    updateGL();
                    break;
                }
                p++;
            }
        }
    }
}

void SkeletonVis::OrderSkeletonsAccordingToTatuDistance(){
    if ( selectedCenterlineTree.first == -1) return;
    similarityDistancesFromSelected.clear();
    similarCenterlines.clear();
    vector<int> numericalAttributes;
    _splomGenerator->GetListOfNumericalAttributes(&numericalAttributes);

    int num = numericalAttributes.size();
    int totalNumerical = num;

    vector< pair<int,float> > tatuLoc;

    float currentTatu = 0;

    CreateTatuMeasures();
    int pos =0;
    for(int i = 0; i < num; i++){
        for(int j = i+1; j < num; j++ ){
            double tatu = tatuValues[pos];

            if (selectedCenterlineTree.first == i && selectedCenterlineTree.second == j)
                currentTatu = tatu;

            tatuLoc.push_back( make_pair(pos, tatu) );
            pos++;
        }
    }

    float maxTatuDif = 0;
    float cv = currentTatu;//(currentTatu - minTatu)/(maxTatu - minTatu);
    pos = 0;
    for(int i = 0; i < totalNumerical; i++){
        for(int j = i+1; j < totalNumerical; j++ ){
            float ov =tatuLoc.at(pos).second; //  (tatuLoc.at(pos).second - minTatu)/(maxTatu-minTatu);
            float v = fabs(cv-ov);
            similarityDistancesFromSelected.push_back(v);
            if (v > maxTatuDif) maxTatuDif = v;



            pos++;
        }
    }

    pos = 0;

    for(int i = 0; i < totalNumerical; i++){
        for(int j = i+1; j < totalNumerical; j++ ){

            float v =  similarityDistancesFromSelected.at(pos)/maxTatuDif;
            pos++;
            if ( v < similarityTreshold){
                similarCenterlines.push_back(make_pair(i, j));
            }

        }
    }

}

void SkeletonVis::OrderSkeletonsAccordingToHausDorffDistance(){
    if ( selectedCenterlineTree.first == -1) return;

    similarityDistancesFromSelected.clear();
    similarCenterlines.clear();
    vector<int> numericalAttributes;
    _splomGenerator->GetListOfNumericalAttributes(&numericalAttributes);

    int totalNumerical = numericalAttributes.size();

    HImageType::Pointer dtimg1 = HImageType::New();
    HImageType::Pointer skimg1 = HImageType::New();

    if (modelComparison){
        CreateITKModelImage(dtimg1,true);
        CreateITKModelImage(skimg1,false);
    }
    else {
        CreateITKImage(numericalAttributes.at(selectedCenterlineTree.first), numericalAttributes.at(selectedCenterlineTree.second), dtimg1,true);
        CreateITKImage(numericalAttributes.at(selectedCenterlineTree.first), numericalAttributes.at(selectedCenterlineTree.second), skimg1, false);
    }


    typedef itk::ImageFileWriter< HImageType >  WriterType;
    WriterType::Pointer writer = WriterType::New();
    string key = "model.vtk";
    writer->SetInput(skimg1);
    writer->SetFileName(key);
    writer->Update();



    vector< float> nonNormalizedDistances;
    vector< pair<int, int> > allPairs;
    float maxHausDorffDistance = 0;


    HImageType::Pointer dtimg2 = HImageType::New();
    HImageType::Pointer skimg2 = HImageType::New();


    for(int i = 0; i < totalNumerical; i++){
        for(int j = i+1; j < totalNumerical; j++ ){

            CreateITKImage(numericalAttributes.at(i), numericalAttributes.at(j), dtimg2, true);
            CreateITKImage(numericalAttributes.at(i), numericalAttributes.at(j), skimg2, false);

            double d1 = maxInS(skimg1,dtimg2);
            double d2 = maxInS(skimg2,dtimg1);

            double d = d1;
            if (d2 > d) d2 = d;

            nonNormalizedDistances.push_back(d);
            if ( d > maxHausDorffDistance) maxHausDorffDistance = d;
            allPairs.push_back( make_pair(i,j));
        }
    }
    //

    for(unsigned int i = 0; i < nonNormalizedDistances.size(); i++){
        float v = nonNormalizedDistances.at(i) /  maxHausDorffDistance;
        similarityDistancesFromSelected.push_back(v);
        if ( v < similarityTreshold ){
            similarCenterlines.push_back(allPairs.at(i));
        }
    }

}

void SkeletonVis::OrderSkeletonsAccordingToFrechetDistance(){
   //
   if ( selectedCenterlineTree.first == -1) return;


   float* curve1;
   int n_size;
   float maxFrechet = 0;
   similarityDistancesFromSelected.clear();
   similarCenterlines.clear();

   vector<int> numericalAttributes;
   _splomGenerator->GetListOfNumericalAttributes(&numericalAttributes);
   int totalNumerical = numericalAttributes.size();

   //***********************************************

   if (!modelComparison){

       vector<LocalPoint> skeleton;
       GetLargestCenterlinePoints(numericalAttributes.at(selectedCenterlineTree.first), numericalAttributes.at(selectedCenterlineTree.second), -1,-1, &skeleton);
       //float curve1[skeleton.size()*3];
       curve1 = (float*) malloc(sizeof(float)*skeleton.size()*3);

       for(unsigned int i = 0; i < skeleton.size(); i++){
           curve1[i*3 + 0] = skeleton.at(i).p[0];
           curve1[i*3 + 1] = skeleton.at(i).p[1];
           curve1[i*3 + 2] = 0;
       }
       n_size = skeleton.size();

   }
   //****************************************

   else {

       if ( useFreeFormModel){

           vector<LocalPoint> freeFormPoints;
           if (useFreeFormModel){
                 freeFromDialog->GetCurve(&freeFormPoints);
           }
           curve1 = (float*) malloc(sizeof(float)*freeFormPoints.size()*3);

           for(unsigned int i = 0; i < freeFormPoints.size(); i++){
               curve1[i*3 + 0] = freeFormPoints.at(i).p[0];
               curve1[i*3 + 1] = freeFormPoints.at(i).p[1];
               curve1[i*3 + 2] = 0;
           }

           n_size = freeFormPoints.size();

       }
       else {
           int size = _splomGenerator->GetSizeTexture();
           curve1 = (float*) malloc(sizeof(float)*size*3);
           for(int i = 0; i < size; i++){

               float x = static_cast<float>(i)/static_cast<float>(size);
               float y = a_model*x*x + b_model*x + c_model;
               curve1[i*3 + 0] = x;
               curve1[i*3 + 1] = y;
               curve1[i*3 + 2] = 0;
           }
           n_size = size;
       }

   }

   for(int i = 0; i < totalNumerical; i++){
       for(int j = i+1; j < totalNumerical; j++ ){

           vector<LocalPoint> other;
           GetLargestCenterlinePoints(numericalAttributes.at(i), numericalAttributes.at(j), -1,-1, &other);

           float curve2[other.size()*3];
           float curve3[other.size()*3];

           for(unsigned int k = 0; k < other.size(); k++){
               curve2[k*3 + 0] = other.at(k).p[0];
               curve2[k*3 + 1] = other.at(k).p[1];
               curve2[k*3 + 2] = 0;

               curve3[k*3 + 0] = other.at(other.size() -k-1).p[0];
               curve3[k*3 + 1] = other.at(other.size() -k-1).p[1];
               curve3[k*3 + 2] = 0;

           }

           int n1 = n_size -1;
           int n2 = other.size() -1;

           float frechet = MathHelper::FrechetDistanceBtwTwoLines(curve1, curve2, n_size, other.size(), n1, n2, 0);
           float frechet2 = MathHelper::FrechetDistanceBtwTwoLines(curve1, curve3, n_size, other.size(), n1, n2, 0);
           double v = min(frechet, frechet2);

           if ( v > maxFrechet)
               maxFrechet = v;

           similarityDistancesFromSelected.push_back(v);
       }
   }

   int pos = 0;
   for(int i = 0; i < totalNumerical; i++){
       for(int j = i+1; j < totalNumerical; j++ ){

           float frechet =  similarityDistancesFromSelected.at(pos)/maxFrechet;
           pos++;

           if ( frechet < similarityTreshold){
               similarCenterlines.push_back(make_pair(i, j));
           }
       }
    }

   free(curve1);
}


void SkeletonVis::ChangeSimilarityThreshold(double value){
    similarityTreshold = value;

    similarCenterlines.clear();

    if (typeOfMeasureMD == SkeletonVis::Frechet){
          OrderSkeletonsAccordingToFrechetDistance();
    }
    else if ( typeOfMeasureMD == SkeletonVis::Tatu){
          OrderSkeletonsAccordingToTatuDistance();
    }
    else if ( typeOfMeasureMD == SkeletonVis::Hausdorff){
        OrderSkeletonsAccordingToHausDorffDistance();
    }

    updateGL();
}


void SkeletonVis::CreateMDSProjection(){


    if (typeOfMeasureMD == SkeletonVis::Frechet){
        CreateMDSProjectionFrechet();
    }
    else if ( typeOfMeasureMD == SkeletonVis::Wilkinson){
        CreateMDSProjectionScagnostics();
    }
    else if ( typeOfMeasureMD == SkeletonVis::Tatu){
        CreateMDSProjectionTatu();
    }
    else if ( typeOfMeasureMD == SkeletonVis::Hausdorff){
        CreateMDSProjectionHausdorff();
    }
}

void SkeletonVis::CreateMDSProjectionHausdorff(){

    vector<int> numericalAttributes;
    _splomGenerator->GetListOfNumericalAttributes(&numericalAttributes);
    int num = numericalAttributes.size();
    int totalImages = ((num)*(num-1))/2; // n*(n-1) /2

    vector< pair<int,int> > pairs;
    for(int i = 0; i < num; i++){
        for(int j = i+1; j < num; j++ ){
            pairs.push_back(make_pair(numericalAttributes.at(i),numericalAttributes.at(j)));
        }
    }

    double* similarityMatrix = (double*) malloc(sizeof(double)*totalImages*totalImages);

    HImageType::Pointer dtimg1 = HImageType::New();
    HImageType::Pointer skimg1 = HImageType::New();

    HImageType::Pointer dtimg2 = HImageType::New();
    HImageType::Pointer skimg2 = HImageType::New();


    for(int i = 0; i < totalImages; i++){
        CreateITKImage(pairs.at(i).first, pairs.at(i).second, dtimg1,true);
        CreateITKImage(pairs.at(i).first, pairs.at(i).second, skimg1, false);

          for(int j = i; j< totalImages; j++){
             CreateITKImage(pairs.at(j).first, pairs.at(j).second, dtimg2, true);
             CreateITKImage(pairs.at(i).first, pairs.at(i).second, skimg2, false);

             double d1 = maxInS(skimg1,dtimg2);
             double d2 = maxInS(skimg2,dtimg1);


             double d = d1;
             if (d2 > d) d2 = d;



             similarityMatrix[i*totalImages +j] = d*d;
             similarityMatrix[j*totalImages +i] = d*d;

        }
    }


    std::cout << "Iters? " << std::endl;
    vector< pair<float,float> >  points;
    Statistics::CalculateMDS(similarityMatrix, totalImages,&points);

    mdsPoints.clear();
    for(int i = 0; i < totalImages; i++){
        LocalPoint p;


        p.p[0] = points.at(i).first;
        p.p[1] = points.at(i).second;
        std::cout << i << " .. " << p.p[0] << " , " << p.p[1] << std::endl;
        mdsPoints.push_back(p);
    }

    free(similarityMatrix);
}

double SkeletonVis::maxInS(HImageType::Pointer skeleton, HImageType::Pointer dt){
   //max_{x in S1}( DT2(x)
   //one simply walk over the points of one of the skeletons and measure the DT
   //of the other skeleton there, to find how close you are to the other skeleton

    double maxVal = 0;

    itk::ImageRegionIterator< HImageType>  it(skeleton, skeleton->GetLargestPossibleRegion());
    itk::ImageRegionIterator< HImageType>  dit(dt, dt->GetLargestPossibleRegion());

    while(!it.IsAtEnd()){

        if(it.Get() != itk::NumericTraits<HImageType::PixelType>::Zero){
           float val = dit.Get();
           if ( val >  maxVal){
               maxVal = val;
           }
        }
        ++it;
        ++dit;
    }

    return maxVal;
}

void SkeletonVis::CreateITKModelImage(HImageType::Pointer img, bool DT){
    int size = _splomGenerator->GetSizeTexture();
    HImageType::RegionType region;
    HImageType::IndexType index;
    index[0] = 0;      index[1] = 0;
    HImageType::SizeType sizeT;
    sizeT[0] = size;
    sizeT[1] = size;
    region.SetIndex(index);
    region.SetSize(sizeT);
    img->SetRegions(region);
    img->Allocate(true);
    img->FillBuffer(itk::NumericTraits< HImageType::PixelType >::Zero);
    img->Update();


    if ( useFreeFormModel){
        vector<LocalPoint> freeFormPoints;
        freeFromDialog->GetCurve(&freeFormPoints);

        for(unsigned int i = 0; i < freeFormPoints.size(); i++){
            float x = freeFormPoints.at(i).p[0];
            float y = freeFormPoints.at(i).p[1];

            HImageType::IndexType idx;
            idx[0] = x*size;
            idx[1] = size - y*size;

            img->SetPixel(idx, itk::NumericTraits< HImageType::PixelType >::One);
        }

    }
    else {
        for(int i = 0; i < size; i++){

            float x = static_cast<float>(i)/static_cast<float>(size);
            float y = a_model*x*x + b_model*x + c_model;

            HImageType::IndexType idx;
            idx[0] = x*size;
            idx[1] = size - y*size;

            img->SetPixel(idx, itk::NumericTraits< HImageType::PixelType >::One);
        }


    }


    img->Update();

    if (DT){
        typedef itk::SignedMaurerDistanceMapImageFilter< HImageType, HImageType >  FilterType;
        typename FilterType::Pointer filter = FilterType::New();

        filter->SetInput( img );
        filter->SetBackgroundValue(itk::NumericTraits< HImageType::PixelType >::Zero);
        filter->SetSquaredDistance(false);
        filter->SetUseImageSpacing(false);
        filter->Update();

        HImageType::Pointer tmp = filter->GetOutput();

        itk::ImageRegionIterator< HImageType>  it(img, img->GetLargestPossibleRegion());
        itk::ImageRegionIterator< HImageType>  dit(tmp, tmp->GetLargestPossibleRegion());

        while(!it.IsAtEnd()){

            it.Set(dit.Get());
            ++it;
            ++dit;
        }

    }

}

void SkeletonVis::CreateITKImage(int attr1, int attr2, HImageType::Pointer img, bool DT){

    if ( _splomGenerator->IsRunning()) return;

    HImageType::RegionType region;
    HImageType::IndexType index;
    index[0] = 0;      index[1] = 0;
    HImageType::SizeType sizeT;
    int size = _splomGenerator->GetSizeTexture();

    sizeT[0] = size;
    sizeT[1] = size;
    region.SetIndex(index);
    region.SetSize(sizeT);
    img->SetRegions(region);
    img->Allocate(true);
    img->FillBuffer(itk::NumericTraits< HImageType::PixelType >::Zero);
    img->Update();

    short singlePlane[size*size*3];

    GenerateSkeleton(attr1, attr2, -1, -1, 1 , singlePlane);

    int k = 0;
    itk::ImageRegionIterator< HImageType>  it(img, img->GetLargestPossibleRegion());
    while(!it.IsAtEnd()){

        if ( singlePlane[k*3 + 0] != 0 && singlePlane[k*3 + 0] == singlePlane[k*3 +1] && singlePlane[k*3 + 0] == singlePlane[k*3 +2] ){
            it.Set(itk::NumericTraits< HImageType::PixelType >::One);
        }
        else {
            it.Set(itk::NumericTraits< HImageType::PixelType >::Zero);
        }
        ++it;
        k++;
    }
    img->Update();

    if (DT){
        typedef itk::SignedMaurerDistanceMapImageFilter< HImageType, HImageType >  FilterType;
        typename FilterType::Pointer filter = FilterType::New();

        filter->SetInput( img );
        filter->SetBackgroundValue(itk::NumericTraits< HImageType::PixelType >::Zero);
        filter->SetSquaredDistance(false);
        filter->SetUseImageSpacing(false);
        filter->Update();

        HImageType::Pointer tmp = filter->GetOutput();

        itk::ImageRegionIterator< HImageType>  it(img, img->GetLargestPossibleRegion());
        itk::ImageRegionIterator< HImageType>  dit(tmp, tmp->GetLargestPossibleRegion());

        while(!it.IsAtEnd()){

            it.Set(dit.Get());
            ++it;
            ++dit;
        }

    }
}

void SkeletonVis::CreateMDSProjectionTatu(){


    vector<int> numericalAttributes;
    _splomGenerator->GetListOfNumericalAttributes(&numericalAttributes);
    int num = numericalAttributes.size();
    int totalImages = ((num)*(num-1))/2; // n*(n-1) /2

    WilkinsonScagnostics scagnosticsMeasures;
    scagnosticsMeasures.SetSplomGenerator(_splomGenerator);


    vector< pair<int,float> > tatuLoc;

    CreateTatuMeasures();


    int pos =0;
    for(int i = 0; i < num; i++){
        for(int j = i+1; j < num; j++ ){
            double tatu = tatuValues[pos];
            tatuLoc.push_back( make_pair(pos, tatu) );
            pos++;
        }
    }

    float maxTatu = tatuLoc.at(0).second;
    float minTatu = tatuLoc.at(0).second;

    for(unsigned int k = 0; k < tatuLoc.size(); k++){
        if ( tatuLoc.at(k).second > maxTatu ) maxTatu = tatuLoc.at(k).second;
        if ( tatuLoc.at(k).second < minTatu ) minTatu = tatuLoc.at(k).second;
    }

    double* similarityMatrix = (double*) malloc(sizeof(double)*totalImages*totalImages);

    for(int i = 0; i< totalImages; i++){
        for(int j = i; j< totalImages; j++){

            // now we need to calculate the distance btw the points
            // we normalize everything except the monotonic
            float v1 = (tatuLoc.at(i).second - minTatu)/(maxTatu -minTatu);
            float v2 = (tatuLoc.at(j).second - minTatu)/(maxTatu -minTatu);
            float d = fabs(v1-v2);

            similarityMatrix[i*totalImages +j] = d*d;
            similarityMatrix[j*totalImages +i] = d*d;
        }
    }

    vector< pair<float,float> > points;
    Statistics::CalculateMDS(similarityMatrix, totalImages,&points);

    mdsPoints.clear();
    for(int i = 0; i < totalImages; i++){
        LocalPoint p;
        p.p[0] = points.at(i).first;
        p.p[1] = points.at(i).second;
        mdsPoints.push_back(p);
    }

   free(similarityMatrix);


}

void SkeletonVis::CreateMDSProjectionScagnostics(){


    CreateScagnosticsMeasures();

    vector<int> numericalAttributes;
    _splomGenerator->GetListOfNumericalAttributes(&numericalAttributes);
    int num = numericalAttributes.size();
    int totalImages = ((num)*(num-1))/2; // n*(n-1) /2

    double measures[totalImages][9];
    double maxs[9];
    double mins[9];

    for(int i = 0; i < 9; i++){
        maxs[i] = scagnosticsMaxsMins[i*2 + 0];
        mins[i] = scagnosticsMaxsMins[i*2 + 1];
    }

    int pos = 0;
    for(int j= 0; j < totalImages; j++){

        for(int k = 0; k < 9; k++){
            measures[j][k] = scagnosticsValues[pos];
            pos++;
        }
    }
    // We have 9 measures for scagnostics and their value for all the scatterplots
    double* similarityMatrix = (double*) malloc(sizeof(double)*totalImages*totalImages);

    for(int i = 0; i< totalImages; i++){
        for(int j = i; j< totalImages; j++){

            // now we need to calculate the distance btw the points
            // we normalize everything except the monotonic
            float sum = 0;
            for(int k = 0; k < 8; k++){
                float v1 = (measures[i][k] - mins[k])/(maxs[k] - mins[k]);

                if (maxs[k] - mins[k] < 0.000001) v1 = 0;
                float v2 = (measures[j][k] - mins[k])/(maxs[k] - mins[k]);
                if (maxs[k] - mins[k] < 0.000001) v2 = 0;

                //normalized value
                sum += (v1 -v2)*(v1-v2);
            }

            sum += (measures[i][8]-measures[j][8])*(measures[i][8]-measures[j][8]); // monotonic

            std::cout << sum << " ... " << std::endl;
            // no need to multiply by itself since we are taking a root before and then powering by 2
            // so it results in the same
            similarityMatrix[i*totalImages +j] = sum;
            similarityMatrix[j*totalImages +i] = sum;
        }
    }

    vector< pair<float,float> > points;
    Statistics::CalculateMDS(similarityMatrix, totalImages,&points);

    mdsPoints.clear();
    for(int i = 0; i < totalImages; i++){
        LocalPoint p;
        p.p[0] = points.at(i).first;
        p.p[1] = points.at(i).second;
        mdsPoints.push_back(p);
    }

   free(similarityMatrix);
}

void SkeletonVis::CreateMDSProjectionFrechet(){
    vector<int> numericalAttributes;
    _splomGenerator->GetListOfNumericalAttributes(&numericalAttributes);
    int num = numericalAttributes.size();
    int totalImages = ((num)*(num-1))/2; // n*(n-1) /2

    vector< pair<int,int> > pairs;
    for(int i = 0; i < num; i++){
        for(int j = i+1; j < num; j++ ){
            pairs.push_back(make_pair(numericalAttributes.at(i),numericalAttributes.at(j)));
        }
    }

    double* similarityMatrix = (double*) malloc(sizeof(double)*totalImages*totalImages);

    for(int i = 0; i< totalImages; i++){
         // i = row

        vector<LocalPoint> c1;
        GetLargestCenterlinePoints( pairs.at(i).first, pairs.at(i).second, -1,-1, &c1);

        float curve1[c1.size()*3];
        for(unsigned int k = 0; k < c1.size(); k++){
            curve1[k*3 + 0] = c1.at(k).p[0];
            curve1[k*3 + 1] = c1.at(k).p[1];
            curve1[k*3 + 2] = 0;
        }


         for(int j = i; j< totalImages; j++){
               //j = col
             vector<LocalPoint> c2;
             GetLargestCenterlinePoints( pairs.at(j).first, pairs.at(j).second, -1,-1, &c2);

             float curve2[c2.size()*3];
             float curve3[c2.size()*3];

             for(unsigned int k = 0; k < c2.size(); k++){
                 curve2[k*3 + 0] = c2.at(k).p[0];
                 curve2[k*3 + 1] = c2.at(k).p[1];
                 curve2[k*3 + 2] = 0;

                 curve3[k*3 + 0] = c2.at(c2.size() -k-1).p[0];
                 curve3[k*3 + 1] = c2.at(c2.size() -k-1).p[1];
                 curve3[k*3 + 2] = 0;

             }
             //------------------------------
             int n1 = c1.size() -1;
             int n2 = c2.size() -1;
             double frechet = MathHelper::FrechetDistanceBtwTwoLines(curve1, curve2, n1 +1, n2 +1, n1, n2, 0);
             double frechet2 = MathHelper::FrechetDistanceBtwTwoLines(curve1, curve3, n1 +1, n2 +1, n1, n2, 0);

             double v = min(frechet, frechet2);

             similarityMatrix[i*totalImages +j] = v*v;
             similarityMatrix[j*totalImages +i] = v*v;
         }
    }

    vector< pair<float,float> > points;
    Statistics::CalculateMDS(similarityMatrix, totalImages,&points);
    mdsPoints.clear();

    for(int i = 0; i < totalImages; i++){
        LocalPoint p;
        p.p[0] = points.at(i).first;
        p.p[1] = points.at(i).second;
        mdsPoints.push_back(p);
    }

   free(similarityMatrix);
}

void SkeletonVis::mouseHistogramDoubleClick(QMouseEvent *event){
    double pX, pY;
    GetWorldCoordinates(event->x(), event->y(), &pX, &pY);
    int cols = _cooccurCalc->GetNumberOfCategoricalDims();
    int sizes[cols];
    int rows = _cooccurCalc->GetTotalSumCategoricalSize(sizes);

    if ( pX >= 0.05 && 0 <= pY && pY <= 0.05){
      // Hovering in the col locations

        float startingColPos = 0;
        for(int col = 0; col < cols ; col++){
            int actualCol = columnOrder.at(col);

            float percent = static_cast<float>(sizes[actualCol]) / static_cast<float>(rows);
            float colSize = percent*0.95;

            float distance = pow(pX - (startingColPos + colSize*0.5 + 0.05), 2.0) + pow(pY -0.02,2.0 );
            distance = sqrt(distance);

            if ( distance < 0.015){
                this->ReorderAccordingToColumn(col);
                return;
            }
            startingColPos += colSize;
        }
    }
}

void SkeletonVis::mouseScagnosticsDoubleClick(QMouseEvent *event){

}

void SkeletonVis::mouseSkeletonMovement(QMouseEvent* event){
    int x = event->x();
    int y = event->y();

    if (isLeftMouseActive) {
        double pX, pY;
        GetWorldCoordinates(event->x(), event->y(), &pX, &pY);

        pY = 1.0 - pY;
        if ( pX > pY && selectingOnSkeletonType == 1){
            MovingOnSkeletonGrid(pX, 1.0 - pY);
        }
        else if (selectingOnSkeletonType == 2){
             MovingOnSkeletonCenterline(pX, 1.0 - pY);
        }
    }
    else if (isRightMouseActive) {

        transX += 2.0 * double(x - oldMouseX) / scale / imgSize;
        transY -= 2.0 * double(y - oldMouseY) / scale / imgSize;
        updateGL();

        if (multiResolution )
             MovementSkeleton();
    }
    oldMouseX = x; oldMouseY = y;
}

void SkeletonVis::SwapRowsTopHistogram(){
    int cols = _cooccurCalc->GetNumberOfCategoricalDims();
    int sizes[cols];
    int rows = _cooccurCalc->GetTotalSumCategoricalSize(sizes);
    if (currentRowSelected == rows-1) return; // no-where else to go
    vector< pair<int,int> > newRowOrder;

    if ( currentRowSelected == -1 ) return;
    int n =  rowOrder.size();
    for(int i = 0; i < n; i++){

        if ( i == currentRowSelected){
            newRowOrder.push_back(rowOrder.at(currentRowSelected+1));

        }
        else if ( i == currentRowSelected + 1){
            newRowOrder.push_back(rowOrder.at(currentRowSelected));

        }
        else {
            newRowOrder.push_back( rowOrder.at(i));

        }
    }
    rowOrder.clear();
    for(int i = 0; i < n; i++){
        rowOrder.push_back(newRowOrder.at(i));
    }

    currentRowSelected += 1;
    hoverRow = currentRowSelected;
    updateGL();


}
void SkeletonVis::SwapRowsBottomHistogram(){

    if ( currentRowSelected == 0 ) return; // no-where else to go
    if ( currentRowSelected == -1 ) return;

    vector< pair<int,int> > newRowOrder;

    int n =  rowOrder.size();
    for(int i = 0; i < n; i++){

        if ( i == currentRowSelected){
            newRowOrder.push_back(rowOrder.at(currentRowSelected-1));
        }
        else if ( i == currentRowSelected -1){
             newRowOrder.push_back(rowOrder.at(currentRowSelected));
        }
        else {
             newRowOrder.push_back( rowOrder.at(i));
        }
    }

    rowOrder.clear();
    for(int i = 0; i < n; i++){
        rowOrder.push_back(newRowOrder.at(i));
    }

    currentRowSelected -= 1;
    hoverRow = currentRowSelected;
    updateGL();

}


void SkeletonVis::SwapColsRightHistogram(){

    if (currentColSelected == -1) return;

    vector< int > newColOrder;
    int cols = _cooccurCalc->GetNumberOfCategoricalDims();
    int sizes[cols];
    int rows = _cooccurCalc->GetTotalSumCategoricalSize(sizes);

    int n = columnOrder.size();
    for(int i = 0; i < n; i++){
        if ( i == currentColSelected){
            newColOrder.push_back(columnOrder.at(currentColSelected+1));
        }
        else if ( i == currentColSelected + 1){
            newColOrder.push_back(columnOrder.at(currentColSelected));
        }
        else {
            newColOrder.push_back( columnOrder.at(i));
        }
    }
    columnOrder.clear();
    for(int i = 0; i < n; i++){
        columnOrder.push_back(newColOrder.at(i));
    }

    currentColSelected += 1;
    hoverColumn= currentColSelected;

    float startingColPos = 0;
    for(int col = 0; col < cols ; col++){
        int actualCol = columnOrder.at(col);
        float percent = static_cast<float>(sizes[actualCol]) / static_cast<float>(rows);
        float colSize = percent*0.95;
        if ( col == currentColSelected)
           hoverSphereLoc[0] = (startingColPos + colSize*0.5 + 0.05);
        startingColPos += colSize;
    }
    updateGL();
}

void SkeletonVis::SwapColsLeftHistogram(){

    if (currentColSelected == -1) return;

    vector< int > newColOrder;
    int cols = _cooccurCalc->GetNumberOfCategoricalDims();
    int sizes[cols];
    int rows = _cooccurCalc->GetTotalSumCategoricalSize(sizes);
    int n = columnOrder.size();

    for(int i = 0; i < n; i++){
        if ( i == currentColSelected){
            newColOrder.push_back(columnOrder.at(currentColSelected-1));
        }
        else if ( i == currentColSelected - 1){
            newColOrder.push_back(columnOrder.at(currentColSelected));
        }
        else {
            newColOrder.push_back( columnOrder.at(i));
        }
    }
    columnOrder.clear();
    for(int i = 0; i < n; i++){
        columnOrder.push_back(newColOrder.at(i));
    }

    currentColSelected -= 1;
    hoverColumn= currentColSelected;

    float startingColPos = 0;
    for(int col = 0; col < cols ; col++){
        int actualCol = columnOrder.at(col);
        float percent = static_cast<float>(sizes[actualCol]) / static_cast<float>(rows);
        float colSize = percent*0.95;
        if ( col == currentColSelected)
           hoverSphereLoc[0] = (startingColPos + colSize*0.5 + 0.05);
        startingColPos += colSize;

    }
    updateGL();
}


void SkeletonVis::mouseHistogramMovement(QMouseEvent *event){

    hoverRow = -1; hoverColumn = -1;

    double pX, pY;
    GetWorldCoordinates(event->x(), event->y(), &pX, &pY);
    int cols = _cooccurCalc->GetNumberOfCategoricalDims();
    if (cols <= 1) return;

    int sizes[cols];
    int rows = _cooccurCalc->GetTotalSumCategoricalSize(sizes);
    float rowSize = 0.95/rows;

    if (isLeftMouseActive && actionOnHistogram != -1){
        // Check whether there is a selection of row/order
        // if not then we are brushing instead
        if (actionOnHistogram == 1){
            std::cout << "We are moving the rows around " << std::endl;
             // if the Y distance is over a rowSize we swap
             float currentY = currentRowSelected*rowSize + rowSize*0.5 + 0.05; // Midpoint
             hoverRow = currentRowSelected;
             if ( fabs(pY -currentY) > rowSize ){
                 // Difference is over a rowSize

                 if ( pY > currentY ) {
                     SwapRowsTopHistogram();
                 }
                 else if ( pY < currentY){
                     SwapRowsBottomHistogram();
                 }
             }
        }

        else if ( actionOnHistogram == 2){
            // We would like to move the column
            hoverColumn = currentColSelected;
            hoverSphereLoc[0] = pX;

            float startingColPos = 0;
            for(int col = 0; col < cols ; col++){
                int actualCol = columnOrder.at(col);
                float percent = static_cast<float>(sizes[actualCol]) / static_cast<float>(rows);
                float colSize = percent*0.95;

                // only check to the left & to the right
                float distance = pow(pX - (startingColPos + colSize*0.5 + 0.05), 2.0) + pow(pY -0.02,2.0 );
                distance = sqrt(distance);
                if ( distance < 0.015){

                    if ( col == currentColSelected -1 ){
                         SwapColsLeftHistogram();
                    }
                    else if ( col == currentColSelected +1){
                          SwapColsRightHistogram();
                    }
                }

                startingColPos += colSize;

            }
        }
        else if (actionOnHistogram == 3){
            // brushing action ...

            // Only add if its at least a radius away...
            LocalPoint current;
            current.p[0] = pX; current.p[1] = pY;

            float minDistance = 999;
            for(unsigned int j = 0; j< brushingArea.size(); j++){

                LocalPoint other = brushingArea.at(j);
                float distance = sqrt(pow(pX - other.p[0], 2.0)+pow(pY - other.p[1],2.0));
                if ( distance < minDistance){
                    minDistance = distance;
                }
            }

            if ( minDistance > brushHistSize*0.75){
                brushingArea.push_back(current);
            }
            updateGL();

        }
    }
    else if (isRightMouseActive){
        // Movement of the image
        int x = event->x();
        int y = event->y();
        transX += 2.0 * double(x - oldMouseX) / scale / imgSize;
        transY -= 2.0 * double(y - oldMouseY) / scale / imgSize;
        updateGL();
        oldMouseX = x; oldMouseY = y;

    }
    else {
        // Hovering about the image...
        // Check for hovering location
        double pX, pY;
        GetWorldCoordinates(event->x(), event->y(), &pX, &pY);
        int cols = _cooccurCalc->GetNumberOfCategoricalDims();
        int sizes[cols];
        int rows = _cooccurCalc->GetTotalSumCategoricalSize(sizes);
        float rowSize = 0.95/rows;

        if ( 0 <= pX && pX <= 0.05 && pY >= 0.05){
          // Hovering in the row locations
          for(int row = 0; row < rows; row++){
              if ( row*rowSize +0.05 <= pY && pY < (row+1)*rowSize +0.05){
                  hoverRow = row;
                  updateGL();
                  return;
              }
          }
        }
        else if ( pX >= 0.05 && 0 <= pY && pY <= 0.05){
          // Hovering in the col locations

            float startingColPos = 0;
            for(int col = 0; col < cols ; col++){
                int actualCol = columnOrder.at(col);

                float percent = static_cast<float>(sizes[actualCol]) / static_cast<float>(rows);
                float colSize = percent*0.95;

                float distance = pow(pX - (startingColPos + colSize*0.5 + 0.05), 2.0) + pow(pY -0.02,2.0 );
                distance = sqrt(distance);

                if ( distance < 0.015){
                    hoverColumn = col;
                    hoverSphereLoc[0] = pX;
                    hoverSphereLoc[1] = 0.02;
                    updateGL();
                    return;
                }
                startingColPos += colSize;
            }
        }
    }

    updateGL();

}
void SkeletonVis::mouseScagnosticsMovement(QMouseEvent *event){
    int x = event->x();
    int y = event->y();

   if (isRightMouseActive) {

        transX += 2.0 * double(x - oldMouseX) / scale / imgSize;
        transY -= 2.0 * double(y - oldMouseY) / scale / imgSize;
        updateGL();
    }
    oldMouseX = x; oldMouseY = y;
}

void SkeletonVis::mouseMoveEvent(QMouseEvent *event){


    switch(currentMultiVisMethod){
        case SkeletonVis::SPLOM:
             mouseSkeletonMovement(event);
        break;
        case SkeletonVis::CategoricalHist:
             mouseHistogramMovement(event);
        break;
        case SkeletonVis::Scagnostics:
             mouseScagnosticsMovement(event);
        break;
        case SkeletonVis::MDS:
             mouseMDSMovement(event);
        break;
    }
}

void SkeletonVis::MovingOnSkeletonGrid(double pX, double pY){
    bool finished = (_splomGenerator != NULL && _splomGenerator->IsFinished());
    if (!finished) return;


    int totalNumerical = _splomGenerator->GetTotalNumericalAttributes();
    float blockSize = 1.0 / (totalNumerical);
    for(int i = 0; i < totalNumerical; i++){
        for(int j = i; j < totalNumerical; j++ ){

            int x = j;
            int y = totalNumerical - i- 1;

            float startX = x*blockSize;
            float endX = (x+1)*blockSize;
            float startY = y*blockSize;
            float endY = (y+1)*blockSize;
            if (InRange(pX,pY, startX,endX, startY,endY)){

                LocalPoint lastPoint = areaSelectionOutlineInGrid.at(areaSelectionOutlineInGrid.size() -1 );
                float newDir[3] = {static_cast<float>(pX - lastPoint.p[0]), static_cast<float>( pY - lastPoint.p[1]) , 0};
                float size = sqrt(pow(newDir[0],2.0) + pow(newDir[1],2.0) );
                newDir[0] /= size;
                newDir[1] /= size;

                bool run = false;


                if ( selectedGrid.first == -1){
                    selectedGrid.first = i;
                    selectedGrid.second = j;
                    run = true;
                }
                else if ( selectedGrid.first == i && selectedGrid.second == j){
                    run = true;
                }

                if (run){

                    if ( currentMovingDirectionInGrid[0] == 0 &&  currentMovingDirectionInGrid[1] == 0 ){
                        // if there has been no movement before, we add put the current direction ...
                        currentMovingDirectionInGrid[0] = newDir[0];
                        currentMovingDirectionInGrid[1] = newDir[1];
                        // and forget about it...
                    }
                    else {
                        float angleChange = MathHelper::AngleBwVectors(currentMovingDirectionInGrid, newDir,true);

                        if  (angleChange > 10){
                            LocalPoint p;
                            p.p[0] = pX;
                            p.p[1] = pY;
                            areaSelectionOutlineInGrid.push_back(p);
                            currentMovingDirectionInGrid[0] = newDir[0];
                            currentMovingDirectionInGrid[1] = newDir[1];
                            updateGL();
                        }

                    }
                }
           }
        }
    }

    //
}

void SkeletonVis::MovingOnSkeletonCenterline(float pX, float pY){

    int totalNumerical = _splomGenerator->GetTotalNumericalAttributes();
    float blockSize = 1.0 / (totalNumerical);

    for(int i = 0; i < totalNumerical; i++){
        // (1,1)  (0,1)  (0,0) ( 1,0)

        for(int j = i; j < totalNumerical; j++ ){
        //for(int j = 0; j < i; j++ ){
            int x = i;
            int y = totalNumerical -j- 1;
            if ( i == j){

            }
            else {

                float startX = x*blockSize;
                float endX = (x+1)*blockSize;
                float startY = y*blockSize;
                float endY = (y+1)*blockSize;
                if (InRange(pX,pY, startX,endX, startY,endY)){

                    if ( i == selectedCenterline.first && j == selectedCenterline.second){
                       LocalPoint p;
                       p.p[0] = pX; p.p[1] = pY;
                       brushingArea.push_back(p);
                       updateGL();
                    }
                }
            }
        }
    }


}

void SkeletonVis::ZoomSkeleton(){
    bool work = (_splomGenerator->HasData() && _splomGenerator->IsFinished());
    if (!work) return;
    if (runningCenterline) return;
    if ( tex3DCombined == NULL ) return;


    double startX, startY, endX, endY;

    glGetIntegerv( GL_VIEWPORT, viewport );

    //******************+
    //
    /*   The screen goes with the coordinates as
     *  what the clicks give me
     *  0,1  -- 1 ,1
     *   |        |
     *  0,0   -- 1,0
     *
     * But in the opengl fashion
     * the 0,0 is the upper corner
     * 0,0  -- ...... windowX, 0       --
     *  |
     * 0, windowY ... windowX, windowY  --
     * /////////////////////*/


    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(0.0, 1.0, 0.0, 1.0);

    // Setup modelview matrix
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glScalef(scale, scale, 1.0);
    glTranslatef(transX, transY, 0.0);

    GetWorldCoordinates(0,viewport[3]-0,&startX, &startY);
    GetWorldCoordinates(viewport[2],viewport[3]-viewport[3],&endX, &endY);

    // if either of the for boxes groups is in the range we are looking into..
    //

    set<pair<int,int> > tmpGrid;
    set<pair<int,int> > tmpCenterline;

    gridsToGen.clear();
    centerlineToGen.clear();
    // Let's calculate which attributes I am watching
    // or which combinations? ....

    QTime t;
    t.start();
   if(_splomGenerator != NULL && _splomGenerator->IsFinished()){
       //When there is nothing left
       //
       int totalNumerical = _splomGenerator->GetTotalNumericalAttributes();
       float blockSize = 1.0 / (totalNumerical);
       vector<int> numericalAttributes;
       _splomGenerator->GetListOfNumericalAttributes(&numericalAttributes);
       int totalInGrid = 0;
       int totalInCenterline = 0;

       // Get the number of different rows and columns...
       set<int> columns;
       set<int> rows;

       // Let's see how many grids...
       for(int i = 0; i < totalNumerical; i++){
           for(int j = i; j < totalNumerical; j++ ){

               int x = j;
               int y = totalNumerical - i- 1;
               int adds[4][2] = { {0,0},{1,0},{0,1},{1,1}};

               if (i==j) continue;
               // Regretfully the relationship is not quite that nice
               // The centerline need the grid ones to be generated
               // but the grid ones have no need for the centerlines...
               // so we need to check that as well...
               for(int k = 0; k < 4; k++)
                   if ( InRange( (x+adds[k][0])*blockSize,(y+adds[k][1])*blockSize ,
                                 startX, endX, startY, endY)){
                       totalInGrid++;

                       pair<int, int> p = make_pair(numericalAttributes.at(i),numericalAttributes.at(j));
                       AddToVector(&gridsToGen, p);
                       tmpGrid.insert(p);
                       rows.insert(i);
                       columns.insert(j);
                       break;
                   }

               x = i;
               y = totalNumerical -j- 1;

               // So total in grid, are set of the ones in grid and the ones in centerlines
               // the ones in centerline are only the ones in centerline ....
               for(int k = 0; k < 4; k++)
                   if ( InRange( (x+adds[k][0])*blockSize,(y+adds[k][1])*blockSize ,
                                 startX, endX, startY, endY)){
                       totalInCenterline++;
                       pair<int, int> p = make_pair(numericalAttributes.at(i),numericalAttributes.at(j));

                       AddToVector(&gridsToGen, p);
                       AddToVector(&centerlineToGen, p);
                       tmpCenterline.insert(p);
                       rows.insert(j);
                       columns.insert(i);
                       break;
                   }
           }
       }

       int maxNum = max(rows.size(), columns.size());

       if ( maxNum == 0){

           return;
       }
       int res = min(viewport[2],viewport[3]);
       int blockRes = res/maxNum;

       int desiredSize = skelft2DSize(blockRes, blockRes);


       if ( desiredSize < 64 ) desiredSize = 64;

       if ( desiredSize !=  latestDesiredSize){
           latestDesiredSize = desiredSize;
           emit ChangeResolution();
           return;
       }
   }
}

void SkeletonVis::mouseSkeletonRelease(QMouseEvent *event){
    if ( event->button() == Qt::RightButton){
         isRightMouseActive = false;
    }
    if ( event->button() == Qt::LeftButton){
        // We have area selection interaction
        isLeftMouseActive = false;
        if ( areaSelectionInGrid){
            areaSelectionInGrid = false;
            GetPointsInGrid();
            updateGL();
        }
        else if (brushingArea.size() > 0) {
            // We have brushing interaction
            // need to draw the centerlines ...
            int totalNumerical = _splomGenerator->GetTotalNumericalAttributes();

            int i = selectedCenterline.first;
            int j = selectedCenterline.second;
            int x = i;
            int y = totalNumerical -j- 1;
            float blockSize = 1.0 / (totalNumerical);
            int size = _splomGenerator->GetSizeTexture();

            float startX = x*blockSize;
            float startY = y*blockSize;

            // This define the area of the current Centerline image...

            vector<LocalPoint> skeleton;
            vector<SkeletonSegment*> currentTree;

            vector<int> numericalAttributes;
            _splomGenerator->GetListOfNumericalAttributes(&numericalAttributes);

            //GetCenterlinePoints(numericalAttributes.at(i), numericalAttributes.at(j), -1,-1, &skeleton);
            GetSimplifiedCenterlinePoints(numericalAttributes.at(i), numericalAttributes.at(j), -1,-1, &skeleton);
            double maxD, minD, maxDistance;
            GetSkeletonSegments(numericalAttributes.at(i), numericalAttributes.at(j), -1,-1, &currentTree, &maxD, &minD, &maxDistance);

            if (currentTree.empty()) return; // There are no elements

            LocalPoint firstTangent = currentTree.at(0)->tangents.at(0);
            LocalPoint firstPoint = currentTree.at(0)->imageCoorPoints.at(0);

            double rightAngle = GetRightAngle(firstTangent, firstPoint); // either -90 or 90

            std::cout << "Released " << numericalAttributes.at(i) << " , " << numericalAttributes.at(j) << ":: " << rightAngle << std::endl;

            selectedCenterlinePoints.clear();
            double axis[3];
            axis[0] = 0.0; axis[1] = 0; axis[2] = 1.0;

            for(unsigned int idx = 0; idx < skeleton.size(); idx++){
                //
                LocalPoint inImage = skeleton.at(idx);
                // The image points are not exactly the same, there is a rotation done to the texture so it works
                //
                double dv[3] = {0,0,0};
                for(unsigned int tree = 0; tree < currentTree.size(); tree++){
                    for(unsigned int h = 0; h < currentTree.at(tree)->pointsInSegment.size(); h++ ){
                        if ( currentTree.at(tree)->pointsInSegment.at(h).second == static_cast<int>(idx)){
                              dv[0] = currentTree.at(tree)->tangents.at(h).p[0];
                              dv[1] = currentTree.at(tree)->tangents.at(h).p[1];
                        }
                    }
                }
                double orthoP[3];
                MathHelper::Normalize(dv);
                MathHelper::RotatePointAroundAxis(90, axis, dv, orthoP);
                MathHelper::Normalize(orthoP);

                float inImageX = inImage.p[1] / size;
                float inImageY = 1.0 - inImage.p[0] / size;

                LocalPoint inGrid;
                inGrid.p[0] = startX + inImageX*blockSize;
                inGrid.p[1] = startY + inImageY*blockSize;

                LocalPoint inGridR;

                inGridR.p[0] = inGrid.p[0] + 0.05*orthoP[0]*blockSize;
                inGridR.p[1] = inGrid.p[1] + 0.05*orthoP[1]*blockSize;

                for(unsigned int k = 0; k < brushingArea.size(); k++){
                    LocalPoint p = brushingArea.at(k);
                    float distance = sqrt(pow(p.p[0]-inGrid.p[0],2.0) + pow(p.p[1]- inGrid.p[1], 2.0));

                    if ( distance <brushSkeletonSize){
                        //selectedCenterlinePoints.push_back(inGrid);
                        selectedCenterlinePoints.push_back(make_pair(inGrid,inGridR));
                    }
                }
            }

            std::cout << ".... " << std::endl;
            updateGL();
        }
    }
}
void SkeletonVis::mouseHistogramRelease(QMouseEvent *event){
    currentRowSelected = -1; currentColSelected = -1;


    if ( actionOnHistogram == 3){
        // The list
        if ( highlightedPoints != NULL )
            free(highlightedPoints);



        int cols = _cooccurCalc->GetNumberOfCategoricalDims();

        if (cols <= 1) return;
        Dataset* data = _splomGenerator->GetData();
        Statistics* stats = _splomGenerator->GetStats();

        highlightedPoints = (bool*) malloc(sizeof(bool)*data->GetTotalNumberOfElements() );
        memset(highlightedPoints, false, sizeof(bool)*data->GetTotalNumberOfElements() ); // initialize everything to false
        vector<int> indices;
        _cooccurCalc->GetCategoricalIndices(&indices);

        int totalPoints = 0;
        int sizes[cols];
        int rows = _cooccurCalc->GetTotalSumCategoricalSize(sizes);
        // For every cell in the histogram, let's check if it falls in any brush circle
        // if so, highlight all the samples that belong there....
        float rowSize = 0.95 / rows;
        float startingColPos = 0;

        for(int col = 0; col < cols ; col++){
            int actualCol = columnOrder.at(col);

            float percent = static_cast<float>(sizes[actualCol]) / static_cast<float>(rows);
            float colSize = percent*0.95;
            float st = startingColPos;
            float end = startingColPos + colSize;

            for(int row = 0; row < rows; row++){
                pair<int,int> actualRow = rowOrder.at(row);
                int rowBelongsTo = rowBelonging.at(actualRow.first).first;
                int idxTo = rowBelonging.at(actualRow.first).second;

                if ( indices.at(actualCol) == rowBelongsTo){

                     float midX = ((st+0.05) + (end+0.05))/2.0;
                     float midY = ((row*rowSize  +0.05) + ((row+1)*rowSize+0.05))/2.0;

                     for(unsigned int k = 0; k < brushingArea.size(); k++){
                         LocalPoint p = brushingArea.at(k);
                         float distance = sqrt(pow(p.p[0]-midX,2.0) + pow(p.p[1]- midY, 2.0));
                         if ( distance < brushHistSize){
                             vector<int> indices;
                             Attribute attr = stats->GetAttr(rowBelongsTo);
                             string valueS = attr.GetDims().at(idxTo).GetValue();
                             QString tmp = QString::fromStdString(valueS);
                             int value = tmp.toInt();
                             data->FilterDatasetByCategoricalValue(rowBelongsTo, value, &indices );
                             for(unsigned int h = 0; h < indices.size(); h++){
                                 highlightedPoints[indices.at(h)] = true;
                             }
                             totalPoints++;
                             break;
                         }
                     }
                }
            }

            startingColPos += colSize;
        }

         if (totalPoints > 0)
            emit SamplesOnHistogramSelected();

         if ( totalPoints == 0)
             brushingArea.clear();
        actionOnHistogram = -1;
        updateGL();
    }

    if ( event->button() == Qt::RightButton){
         isRightMouseActive = false;
    }
    if ( event->button() == Qt::LeftButton){
        // We have area selection interaction
        isLeftMouseActive = false;
    }
}

void SkeletonVis::mouseScagnosticsRelease(QMouseEvent *event){
    isRightMouseActive = false;

}


void SkeletonVis::MovementSkeleton(){
    bool work = (_splomGenerator->HasData() && _splomGenerator->IsFinished());
    if (!work) return;
    if (runningCenterline) return;
    if ( tex3DCombined == NULL ) return;
    double startX, startY, endX, endY;

    glGetIntegerv( GL_VIEWPORT, viewport );

    //******************+
    //
    /*   The screen goes with the coordinates as
     *  what the clicks give me
     *  0,1  -- 1 ,1
     *   |        |
     *  0,0   -- 1,0
     *
     * But in the opengl fashion
     * the 0,0 is the upper corner
     * 0,0  -- ...... windowX, 0       --
     *  |
     * 0, windowY ... windowX, windowY  --
     * /////////////////////*/


    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(0.0, 1.0, 0.0, 1.0);

    // Setup modelview matrix
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glScalef(scale, scale, 1.0);
    glTranslatef(transX, transY, 0.0);

    GetWorldCoordinates(0,viewport[3]-0,&startX, &startY);
    GetWorldCoordinates(viewport[2],viewport[3]-viewport[3],&endX, &endY);

    // if either of the for boxes groups is in the range we are looking into..
    //

    set<pair<int,int> > tmpGrid;
    set<pair<int,int> > tmpCenterline;


    set<pair<int,int> > previouslyGenGrid;
    set<pair<int,int> > previouslyGenCenterlines;

    vector<pair<int,int> >::iterator it;
    for( it = centerlineToGen.begin(); it != centerlineToGen.end()  ; ++it ){
        previouslyGenCenterlines.insert(*it);
    }

    for( it = gridsToGen.begin(); it != gridsToGen.end()  ; it++ ){
        previouslyGenGrid.insert(*it);
    }

    gridsToGen.clear();
    centerlineToGen.clear();
    // Let's calculate which attributes I am watching
    // or which combinations? ....

    QTime t;
    t.start();
   if(_splomGenerator != NULL && _splomGenerator->IsFinished()){
       //When there is nothing left
       //
       int totalNumerical = _splomGenerator->GetTotalNumericalAttributes();
       float blockSize = 1.0 / (totalNumerical);
       vector<int> numericalAttributes;
       _splomGenerator->GetListOfNumericalAttributes(&numericalAttributes);
       int totalInGrid = 0;
       int totalInCenterline = 0;

       // Get the number of different rows and columns...
       set<int> columns;
       set<int> rows;

       // Let's see how many grids...
       for(int i = 0; i < totalNumerical; i++){
           for(int j = i; j < totalNumerical; j++ ){

               int x = j;
               int y = totalNumerical - i- 1;
               int adds[4][2] = { {0,0},{1,0},{0,1},{1,1}};

               if (i==j) continue;
               // Regretfully the relationship is not quite that nice
               // The centerline need the grid ones to be generated
               // but the grid ones have no need for the centerlines...
               // so we need to check that as well...
               for(int k = 0; k < 4; k++)
                   if ( InRange( (x+adds[k][0])*blockSize,(y+adds[k][1])*blockSize ,
                                 startX, endX, startY, endY)){
                       totalInGrid++;

                       pair<int, int> p = make_pair(numericalAttributes.at(i),numericalAttributes.at(j));
                       AddToVector(&gridsToGen, p);
                       tmpGrid.insert(p);
                       rows.insert(i);
                       columns.insert(j);
                       break;
                   }

               x = i;
               y = totalNumerical -j- 1;

               // So total in grid, are set of the ones in grid and the ones in centerlines
               // the ones in centerline are only the ones in centerline ....
               for(int k = 0; k < 4; k++)
                   if ( InRange( (x+adds[k][0])*blockSize,(y+adds[k][1])*blockSize ,
                                 startX, endX, startY, endY)){
                       totalInCenterline++;
                       pair<int, int> p = make_pair(numericalAttributes.at(i),numericalAttributes.at(j));

                       AddToVector(&gridsToGen, p);
                       AddToVector(&centerlineToGen, p);
                       tmpCenterline.insert(p);
                       rows.insert(j);
                       columns.insert(i);
                       break;
                   }
           }
       }
       int maxNum = max(rows.size(), columns.size());
       if ( maxNum == 0){     return;     }
       //*************************************************//


       std::vector< pair<int,int> > diff;
       std::set_difference( previouslyGenCenterlines.begin(),  previouslyGenCenterlines.end(),centerlinesGen.begin(), centerlinesGen.end(),
                            std::inserter(diff,diff.begin()));

       std::vector< pair<int,int> > diff2;
       std::set_difference( centerlinesGen.begin(), centerlinesGen.end(),previouslyGenCenterlines.begin(),  previouslyGenCenterlines.end(),
                            std::inserter(diff2,diff2.begin()));

       if ( diff.size() > 0){
           emit ChangeResolution();
       }
   }

}

void SkeletonVis::PassAsGeneratedGrid(){
    gridsGen.clear();
    for(unsigned int i = 0; i < gridsToGen.size();i++)
        gridsGen.insert(gridsToGen.at(i));

}

void SkeletonVis::AddToVector( vector< pair<int, int> >* vec, pair<int, int> loc){
     bool inside = false;

     for(unsigned int i = 0; i < vec->size(); i++){
         if ( vec->at(i).first == loc.first && vec->at(i).second == loc.second )
             inside = true;
     }

     if (!inside){
         vec->push_back(loc);
     }
}

bool SkeletonVis::InRange(float px, float py, float startX, float endX, float startY, float endY){
   bool inside = false;
   if (  px >= startX && px <= endX){
        if ( py >= startY && py <= endY ){
             inside = true;
        }
  }
  return inside;
}

void SkeletonVis::mouseReleaseEvent(QMouseEvent *event){

    switch(currentMultiVisMethod){
        case SkeletonVis::SPLOM:
             mouseSkeletonRelease(event);
        break;
        case SkeletonVis::CategoricalHist:
             mouseHistogramRelease(event);
        break;
        case SkeletonVis::Scagnostics:
             mouseScagnosticsRelease(event);
        break;
    case SkeletonVis::MDS:
          mouseMDSRelease(event);
        break;
    }
}

void SkeletonVis::GetPointsInGrid(){
   std::cout << "Currently selected grid was " << selectedGrid.first << " , " << selectedGrid.second << std::endl;


   if ( areaSelectionOutlineInGrid.size() < 3) return; // Not enough points for a triangle...
   if ( selectedGrid.first == -1) return;


   // We try to figure out the points.. that were selected ....
   int totalNumerical = _splomGenerator->GetTotalNumericalAttributes();

   float blockSize = 1.0 / (totalNumerical);

   int x = selectedGrid.second;
   int y = totalNumerical - selectedGrid.first - 1;

   /* This acts as the region used  */

   float startX = x*blockSize;
   float endX = (x+1)*blockSize;
   float startY = y*blockSize;
   float endY = (y+1)*blockSize;

   vector<int> list;
   _splomGenerator->GetListOfNumericalAttributes(&list);

   // The list
   if ( highlightedPoints != NULL )
       free(highlightedPoints);

   Dataset* data = _splomGenerator->GetData();
   Statistics* stats = _splomGenerator->GetStats();

   highlightedPoints = (bool*) malloc(sizeof(bool)*data->GetTotalNumberOfElements() );
   memset(highlightedPoints, false, sizeof(bool)*data->GetTotalNumberOfElements() ); // initialize everything to false
   // Now we need to check the number of points selected ...

   int totalPoints = 0;

   int attr1 = list.at( selectedGrid.first);
   int attr2 = list.at( selectedGrid.second);

   float minAttr1 = stats->GetMinimumInAttribute(attr1), maxAttr1 = stats->GetMaximumInAttribute(attr1);
   float range1 = maxAttr1 - minAttr1;

   float minAttr2 = stats->GetMinimumInAttribute(attr2), maxAttr2 = stats->GetMaximumInAttribute(attr2);
   float range2 = maxAttr2 - minAttr2;


   for(int k = 0; k < data->GetTotalNumberOfElements(); k++){

       float v1 = data->GetElementValue(k, attr1);
       float v2 = data->GetElementValue(k, attr2);

       if ( v1 == -99999 || v2 == -99999)
       {
            continue;
       }

       if ( _splomGenerator->GetAttributeToCluster() != -1){
           //TODO if the clustering attribute is also negative...
           int value = data->GetElementValue(k, _splomGenerator->GetAttributeToCluster() );
           if ( value == -99999) continue;
       }

       float mx = (data->GetElementValue(k, attr1) -minAttr1)/range1;
       float my = (data->GetElementValue(k, attr2) -minAttr2)/range2;

       //my = 1.0 - my;

       double p[3];
       p[0] = (endX-startX)*mx + startX;// mx is in range 0..1, we re-range it from startX and endX
       p[1] = (endY-startY)*my + startY;
       p[2] = 0;
       //mx and my are in range 0...1
       for(unsigned int j = 2; j < areaSelectionOutlineInGrid.size() ; j++){
            LocalPoint p1 = areaSelectionOutlineInGrid.at(0);
            LocalPoint p2 = areaSelectionOutlineInGrid.at(j-1);
            LocalPoint p3 = areaSelectionOutlineInGrid.at(j);

            double a[3] = {p1.p[0], p1.p[1], 0};
            double b[3] = {p2.p[0], p2.p[1], 0};
            double c[3] = {p3.p[0], p3.p[1], 0};

            if ( MathHelper::IsInTriangle(p, a,b,c)){
                totalPoints++;
                highlightedPoints[k] = true;

                break;
            }
       }
   }
   std::cout << "Total points selected are " << totalPoints << std::endl;
   if (totalPoints > 0)
        emit SamplesOnGridSelected();

}

void SkeletonVis::Redraw(){
    updateGL();
}

void  SkeletonVis::ToggleDensityContour(bool densityBool){
    useDensityContour = densityBool;
}

void  SkeletonVis::SetDensityContourValue(float value){
    densityContourPercentage = value;
}



void SkeletonVis::GenerateDensityContours(int clusteringVariable, int filterVariable ){
    int size = _splomGenerator->GetSizeTexture();
    short singlePlane[size*size];
    vector<int> numericalAttributes;
    _splomGenerator->GetListOfNumericalAttributes(&numericalAttributes);
    int num = numericalAttributes.size();

    DensityImageType::RegionType region;
    DensityImageType::IndexType start;
    start[0] = 0; start[1] = 0;
    DensityImageType::SizeType sizeI;
    sizeI[0] = size; sizeI[1] = size;
    region.SetIndex(start);
    region.SetSize(sizeI);


    float spacing[2] = {1.0f/static_cast<float>(size),1.0f/static_cast<float>(size) };
    float origin[2] = {0,0};
    DensityImageType::Pointer densityImage = DensityImageType::New();
    densityImage->SetRegions(region);
    densityImage->SetOrigin(origin);
    densityImage->SetSpacing(spacing);
    densityImage->Allocate();

    typedef itk::ImageToVTKImageFilter<DensityImageType> ProbabilityConnector;
    ProbabilityConnector::Pointer conv = ProbabilityConnector::New();
    conv->SetInput(densityImage);
    conv->Update();

    vtkSmartPointer<vtkImageData> imageData = conv->GetOutput();

    contours.clear();

    for(int i = 0; i < num;i++){
        for(int j = i+1; j < num; j++){
            _splomGenerator->GetDensityImage(numericalAttributes.at(i), numericalAttributes.at(j), clusteringVariable,filterVariable, singlePlane);
            //
            int maxDensity = 0;
            short* actualData = (short*)imageData->GetScalarPointer();
            for(int k = 0; k < size*size; k++){
                if ( singlePlane[k] > maxDensity) maxDensity = singlePlane[k];
                actualData[k]= singlePlane[k];
            }
            imageData->Modified();
            //

            float valueToContour = densityContourPercentage*maxDensity;

            vtkSmartPointer<vtkContourFilter> contourFilter = vtkSmartPointer<vtkContourFilter>::New();
            contourFilter->SetInputData(imageData);
            contourFilter->SetNumberOfContours(1);
            contourFilter->SetValue(0, valueToContour);
            contourFilter->UseScalarTreeOn();
            contourFilter->ComputeScalarsOn();
            contourFilter->Update();

            vtkSmartPointer<vtkPolyData> newContour = contourFilter->GetOutput();
            contours.push_back(newContour);
        }
    }
    updateGL();
}

void SkeletonVis::Generate3DTexture(int clusteringVariable, int filterVariable, bool blendingChange){

    int size = _splomGenerator->GetSizeTexture();
    short singlePlane[size*size*3];

    vector<int> numericalAttributes;
    _splomGenerator->GetListOfNumericalAttributes(&numericalAttributes);
    int num = numericalAttributes.size();
    int totalImages = ((num)*(num-1))/2; // n*(n-1) /2

    changeInTexture = true;
    // if we are doing blending change, then the size won't change...

    bool released = true;

    if ( tex3DCombined != NULL){
        glBindTexture(GL_TEXTURE_3D, 0);
        glDeleteTextures(1, &texture3DCombined);
        if (!blendingChange){
            free(tex3DCombined);
        }
        else
            released = false;
    }

    if ( released)
        tex3DCombined = (GLubyte*) malloc(sizeof(GLubyte)*size*size*(totalImages+1)*4*2); // Double for the regular texture, and the centerlines of the same size
    //
    int pos = 0;
    int currentImage = 0;


    //std::cout << "Calling to generate the image? " << clusteringVariable << " .. " << filterVariable  << std::endl;
    //float color[3];
    for(int i = 0; i < num;i++){
        //SPLOMThread::GetColor(i, num, 1.0, color);
        for(int j = i+1; j < num; j++){


            if (_splomGenerator->GetTypeBlend() == 4)
               _splomGenerator->GetImage(numericalAttributes.at(i), numericalAttributes.at(j), clusteringVariable,filterVariable, singlePlane);
            else
                _splomDrawer->GetImage(numericalAttributes.at(i), numericalAttributes.at(j), clusteringVariable,filterVariable, singlePlane);

           // Let's change the value, testing the location of the 3D texture coloring...
            for(int k = 0; k < size*size;k++){
                if ( singlePlane[k*3 + 0] == 0 && singlePlane[k*3 + 0] == singlePlane[k*3 +1] && singlePlane[k*3 + 0] == singlePlane[k*3 +2] )
                {
                    singlePlane[k*3 + 0] = 255*0.95;
                    singlePlane[k*3 + 1] = 255*0.95;
                    singlePlane[k*3 + 2] = 255*0.95;

                }
                else if ( singlePlane[k*3 + 0] == 255 && clusteringVariable == -1 && singlePlane[k*3 + 0] == singlePlane[k*3 +1] && singlePlane[k*3 + 0] == singlePlane[k*3 +2] )
                {
                    singlePlane[k*3 + 0] = 0;
                    singlePlane[k*3 + 1] = 0;
                    singlePlane[k*3 + 2] = 0;
                }


                tex3DCombined[pos] = singlePlane[k*3 + 0]; pos++;
                tex3DCombined[pos] = singlePlane[k*3 + 1]; pos++;
                tex3DCombined[pos] = singlePlane[k*3 + 2]; pos++;
                tex3DCombined[pos] = 255; pos++;

             }
           currentImage++;
        }
    }

    for(int k = 0; k < size*size; k++){
        tex3DCombined[pos] = 0; pos++;
        tex3DCombined[pos] = 0; pos++;
        tex3DCombined[pos] = 0; pos++;
        tex3DCombined[pos] = 255; pos++;
    }


    if ( released){
        // if generating new memory
        memset(&(tex3DCombined[pos]), 255*0.95, size*size*sizeof(GLubyte)*4*(totalImages+1));
    }

}

void  SkeletonVis::wheelEvent(QWheelEvent *event){
    int delta = event->delta();
    if ( delta > 0)
        scale = scale * 1.10;
    else
        scale = scale / 1.10;
    updateGL();
    if ( currentMultiVisMethod == SkeletonVis::SPLOM && multiResolution)
        ZoomSkeleton();
}

void SkeletonVis::GetWorldCoordinates(double x, double y, double *xAtPress, double *yAtPress){
    double x_loc = x;
    double y_loc = y;
    double z;

    GLint viewport[4];
    GLdouble modelview[16];
    GLdouble projection[16];
    GLfloat winX, winY;

    glGetDoublev( GL_MODELVIEW_MATRIX, modelview );
    glGetDoublev( GL_PROJECTION_MATRIX, projection );
    glGetIntegerv( GL_VIEWPORT, viewport );

    double diff = height() - viewport[3];

    y_loc -= diff;

    winX = (float)x_loc;
    winY = (float)viewport[3] - (float)y_loc;

    gluUnProject( winX, winY, 0, modelview, projection, viewport, xAtPress, yAtPress, &z);

}

void SkeletonVis::GenerateSkeleton(int attr1, int attr2, int clusteringAttribute, int valuetoFilter, int type, short texture[]){
    // We need to generate the skeleton images
    int size = _splomGenerator->GetSizeTexture();
    short singlePlane[size*size*3];

    short densityPlane[size*size];

    // check whether it is in the memory... if its not then let's generate it, otherwise just return it
    bool found = false;
    int index = -1;

    stringstream ss;
    ss << (attr1*100 + attr2);
    string key = ss.str();


    float thresInfo[5] = {this->thresStart, this->thresStride, this->thresLimit, static_cast<float>(useDensityContour), densityContourPercentage};
    for(int k = 0; k < (int)centerlineMemory.size(); k++){
        if ( centerlineMemory.at(k)->IsSame(attr1, attr2, clusteringAttribute, valuetoFilter,size, thresInfo)){
            index = k;  found = true;   break;
        }
    }


    if (!found){

        QTime t1;
        t1.start();
        vector<LocalPoint> centerlinePoints;
        _splomGenerator->GetImage(attr1, attr2, clusteringAttribute,valuetoFilter, singlePlane);

        if (useDensityContour){

            _splomGenerator->GetDensityImage(attr1, attr2, clusteringAttribute, valuetoFilter, densityPlane);

            float maxDensity = 0;
            for(int i = 0; i < size*size; i++)
                if ( densityPlane[i] > maxDensity) maxDensity = densityPlane[i];

            for(int i = 0; i < size*size; i++)
            {
                if ( densityPlane[i] < maxDensity*densityContourPercentage){
                    singlePlane[i*3 + 0 ] = 0;
                    singlePlane[i*3 + 1 ] = 0;
                    singlePlane[i*3 + 2 ] = 0;
                }
                else {
                    singlePlane[i*3 + 0 ] = 255;
                    singlePlane[i*3 + 1 ] = 255;
                    singlePlane[i*3 + 2 ] = 255;
                }
            }


        }
        //std::cout << "Image " << t1.elapsed() << std::endl;
        generator.GetSkeleton(singlePlane, size, texture, type, &centerlinePoints,key );
        centerlineMemory.push_back( new Memoization(attr1, attr2, clusteringAttribute, valuetoFilter, size, texture, centerlinePoints,thresInfo));

        vector<SkeletonSegment*> segments;
        vector<LocalPoint> simplifiedCenterlines;
        vector<int> largestCenterline;
        double maxD, minD;
        double maxDistance = 0;

        CreateDensityOnCenterline(attr1, attr2, clusteringAttribute, valuetoFilter, &simplifiedCenterlines,&segments, &maxD, &minD, &maxDistance, &largestCenterline);

        centerlineMemory.at( centerlineMemory.size() -1)->SetSimplification(&segments, &simplifiedCenterlines, &largestCenterline, maxD, minD, maxDistance);
    }
    else {
        centerlineMemory.at(index)->GetSkeleton(texture);
    }
}
bool SkeletonVis::GetCenterlinePoints(int attr1, int attr2, int clusteringAttribute, int valuetoFilter, vector<LocalPoint>* points){
    int size = _splomGenerator->GetSizeTexture();
    bool found = false;
    int index = -1;
    float thresInfo[5] = {this->thresStart, this->thresStride, this->thresLimit, static_cast<float>(useDensityContour), densityContourPercentage};

    for(int k = 0; k < (int)centerlineMemory.size(); k++){
        if ( centerlineMemory.at(k)->IsSame(attr1, attr2, clusteringAttribute, valuetoFilter,size, thresInfo)){
            index = k;
            found = true;
            break;
        }
    }

    if (found){
         centerlineMemory.at(index)->GetCenterlinePoints(points);
    }
    // If not found, then call Generate Skeleton...
    return found;
}
bool SkeletonVis::GetSimplifiedCenterlinePoints(int attr1, int attr2, int clusteringAttribute, int valuetoFilter, vector<LocalPoint>* points){
    int size = _splomGenerator->GetSizeTexture();
    bool found = false;
    int index = -1;

    float thresInfo[5] = {this->thresStart, this->thresStride, this->thresLimit,static_cast<float>(useDensityContour), densityContourPercentage};

    for(int k = 0; k < (int)centerlineMemory.size(); k++){
        if ( centerlineMemory.at(k)->IsSame(attr1, attr2, clusteringAttribute, valuetoFilter,size, thresInfo)){
            index = k;
            found = true;
            break;
        }
    }

    if (found){
         centerlineMemory.at(index)->GetSimplifiedCenterlinePoints(points);
    }
    // If not found, then call Generate Skeleton...
    return found;
}


bool SkeletonVis::GetLargestCenterlinePoints(int attr1, int attr2, int clusteringAttribute, int valuetoFilter, vector<LocalPoint>* points){
    int size = _splomGenerator->GetSizeTexture();
    bool found = false;
    int index = -1;

    float thresInfo[5] = {this->thresStart, this->thresStride, this->thresLimit, static_cast<float>(useDensityContour), densityContourPercentage};

    for(int k = 0; k < (int)centerlineMemory.size(); k++){
        if ( centerlineMemory.at(k)->IsSame(attr1, attr2, clusteringAttribute, valuetoFilter,size, thresInfo)){
            index = k;
            found = true;
            break;
        }
    }

    if (found){
          centerlineMemory.at(index)->GetLargestCenterlinePoints(points);
        // centerlineMemory.at(index)->GetSimplifiedCenterlinePoints(points);
    }
    // If not found, then call Generate Skeleton...
    return found;
}

bool SkeletonVis::GetSkeletonSegments(int attr1, int attr2, int clusteringAttribute, int valuetoFilter, vector<SkeletonSegment*>* segments, double* maxD, double* minD, double *maxDistance){
    int size = _splomGenerator->GetSizeTexture();
    bool found = false;
    int index = -1;

    float thresInfo[5] = {this->thresStart, this->thresStride, this->thresLimit, static_cast<float>(useDensityContour), densityContourPercentage};

    for(int k = 0; k < (int)centerlineMemory.size(); k++){
        if ( centerlineMemory.at(k)->IsSame(attr1, attr2, clusteringAttribute, valuetoFilter,size, thresInfo )){
            index = k;
            found = true;
            break;
        }
    }

    if (found){
         centerlineMemory.at(index)->GetSkeletonSegments(segments, maxD, minD, maxDistance);
    }
    // If not found, then call Generate Skeleton...
    return found;
}

void SkeletonVis::CreateOnePixelWidthSkeleton(vector<LocalPoint>* inputSkeleton, vector<LocalPoint>* outputSkeleton){

    // Create an empty image
    // Create an empty image
    typedef itk::Image<unsigned char, 2> ImageType;
    ImageType::Pointer image = ImageType::New();
    ImageType::IndexType start;
    start.Fill(0);

    ImageType::SizeType size;
    size.Fill(_splomGenerator->GetSizeTexture());

    ImageType::RegionType region(start, size);
    image->SetRegions(region);
    image->Allocate();
    image->FillBuffer(0);
    // load the points
    for(unsigned int i = 0; i < inputSkeleton->size(); i++){

        ImageType::IndexType idx;
        idx[0] = inputSkeleton->at(i).p[0];
        idx[1] = inputSkeleton->at(i).p[1];
        image->SetPixel(idx, itk::NumericTraits<SPLOMImageType::PixelType>::OneValue() );
    }

    //****************************
    // apply filter
    /*typedef itk::BinaryThinningImageFilter < ImageType, ImageType> BinaryThinningImageFilterType;
    BinaryThinningImageFilterType::Pointer binaryFilter = BinaryThinningImageFilterType::New();
    binaryFilter->SetInput(image);
    binaryFilter->Update();
   // get the points

    ImageType::Pointer thinned = binaryFilter->GetOutput(); //laplacianSharpeningImageFilter->GetOutput();*/

    // The image is thinned, but there are still issues with it...
    // There are certain patterns, corners wise to simplify the issues
    ImageType::SizeType radius;
    radius[0] = 1;
    radius[1] = 1;

    // Possible cases, 1 is filled, 0 is unfilled, and 2 is doesn't matter
    int cases[20][9] = {{0,0,0,1,1,0,0,1,2},   {0,0,0,0,1,1,2,1,0},
                        {0,1,2,1,1,0,0,0,0},   {2,1,0,0,1,1,0,0,0},
                        {2,0,2,1,1,1,0,1,0},   {2,1,0,0,1,1,2,1,0},
                        {0,1,0,1,1,1,2,0,2},   {0,1,2,1,1,0,0,1,2},
                        {0,0,0,0,1,1,0,1,0},   {0,1,0,1,1,0,0,0,0},
                        {0,1,0,0,1,1,0,0,0},   {0,0,0,1,1,0,0,1,0},
                        {0,0,1,0,1,1,0,1,0},   {0,0,0,1,1,0,0,1,1},
                        {0,1,0,1,1,0,1,0,0},   {1,1,0,0,1,1,0,0,0},
                        {0,0,0,0,1,1,1,1,0},   {1,0,0,1,1,0,0,1,0},
                        {0,1,0,0,1,1,0,0,1},   {0,1,1,1,1,0,0,0,0}};

    int specialCases[4][9] = { {1,1,1,1,1,0,1,0,0},
                               {1,1,1,0,1,1,0,0,1},
                               {1,0,0,1,1,0,1,1,1},
                               {0,0,1,0,1,1,1,1,1}};

    int fillInCase[4][2] = { {1,3},{1,5} , {3,7},{5,7}};

    itk::NeighborhoodIterator<ImageType> niterator(radius, image, image->GetLargestPossibleRegion());

    int changes = 0;

    while(!niterator.IsAtEnd()){

        if ( niterator.GetCenterPixel() != 0){
           // Center pixel as something...

            for(int tCase = 0; tCase< 20; tCase++){
                int fulfillsProperty = 0;

                for(int i = 0; i < 9; i++){

                    if ( cases[tCase][i] == 2) {fulfillsProperty++; continue; }

                    bool IsInBounds;
                    int val = niterator.GetPixel(i, IsInBounds);
                    if ( IsInBounds){
                        if ( cases[tCase][i] == val)
                            fulfillsProperty++;
                    }
                }

                if ( fulfillsProperty == 9){
                    changes++;
                    niterator.SetCenterPixel(0);
                }
            }


            for(int sCase = 0; sCase< 4; sCase++){
                int fulfillsProperty = 0;

                for(int i = 0; i < 9; i++){

                    if ( cases[sCase][i] == 2) {fulfillsProperty++; continue; }

                    bool IsInBounds;
                    int val = niterator.GetPixel(i, IsInBounds);
                    if ( IsInBounds){
                        if ( specialCases[sCase][i] == val)
                            fulfillsProperty++;
                    }
                }

                if ( fulfillsProperty == 9){
                    changes++;

                    niterator.SetPixel(fillInCase[sCase][0],0);
                    niterator.SetPixel(fillInCase[sCase][1],0);
                    //niterator.SetCenterPixel(0);
                }
            }

        }
        ++niterator;
    }

    //
    itk::ImageRegionIteratorWithIndex< ImageType>  it(image, image->GetLargestPossibleRegion());
    while(!it.IsAtEnd()){
        ImageType::IndexType idx = it.GetIndex();

        if ( it.Get() != 0){
           LocalPoint p;
           p.p[0] = idx[0];
           p.p[1] = idx[1];
           outputSkeleton->push_back(p);
        }
        ++it;
    }

}


void SkeletonVis::MapPoints(vector<LocalPoint>* mappedPoints,int attr1, int attr2, int clusteringAttribute, int valueToFilter){
    Dataset* data = _splomGenerator->GetData();
    Statistics* stats = _splomGenerator->GetStats();

    float minAttr1 = stats->GetMinimumInAttribute(attr1), maxAttr1 = stats->GetMaximumInAttribute(attr1);
    float range1 = maxAttr1 - minAttr1;

    float minAttr2 = stats->GetMinimumInAttribute(attr2), maxAttr2 = stats->GetMaximumInAttribute(attr2);
    float range2 = maxAttr2 - minAttr2;

    int totalElements = data->GetTotalNumberOfElements();
    for(int k = 0; k < totalElements; k++){
        float v1 = data->GetElementValue(k, attr1);
        float v2 = data->GetElementValue(k, attr2);

        if ( v1 == -99999 || v2 == -99999)
        {   continue; }

        if ( clusteringAttribute != -1){
            //TODO if the clustering attribute is also negative...
            int value = data->GetElementValue(k, clusteringAttribute);
            if ( value == -99999) continue;
        }

        if (clusteringAttribute != -1 && data->GetElementValue(k,clusteringAttribute) !=  valueToFilter) continue;

        LocalPoint m;

        double val1 = (data->GetElementValue(k, attr1) -minAttr1)/range1;
        double val2 = (data->GetElementValue(k, attr2) -minAttr2)/range2;

        m.p[0] = val1;
        m.p[1] = val2;
        mappedPoints->push_back(m);
    }
}

void SkeletonVis::GetEndPointsAndBifurcations(vector<LocalPoint>* centerlinePoints, vector<int>* endPoints, vector<int>* bifurcations){

    for(unsigned int i = 0; i < centerlinePoints->size(); i++){

        LocalPoint c = centerlinePoints->at(i);
        int numNeighbors =0;
        for(unsigned int j = 0; j < centerlinePoints->size(); j++){
            if (i ==j) continue;

            int xDistance = abs(c.p[0]- centerlinePoints->at(j).p[0]);
            int yDistance = abs(c.p[1]- centerlinePoints->at(j).p[1]);

            if ( xDistance <= 1 && yDistance <= 1){      numNeighbors++;   }
            if (numNeighbors > 1) break;
        }
        if (numNeighbors == 1)
          endPoints->push_back(i);
        if (numNeighbors > 2)
            bifurcations->push_back(i);
    }

}

void SkeletonVis::CreateParentsAndTangents(vector<SkeletonSegment*>* currentTree, int* totalCenterlinePoints){
    for(unsigned int i = 0; i < currentTree->size(); i++){
             currentTree->at(i)->indexOfSegment = i;
             (*totalCenterlinePoints) += currentTree->at(i)->pointsInSegment.size();
             if ( i == 0)
                 currentTree->at(i)->indexOfParentSegment = -1;
             else {
                 int bifurcationIndex = currentTree->at(i)->parentIdx.at(0);
                 // now I have the id in the centerline of the parent...
                 // check all the previous ones, I should have it already generated
                 for(unsigned int j = 0; j < i; j++){
                     int n = currentTree->at(j)->pointsInSegment.size();
                     if ( currentTree->at(j)->pointsInSegment.at(n-1).second == bifurcationIndex){
                         currentTree->at(i)->indexOfParentSegment = j;
                         break;
                     }
                 }
             }
             // The tangents that we are going to calculate are in the image space ...
             vector<LocalPoint>* imageCoors = &(currentTree->at(i)->imageCoorPoints);
             for(unsigned int j = 0; j < imageCoors->size();j++){
                 double newTangent[3] = {0,0,0};
                 MathHelper::GetTangent(imageCoors, newTangent, j);

                 // if tangent is 0, get the parent, and calculate the tangent according to that
                 if ( MathHelper::Norm(newTangent, 3) == 0){
                     LocalPoint parent;
                     // Get parent
                     int bifurcationIndex = currentTree->at(i)->parentIdx.at(0);
                     for(unsigned int h = 0; h < i; h++){
                         int n = currentTree->at(h)->pointsInSegment.size();
                         if ( currentTree->at(h)->pointsInSegment.at(n-1).second == bifurcationIndex)
                         {
                             parent = currentTree->at(h)->imageCoorPoints.at(n-1);
                             break;
                         }
                     }
                     // Get current
                     LocalPoint currentPoint = imageCoors->at(j);
                     newTangent[0] = currentPoint.p[0] - parent.p[0];
                     newTangent[1] = currentPoint.p[1] - parent.p[1];
                     newTangent[2] = 0;
                     MathHelper::Normalize(newTangent,3);
                 }
                 LocalPoint tangent;
                 tangent.p[0] = newTangent[0];
                 tangent.p[1] = newTangent[1];
                 currentTree->at(i)->tangents.push_back(tangent);
             }
     }

}

void SkeletonVis::CreateTreeDensity(vector<SkeletonSegment*>* currentTree ,double* maxDensity, double* minDensity, double* maxDistance,
                                    int totalCenterlinePoints, int attr1, int attr2, int clusteringAttribute, int valueToFilter){

    // Find which is the appropiate angle for querying...
    LocalPoint firstTangent = currentTree->at(0)->tangents.at(0);
    LocalPoint firstPoint = currentTree->at(0)->imageCoorPoints.at(0);
    double rightAngle = GetRightAngle(firstTangent, firstPoint); //

    // Now we can calculate the density in both directions...
    // For each point in the centerline, we calculate every sampling distance of (?)
    // a gaussian radial kernel which we aggregate points according to their distance
    // for left and right, and then substract the values right from left
    // and color code according to that...

    *maxDensity = -999999;
    *minDensity = 999999;

    //std::cout << "Total centerline points " << totalCenterlinePoints << std::endl;
    double axis[3], orthoP[3];
    axis[0] = 0.0; axis[1] = 0; axis[2] = 1.0;
    double dv[3] = {0,0,0};

    int numDebugPoints = 20;

    float* centerlinePointsV = (float*)malloc(sizeof(float)*totalCenterlinePoints*2);
    float* dir = (float*)malloc(sizeof(float)*totalCenterlinePoints*2);
    float* density = (float*)malloc(sizeof(float)*totalCenterlinePoints);
    float* distanceToBorder = (float*) malloc( sizeof(float)* totalCenterlinePoints);

    float* dbgPoints = (float*)malloc(sizeof(float)*numDebugPoints*2);

    for(int k = 0; k < numDebugPoints*2;k++)
        dbgPoints[k] = -1;

    int dimT = _splomGenerator->GetSizeTexture();
    short singlePlane[dimT*dimT];

    _splomGenerator->GetDensityImage(attr1, attr2, clusteringAttribute,valueToFilter, singlePlane);

    int pos = 0;
    for(unsigned int i = 0; i < currentTree->size(); i++){
        vector<LocalPoint>* imageCoors = &(currentTree->at(i)->imageCoorPoints);
        for(unsigned int j = 0; j < imageCoors->size();j++){
            LocalPoint tangent = currentTree->at(i)->tangents.at(j);
            dv[0] = tangent.p[0];  dv[1] = tangent.p[1];
            MathHelper::RotatePointAroundAxis(rightAngle, axis, dv, orthoP);

            centerlinePointsV[pos] = imageCoors->at(j).p[0];
            dir[pos] = orthoP[0];
            pos++;
            centerlinePointsV[pos] = imageCoors->at(j).p[1];
            dir[pos] = orthoP[1];
            pos++;
        }
    }



    InterpolatedPropertyMapping(dimT, centerlinePointsV, dir, singlePlane, 0.005, totalCenterlinePoints, &density,
                                0, &distanceToBorder, &dbgPoints, numDebugPoints);
    //LinearInterpolatedDensityMapping(dimT, centerlinePointsV, dir, singlePlane, 0.01, totalCenterlinePoints, &density);
    //NonLinearInterpolatedDensityMapping(dimT, centerlinePointsV, dir, singlePlane, 0.03, totalCenterlinePoints, &density);
    pos = 0;
    *maxDistance = 0;

    for(unsigned int i = 0; i < currentTree->size(); i++){
        for(unsigned int j = 0; j < currentTree->at(i)->imageCoorPoints.size();j++){
            currentTree->at(i)->densities.push_back( density[pos]);
            currentTree->at(i)->distances.push_back( distanceToBorder[pos]);
            if (density[pos] > *maxDensity) *maxDensity = density[pos];
            if (density[pos] < *minDensity) *minDensity = density[pos];

            if (distanceToBorder[pos] > *maxDistance) *maxDistance = distanceToBorder[pos];

            pos++;
        }
        for(int k = 0; k < numDebugPoints; k++){

            LocalPoint t;
            t.p[0] = dbgPoints[k*2 + 0 ];
            t.p[1] = dbgPoints[k*2 + 1 ];

            currentTree->at(0)->debugQueryPoints.push_back(t);
        }
    }

    //std::cout << "Max distance? " << *maxDistance << std::endl;
    free(centerlinePointsV);
    free(dir);
    free(density);
    free(dbgPoints);
    free(distanceToBorder);
}



void SkeletonVis::GetLargestContinuousCenterline( vector<SkeletonSegment *> *currentTree, vector<int>* largestCenterline){

    // Let's do a greedy approach,
    // First we get all the end points in the current tree...
    // now all the end points should be "technically" connected to one other
    vector< pair<LocalPoint,int> > points;

    for(unsigned int h = 0; h < currentTree->size(); h++){
        SkeletonSegment* currentSegment = currentTree->at(h);
        for(unsigned int k= 0; k < currentSegment->pointsInSegment.size(); k++){
            pair<LocalPoint, int> p = currentSegment->pointsInSegment.at(k);
            points.push_back(p);
        }
    }
    vector<int> endpoints;



    vtkSmartPointer<vtkMutableDirectedGraph> graph =
         vtkSmartPointer<vtkMutableDirectedGraph>::New();

    vector< vector<int> > allneighbors;
    vtkSmartPointer<vtkPoints> vtkpoints =
       vtkSmartPointer<vtkPoints>::New();

    vtkIdType idsInGraph[points.size()];
    for(unsigned int i = 0; i < points.size(); i++){
      vtkIdType v = graph->AddVertex();
      idsInGraph[i] = v;
      LocalPoint c = points.at(i).first;
      vtkpoints->InsertNextPoint(c.p[0], c.p[1], 0.0);
    }



    for(unsigned int i = 0; i < points.size(); i++){

        int numNeighbors = 0;

        LocalPoint c = points.at(i).first;

        vector<int> neighbors;

        for(unsigned int j = 0; j < points.size(); j++){
            if (i == j) continue;

            LocalPoint o = points.at(j).first;
            int xDistance = abs(c.p[0]- o.p[0]);
            int yDistance = abs(c.p[1]- o.p[1]);

            if ( xDistance <= 1 && yDistance <= 1){
                numNeighbors++;

                graph->AddEdge(  idsInGraph[i], idsInGraph[j] );
                graph->AddEdge(  idsInGraph[j], idsInGraph[i] );

                neighbors.push_back(j);
            }
        }

        allneighbors.push_back(neighbors);
        if ( numNeighbors == 1)
            endpoints.push_back(i);
    }
    //
    graph->SetPoints(vtkpoints);

    float largestSize = -1;

    LocalPoint st(0,0);
    LocalPoint end;


    vtkSmartPointer<vtkGraphToPolyData> graphToPolyData =
       vtkSmartPointer<vtkGraphToPolyData>::New();
     graphToPolyData->SetInputData(graph);
     graphToPolyData->Update();


    for(unsigned int i = 0; i < endpoints.size(); i++){
        for(unsigned int j = i+1; j < endpoints.size(); j++){
            // for all combinations of end points, we start from endpoint i
            // and we follow the neighbors
            //for(int k = 0; k < points.size(); k++) seen[k] = false; // initialize  to false, memset would be faster
            //seen[endpoints.at(i)] = true;
            vector<int> currentCenterline;

            vtkSmartPointer<vtkDijkstraGraphGeodesicPath> dijkstra =
                vtkSmartPointer<vtkDijkstraGraphGeodesicPath>::New();
              dijkstra->SetInputConnection(graphToPolyData->GetOutputPort());
              dijkstra->SetStartVertex( idsInGraph[endpoints.at(i)]);
              dijkstra->SetEndVertex(idsInGraph[endpoints.at(j)]);
              dijkstra->Update();

            vtkSmartPointer<vtkPolyData> res = dijkstra->GetOutput();
            vtkSmartPointer<vtkCellArray> cellArray = vtkSmartPointer<vtkCellArray>::New();
            vtkSmartPointer<vtkPoints> dijkpoints = res->GetPoints();
            vtkSmartPointer<vtkIdList> list = vtkSmartPointer<vtkIdList>::New();

            cellArray = res->GetLines();
            cellArray->GetNextCell(list);


            currentCenterline.clear();

            for(int k = 0; k < list->GetNumberOfIds(); k++){
                double p[2];
                dijkpoints->GetPoint(list->GetId(k),p);


                for(unsigned int h = 0; h < points.size(); h++){
                    LocalPoint o = points.at(h).first;

                    float d = sqrt( pow(p[0] - o.p[0],2.0) +  pow(p[1] - o.p[1],2.0) );
                    if ( d < 0.5){
                        currentCenterline.push_back(h);
                        break;
                    }

                }
            }


            float n = GetGeodesicDistance(&points, &currentCenterline);

            if ( n > largestSize ){
                largestSize = n;
                largestCenterline->clear();
                for(unsigned int k = 0; k < currentCenterline.size(); k++){
                     largestCenterline->push_back(   points.at(currentCenterline.at(k)).second);
                }

                st = points.at(currentCenterline.at(0)).first;
                end = points.at(currentCenterline.back()).first;
            }
        }

    }

    // actually the closest to the origin

    float dSt = sqrt(st.p[0]*st.p[0] + st.p[1]*st.p[1]);
    float dE = sqrt(end.p[0]*end.p[0] + end.p[1]*end.p[1]);
    if ( dE < dSt ){
       reverse(largestCenterline->begin(), largestCenterline->end());
    }


}


double SkeletonVis::GetGeodesicDistance( vector< pair<LocalPoint,int> >* points, vector<int>* ids ){
    vector<LocalPoint> actualPoints;

    for(unsigned int i = 0; i < ids->size(); i++){
        actualPoints.push_back( points->at(ids->at(i)).first);
    }
  float totalDistance = 0;


    for(unsigned int j = 1; j < actualPoints.size(); j++){
         LocalPoint a = actualPoints.at(j-1);
         LocalPoint b = actualPoints.at(j);

         float d = sqrt( pow(a.p[0] -b.p[0],2.0) + pow(a.p[1] - b.p[1],2.0));
         totalDistance += d;
    }


   /* LocalPoint a = actualPoints.at(0);
    LocalPoint b = actualPoints.at(actualPoints.size()-1 );

    float d = sqrt( pow(a.p[0] -b.p[0],2.0) + pow(a.p[1] - b.p[1],2.0));
    totalDistance += d;
   */

    return totalDistance;
}

void SkeletonVis::CreateDensityOnCenterline(int attr1, int attr2, int clusteringAttribute, int valueToFilter, vector<LocalPoint>* simplifiedCenterline,
                                              vector<SkeletonSegment*>* currentTree, double* maxDensity, double* minDensity, double* maxDistance, vector<int>* largestCenterline){
    vector<LocalPoint> tmpCenterlinePoints, centerlinePoints;

    bool found = GetCenterlinePoints(attr1, attr2, clusteringAttribute, valueToFilter,&tmpCenterlinePoints);
    //bool found = GetCenterlinePoints(attr, &tmpCenterlinePoints);

    if (!found) return;

    QTime t1;
    t1.start();

    //std::cout << std::endl;
    CreateOnePixelWidthSkeleton(&tmpCenterlinePoints, &centerlinePoints);
    // We have the centerline points that we want...
    // now we get the points
    // We assume a single centerline i.e. that it is not disconnected
    // Find the lowest centerline point,
    // first we need to find the end points... i.e. find those that have only one neighbor
    vector<int> endPoints, bifurcations;
    bool visited[centerlinePoints.size()]; // For generating the segments

    GetEndPointsAndBifurcations(&centerlinePoints, &endPoints, &bifurcations);
    bool seenEndPoints[endPoints.size()];
    int totalSeenEndPoints = 0;

    for(unsigned int k = 0; k < endPoints.size(); k++){
        seenEndPoints[k] = false;
    }
    for(unsigned int i = 0; i < centerlinePoints.size(); i++){visited[i] = false;         simplifiedCenterline->push_back(centerlinePoints.at(i));}

    int startingPoint = GetStartingPoint(&centerlinePoints, &endPoints, seenEndPoints);
    if ( startingPoint == -1){  return ;   }

    //*****************************************
    vector< vector<SkeletonSegment*>> possibleTrees;
    int idxOfGreatestTree = -1;
    int sizeOfGreatestTree = -1;

    int totEndPoints = endPoints.size();
    while(totalSeenEndPoints != totEndPoints){

        vector<SkeletonSegment*> tmpTree;
        //int sizeofTree = CreateTree(startingPoint, &tmpTree, seenEndPoints, centerlinePoints,endPoints, visited);
        int sizeofTree = CreateTree2(startingPoint, &tmpTree, seenEndPoints, centerlinePoints,endPoints, visited);



        possibleTrees.push_back(tmpTree); // For memory management, we need to delete the skeleton segments after we find
        // the one we keep
        if (sizeofTree > sizeOfGreatestTree ){
            idxOfGreatestTree = possibleTrees.size() -1;
            sizeOfGreatestTree = sizeofTree;
        }
        totalSeenEndPoints = 0;
        for(unsigned int k = 0; k < endPoints.size(); k++){
           if( seenEndPoints[k] ) totalSeenEndPoints++;
        }

        if (totalSeenEndPoints != totEndPoints){
            for(int k = 0; k < totEndPoints; k++){
               if(!seenEndPoints[k]){
                   startingPoint = endPoints.at(k);
                   break;
               }
            }
        }

    }

    //**********************************************************

    for(unsigned int i = 0; i <  possibleTrees.at(idxOfGreatestTree).size(); i++){
        currentTree->push_back(possibleTrees.at(idxOfGreatestTree).at(i));
    }

    // now we have the tree itself, the zeroth segment, at zeroth position
    // let's set the index values, for which are the parents index
    int totalCenterlinePoints = 0;
    CreateParentsAndTangents(currentTree,&totalCenterlinePoints);

    CreateTreeDensity(currentTree, maxDensity, minDensity, maxDistance, totalCenterlinePoints, attr1, attr2, clusteringAttribute, valueToFilter);

    // std::cout << "Attrs " << attr1 << " , " << attr2 << std::endl;
    //std::cout << simplifiedCenterline->size() <<" , " << totalCenterlinePoints << std::endl;

    /*
    GetLargestContinuousCenterline(currentTree, totalCenterlinePoints, largestCenterline );



    vector<int> largestCenterline2;*/
    GetLargestContinuousCenterline( currentTree, largestCenterline);
}

void SkeletonVis::ChangeMultiVisMeasure(int idx){
    typeOfMeasureMD = idx;
    mdsPoints.clear();
    ChangeMultiVisMethod(currentMultiVisMethod);
}

double SkeletonVis::Gaussian2D(LocalPoint center, LocalPoint point, double amplitude, const double stddev2){

    double t1 =  pow(point.p[0] - center.p[0],2)/ (stddev2);
    double t2 =  pow(point.p[1] - center.p[1],2)/ (stddev2);
    double gaussianValue = amplitude*exp(-(t1+t2));
    if ( isnan(gaussianValue)){
        //return 0;
    }
    return gaussianValue;
}

int SkeletonVis::GetStartingPoint(vector<LocalPoint>* centerlinePoints, vector<int>* endPoints, bool seen[]){
    int startingPoint = -1;
    float minDistance = 9999;
    // and then we find the point with the minimum distance
    int size = _splomGenerator->GetSizeTexture();

    for(unsigned int i = 0; i < endPoints->size();i++){

        if (seen[i]) continue; //
        LocalPoint inImage = centerlinePoints->at(endPoints->at(i));
        float inImageX = inImage.p[1] / size;
        // float inImageY = 1.0 - inImage.p[0] / size;
        float distance = inImageX; //sqrt(inImageX*inImageX + inImageY*inImageY);
        if ( distance < minDistance)
        {
            minDistance = distance;
            startingPoint = endPoints->at(i);
        }
    }
    // now we've got the endpoints, & the first starting centerline location...
    return startingPoint;
}

double SkeletonVis::GetRightAngle(LocalPoint firstTangent, LocalPoint firstPoint){
    double axis[3];
    axis[0] = 0.0; axis[1] = 0; axis[2] = 1.0;
    double dv[3] = {firstTangent.p[0], firstTangent.p[1],0};

    double orthoP[3];
    double orthoN[3];

    MathHelper::Normalize(dv);
    MathHelper::RotatePointAroundAxis(90, axis, dv, orthoP);
    MathHelper::RotatePointAroundAxis(-90, axis, dv, orthoN);

    // now orthoP is a orthogonal direction to the vector, positive 90 degrees
    MathHelper::Normalize(orthoP);
    MathHelper::Normalize(orthoN);

    // Need to figure out which one is left & which one is right...
    // We define right as :
    //    if the orthogonal direction is parallel to x axis, y = 0 then the direction that is positive x
    if ( fabs(orthoN[1]) < 0.0001 ){ //

        if (orthoN[0] > 0 ) return -90;
        else if (orthoP[0] > 0 ) return 90;
    }
    //    if the direction is not parallel, then the one that with intersection to the x is positive
    // x = 0
    double t = -firstPoint.p[1]/ orthoP[1];
    if (t > 0 ) return 90;
    else return -90;
}




int SkeletonVis::CreateTree2(int startingPoint, vector<SkeletonSegment*>* currentTree, bool seenEndPoints[],
                             vector<LocalPoint> centerlinePoints, vector<int> endPointsIndices,bool visited[], bool debugFunction){

    int size = _splomGenerator->GetSizeTexture();

    stack< pair<int, int> > queryPoints;
    queryPoints.push( make_pair(startingPoint,-1));

    SkeletonSegment* segment = new SkeletonSegment();

    set<int> createNewSegmentIds;

    while(!queryPoints.empty()){
         // Get the top of the queue
         pair<int, int> cInfo =queryPoints.top();
         int cIdx = cInfo.first;

         for(unsigned int h = 0; h < endPointsIndices.size(); h++){
             if (cIdx == endPointsIndices.at(h)){
                 seenEndPoints[h] = true;
             }
         }
         LocalPoint c = centerlinePoints.at(cIdx);
         //** if the current id is in create new segment ids, then we create a new segment
         //bool afterBifur = false;
         //Ask if I should create a new segment

         if ( createNewSegmentIds.count(cIdx) == 1){
             if (debugFunction) std::cout << "The id of point " << c.p[0] << " , " << c.p[1] << " is for a new segment " << std::endl;

             currentTree->push_back(segment); // save the current segment....
             segment = new SkeletonSegment();
             //afterBifur = true;
         }
         segment->pointsInSegment.push_back(make_pair(c, cIdx));

         LocalPoint imageP;
         imageP.p[0] = c.p[1] / size;
         imageP.p[1] = 1.0 - c.p[0] / size;
         segment->imageCoorPoints.push_back(imageP);

         segment->parentIdx.push_back(cInfo.second);
         visited[cIdx] = true;
         //
         queryPoints.pop();
         // find the neighbors of the current top
         // the non-visited ones..
         vector<int> neighbors;

         for(unsigned int j = 0; j < centerlinePoints.size(); j++){
             if (cIdx == static_cast<int>(j)) continue;

             int xDistance = abs(c.p[0]- centerlinePoints.at(j).p[0]);
             int yDistance = abs(c.p[1]- centerlinePoints.at(j).p[1]);

             if ( xDistance <= 1 && yDistance <= 1 && !visited[j]){
                neighbors.push_back(j);
             }
         }
         if ( neighbors.size() == 1){ // simple non-bifurcation
               queryPoints.push( make_pair( neighbors.at(0), cIdx) );
         }
         else if ( !neighbors.empty()){ // bifurcation
                for(unsigned int  j =0 ; j< neighbors.size(); j++){
                    queryPoints.push( make_pair(neighbors.at(j), cIdx));
                    createNewSegmentIds.insert(neighbors.at(j));
                    visited[neighbors.at(j)] = true;
                }
         }
    }
    currentTree->push_back(segment); // save the current segment....

    int totalSize =0;
    for(unsigned int i = 0; i < currentTree->size(); i++){
          SkeletonSegment* s = currentTree->at(i);
          totalSize += s->pointsInSegment.size();
    }
    return totalSize;
}

void SkeletonVis::ChangeSPLOMOrderingTatu(){

    // We have the selected tree
    // i, j
    similarityDistancesFromSelected.clear();
    similarCenterlines.clear();

    orderAccordingToSimilarity.clear();
    orderAccordingToSimilarity.push_back(selectedCenterlineTree.first);
    orderAccordingToSimilarity.push_back(selectedCenterlineTree.second);

    int num = _splomGenerator->GetTotalNumericalAttributes();

    vector<int> numericalAttributes;
    _splomGenerator->GetListOfNumericalAttributes(&numericalAttributes);
    CreateTatuMeasures();

    float currentTatu = 0;
    int pos =0;
    vector< pair<int,float> > tatuLoc;
    vector< pair<int, int> > originalPairs;

    for(int i = 0; i < num; i++){
        for(int j = i+1; j < num; j++ ){
            double tatu = tatuValues[pos];

            if (selectedCenterlineTree.first == i && selectedCenterlineTree.second == j)
                currentTatu = tatu;

            tatuLoc.push_back( make_pair(pos, tatu) );
            originalPairs.push_back( make_pair(i,j));
            pos++;
        }
    }

    std::cout << "Going to generate SPLOM Ordering based on Tatu " << std::endl;

    while ( static_cast<int>(orderAccordingToSimilarity.size()) != num){

        float minAvgSimilarity = 9999999;
        int newItem = -1;

        for(int k = 0; k < num; k++){
            if ( count(orderAccordingToSimilarity.begin(),orderAccordingToSimilarity.end(), k )) continue;
            // if we already added it, then we have to reason to use it

            // calculate the avg similarity with the already sorted
            float avgSim = 0;

            for(unsigned  int h = 0; h < orderAccordingToSimilarity.size(); h++){
                int attr1 = k;
                int attr2 = orderAccordingToSimilarity.at(h);

                if ( attr1 > attr2){
                    attr1 = orderAccordingToSimilarity.at(h);
                    attr2 = k;
                }
                // ***************+
                // Since we are using Tatu we can search in the original pairs

                for(unsigned int pairs = 0; pairs < originalPairs.size(); pairs++){

                    if ( originalPairs.at(pairs).first == attr1 && originalPairs.at(pairs).second){
                        avgSim += fabs(tatuLoc.at(pairs).second- currentTatu);
                    }
                }
            }
            avgSim /= originalPairs.size();
            if ( avgSim < minAvgSimilarity){
                minAvgSimilarity = avgSim;
                newItem = k;
            }
        }
        orderAccordingToSimilarity.push_back(newItem);
    }


    reOrderedSplom.clear();



    for(unsigned int i = 0; i < orderAccordingToSimilarity.size(); i++){
        for(unsigned int j = i+1; j < orderAccordingToSimilarity.size(); j++){

             int attr1 = orderAccordingToSimilarity.at(i);
             int attr2 = orderAccordingToSimilarity.at(j);

             if ( attr1 > attr2){
                attr2 = orderAccordingToSimilarity.at(i);
                attr1 = orderAccordingToSimilarity.at(j);
             }

             for(unsigned int k = 0; k < originalPairs.size(); k++){
                 if ( originalPairs.at(k).first == attr1 && originalPairs.at(k).second == attr2){
                     reOrderedSplom.push_back(k);
                     break;
                 }
             }
        }
    }



}

void SkeletonVis::ChangeSPLOMOrderingHausdorff(){
    // We have the selected tree
    // i, j
    similarityDistancesFromSelected.clear();
    similarCenterlines.clear();

    orderAccordingToSimilarity.clear();
    orderAccordingToSimilarity.push_back(selectedCenterlineTree.first);
    orderAccordingToSimilarity.push_back(selectedCenterlineTree.second);


    vector<int> numericalAttributes;
    _splomGenerator->GetListOfNumericalAttributes(&numericalAttributes);

    HImageType::Pointer dtimg1 = HImageType::New();
    HImageType::Pointer skimg1 = HImageType::New();

    CreateITKImage(numericalAttributes.at(selectedCenterlineTree.first), numericalAttributes.at(selectedCenterlineTree.second), dtimg1,true);
    CreateITKImage(numericalAttributes.at(selectedCenterlineTree.first), numericalAttributes.at(selectedCenterlineTree.second), skimg1, false);

    int num = _splomGenerator->GetTotalNumericalAttributes();

    int pos =0;
    vector< pair<int,float> > tatuLoc;
    vector< pair<int, int> > originalPairs;

    for(int i = 0; i < num; i++){
        for(int j = i+1; j < num; j++ ){

            HImageType::Pointer dtimg2 = HImageType::New();
            HImageType::Pointer skimg2 = HImageType::New();

            CreateITKImage(numericalAttributes.at(i), numericalAttributes.at(j), dtimg2, true);
            CreateITKImage(numericalAttributes.at(i), numericalAttributes.at(j), skimg2, false);

            double d1 = maxInS(skimg1,dtimg2);
            double d2 = maxInS(skimg2,dtimg1);

            double d = d1;
            if (d2 > d) d2 = d;


            tatuLoc.push_back( make_pair(pos, d) );
            originalPairs.push_back( make_pair(i,j));
            pos++;
        }
    }

    std::cout << "Going to generate SPLOM Ordering based on Hausdorff " << std::endl;

    while (static_cast<int>(orderAccordingToSimilarity.size()) != num){

        float minAvgSimilarity = 9999999;
        int newItem = -1;

        for(int k = 0; k < num; k++){
            if ( count(orderAccordingToSimilarity.begin(),orderAccordingToSimilarity.end(), k )) continue;
            // if we already added it, then we have to reason to use it

            // calculate the avg similarity with the already sorted
            float avgSim = 0;

            for(unsigned int h = 0; h < orderAccordingToSimilarity.size(); h++){
                int attr1 = k;
                int attr2 = orderAccordingToSimilarity.at(h);

                if ( attr1 > attr2){
                    attr1 = orderAccordingToSimilarity.at(h);
                    attr2 = k;
                }
                // ***************+
                // Since we are using Tatu we can search in the original pairs

                for(unsigned int pairs = 0; pairs < originalPairs.size(); pairs++){

                    if ( originalPairs.at(pairs).first == attr1 && originalPairs.at(pairs).second){
                        avgSim += fabs(tatuLoc.at(pairs).second);
                    }
                }
            }
            avgSim /= originalPairs.size();
            if ( avgSim < minAvgSimilarity){
                minAvgSimilarity = avgSim;
                newItem = k;
            }
        }
        orderAccordingToSimilarity.push_back(newItem);
    }


    reOrderedSplom.clear();



    for(unsigned int i = 0; i < orderAccordingToSimilarity.size(); i++){
        for(unsigned int j = i+1; j < orderAccordingToSimilarity.size(); j++){

             int attr1 = orderAccordingToSimilarity.at(i);
             int attr2 = orderAccordingToSimilarity.at(j);

             if ( attr1 > attr2){
                attr2 = orderAccordingToSimilarity.at(i);
                attr1 = orderAccordingToSimilarity.at(j);
             }

             for(unsigned int k = 0; k < originalPairs.size(); k++){
                 if ( originalPairs.at(k).first == attr1 && originalPairs.at(k).second == attr2){
                     reOrderedSplom.push_back(k);
                     break;
                 }
             }
        }
    }

}
void SkeletonVis::ChangeSPLOMOrderingFrechet(){

    similarityDistancesFromSelected.clear();
    similarCenterlines.clear();

    orderAccordingToSimilarity.clear();
    orderAccordingToSimilarity.push_back(selectedCenterlineTree.first);
    orderAccordingToSimilarity.push_back(selectedCenterlineTree.second);


    vector<int> numericalAttributes;
    _splomGenerator->GetListOfNumericalAttributes(&numericalAttributes);


    float* curve1;
    int n_size;
    float maxFrechet = 0;
    similarityDistancesFromSelected.clear();
    similarCenterlines.clear();


    int totalNumerical = numericalAttributes.size();
    int num = totalNumerical;
    //***********************************************

    if (!modelComparison){

        vector<LocalPoint> skeleton;
        GetLargestCenterlinePoints(numericalAttributes.at(selectedCenterlineTree.first), numericalAttributes.at(selectedCenterlineTree.second), -1,-1, &skeleton);
        //float curve1[skeleton.size()*3];
        curve1 = (float*) malloc(sizeof(float)*skeleton.size()*3);

        for(unsigned int i = 0; i < skeleton.size(); i++){
            curve1[i*3 + 0] = skeleton.at(i).p[0];
            curve1[i*3 + 1] = skeleton.at(i).p[1];
            curve1[i*3 + 2] = 0;
        }
        n_size = skeleton.size();

    }
    //****************************************

    else {

        if ( useFreeFormModel){

            vector<LocalPoint> freeFormPoints;
            if (useFreeFormModel){
                  freeFromDialog->GetCurve(&freeFormPoints);
            }
            curve1 = (float*) malloc(sizeof(float)*freeFormPoints.size()*3);

            for(unsigned int i = 0; i < freeFormPoints.size(); i++){
                curve1[i*3 + 0] = freeFormPoints.at(i).p[0];
                curve1[i*3 + 1] = freeFormPoints.at(i).p[1];
                curve1[i*3 + 2] = 0;
            }

            n_size = freeFormPoints.size();

        }
        else {
            int size = _splomGenerator->GetSizeTexture();
            curve1 = (float*) malloc(sizeof(float)*size*3);
            for(int i = 0; i < size; i++){

                float x = static_cast<float>(i)/static_cast<float>(size);
                float y = a_model*x*x + b_model*x + c_model;
                curve1[i*3 + 0] = x;
                curve1[i*3 + 1] = y;
                curve1[i*3 + 2] = 0;
            }
            n_size = size;
        }

    }




    int pos =0;
    vector< pair<int,float> > tatuLoc;
    vector< pair<int, int> > originalPairs;

    for(int i = 0; i < num; i++){
        for(int j = i+1; j < num; j++ ){

            vector<LocalPoint> other;
            GetLargestCenterlinePoints(numericalAttributes.at(i), numericalAttributes.at(j), -1,-1, &other);

            float curve2[other.size()*3];
            float curve3[other.size()*3];

            for(unsigned int k = 0; k < other.size(); k++){
                curve2[k*3 + 0] = other.at(k).p[0];
                curve2[k*3 + 1] = other.at(k).p[1];
                curve2[k*3 + 2] = 0;

                curve3[k*3 + 0] = other.at(other.size() -k-1).p[0];
                curve3[k*3 + 1] = other.at(other.size() -k-1).p[1];
                curve3[k*3 + 2] = 0;

            }

            int n1 = n_size -1;
            int n2 = other.size() -1;

            float frechet = MathHelper::FrechetDistanceBtwTwoLines(curve1, curve2, n_size, other.size(), n1, n2, 0);
            float frechet2 = MathHelper::FrechetDistanceBtwTwoLines(curve1, curve3, n_size, other.size(), n1, n2, 0);
            double v = min(frechet, frechet2);

            if ( v > maxFrechet)
                maxFrechet = v;


            tatuLoc.push_back( make_pair(pos, v) );
            originalPairs.push_back( make_pair(i,j));
            pos++;
        }
    }

    std::cout << "Going to generate SPLOM Ordering based on Hausdorff " << std::endl;

    while ( static_cast<int>(orderAccordingToSimilarity.size()) != num){

        float minAvgSimilarity = 9999999;
        int newItem = -1;

        for(int k = 0; k < num; k++){
            if ( count(orderAccordingToSimilarity.begin(),orderAccordingToSimilarity.end(), k )) continue;
            // if we already added it, then we have to reason to use it

            // calculate the avg similarity with the already sorted
            float avgSim = 0;

            for( unsigned int h = 0; h < orderAccordingToSimilarity.size(); h++){
                int attr1 = k;
                int attr2 = orderAccordingToSimilarity.at(h);

                if ( attr1 > attr2){
                    attr1 = orderAccordingToSimilarity.at(h);
                    attr2 = k;
                }
                // ***************+
                // Since we are using Tatu we can search in the original pairs

                for(unsigned int pairs = 0; pairs < originalPairs.size(); pairs++){

                    if ( originalPairs.at(pairs).first == attr1 && originalPairs.at(pairs).second){
                        avgSim += fabs(tatuLoc.at(pairs).second);
                    }
                }
            }
            avgSim /= originalPairs.size();
            if ( avgSim < minAvgSimilarity){
                minAvgSimilarity = avgSim;
                newItem = k;
            }
        }
        orderAccordingToSimilarity.push_back(newItem);
    }


    reOrderedSplom.clear();



    for(unsigned int i = 0; i < orderAccordingToSimilarity.size(); i++){
        for(unsigned int j = i+1; j < orderAccordingToSimilarity.size(); j++){

             int attr1 = orderAccordingToSimilarity.at(i);
             int attr2 = orderAccordingToSimilarity.at(j);

             if ( attr1 > attr2){
                attr2 = orderAccordingToSimilarity.at(i);
                attr1 = orderAccordingToSimilarity.at(j);
             }

             for(unsigned int k = 0; k < originalPairs.size(); k++){
                 if ( originalPairs.at(k).first == attr1 && originalPairs.at(k).second == attr2){
                     reOrderedSplom.push_back(k);
                     break;
                 }
             }
        }
    }









}

void SkeletonVis::ToggleFreeForm(){
   std::cout <<" Going to toggle free form " << std::endl;

   if ( modelComparison){

       useFreeFormModel = !useFreeFormModel;
       std::cout << "Use Free Form model? " << useFreeFormModel << std::endl;

       if ( useFreeFormModel){
           freeFromDialog = new FreeFormWidget();
           freeFromDialog->show();
           connect(freeFromDialog,SIGNAL(SomethingDrawn()), this, SLOT(NewFreeFormModel()));
       }
   }
}

void SkeletonVis::NewFreeFormModel(){

    std::cout << "New model drawn? " << std::endl;
    Redraw();
}
