#include "freeformwidget.h"
#include "ui_freeformwidget.h"

FreeFormWidget::FreeFormWidget(QWidget *parent) :
    QGLWidget( QGLFormat(QGL::SampleBuffers), parent),
    ui(new Ui::FreeFormWidget)
{
    ui->setupUi(this);
}

FreeFormWidget::~FreeFormWidget()
{
    delete ui;
}

void FreeFormWidget::Init(){
    setFocusPolicy(Qt::StrongFocus);
}


void FreeFormWidget::paintGL(){
    glEnable( GL_LINE_SMOOTH );
    glEnable( GL_POLYGON_SMOOTH );
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    glViewport(0, 0, (GLsizei) viewport[2], (GLsizei) viewport[3]);
    glClearColor(0.95, 0.95, 0.95, 1.0);
    glClear(GL_COLOR_BUFFER_BIT);

    // Setup projection matrix
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(0.0, 1.0, 0.0, 1.0);

    std::cout << "Drawing something? " << std::endl;

    glColor3f(0,0,0);
    glLineWidth(3.0);

    glBegin(GL_LINE_STRIP);
        for(int i = 0; i < this->pointsDrawn.size(); i++){
            glVertex2f( pointsDrawn.at(i).p[0],pointsDrawn.at(i).p[1]);
        }
    glEnd();
}

void FreeFormWidget::initializeGL(){
    glutInitDisplayMode(GLUT_RGBA | GLUT_ALPHA);
    glViewport(0,0, width(), height());
    glGetIntegerv( GL_VIEWPORT, viewport );
    std::cout << "width & height "<< width() << " , " << height() << std::endl;
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(0.0, 1.0, 0.0, 1.0);
    glMatrixMode(GL_MODELVIEW);
    glGetIntegerv( GL_VIEWPORT, viewport );
}
void FreeFormWidget::resizeGL(int width, int height){
    int side = qMin(width, height);
    //std::cout << "resize gl in free form?? "<< width << " , "<< height <<std::endl;

     glViewport((width - side) / 2, (height - side) / 2, side, side);
     glGetIntegerv( GL_VIEWPORT, viewport );

    //float blockSize = 2.0;
     glMatrixMode(GL_PROJECTION);
     glLoadIdentity();
     gluOrtho2D(0.0, 1.0, 0.0, 1.0);
     glMatrixMode(GL_MODELVIEW);

     //std::cout << viewport[0] << " .. "<< viewport[1] << "..." << viewport[2] << " .." << viewport[3] <<std::endl;
     update();
     updateGL();
}
void FreeFormWidget::mousePressEvent(QMouseEvent* event){


    std::cout << "Mouse is pressed " << std::endl;
    isMousePressed = true;
    pointsDrawn.clear();

    LocalPoint p;

    double x;
    double y;

    GetWorldCoordinates(event->x(), event->y(), &x, &y);
    p.p[0] = x;
    p.p[1] = y;
    std::cout << x << " , " << y << std::endl;
    pointsDrawn.push_back(p);
    updateGL();

}
void FreeFormWidget::mouseDoubleClickEvent(QMouseEvent* event){
    pointsDrawn.clear();
    updateGL();

}
void FreeFormWidget::mouseMoveEvent(QMouseEvent *event){
    std::cout << "Mouse is moved " << std::endl;

    if (isMousePressed){
        LocalPoint p;

        double x;
        double y;

        GetWorldCoordinates(event->x(), event->y(), &x, &y);
        p.p[0] = x;
        p.p[1] = y;
        std::cout << x << " , " << y << std::endl;
        pointsDrawn.push_back(p);
    }
    updateGL();


}
void FreeFormWidget::mouseReleaseEvent(QMouseEvent *event){
    std::cout << "Mouse is released" << std::endl;

    isMousePressed = false;
    emit SomethingDrawn();
    updateGL();

}

void FreeFormWidget::GetCurve(vector<LocalPoint>* pts){

    for(unsigned int i = 0; i < pointsDrawn.size(); i++){
        pts->push_back(pointsDrawn.at(i));
    }
}

void FreeFormWidget::GetImage(HImageType::Pointer img, int size){

}

void FreeFormWidget::wheelEvent(QWheelEvent *event){

}

void FreeFormWidget::GetWorldCoordinates(double x, double y, double *xAtPress, double *yAtPress){
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
