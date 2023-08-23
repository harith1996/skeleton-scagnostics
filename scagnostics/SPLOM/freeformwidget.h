#ifndef FREEFORMWIDGET_H
#define FREEFORMWIDGET_H


/*** Free form widget draws
 *
 *  A model to compare against.
 *
 *  Basically we store the points
 *  and then we either give a curve
 *  or an itk image to compare against...
 *
 * */


#include <QWidget>
#include <QGLWidget>
#include <GL/glut.h>
#include <QMouseEvent>

#include "../dataset.h" // Just for LocalPoint
#include "itkNumericTraits.h"

namespace Ui {
class FreeFormWidget;
}

class FreeFormWidget : public QGLWidget
{
    Q_OBJECT

public:
    explicit FreeFormWidget(QWidget *parent = 0);
    ~FreeFormWidget();

    void Init();

    void GetCurve(vector<LocalPoint>* pts);
    void GetImage(HImageType::Pointer img, int size);
    void GetWorldCoordinates(double x, double y, double *xAtPress, double *yAtPress);

  /*  void Press(QMouseEvent* event){ mousePressEvent(event);}
    void Moved(QMouseEvent* event){ mouseMoveEvent(event);}
    void Release(QMouseEvent* event){ mouseReleaseEvent(event);}
    void DoubleClick(QMouseEvent* event){ mouseDoubleClickEvent(event);}
*/
protected:
   void paintGL();
   void initializeGL();
   void resizeGL(int width, int height);
   virtual void mousePressEvent(QMouseEvent* event);
   virtual void mouseDoubleClickEvent(QMouseEvent* event);
   virtual void mouseMoveEvent(QMouseEvent *event);
   virtual void mouseReleaseEvent(QMouseEvent *event);
   virtual void wheelEvent(QWheelEvent *event);
signals:
   void SomethingDrawn();
private:
    Ui::FreeFormWidget *ui;
    GLint viewport[4];

    vector<LocalPoint> pointsDrawn;
    bool isMousePressed;
};

#endif // FREEFORMWIDGET_H
