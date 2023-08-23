#ifndef GUI_IMAGECONTROLLER_HPP
#define GUI_IMAGECONTROLLER_HPP

#include <gtkmm.h>
#include <functional>

struct ImageController {
    Gtk::Image *image;
    Gtk::Viewport *viewport;
    Gtk::ScrolledWindow *viewportScrolledWindow;
    Glib::RefPtr<Gtk::Adjustment> adjustmentZoom;
    std::function<void(ImageController *imgc)> renderCallback;
    std::function<void(float x, float y)> clickCallback;

    void init();

    void setImage(cairo_surface_t *surface);

    double getZoom();

    private:

    void renderImage();
    void onViewportRealize();
    bool onMouseMove(GdkEventMotion *ev);
    bool onMouseButton(GdkEventButton *ev);
    bool onMouseWheel(GdkEventScroll *ev);
    void setViewportCursor();
    bool onImageMouseButton(GdkEventButton *ev);

    bool dragging = false;
    double dragLastX, dragLastY;
    int imageHeight = 1, imageWidth = 1;
};

#endif

