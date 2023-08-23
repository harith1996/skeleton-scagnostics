#include "imagecontroller.hpp"

void ImageController::init() {
    adjustmentZoom->signal_value_changed().connect(sigc::mem_fun(*this, &ImageController::renderImage));
    if(!viewport->get_realized()) {
        viewport->signal_realize().connect(sigc::mem_fun(*this, &ImageController::onViewportRealize));
    } else {
        onViewportRealize();
    }
}

void ImageController::onViewportRealize() {
    auto view_window = viewport->get_view_window();
    view_window->set_events(view_window->get_events() | Gdk::BUTTON_PRESS_MASK | Gdk::BUTTON_RELEASE_MASK | Gdk::SCROLL_MASK | Gdk::POINTER_MOTION_MASK);

    viewport->signal_button_press_event().connect(sigc::mem_fun(*this, &ImageController::onMouseButton));
    viewport->signal_button_release_event().connect(sigc::mem_fun(*this, &ImageController::onMouseButton));
    viewport->signal_motion_notify_event().connect(sigc::mem_fun(*this, &ImageController::onMouseMove));
    viewport->signal_scroll_event().connect(sigc::mem_fun(*this, &ImageController::onMouseWheel));

    view_window = viewport->get_bin_window();
    view_window->set_events(view_window->get_events() | Gdk::BUTTON_RELEASE_MASK);
    viewport->signal_button_release_event().connect(sigc::mem_fun(*this, &ImageController::onImageMouseButton));
}

void ImageController::setImage(cairo_surface_t *surface_c) {
    int imgWidth = cairo_image_surface_get_width(surface_c);
    int imgHeight = cairo_image_surface_get_height(surface_c);

    imageWidth = imgWidth;
    imageHeight = imgHeight;

    float relX = 0.5;
    float relY = 0.5;
    auto ha = viewportScrolledWindow->get_hadjustment();
    auto va = viewportScrolledWindow->get_vadjustment();

    if(ha->get_upper() != ha->get_lower()) {
        relX = (ha->get_value() + ha->get_page_size()/2 - ha->get_lower()) / (ha->get_upper() - ha->get_lower());
    }
    if(va->get_upper() != va->get_lower()) {
        relY = (va->get_value() + va->get_page_size()/2 - va->get_lower()) / (va->get_upper() - va->get_lower());
    }

    auto surface = Cairo::RefPtr<Cairo::Surface>(new Cairo::Surface(surface_c, false));

    Glib::RefPtr<Gdk::Pixbuf> pixbuf { Gdk::Pixbuf::create(surface, 0, 0, imgWidth, imgHeight) };
    image->set(pixbuf);

    viewport->check_resize();
    viewportScrolledWindow->check_resize();

    ha->set_value(ha->get_lower() + relX * (ha->get_upper() - ha->get_lower()) - ha->get_page_size()/2);
    va->set_value(va->get_lower() + relY * (va->get_upper() - va->get_lower()) - va->get_page_size()/2);
}

double ImageController::getZoom() {
    return adjustmentZoom->get_value();
}

void ImageController::renderImage() {
    if(renderCallback) {
        renderCallback(this);
    }
}

void ImageController::setViewportCursor() {
    viewport->get_view_window()->set_cursor(Glib::wrap(gdk_cursor_new_from_name(viewport->get_view_window()->get_display()->gobj(),
                    dragging ? "grabbing" : "default")));
}

bool ImageController::onMouseButton(GdkEventButton *ev) {
    if(ev->type == GDK_BUTTON_PRESS && ev->button == 1) {
        dragging = true;
        dragLastX = ev->x;
        dragLastY = ev->y;
        setViewportCursor();

        return true;
    } else if(ev->type == GDK_BUTTON_RELEASE && ev->button == 1) {
        dragging = false;
        setViewportCursor();

        return true;
    }
    return false;
}

bool ImageController::onMouseMove(GdkEventMotion *ev) {
    if(!dragging) return false;

    auto ha = viewportScrolledWindow->get_hadjustment();
    auto va = viewportScrolledWindow->get_vadjustment();

    ha->set_value(ha->get_value() - ev->x + dragLastX);
    va->set_value(va->get_value() - ev->y + dragLastY);

    dragLastX = ev->x;
    dragLastY = ev->y;

    return true;
}

bool ImageController::onMouseWheel(GdkEventScroll *ev) {
    if(ev->direction == GDK_SCROLL_UP) {
        adjustmentZoom->set_value(adjustmentZoom->get_value() + 2 * adjustmentZoom->get_step_increment());
        return true;
    } else if(ev->direction == GDK_SCROLL_DOWN) {
        adjustmentZoom->set_value(adjustmentZoom->get_value() - 2 * adjustmentZoom->get_step_increment());
        return true;
    }
}

bool ImageController::onImageMouseButton(GdkEventButton *ev) {
    if(clickCallback && ev->type == GDK_BUTTON_RELEASE && ev->button == 3) {
        int dx = ev->x, dy = ev->y;
        viewportScrolledWindow->translate_coordinates(*image, dx, dy, dx, dy);

        //actual image is centered in image widget
        if(image->get_allocated_width() > imageWidth) {
            dx -= (image->get_allocated_width() - imageWidth) / 2;
        }
        if(image->get_allocated_height() > imageHeight) {
            dy -= (image->get_allocated_height() - imageHeight) / 2;
        }

        clickCallback(dx * 1.0 / imageWidth, dy * 1.0 / imageHeight);
        return true;
    }
    return false;
}
