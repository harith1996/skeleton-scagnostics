#ifndef GUI_PLIST_HPP
#define GUI_PLIST_HPP

#include "imagecontroller.hpp"

void projectionListReset();

void projectionListExtend();

void projectionListMap();

void projectionListMapTSNE();

void projectionListMapDrawOn(ImageController *imgc);

void projectionListMapClick(float x, float y);

#endif
