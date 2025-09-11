#include "FinalPIDSvc/WorkingPoint.h"

WP::WP() {}
WP::~WP() {}

const std::map<int, const char*> WP::WPs = {
    {0, "WP98"},
    {1, "WP90"},
    {2, "WP70"},
    {3, "WP50"},
    {4, "noLep"},
    {5, "Best"}
};