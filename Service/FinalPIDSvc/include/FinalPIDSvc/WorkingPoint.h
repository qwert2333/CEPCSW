#include <stdlib.h>
#include <iostream>
#include <map>

class WP {
    public:
        WP();
        ~WP();

        enum _wp
        {
            is98 = 0,
            is90 = 1,
            is70 = 2,
            is50 = 3,
            noLep = 4,
            Best = 5,
            MAX_NUM_OF_WPS
        };

        static const int num_of_WPs = MAX_NUM_OF_WPS;

        static const std::map<int, const char*> WPs;
};