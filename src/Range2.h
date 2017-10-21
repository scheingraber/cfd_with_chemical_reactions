#ifndef RANGE2_HNTI3559
#define RANGE2_HNTI3559

#include "Range.h"

class Range2 {
    public:
        Range2 ();
        Range2 (int ilow, int ihigh, int jlow, int jhigh);
        Range2 (Range i, Range j);

        // BEWARE: public, direct access possible
        Range i, j;
};


#endif /* end of include guard: RANGE2_HNTI3559 */
