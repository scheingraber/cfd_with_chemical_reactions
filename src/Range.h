#ifndef RANGE_A37M4G5V
#define RANGE_A37M4G5V

class Range {
    public:
        Range ();
        Range (int low, int high);

        // BEWARE: public, direct access possible
        int low, high;
};

#endif /* end of include guard: RANGE_A37M4G5V */
