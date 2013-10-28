/* OPTICS includes */
#ifndef _GLOBALS
#define _GLOBALS
class globalData {
public:
   globalData(void) { dim=-1; }
   int dim;         // dimension of the data
   double eps;      // eps for the clustering algorithm
   int minPts;      // minPts value for the clustering algorithm
   int dtype;       // way to calculate distances; should use enum but...
                    // dtype=0 is euclidean
                    // dtype=1 is 3-d euclidean with cos(dist)
};
#endif
