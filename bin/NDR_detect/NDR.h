#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <set>
#include <cmath>
#include <zlib.h>
#include <pthread.h>
#include <inttypes.h>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/math/distributions/chi_squared.hpp>

#include "gzstream.h"

class NDR{
    public:
        NDR();
        ~NDR();
        void append_NDR();
        void report_NDR();
        void chisq_NDR();
        
}