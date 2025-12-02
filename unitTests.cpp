#include <cassert>
#include <vector>
#include <cmath>
#include <iostream>
#include "TimeFrequency.h"

int main() {
    TimeFrequency tf(5, 10, 1700, 2000, 0.1f, 50e-6);

    float data[] = { 1, 7, 0, 1, 2 };
    // 1. Test computeChannelMedian
    {
        float median = 0.f;
        tf.computeChannelMedian(data, 4, median, 0);
        assert(median == 1.f); //0,1,*1*,2,7
    }

    // 2. Test computeStdDev
    {
        float median = 1.f;
        double stdDev = -1;
        double c = 1e-6;
        double expected = 2.756809750418044;
        tf.computeStdDev(data, 5, median, stdDev);
        assert(std::fabs(stdDev - expected) < c);
    }

    // 3. Test clean data
    {
        float data[] = { 1, 1000, 2 };
        float median = 2.f;
        float threshold = 10.f;
        unsigned flag = tf.cleanData(data, 3, median, threshold);
        assert(flag == 1);
        assert(data[1] == median);
    }

    // Test caching
    {
        tf.setCacheValid(true);
        float median = 0.f;
        tf.computeChannelMedian(nullptr, 5, median, 0);
    }

    std::cout << "All tests passed" << std::endl;
    return 0;
}
