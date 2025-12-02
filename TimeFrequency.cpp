#include "TimeFrequency.h"
#include <algorithm>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <iostream>

using namespace std;

// Float constructor
TimeFrequency::TimeFrequency(
    vector<float> spectra,
    unsigned numChannels,
    size_t numSpectra,
    unsigned startFq,
    unsigned endFq,
    float chWidth,
    double samplingT)
    : data_(std::move(spectra)),
    numChannels_(numChannels),
    numSpectra_(numSpectra),
    startFq_(startFq),
    endFq_(endFq),
    chWidth_(chWidth),
    samplingT_(samplingT),
    cacheValid_(false)
{
    cachedMedians_.resize(numChannels_, 0.0f);
} // Resizes cachedMedians_ vector to hold one median per channel

// 8-bit constructor
TimeFrequency::TimeFrequency(
    const std::vector<uint8_t>& data8,
    unsigned numChannels,
    size_t numSpectra,
    unsigned startFq,
    unsigned endFq,
    float chWidth,
    double samplingT)
    : numChannels_(numChannels),
    numSpectra_(numSpectra),
    startFq_(startFq),
    endFq_(endFq),
    chWidth_(chWidth),
    samplingT_(samplingT),
    cacheValid_(false)
{
    loadFromIntegerData(data8);
    cachedMedians_.resize(numChannels_, 0.0f);
}

// 16-bit constructor
TimeFrequency::TimeFrequency(
    const std::vector<uint16_t>& data16,
    unsigned numChannels,
    size_t numSpectra,
    unsigned startFq,
    unsigned endFq,
    float chWidth,
    double samplingT)
    : numChannels_(numChannels),
    numSpectra_(numSpectra),
    startFq_(startFq),
    endFq_(endFq),
    chWidth_(chWidth),
    samplingT_(samplingT),
    cacheValid_(false)
{
    loadFromIntegerData(data16);
    cachedMedians_.resize(numChannels_, 0.0f);
}

// Compute median of a channel's data
// If cache is valid, return cached value
// For even-sized vectors, median is average of middle element and max of lower half
// Throws exceptions for invalid pointers, sizes, or channel index
void TimeFrequency::computeChannelMedian(const float* channelPtr,
    size_t vectorSize, 
    float& median,
    unsigned channel) 
{
    if (cacheValid_) {
        median = cachedMedians_[channel];
        return;
    }
    if (!channelPtr) throw invalid_argument("Channel pointer is null");
    if (vectorSize == 0) throw invalid_argument("Vector size cannot be 0");
    if (channel >= numChannels_) throw out_of_range("Channel index exceeds the number of channels");

    vector<float> temp(channelPtr, channelPtr + vectorSize);
    size_t n = vectorSize / 2;
    nth_element(temp.begin(), temp.begin() + n, temp.end());
    if (vectorSize % 2 == 1) {
        median = temp[n];
    }
    else {
        median = 0.5f * (temp[n] + *std::max_element(temp.begin(), temp.begin() + n));
    }
    cachedMedians_[channel] = median;
}

// Sets the cacheValid_ flag
// When true, computeChannelMedian returns cached medians instead of recomputing
void TimeFrequency::setCacheValid(bool valid) {
    cacheValid_ = valid;
}

// Compute standard deviation given a channel's data and median
// Throws if inputs are invalid
void TimeFrequency::computeStdDev(const float* channelPtr, 
    size_t vectorSize, 
    float median, 
    double& stdDeviation)
{
    if (!channelPtr) throw invalid_argument("Channel pointer is null");
    if (vectorSize == 0) throw invalid_argument("Vector size cannot be 0");
    if (!isfinite(median)) throw invalid_argument("Value provided for median is invalid.");

    auto begin = channelPtr;
    auto end = begin + vectorSize;
    float sumSq = 0.0;
    for (auto it = begin; it != end; it++) {
        float delta = *it - median;
        sumSq += delta * delta;
    }
    stdDeviation = sqrt(sumSq / vectorSize);
}

// Replaces values in a channel that exceed the threshold with the median
// Returns 1 if the channel was flagged and modified, 0 otherwise
// Throws exceptions for invalid inputs
unsigned TimeFrequency::cleanData(float* channelPtr,
    size_t vectorSize, 
    float median, 
    float threshold) 
{
    if (!channelPtr) throw invalid_argument("Channel pointer is null");
    if (vectorSize == 0) throw invalid_argument("Vector size cannot be 0");
    if (!isfinite(median)) throw invalid_argument("Value provided for median is invalid.");

    auto begin = channelPtr;
    auto end = begin + vectorSize;

    int channelFlagged = 0;
    for (auto it = begin; it != end; ++it) {
        if (*it > threshold) {
            channelFlagged = 1;
            break;
        }
    }

    if (channelFlagged) {
        for (auto it = begin; it != end; ++it) {
            *it = median;
        }
    }

    return channelFlagged;
}

// Alternative for computing median: approximate median using random sampling (tested but not implemented)
//float approximate_median(float* v, size_t N, int num_samples = 4) {
//    static thread_local std::mt19937 gen(std::random_device{}());
//    std::uniform_int_distribution<> dis(0, N - 1);
//
//    size_t subset_size = N / num_samples;
//    std::vector<float> subset(subset_size);
//    double sum_of_medians = 0.0;
//
//    for (int s = 0; s < num_samples; ++s) {
//        // fill subset
//        for (size_t i = 0; i < subset_size; ++i)
//            subset[i] = v[dis(gen)];
//
//        // compute median
//        size_t n = subset_size / 2;
//        std::nth_element(subset.begin(), subset.begin() + n, subset.end());
//        //sort(subset.begin(), subset.end());
//        double median;
//        if (subset_size % 2 == 1) {
//            median = subset[n];
//        }
//        else {
//            //std::nth_element(subset.begin(), subset.begin() + (n - 1), subset.end());
//            median = 0.5 * (subset[n - 1] + subset[n]);
//        }
//        sum_of_medians += median;
//    }
//
//    return static_cast<float>(sum_of_medians / num_samples);
//}
