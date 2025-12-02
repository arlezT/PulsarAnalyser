#pragma once
#include <vector>
#include <cstddef>
#include <cstdint>

/*
 * TimeFrequency:
 * Handles per-channel time-frequency data processing.
 * Provides functions to compute median, standard deviation, and clean outliers.
 * Supports caching of channel medians for faster repeated computations.
 */

class TimeFrequency {
public:

    // Constructor for float input data
    TimeFrequency(
        std::vector<float> data,
        unsigned numChannels,
        size_t numSpectra,
        unsigned startFq,
        unsigned endFq,
        float chWidth,
        double samplingT);

    // Constructor for 8-bit unsigned integer input data
    TimeFrequency(
        const std::vector<uint8_t>& data8,
        unsigned numChannels,
        size_t numSpectra,
        unsigned startFq,
        unsigned endFq,
        float chWidth,
        double samplingT);

    // Constructor for 16-bit unsigned integer input data
    TimeFrequency(
        const std::vector<uint16_t>& data16,
        unsigned numChannels,
        size_t numSpectra,
        unsigned startFq,
        unsigned endFq,
        float chWidth,
        double samplingT);

    // Accessors for the underlying data and metadata.
    // data(): returns a reference to the internal data buffer
    // updateData(): replaces the stored time-frequency data
    // getNumChannels() / getNumSpectra(): metadata describing the layout
    // getChannelPtr(): returns a pointer to the first element of a given channel
    const std::vector<float>& data() const { return data_; }
    std::vector<float>& data() { return data_; }
    void updateData(const std::vector<float>& newData) { data_ = newData; }
    unsigned getNumChannels() const { return numChannels_; }
    size_t getNumSpectra() const { return numSpectra_; }
    float* getChannelPtr(unsigned ch) {return &data_[ch * numSpectra_];}
    const float* getChannelPtr(unsigned ch) const {return &data_[ch * numSpectra_];}


    // Computes median of a channel's data.
    // If cache is valid, returns cached value. Otherwise, computes and updates cache.
    void computeChannelMedian(const float* channelPtr, size_t vectorSize, float& median, unsigned channel);

    // Computes standard deviation of a channel given its median.
    void computeStdDev(const float* channelPtr, size_t vectorSize, float median, double& stdDeviation);
    
    // Replaces values exceeding threshold with median. Returns 1 if channel was modified, 0 otherwise.
    unsigned cleanData(float* channelPtr, size_t vectorSize, float median, float threshold);
    
    // Sets whether cached medians are used in computeChannelMedian
    void setCacheValid(bool valid); //Sets whether cached medians are used in computeChannelMedian

private: 

    // Helper for converting integer input data to float
    template<typename T>
    void loadFromIntegerData(const std::vector<T>& input) 
    {
        data_.resize(input.size());
        for (size_t i = 0; i < input.size(); ++i){
            data_[i] = static_cast<float>(input[i]);
    }}

    std::vector<float> data_; // Data is stored in the class
    std::vector<float> cachedMedians_; //Stores cached medians per channel
    bool cacheValid_;

	// Data attributes
    unsigned numChannels_; //4096;
    size_t numSpectra_;  //10000;
    unsigned startFq_; //1700;
    unsigned endFq_; //2000;
    float chWidth_; //(last_channel_fq - first_channel_fq) / num_channels;
    double samplingT_; //50e-6;
};