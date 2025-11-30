#pragma once
#include <vector>
class TimeFrequency {
public:
	TimeFrequency(
        unsigned numChannels_, 
        unsigned numSpectra_, 
        unsigned startFq_, 
        unsigned endFq_, 
        float chWidth_,
        double samplingT_);
private: 
    std::vector<float> broadbandData_;
	// Data attributes
    unsigned numChannels_; //4096;
    unsigned numSpectra_;  //10000;
    unsigned startFq_; //1700;
    unsigned endFq_; //2000;
    float chWidth_; //(last_channel_fq - first_channel_fq) / num_channels;
    double samplingT_; //50e-6;
};