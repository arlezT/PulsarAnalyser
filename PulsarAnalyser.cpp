#include <cstdio>
#include <iostream>
#include <vector>
#include <random>
#include <fstream>
#include <filesystem>
#include <chrono> 

using namespace std;

#define NUMBER_OF_CHANNELS 4096   // number of channels per spectra
#define NUMBER_OF_SPECTRA 1e4      // number of spectra in half second data (0.5/50us) 50 micro seconds is sampling time
#define MAX_FLOAT_VALUE 1e7 // not ideal but we dont expect the value to be lrger than this.
#define THRESHOLD 4.5 // threshold to find the outliers

/* A central component of RFI mitigation is the calculation of reliable channel statistics across the
time–frequency data. For each frequency channel, we examine how the signal varies over time and compute
a robust estimate of its typical behaviour. This method estimates median value per frequency channel.
This involves taking all time samples corresponding to a single channel, forming a one-dimensional
array of power values, and computing the median of that distribution. The resulting median vector will
give the baseline power estimate. */

void calculate_median(std::vector<float> chunk, std::vector<float>& median)
{
    /* The provided implementation is deliberately very rudimentary and serves only as a starting point.
    You are encouraged to design and implement your own approach, improve the structure, and introduce
    more efficient or more robust methods wherever appropriate. We also encourage multiple implementations
    and a provision to select algorithm.
    */

    for (unsigned channel = 0; channel < NUMBER_OF_CHANNELS; ++channel)
    {
        std::vector<float> temp(chunk.size() / NUMBER_OF_CHANNELS);

        unsigned sample = 0;
        //cout<<"calculating median for channel: "<<channel<<" \r";
        cout << "calculating median for channel: " << channel << endl;
        for (unsigned iter = 0; iter < temp.size() / 2 + 1; ++iter)
        {

            float min = MAX_FLOAT_VALUE;
            unsigned arg = 0;
            for (unsigned i = 0; i < temp.size(); ++i)
            {
                //std::cout << "min at iteration " << i << " = " << min << endl;
                if (min > chunk[temp.size() * channel + i])
                {
                    min = chunk[temp.size() * channel + i];
                    arg = i;
                }
            }
            temp[iter] = min;
            chunk[temp.size() * channel + arg] = MAX_FLOAT_VALUE;
        }
        median[channel] = (temp[temp.size() / 2] + temp[temp.size() / 2 - 1]) / 2.0;
    }
    std::cout << "\n done estimating median for the chunk \n";
}

/*
This method will compute a per-channel standard deviation to quantify the typical fluctuations around the median value.
For each frequency channel, subtract the channel’s median from every corresponding time sample, and then evaluate
the standard deviation of these residuals. This provides a measure of the intrinsic variability in that channel and
helps identify channels with unusually high scatter, which are often indicative of intermittent RFI. The combination
of the median spectrum and the per-channel standard deviation forms a robust statistical baseline for subsequent
RFI detection and cleaning operations.
*/
void calculate_std(std::vector<float> chunk, std::vector<float> median, std::vector<float>& std_dev)
{
    double sum;
    for (unsigned int channel = 0; channel < NUMBER_OF_CHANNELS; ++channel)
    {
        sum = 0.0;

        for (unsigned sample = 0; sample < chunk.size() / NUMBER_OF_CHANNELS; ++sample)
        {
            sum += std::pow(chunk[channel * chunk.size() / NUMBER_OF_CHANNELS + sample] - median[channel], 2);
        }
        std_dev[channel] = std::sqrt(sum / (chunk.size() / NUMBER_OF_CHANNELS));
    }
}

/*
This method cleans the TF data Using the median spectrum and per-channel standard deviations, you will identify outliers within each frequency channel.
For every time sample, compute the deviation from the channel’s median and flag samples whose absolute deviation exceeds
a chosen threshold (for example, several times the per-channel standard deviation). These flagged values are considered
RFI-affected outliers. Once identified, replace each outlier with the corresponding channel median (here we replace whole channel with median value flagging the whole channel).
This produces a cleaned time–frequency dataset in which strong, isolated RFI features are suppressed while the underlying astrophysical
signal is preserved as far as possible.
*/
void clean_data(std::vector<float>& chunk, std::vector<float> median, std::vector<float> std_dev)
{
    unsigned number_of_channels_flagged = 0;
    for (unsigned channel = 0; channel < NUMBER_OF_CHANNELS; ++channel)
    {
        int flag = 0;
        for (unsigned int sample = 0; sample < chunk.size() / NUMBER_OF_CHANNELS; ++sample)
        {
            if (chunk[channel * chunk.size() / NUMBER_OF_CHANNELS + sample] > (THRESHOLD * std_dev[channel] + median[channel])) // I have used threshold of 4.5 (this must be option in your code.)
            {
                flag = 1;
            }
        }
        if (flag == 1)
        {
            number_of_channels_flagged++;
            for (unsigned int sample = 0; sample < chunk.size() / NUMBER_OF_CHANNELS; ++sample) chunk[channel * chunk.size() / NUMBER_OF_CHANNELS + sample] = median[channel];
        }

    }
    std::cout << "channel flagged: " << number_of_channels_flagged << "\n";
}

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
//
//        sum_of_medians += median;
//    }
//
//    return static_cast<float>(sum_of_medians / num_samples);
//}

void computeMedian(vector<float>& data, size_t refIdx, size_t vectorSize,float& median) {
    float* begin = &data[refIdx];
    float* end = &data[refIdx + vectorSize];
    size_t  n = vectorSize / 2;
    nth_element(begin, begin + n, end);
    if (vectorSize % 2 == 1) {
        //medians.push_back(begin[n]); 
        median = begin[n];
    }
    else {
        float upperMedian = begin[n];
        float lowerMedian = *std::max_element(begin, begin + n);
        median=(upperMedian + lowerMedian)*0.5f;
    }
}


int main(int argc, char* argv[])
{
    /*
    The TF data in the data.bin is read in chunks of size of half a second. This is the typical size we use to proceess the data stream for RFI. Hence for 1 second
    data we read the data.bin in two iterations. For each chunk we obtain channel staistics (median and standard deviation). We use these robost statistics to find the
    outliers and replace the effected channels with median value. The cleaned TF is written out in clean_data.bin file.
    */

    /*
    Time-frequency data parameter;

    number of frequency channels = 4096;
    number of spectra = 10000 equivalent to half second.
    frequency of first channel = 1700 MHz
    frequency of last channel = 2000 MHz
    channel width = (2000-1700)/4096 MHz
    sampling time = 50e-6 (50 micro seconds)

    hint: create a class TimeFrequency for example look (https://gitlab.com/ska-telescope/pss/ska-pss-cheetah/-/blob/main/cheetah/data/TimeFrequency.h)
    this is just an example dont bother about the architeture. please feel free to define a simpler version. Make sure that you use the datatype template parameter.
    The datatype of the data.bin is in floats but we need TimeFrequency to be able to deal with 8bit unsinged integers and 16bit unsigned integers and have all the
    above values as private members.

    */

    std::vector<float> chunk(NUMBER_OF_SPECTRA * NUMBER_OF_CHANNELS);
    std::ifstream in("data.bin", std::ios::binary);
    std::ofstream out("cleaned_data.bin", std::ios::binary);

    //std::vector<float> median(NUMBER_OF_CHANNELS);
    //std::vector<float> std_dev(NUMBER_OF_CHANNELS);

    if (!in.is_open()) {
        std::cerr << "ERROR: could not open data.bin" << std::endl;
        return 1;
    }

    else {

        for (unsigned i = 0; i < 2; ++i)
        {
            auto t1 = chrono::high_resolution_clock::now();

            in.read(reinterpret_cast<char*>(chunk.data()), chunk.size() * sizeof(float));

            vector<float> modifiedData;
            modifiedData = chunk;
            
            vector<float> medians;
            medians.reserve(NUMBER_OF_CHANNELS);

            size_t chunkSize = chunk.size() / NUMBER_OF_CHANNELS;
            float median = 0.0f;
            float stdDev;
            unsigned idx = 0;
            unsigned flags = 0;

            for (std::size_t fC = 0; fC < modifiedData.size(); fC += chunkSize) {
                computeMedian(modifiedData, fC, chunkSize, median);
                auto begin = modifiedData.begin() + fC;
                auto end = begin + chunkSize;
                float sumSq = 0.0;
                for (auto it = begin; it != end; it++) {
                    float delta = *it - median;
                    sumSq += delta * delta;
                }
                stdDev = sqrt(sumSq / (chunkSize - 1));
                float adjustedThreshold = THRESHOLD * stdDev + median;
                for (auto it = begin; it != end; it++) {
                    if (*it > adjustedThreshold) {
                        flags++;
                        *it = median;
                    }
                }
                idx++;
            }

            auto t2 = chrono::high_resolution_clock::now();
            auto ms_int = chrono::duration_cast<chrono::milliseconds>(t2 - t1);
            std::cout << ms_int.count() << "ms\n";

            //out.write(reinterpret_cast<char*>(chunk.data()), chunk.size()*sizeof(float));
        }
        out.close();
        in.close();
    }
}

////////////////////////////////////////////////////////////////////////////////////