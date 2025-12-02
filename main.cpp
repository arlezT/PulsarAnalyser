#include <cstdio>
#include <iostream>
#include <vector>
#include <random>
#include <fstream>
#include <filesystem>
#include <chrono> 
#include <thread>

#include "TimeFrequency.h"

using namespace std;

#define NUMBER_OF_CHANNELS 4096   // number of channels per spectra
#define NUMBER_OF_SPECTRA 10000      // number of spectra in half second data (0.5/50us) 50 micro seconds is sampling time
#define MAX_FLOAT_VALUE 1e7 // not ideal but we dont expect the value to be lrger than this.
#define THRESHOLD 4.5 // threshold to find the outliers
#define START_FREQUENCY 1700
#define END_FREQUENCY 2000
#define SAMPLING_TIME 50e-6

//struct RFIModule 
//{
//    TimeFrequency process(const TimeFrequency& input)
//    {
//        TimeFrequency tf = input;
//        // Initialise variables for the channel median and number of channels replaced 
//        unsigned channels_flagged = 0;
//
//        // Determine number of available CPU cores for parallel processing
//        // and calculate how many channels each thread should process
//        unsigned numThreads = thread::hardware_concurrency();
//        if (numThreads == 0) numThreads = 1;
//
//        vector<thread> threads;
//        threads.reserve(numThreads);
//        size_t chPerThread = tf.getNumChannels() / numThreads;
//
//        // Worker lambda: processes a batch of channels
//        // Computes median and standard deviation per channel
//        // Cleans outliers based on an adjusted threshold
//        // Previous median per channel is cached for speed and assumed statistically valid on large datasets
//        // Standard deviation is recomputed for each chunk because it is more sensitive to variation
//
//        auto worker = [&](unsigned startCh, unsigned endCh)
//            {
//                for (unsigned ch = startCh; ch < endCh; ++ch)
//                {
//                    float median;
//                    double stdDev;
//                    float* chP = tf.getChannelPtr(ch);
//                    tf.computeChannelMedian(chP, tf.getNumSpectra(), median, ch);
//                    tf.computeStdDev(chP, tf.getNumSpectra(), median, stdDev);
//
//                    float adjustedThreshold = THRESHOLD * stdDev + median;
//                    unsigned flag = tf.cleanData(chP, tf.getNumSpectra(), median, adjustedThreshold);
//                    if (flag > 0) channels_flagged++;
//                }
//            };
struct RFIModule
{
    void process(TimeFrequency& tf)
    {
        unsigned channels_flagged = 0;

        unsigned numThreads = thread::hardware_concurrency();
        if (numThreads == 0) numThreads = 1;

        vector<thread> threads;
        threads.reserve(numThreads);
        size_t chPerThread = tf.getNumChannels() / numThreads;

        auto worker = [&](unsigned startCh, unsigned endCh)
            {
                for (unsigned ch = startCh; ch < endCh; ++ch)
                {
                    float median;
                    double stdDev;
                    float* chP = tf.getChannelPtr(ch);

                    tf.computeChannelMedian(chP, tf.getNumSpectra(), median, ch);
                    tf.computeStdDev(chP, tf.getNumSpectra(), median, stdDev);

                    float adjustedThreshold = THRESHOLD * stdDev + median;
                    if (tf.cleanData(chP, tf.getNumSpectra(), median, adjustedThreshold))
                        channels_flagged++;
                }
            };

            unsigned start = 0;
            for (unsigned t = 0; t < numThreads; ++t) {
                unsigned end = (t == numThreads - 1) ?
                    tf.getNumChannels() : start + chPerThread;

                threads.emplace_back(worker, start, end);
                start = end;
            }
            for (auto& th : threads) th.join();
            std::cout << "Channels flagged: " << channels_flagged << "\n";
            //return tf;
    }
};

int main(int argc, char* argv[])
{

    std::ifstream in("data.bin", std::ios::binary);
    std::ofstream out("cleaned_data.bin", std::ios::binary);

    if (!in.is_open()) {
        std::cerr << "ERROR: could not open data.bin" << std::endl;
        return 1;
    }

    else {
        std::vector<float> chunk(NUMBER_OF_SPECTRA * NUMBER_OF_CHANNELS);

        TimeFrequency tf(
            chunk,
            NUMBER_OF_CHANNELS,
            NUMBER_OF_SPECTRA,
            START_FREQUENCY,
            END_FREQUENCY,
            (END_FREQUENCY - START_FREQUENCY) / NUMBER_OF_CHANNELS,
            SAMPLING_TIME
        );

        for (unsigned i = 0; i < 2; ++i)
        {
            // Get current timestamp in order to accurately measure processing per chunk of data processed 
            auto t1 = chrono::high_resolution_clock::now();

            // Read data from buffer into vector container
            in.read(reinterpret_cast<char*>(chunk.data()), chunk.size() * sizeof(float));

            // Update the TF object with the data for the run. 
            // The data is reused so that the cache flag persists between runs.
            tf.updateData(chunk);

            // Run RFIM on the TF object.
            RFIModule rfi;
            rfi.process(tf);

            // Cache medians from first run to cut processing time for subsequent runs.
            if (i == 0) tf.setCacheValid(true);

            // Measure and print the elapsed time for processing the chunk
            auto t2 = chrono::high_resolution_clock::now();
            auto ms_int = chrono::duration_cast<chrono::milliseconds>(t2 - t1);
            std::cout << ms_int.count() << "ms\n";
            
            // Write cleaned binary data to output file
            out.write(reinterpret_cast<char*>(tf.data().data()), tf.data().size() * sizeof(float));

        }
        out.close();
        in.close();
    }
}