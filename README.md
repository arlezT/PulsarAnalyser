# PulsarAnalyser

PulsarAnalyser is a C++ application designed for time-frequency analysis of radio-frequency signal.  
It reads raw binary input, computes per-channel statistics, removes outliers, and writes a
cleaned dataset to disk. The project was created as part of a technical assignment for 
the role of Data Processing Software Engineer.

## Features

- Per-channel median and standard deviation calculation  
- Outlier detection and replacement using median-based thresholds  
- Optional caching of channel medians to accelerate subsequent iterations  
- Multi-threaded processing (parallel per-channel cleaning)  
- Cleaned data written back to binary format  
- Clear separation between main application logic and processing routines (`TimeFrequency` class)

## Project Structure
```
PulsarAnalyser/
│
├── main.cpp                # Entry point and main processing loop
├── TimeFrequency.h/.cpp    # Per-channel statistical processing
├── unitTests.cpp           # Unit test scaffolding (manual tests)
├── CMakeLists.txt          # Build configuration
└── README.md               # Project documentation
```

## Building
The project uses CMake and can be built on both Linux and Windows.


## Running
The program expects a binary input file named data.bin, placed in the same working directory as the executable.

