# Aero_ML: An aeroacoustic Machine Learning Repository

![example](example.png)


This machine learning project's goal is to take raw audio data from Cal State Fullerton's anechoic Wind Tunnel, training a machine learning model with 2 microphone inputs and output the performance of 7 beamformed microphones (see image above). This project repo is in MATLAB using the following dependencies:

- Audio Toolbox
- Deep Learning Toolbox
- MATLAB
- Parallel Computing Toolbox
- Signal Processing Toolbox

Future work on this project will include the code transition into Python using Keras with a Tensorflow backend, and further development to increase the accuracy of the model(s).

## Project Overview

Overall this graduate project was a compilation of 3 items. The first was the renovation of the test section and addition of a permanent datum in the form of a dowel pin which also serves as an alignment feature for top and bottom panels in the wind tunnel (WT) test section. The second is the enhancement of the anechoic setup. Lastly and most importantly (required by the 2 previous items) collect a validation set of testing data and train an LSTM or BiLSTM to have a prediction outputted from data that hasn't been seen before by the model.

This git will have a small sample of test files and couple .mat files with most of the data set read into a .mat file structure, however this is just to ensure all the scripts run without error.

### File Setup
ML_data
├── Test data 1
│   ├── MWS200V40P68R0.dat
│   ├── MWS200V40P68R1.dat
│   ├── MWS200V40P68R2.dat
│   ├── MWS200V40P68R3.dat
│   └── MWS200V40P68R4.dat
└── Test data 2
    ├── MLD200V40P68R0.dat
    ├── MLD200V40P68R1.dat
    ├── MLD200V40P68R2.dat
    ├── MLD200V40P68R3.dat
    └── MLD200V40P68R4.dat

Created by: Andrew Bartels, Graduate student Cal State, Fullerton
