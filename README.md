# Aero_ML: An aeroacoustic Machine Learning Repository

![example](example.png)


## Project Overview
This machine learning project's goal is to take raw audio data from Cal State Fullerton's anechoic Wind Tunnel, training a machine learning model with 2 microphone inputs and output the performance of 7 beamformed microphones (see image above). 

Overall this graduate project was a compilation of 3 items. 

1) Renovation of the test section and addition of a permanent datum in the form of a dowel pin which also serves as an alignment feature for top and bottom panels in the wind tunnel (WT) test section. 

2) Enhancement of the anechoic setup. 

3) Collect a validation set of testing data and train an LSTM or BiLSTM to have a prediction outputted from data that hasn't been seen before by the model. (Items 1 & 2 being dependancies))

This git will have a small sample of test files and couple .mat files with most of the data set read into a .mat file structure, however this is just to ensure all the scripts run without error.

# Data Workflow
The workflow the data is as follows:

#### Path Setup and Harness
1) Data Acq in LabVIEW generate files in the `RAW` folder under different sections. See [exp_list.xlsx](./RAW/exp_list.xlsx)

2) Data files are read into Matlab via the `LSTM_harness.m` where all file paths are specified and is a harness that runs the other files needed to preprocess and run the data.

#### Preprocessing and cleaning methods
3) Once the datafiles (`.dat`) are read in, they are beamformed for 2 and 7 microphones (the 7 columns in the `.dat` files). 

4) Once the signals are calibrated for the response of each microphone, mel-spectrogram features are generated for the LSTM to key off the coefficients

5) The model is called and trained with the corresponding layers in the `LSTM_model.m` file

#### Output
6) This is a research project, so the logging is pretty verbose, a user will see multiple plots and files if the training harness is run locally.


For further details on the project see [Final Defense Paper](Defense_paper/Bartels_Project_FINAL.pdf)


### File Setup
```console
RAW
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
```



# Dependancies
This project repo is in MATLAB using the following dependencies:

- Audio Toolbox
- Deep Learning Toolbox
- MATLAB
- Parallel Computing Toolbox
- Signal Processing Toolbox

Future work on this project will include the code transition into Python using Keras with a Tensorflow backend, and further development to increase the accuracy of the model(s).


# Author
- Andrew Bartels, Graduate student Cal State, Fullerton 2019



# License

Copyright © 2019 [Andrew Bartels](https://github.com/kefranabg).<br />
This project is [MIT](https://https://github.com/andrewbartels1/Aero_ML/LICENSE) licensed.