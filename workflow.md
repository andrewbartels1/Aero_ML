# Workflow thoughts
How to generate cleaned data from the raw inputs, fed through a Neural Network

 Description: This is a large scale overview of how the workflow process will
 be though out
##  PreProcessing Class Script

 Dir full of auido raw .dat files

# Even more in the weeds:
 # book (list of dicts) where the layout is:
CRTL + ALT + R


    RAWDATA dir
    ├── TSC (Test Section Calibration)
    │   └── TSC[RUN]V[Velocity]P[Position]R[rep_num]*.dat
    └── MWS (momentum wake sheilding)
        └── MWS[RUN]V[Velocity]P[Position]R[rep_num]*.dat
# Take files and sort into something logical

# Turn all the files into a series of Time-Frequency Representations

# Preprocess and sort for loading into model

 #TODO: use the expandauser to make everything flag inputs into the preprocessing
# claa
# =============================================================================
#  Model Script
# =============================================================================

# Some sort of CNN/Deep Learning network (AlexNet replica?)

"""Possible inputs that could work into the model:
    1) souce signal algorithm matrix input
    2) Time-Freq representations
    3) Any other ideas... predict on raw audio data?

"""
# Make sure learning gradient is recorded etc. (Maybe use Keras?)

# Have the model output a classification (maybe regression) of what is wind tunnel
# noise and what isn't
