# DiDo

A machine learning algorithm to assess the brain injury severity of newborns in the NICU in real time, developed by Danai Bili in 2021-2022. 

# Overview

This project is based on 2 axes:

1. data analysis
2. real time code development

--> data analysis: run feature extraction and clustering separately, to save time, especially when using files from multiple patients
--> real time: run dido.py

# Central Scripts

convertdata.m : MATLAB script to turn bNIRS and DCS processed data into a readable format

datapanda.py : feature extraction. Use 'florence' when you have bNIRS and DCS data, 'cyril' when you have only bNIRS data and 'live' when running the real time script

pipeline_core.py : run, tune and investigate patient clustering, ideally as part of the data analysis methodology

dido.py: run real time algorithm


# Other scripts

All scripts developed during the project are deposited here for completeness.

# Documentation

dido_manual.pdf : instructions on how to use dido.py in the clinic

thesis.pdf : available upon request at danaibili2@gmail.com, complete description of experiment
