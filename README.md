Source code for the implementation of the model described in

# Stochastic parabolic growth promotes coexistence and a relaxed error threshold in RNA-like replicator populations

Mátyás Paczkó<sup>1,2</sup>, Eörs Szathmáry<sup>1,3,4,*</sup>, András Szilágyi<sup>1</sup>

<sup>1</sup>HUN-REN Centre for Ecological Research, Institute of Evolution, 1121 Budapest, Hungary

<sup>2</sup>Doctoral School of Biology, Institute of Biology, ELTE Eötvös Loránd University, 1117
Budapest, Hungary

<sup>3</sup>Center for the Conceptual Foundations of Science, Parmenides Foundation, 82343 Pöcking,
Germany

<sup>4</sup>Department of Plant Systematics, Ecology and Theoretical Biology, Eötvös Loránd
University, 1117 Budapest, Hungary

<sup>*</sup>Correspondence: szathmary.eors@ecolres.hu

Note: The code available from this repository is provided without any warranty or liability.

## Overview

Individual-based, stochastic model framework of parabolic replication for investigating the diversity maintaining ability of sub-exponentially growing replicator systems and the potentially relaxed copying error threshold that characterizes such dynamics.

## Prerequisites

The code available from this repository has been tested on the OS "Linux Mint 18.3" using the "4.15.0-142-generic" kernel with GCC version 5.4.0.

## Basic usage

1. Move to a directory of your choice (in console):
   
   `cd /home/user/my_directory`
   
3. Compile C file using GCC optimizer to yield faster<sup>*</sup> simulation:
   
   `gcc filename.c -lm -Wall -Ofast -o optional_executable_file_name`
   
5. Run the resulting executable file (using nohup to run the process in the background):
   
   `nohup ./optional_executable_file_name &`

<sup>*</sup>Optimizer usage is recommended to reduce execution time, because large population sizes and replication event numbers - as termination condition (corresponding to the combination of *N*=10<sup>5</sup> and *MAX_REPLICATIONS*=10<sup>7</sup> parameter values) - result in considerably increased execution times.
