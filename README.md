# Algamal2022_analysis_w_OASIS
This is how I used OASIS (https://github.com/zhoupc/OASIS_matlab) to analyze the data reported in Algamal et al. 2022, doi: https://doi.org/10.1101/2022.04.27.489759

To install

1- Download OASIS_matlab (https://github.com/zhoupc/OASIS_matlab)

2- Download the files and folders included in this repository (Algamal2022_analysis_w_OASIS)

3- Copy all the files and folders in Algamal2022_analysis_w_OASIS, then paste them inside the folder (OASIS_matlab-master)

4- You can then use the main file (CallCa) to analyze your data, make sure you cahnge the path to the path where OASIS_matlab-master is located (line 1 in CallCa)

To analyze data imported from Suite2p, do the following.

1- Save the Suit2p data to Matlab. This will create a Matlab file named Fall.mat

2- Rename Fall.mat to your desired name for example, exp_1.mat

3- Move the file to the 'data' folder

4- Open CallCa file in Matlab

5- Make sure you cahnge the path to the path where OASIS_matlab-master is located (line 1 in CallCa)

6- Type in your file name next to name in line 2. For example name = 'exp_1';

7- Choose the 'smin' value, minimum event size relevant to noise standard deviation, default is -2.5, which means two and half times SD

8- Optional; change lambda value if known

8-Click Run

9- Event rates will be exported to the analyzed folder as an excel sheet

10- Pearson correlation matices will be stored in the correlations folder



Alternatively, you can store your time series data in a Matlab file (e.g exp_1. mat) that contains the two following variables

1; 'data' ---> fluorescence values over time

2; 'dt' --> frame time (e.g 0.1 s)

Then follow steps 3-8
