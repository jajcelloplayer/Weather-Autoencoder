# Weather-Autoencoder
Code and Data for a spatially varrying Bayesian autoencoder. 
This autoencoder is used to combine and meld multiple climate data products, such as those provided in the data folder.

## Code Files
#### AMCMCUpdate.R
RScript containing the adaptive MCMC function used to train the autoencoder. 

#### AutoencoderScript.R
RScript containing code to run and train the Melding Autoencoder. 

Options to configure the Neural Network are included at the top of the document.
Next, there is an area to specify the prior parameters for the spatial basis function weights (Normal) and the common recreation error (Inverse Gamma).
At the end of the script, there is code to create plots of values of interest. These plots can be customized, changed, or removed without changing the function of the autoencoder itself. 

#### month_key.txt
Table showing which data products are available for each month from January 1998 to September 2018. See Data Files section below for details on each of the data products included.

## Data Files
Rdata files containing the data that motivated the development of this autoencoder. 
Each of the 4 files includes gridded monthly precipitation data for the High Mountain Asia region. Each file and its respective data source are explained below. 

#### aligned_era5.Rdata
Data from the ERA5 reanalysis. Data is included for April-September from 1998 to 2018.

#### aligned_merra2.Rdata
Data from the MERRA2 reanalysis. Data is included for April-September from 1998 to 2017.

#### aphro.Rdata
Data from the APHRODITE reanalysis. Data is included for all months from January 1998 to December 2015.

#### trmm.Rdata
Data from the Tropical Rainfall Measuring Mission. Data is included for all months from January 1998 to December 2017.

## Additional Files
Supplimentary files supporting the work presented in the paper detailing the autoencoder. 

#### ConsensusMovie.mp4
A 'time lapse' showing the autoencoder consensus for all months from January 1998 to December 2017.
