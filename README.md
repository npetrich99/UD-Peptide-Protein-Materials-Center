# Automated Circular Dichroism (CD) Data Analysis in R:

The data processing of circular dichroism (CD) data can be a time-consuming process. These R scripts have been written to aid in this process. To use these R scripts properly, follow these instructions for naming and exporting data, as well as how to download and use the scripts.

## Processing and Exporting Data:

Both scripts ask the user to identify a folder where all the CD data is contained. From there, it uses the names of the files to identify CSV data files. When setting the CD parameters, as well as when exporting the data, ensure the following conditions are met.

1. When running Spectra Measurement for blank samples, in the Information tab, provide the Sample Name in the associated column. In the Comment column, type exactly as written “Blank1”, “Blank2”, etc. with no space between “Blank” and the cell number. Under the Data tab, select Sample-Comment-No. data Format.
   
2. If running the Spectra Measurement for a sample, in the Information tab, provide the Sample Name in the associated column. In the Comment column, type exactly as written "Sample1", "Sample2", etc. with no space between "Sample" and the cell number. Under the Data tab, select Sample-Comment-No. data Format.
   
3. If running the Temperature/Wavelength Scan, in the Information tab, provide the Sample Name in the associated column. Under the Data tab, select the Sample-No. data Format.

## Downloading and Using R:

For those that have never used R, the following instructions provide enough information to allow users to download and implement the CD data analysis scripts. R is an open-source software, which means it is free and accessible to anyone.

1. To download R, access the following website:
      <ol type="a">
      <li>https://cran.r-project.org</li>
      <li>Follow the instructions on that website. Users can download R for Linux, macOS, or Windows at the top of the screen.</li>
      </ol>

2. Next, users should download RStudio. Initially downloading R provides the programming language; however, RStudio is needed as an interface to run and write scripts. RStudio can be downloaded by accessing the following website:

      <ol type="a">
      <li>https://posit.co/download/rstudio-desktop/#download</li>
      <li>For Windows users, the blue button on the right underneath “Install RStudio” can be clicked to start the download.</li>
      <li>For macOS users, scroll down the page to find the proper file to download.</li>
      </ol>

3. Once these two steps have been completed, users should be able to open the R scripts in RStudio.

## Using the Scripts in RStudio:

The scripts provide general instructions for running the data analysis (either for single CD spectra in Spectra Measurement or for multiple CD spectra at different temperatures and melting curves in Temperature/Wavelength Scans); however, users will encounter an error when running these scripts for the first time if they haven’t used R before. R uses libraries of code for different tasks. Because of this, these packages often need to be installed. This should be done the first time before running the code, then it does not need to be performed again. This can simply be done by running the script titled "Initial R Package Installation.R". Otherwise, the following instructions can be followed manually in the RStudio Console:

1.	In the Console in RStudio, type the following line:

      <ol type="a">
  	   <li>install.packages("readr")</li>
      </ol>
        
3.	Press Ctrl + Enter

      <ol type="a">
      <li>This will install the “readr” package. Once this is completed, the console will display “package ‘readr’ successfully unpacked and MD5 sums checked”.</li>
      </ol>
        
3.	Repeat this process for the following other packages:

      <ol type="i">
      <li>tidyverse</li>
      <li>ggplot2</li>
      <li>scales</li>
      <li>extrafont</li>
      <li>ggstar</li>
      <li>shiny</li>
      <li>shinyjs</li>
      <li>svglite</li>
      </ol>


4.	The script is now ready to be run and this process does not need to be repeated.

Once you have completed the package installation once, to use the R scripts, users just need to open the file in RStudio, ensure that all their CSV files are properly named in a folder, and press Code > Run Region > Run All (or hold Ctrl + Alt + R). The code will prompt users to input their cell parameters and select the folder with their files in it. Then the rest of the code will be run to calculate MRE and save plots of the data. For those that want to also plot in Origin or some other plotting software, the calculated MRE is saved in separate CSV files for each cell.
