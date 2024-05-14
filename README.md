The data processing of circular dichroism (CD) data can be a time-consuming process. The following R script has been written to aid in this process. To use this R script properly, follow these instructions for naming and exporting data, as well as how to download and use the script.
Processing and Exporting Data:
This script asks the user to identify a folder where all the CD data is contained. From there, it uses the names of the files to identify CSV data files. When setting the CD parameters, as well as when exporting the data, ensure the following conditions are met.
    1.	When running Spectra Measurement for blank samples, in the Information tab, provide the sample name in the associated column. In the comment column, type exactly as written “Blank1”, “Blank2”, etc. with no space between “Blank” and the number. Under the Data tab, select Sample-Comment-No. data Format.
    2.	When running the Temperature/Wavelength Scan, in the Information tab, provide the sample name in the associated column. Under the Data tab, select the Sample-No. data Format.

Downloading and Using R:
For those that have never used R, the following instructions provide enough information to allow users to download and implement the CD data analysis script. R is an open-source software, which means it is free and accessible to anyone.
    1.	To download R, access the following website:
        a.	https://cran.r-project.org
        b.	Follow the instructions on that website. Users can download R for Linux, macOS, or Windows at the top of the screen.
    2.	Next, users should download RStudio. Initially downloading R provides the programming language; however, RStudio is needed as an interface to run and write scripts. RStudio can be downloaded by accessing the following website:
        a.	https://posit.co/download/rstudio-desktop/#download
        b.	For Windows users, the blue button on the right underneath “Install RStudio” can be clicked to start the download.
        c.	For macOS users, scroll down the page to find the proper file to download.
    3.	Once these two steps have been completed, users should be able to open the R script in RStudio.

Using the Script in R:
The script provides general instructions for running the data analysis; however, users will encounter an error when running the code for the first time if they haven’t used R before. R uses libraries of code for different tasks. Because of this, these packages often need to be installed. This should be done the first time before running the code, then it does not need to be performed again:
    1.	In the Console in RStudio, type the following line:
        a.	install.packages("readr")
    2.	Press Ctrl + Enter
        a.	This will install the “readr” package. Once this is completed, the console will display “package ‘readr’ successfully unpacked and MD5 sums checked”.
    3.	Repeat this process for the following other packages:
        a.	tidyverse
        b.	ggplot2
        c.	scales
        d.	extrafont
        e.	ggstar
        f.	shiny
        g.	shinyjs
        h.	svglite
    4.	The script is now ready to be run and this process does not need to be repeated.

Once you have completed the package installation once, to use the R script, users just need to open the file in RStudio, ensure that all their CSV files are properly named in a folder, and press Code > Run Region > Run All (or hold Ctrl + Alt + R). The code will prompt you to input your cell parameters and select the folder with your files in it. Then the rest of the code will be run to calculate MRE and save plots of the data. For those that want to also plot in Origin or some other plotting software, the calculated MRE is saved in separate CSV files for each cell.
