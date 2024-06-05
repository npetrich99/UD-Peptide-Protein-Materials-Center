
# Circular Dichroism (CD) Spectroscopy Temperature/Wavelength Scan Data Analysis

# Code written by Nolan Petrich

# This script takes raw CD data from the Jasco J-1500 and plots both forward and
# reverse CD spectra for 10 degree increments from 10 degrees Celsius to 90 
# degrees Celsius. Additionally, this code plots melting curves for that data,
# at one degree increments, at the 222 nm wavelength. This includes background 
# subtraction of blank samples. Due to the naming convention of files, the name 
# of the blank sample CSV files must be saved of the form Sample-Comment-No. The
# comment must be labeled as "Blank1", "Blank2", etc. with no spaces between
# the word and number. If this convention is not used, then the name of the 
# blank file must be edited, starting at line 235. When running the code, users 
# will be prompted to select the folder that contains the CSV data files, as 
# well as specify the sample conditions. When providing the wavelength range,
# ensure that the values are those that were measured in the CD data. There are
# also options to edit the lower and upper bounds for the x and y axis for the
# CD spectra and melting curve plots. Default values for all inputs are 
# provided.

# Importing Libraries:
library(readr)
library(tidyverse)
library(ggplot2)
library(scales)
library(extrafont)
library(ggstar)
library(shiny)
library(shinyjs)

# User Cell Information Input:
# Initialize Default Values for Each Cell:
for (i in 1:6) {
  assign(paste0("C", i, "N"), paste("Sample", i), envir = .GlobalEnv)
  assign(paste0("C", i, "C"), 0.1, envir = .GlobalEnv)
  assign(paste0("C", i, "P"), 0.1, envir = .GlobalEnv)
  assign(paste0("C", i, "R"), 29, envir = .GlobalEnv)
}

ui <- fluidPage(
  titlePanel("Enter Cell Parameter Information"),
  useShinyjs(),  # Initialize shinyjs
  sidebarLayout(
    sidebarPanel(
      numericInput("starting_wavelength", label = "Starting Wavelength (nm):",
                   value = 185, min = 150, max = 350),
      numericInput("ending_wavelength", label = "Ending Wavelength (nm):",
                   value = 250, min = 150, max = 350),
      # Dynamic UI for All Cells:
      uiOutput("cell_ui"),
      numericInput("cd_spectra_x_lower", label = "CD Spectra X Lower Bound:",
                   value = 190, min = 185, max = 350),
      numericInput("cd_spectra_x_upper", label = "CD Spectra X Upper Bound:",
                   value = 250, min = 190, max = 350),
      numericInput("cd_spectra_y_lower", label = "CD Spectra Y Lower Bound:",
                   value = -30, min = -60, max = 60),
      numericInput("cd_spectra_y_upper", label = "CD Spectra Y Upper Bound:",
                   value = 60, min = -25, max = 100),
      numericInput("melting_curve_x_lower", 
                   label = "Melting Curve X Lower Bound:",
                   value = 10, min = -1000, max = 1000),
      numericInput("melting_curve_x_upper", 
                   label = "Melting Curve X Upper Bound:",
                   value = 90, min = -1000, max = 1000),
      numericInput("melting_curve_y_lower", 
                   label = "Melting Curve Y Lower Bound:",
                   value = -25, min = -1000, max = 1000),
      numericInput("melting_curve_y_upper", 
                   label = "Melting Curve Y Upper Bound:",
                   value = 0, min = -1000, max = 1000),
      actionButton("submit", "Submit")
    ),
    mainPanel(
      textOutput("values")
    )
  )
)

server <- function(input, output, session) {
  # Number of Cells:
  num_cells <- 6
  
  # Reactive Values for Starting and Ending Wavelengths and Bounds:
  wavelengths <- reactiveValues(
    starting = 185,
    ending = 250
  )
  
  cd_spectra_bounds <- reactiveValues(
    x_lower = 190,
    x_upper = 250,
    y_lower = -30,
    y_upper = 60
  )
  
  melting_curve_bounds <- reactiveValues(
    x_lower = 10,
    x_upper = 90,
    y_lower = -25,
    y_upper = 0
  )
  
  # Generate Dynamic UI for All Cells:
  output$cell_ui <- renderUI({
    cell_ui <- lapply(1:num_cells, function(i) {
      tagList(
        textInput(paste0("sample_name_", i), label = paste("Cell", i,
                  "Sample Name:"), value = get(paste0("C", i, "N"))),
        numericInput(paste0("concentration_", i), label = paste("Cell",
                  i, "Concentration (mM):"), value = get(paste0("C", i, "C"))),
        numericInput(paste0("pathlength_", i), label = paste("Cell", i,
                  "Pathlength (cm):"), value = get(paste0("C", i, "P"))),
        numericInput(paste0("residues_", i), label = paste("Cell", i,
                  "Number of Residues:"), value = get(paste0("C", i, "R")))
      )
    })
    do.call(tagList, cell_ui)
  })
  
  # Store Values in Reactive Values and Global Environment Upon Submission:
  observeEvent(input$submit, {
    for (i in 1:num_cells) {
      values <- reactiveValues(
        C_N = paste("Sample", i),
        C_C = 0.1,
        C_P = 0.1,
        C_R = 29
      )
      
      values$C_N <- input[[paste0("sample_name_", i)]]
      values$C_C <- input[[paste0("concentration_", i)]]
      values$C_P <- input[[paste0("pathlength_", i)]]
      values$C_R <- input[[paste0("residues_", i)]]
      
      assign(paste0("C", i, "N"), values$C_N, envir = .GlobalEnv)
      assign(paste0("C", i, "C"), as.numeric(values$C_C), envir = .GlobalEnv)
      assign(paste0("C", i, "P"), as.numeric(values$C_P), envir = .GlobalEnv)
      assign(paste0("C", i, "R"), as.numeric(values$C_R), envir = .GlobalEnv)
    }
    
    # Store Starting and Ending Wavelengths:
    wavelengths$starting <- as.numeric(input$starting_wavelength)
    wavelengths$ending <- as.numeric(input$ending_wavelength)
    
    # Store CD Spectra Bounds:
    cd_spectra_bounds$x_lower <- as.numeric(input$cd_spectra_x_lower)
    cd_spectra_bounds$x_upper <- as.numeric(input$cd_spectra_x_upper)
    cd_spectra_bounds$y_lower <- as.numeric(input$cd_spectra_y_lower)
    cd_spectra_bounds$y_upper <- as.numeric(input$cd_spectra_y_upper)
    
    # Store Melting Curve Bounds:
    melting_curve_bounds$x_lower <- as.numeric(input$melting_curve_x_lower)
    melting_curve_bounds$x_upper <- as.numeric(input$melting_curve_x_upper)
    melting_curve_bounds$y_lower <- as.numeric(input$melting_curve_y_lower)
    melting_curve_bounds$y_upper <- as.numeric(input$melting_curve_y_upper)
    
    # Save Starting and Ending Wavelengths to Variables in Global Environment:
    assign("starting_wavelength", wavelengths$starting, envir = .GlobalEnv)
    assign("ending_wavelength", wavelengths$ending, envir = .GlobalEnv)
    
    # Save CD Spectra Bounds to Variables in Global Environment:
    assign("cd_spectra_x_lower", cd_spectra_bounds$x_lower, envir = .GlobalEnv)
    assign("cd_spectra_x_upper", cd_spectra_bounds$x_upper, envir = .GlobalEnv)
    assign("cd_spectra_y_lower", cd_spectra_bounds$y_lower, envir = .GlobalEnv)
    assign("cd_spectra_y_upper", cd_spectra_bounds$y_upper, envir = .GlobalEnv)
    
    # Save Melting Curve Bounds to Variables in Global Environment:
    assign("melting_curve_x_lower", melting_curve_bounds$x_lower, 
           envir = .GlobalEnv)
    assign("melting_curve_x_upper", melting_curve_bounds$x_upper, 
           envir = .GlobalEnv)
    assign("melting_curve_y_lower", melting_curve_bounds$y_lower, 
           envir = .GlobalEnv)
    assign("melting_curve_y_upper", melting_curve_bounds$y_upper, 
           envir = .GlobalEnv)
    
    # Close App:
    stopApp()
  })
  
  # Render Text Output:
  output$values <- renderText({
    # Check if Values Have Been Submitted:
    if (input$submit > 0) {
      paste(
        sapply(1:num_cells, function(i) {
          paste(
            "Cell", i, "Sample Name:", get(paste0("C", i, "N")), "\n",
            "Cell", i, "Concentration (mM):", get(paste0("C", i, "C")), "\n",
            "Cell", i, "Pathlength (cm):", get(paste0("C", i, "P")), "\n",
            "Cell", i, "Number of Residues:", get(paste0("C", i, "R")), "\n\n"
          )
        }),
        "Starting Wavelength:", starting_wavelength, "nm\n",
        "Ending Wavelength:", ending_wavelength, "nm\n",
        "CD Spectra X Lower Bound:", cd_spectra_x_lower, "\n",
        "CD Spectra X Upper Bound:", cd_spectra_x_upper, "\n",
        "CD Spectra Y Lower Bound:", cd_spectra_y_lower, "\n",
        "CD Spectra Y Upper Bound:", cd_spectra_y_upper, "\n",
        "Melting Curve X Lower Bound:", melting_curve_x_lower, "\n",
        "Melting Curve X Upper Bound:", melting_curve_x_upper, "\n",
        "Melting Curve Y Lower Bound:", melting_curve_y_lower, "\n",
        "Melting Curve Y Upper Bound:", melting_curve_y_upper
      )
    }
  })
}

shinyApp(ui = ui, server = server)

################################################################################

# Prompt User to Select Directory:
directory <- choose.dir()

# Check if Directory is Selected:
if (!is.null(directory)) {
  # Set the Selected Directory as Working Directory:
  setwd(directory)
  cat("Selected directory set as the working directory:",
      directory, "\n")
} else {
  cat("No directory selected. Exiting.\n")
  stop("No directory selected. Exiting.")
}

# Determine Maximum Rows:
final_row <- ending_wavelength - starting_wavelength + 1
HT_start <- final_row + 22
HT_end <- final_row

# Blank Data:
  # Locating Cell 1 Blank Data:
    blank_file <- list.files(".", "Blank1-1.csv")
    if (length(blank_file) == 0) {
      # File Doesn't Exist, Set CD_B1 to 0:
      CD_B1 <- data.frame(Wavelength_B1 = numeric(final_row), 
                          CD_B1 = numeric(final_row),
                          HT_B1 = numeric(final_row))
    } else {
      # File Exists, Read the Data:
      C1B <- read.csv(blank_file, skip = 20, header = FALSE)
      C1B <- C1B[1:final_row, 1:3]
      header_C1B <- c("Wavelength_B1", "CD_B1", "HT_B1")
      names(C1B) <- header_C1B
  
      # Converting Blank CD Data Columns to Numeric:
      C1B <- data.frame(lapply(C1B, function(x) {
        as.numeric(as.character(x))
      }))
    }
    
  # Locating Cell 2 Blank Data:
    blank_file <- list.files(".", "Blank2-1.csv")
    if (length(blank_file) == 0) {
      # File Doesn't Exist, Set CD_B2 to 0:
      CD_B2 <- 0
    } else {
      # File Exists, Read the Data:
      C2B <- read.csv(blank_file, skip = 20, header = FALSE)
      C2B <- C2B[1:final_row, 1:3]
      header_C2B <- c("Wavelength_B2", "CD_B2", "HT_B2")
      names(C2B) <- header_C2B
      
      # Converting Blank CD Data Columns to Numeric:
      C2B <- data.frame(lapply(C2B, function(x) {
        as.numeric(as.character(x))
      }))
    }    
    
  # Locating Cell 3 Blank Data:
    blank_file <- list.files(".", "Blank3-1.csv")
    if (length(blank_file) == 0) {
      # File Doesn't Exist, Set CD_B3 to 0:
      CD_B3 <- 0
    } else {
      # File Exists, Read the Data:
      C3B <- read.csv(blank_file, skip = 20, header = FALSE)
      C3B <- C3B[1:final_row, 1:3]
      header_C3B <- c("Wavelength_B3", "CD_B3", "HT_B3")
      names(C3B) <- header_C3B
      
      # Converting Blank CD Data Columns to Numeric:
      C3B <- data.frame(lapply(C3B, function(x) {
        as.numeric(as.character(x))
      }))
    }

  # Locating Cell 4 Blank Data:
    blank_file <- list.files(".", "Blank4-1.csv")
    if (length(blank_file) == 0) {
      # File Doesn't Exist, Set CD_B4 to 0:
      CD_B4 <- 0
    } else {
      # File Exists, Read the Data:
      C4B <- read.csv(blank_file, skip = 20, header = FALSE)
      C4B <- C4B[1:final_row, 1:3]
      header_C4B <- c("Wavelength_B4", "CD_B4", "HT_B4")
      names(C4B) <- header_C4B
      
      # Converting Blank CD Data Columns to Numeric:
      C4B <- data.frame(lapply(C4B, function(x) {
        as.numeric(as.character(x))
      }))
    }
  
    # Locating Cell 5 Blank Data:
    blank_file <- list.files(".", "Blank5-1.csv")
    if (length(blank_file) == 0) {
      # File Doesn't Exist, Set CD_B5 to 0:
      CD_B5 <- 0
    } else {
      # File Exists, Read the Data:
      C5B <- read.csv(blank_file, skip = 20, header = FALSE)
      C5B <- C5B[1:final_row, 1:3]
      header_C5B <- c("Wavelength_B5", "CD_B5", "HT_B5")
      names(C5B) <- header_C5B
      
      # Converting Blank CD Data Columns to Numeric:
      C5B <- data.frame(lapply(C5B, function(x) {
        as.numeric(as.character(x))
      }))
    }
  
    # Locating Cell 6 Blank Data:
    blank_file <- list.files(".", "Blank6-1.csv")
    if (length(blank_file) == 0) {
      # File Doesn't Exist, Set CD_B6 to 0:
      CD_B6 <- 0
    } else {
      # File Exists, Read the Data:
      C6B <- read.csv(blank_file, skip = 20, header = FALSE)
      C6B <- C6B[1:final_row, 1:3]
      header_C6B <- c("Wavelength_B6", "CD_B6", "HT_B6")
      names(C6B) <- header_C6B
      
      # Converting Blank CD Data Columns to Numeric:
      C6B <- data.frame(lapply(C6B, function(x) {
        as.numeric(as.character(x))
      }))
    }

################################################################################

# Plotting Circular Dichroism Spectra - Increasing:
  process_cell_f <- function(cell_number) {
    # Getting Ready for Exporting:
    cd_spectra <- "CD Spectra - Increasing"
    
    # Initialize p:
    p <- NULL
    
    # Check if Folder Already Exists:
    if (!file.exists(file.path(directory, cd_spectra))) {
      # If the Folder Does Not Exist, Create It and Set Output Directory:
      dir.create(file.path(directory, cd_spectra), recursive = TRUE)
      output_directory <- paste(directory, cd_spectra, sep = "/")
    } else {
      # If Folder Already Exists, Set Output Directory:
      output_directory <- paste(directory, cd_spectra, sep = "/")
    }
    
    # Locating CD Data File:
    cd_file <- list.files(".", paste0(")-Cell ", cell_number, ".csv"))
    
    # If File Exists:
    if (length(cd_file) > 0) {
      # Reading CD Data File:
      cd_f_data <- read.csv(cd_file, skip = 20, header = FALSE)
      cd_f_data <- cd_f_data[1:final_row, 1:10]
      
      # Naming Columns:
      header_cd_f <- paste0("Wavelength_F", cell_number)
      header_cd_f <- c(header_cd_f, paste0("CD_F", cell_number, "_", c("10C",
                      "20C", "30C", "40C", "50C", "60C", "70C", "80C", "90C")))
      names(cd_f_data) <- header_cd_f
      
      # Combining Raw CD - Forward Data Row-Wise:
      cd_f_combined <- cbind(cd_f_data)
      
      # Converting Data Columns to Numeric:
      cd_f_combined <- as.data.frame(lapply(cd_f_combined, as.numeric))
      
      # Subtracting Blank and Calculating MRE for the CD Data Columns:
      if (exists(paste0("C", cell_number, "B"))) {
        blank_data <- get(paste0("C", cell_number, "B"))
        cd_f_combined[, -1] <- (cd_f_combined[, -1] - 
            blank_data[[paste0("CD_B", cell_number)]]) / 
            (get(paste0("C", cell_number, "P")) * get(paste0("C", 
            cell_number, "C")) * get(paste0("C", cell_number, "R"))) * 0.1
      }
      
      # Read High Tension Voltage:
      ht_f_data <- read.csv(cd_file, skip = HT_start, header = FALSE)
      ht_f_data <- ht_f_data[1:HT_end, 1:10]
      
      # Naming Columns:
      header_ht_f <- paste0("Wavelength_F", cell_number)
      header_ht_f <- c(header_ht_f, paste0("HT_F", cell_number, "_", c("10C",
                     "20C", "30C", "40C", "50C", "60C", "70C", "80C", "90C")))
      names(ht_f_data) <- header_ht_f
      
      # Combining Raw HT - Forward Data Row-Wise:
      ht_f_combined <- cbind(ht_f_data)
      
      # Converting Data Columns to Numeric:
      ht_f_combined <- as.data.frame(lapply(ht_f_combined, as.numeric))
      
      # Filter CD Data Based on High Tension Voltage Greater Than 600 V:
      cd_f_combined <- cd_f_combined[ht_f_combined[[paste0("HT_F", cell_number,
                      "_10C")]] < 600 | is.na(ht_f_combined[[paste0("HT_F",
                      cell_number, "_10C")]]), ]
      
      # Save Data Frame as CSV:
      write.csv(cd_f_combined[, ], 
                file = paste0(output_directory, "/Cell_", cell_number,
                              "_Spectra_Forward.csv"), 
                row.names = FALSE)
      
      # Set Theme:
      theme_set(theme_bw())
      theme_update(
        text = element_text(family = "Gill Sans MT", size = 14, 
                color = "black"), # Set Font
        axis.title.x = element_text(face = "bold"), # Bold X Axis Title
        axis.title.y = element_text(face = "bold", angle = 90), 
                      # Bold Y Axis Title
        axis.text.x = element_text(face = "bold", size = 12, color = "black",
                      margin = margin(t = 5)), #Bold X Axis Text
        axis.text.y = element_text(face = "bold", size = 12, color = "black",
                      margin = margin(r = 5)), # Bold Y Axis Text
        panel.grid.major = element_blank(),  # Remove Major Gridlines
        panel.grid.minor = element_blank(),   # Remove Minor Gridlines
        panel.border = element_rect(color = "black", linewidth = 3), 
                      # Set Border to Black
        plot.background = element_rect(fill = "transparent", color = NA), 
                          # Set Plot Background to Transparent
        panel.background = element_rect(fill = "transparent", color = NA), 
                          # Set Panel Background to Transparent
        axis.ticks.length = unit(c(0.1, 0.05),
                            "inches"), 
                          # Set Lengths for Major and Minor Tick Marks
        axis.ticks = element_line(linewidth = c(1,
                    0.5), color = "black"), 
                    # Set Thicknesses for Major and Minor Tick Marks
        aspect.ratio = 1, # Make Plot Square
        legend.background = element_blank(), # Remove Legend Background
        legend.position = c(0.85, 0.72), # Position Legend
        legend.box.background = element_rect(color = "black", size = 1), 
                                # Introduce Legend Panel
        legend.key.size = unit(0.8, "lines"),  # Adjust Legend Size
        legend.text = element_text(size = 10) # Adjust Legend Text Size
      )
      
      # Plotting:
        p <- ggplot(cd_f_combined) +
          geom_hline(yintercept = 0, linetype = "solid", color = "black", 
                        size = 1) +  # Add Horizontal Line at y = 0
          geom_star(aes(x = get(paste0("Wavelength_F", cell_number)), 
                        y = get(paste0("CD_F", cell_number, "_", "10C")), 
                        starshape = "10C"), color = "#332288", fill = "#332288",
                        size = 0.5) +
          geom_line(aes(x = get(paste0("Wavelength_F", cell_number)),
                        y = get(paste0("CD_F", cell_number, "_", "10C"))),
                        color = "#332288", linetype = "solid", size = 0.5) +
          geom_star(aes(x = get(paste0("Wavelength_F", cell_number)), 
                        y = get(paste0("CD_F", cell_number, "_", "20C")), 
                        starshape = "20C"), color = "#88CCEE", fill = "#88CCEE",
                        size = 0.5) +
          geom_line(aes(x = get(paste0("Wavelength_F", cell_number)),
                        y = get(paste0("CD_F", cell_number, "_", "20C"))),
                        color = "#88CCEE", linetype = "solid", size = 0.5) +
          geom_star(aes(x = get(paste0("Wavelength_F", cell_number)), 
                        y = get(paste0("CD_F", cell_number, "_", "30C")), 
                        starshape = "30C"), color = "#44AA99", fill = "#44AA99",
                        size = 0.5) +
          geom_line(aes(x = get(paste0("Wavelength_F", cell_number)),
                        y = get(paste0("CD_F", cell_number, "_", "30C"))),
                        color = "#44AA99", linetype = "solid", size = 0.5) +
          geom_star(aes(x = get(paste0("Wavelength_F", cell_number)), 
                        y = get(paste0("CD_F", cell_number, "_", "40C")), 
                        starshape = "40C"), color = "#117733", fill = "#117733",
                        size = 0.5) +
          geom_line(aes(x = get(paste0("Wavelength_F", cell_number)),
                        y = get(paste0("CD_F", cell_number, "_", "40C"))),
                        color = "#117733", linetype = "solid", size = 0.5) +
          geom_star(aes(x = get(paste0("Wavelength_F", cell_number)), 
                        y = get(paste0("CD_F", cell_number, "_", "50C")), 
                        starshape = "50C"), color = "#999933", fill = "#999933",
                        size = 0.5) +
          geom_line(aes(x = get(paste0("Wavelength_F", cell_number)),
                        y = get(paste0("CD_F", cell_number, "_", "50C"))),
                        color = "#999933", linetype = "solid", size = 0.5) +
          geom_star(aes(x = get(paste0("Wavelength_F", cell_number)), 
                        y = get(paste0("CD_F", cell_number, "_", "60C")), 
                        starshape = "60C"), color = "#DDCC77", fill = "#DDCC77",
                        size = 0.5) +
          geom_line(aes(x = get(paste0("Wavelength_F", cell_number)),
                        y = get(paste0("CD_F", cell_number, "_", "60C"))),
                        color = "#DDCC77", linetype = "solid", size = 0.5) +
          geom_star(aes(x = get(paste0("Wavelength_F", cell_number)), 
                        y = get(paste0("CD_F", cell_number, "_", "70C")), 
                        starshape = "70C"), color = "#CC6677", fill = "#CC6677",
                        size = 0.5) +
          geom_line(aes(x = get(paste0("Wavelength_F", cell_number)),
                        y = get(paste0("CD_F", cell_number, "_", "70C"))),
                        color = "#CC6677", linetype = "solid", size = 0.5) +
          geom_star(aes(x = get(paste0("Wavelength_F", cell_number)), 
                        y = get(paste0("CD_F", cell_number, "_", "80C")), 
                        starshape = "80C"), color = "#882255", fill = "#882255",
                        size = 0.5) +
          geom_line(aes(x = get(paste0("Wavelength_F", cell_number)),
                        y = get(paste0("CD_F", cell_number, "_", "80C"))),
                        color = "#882255", linetype = "solid", size = 0.5) +
          geom_star(aes(x = get(paste0("Wavelength_F", cell_number)), 
                        y = get(paste0("CD_F", cell_number, "_", "90C")), 
                        starshape = "90C"), color = "#AA4499", fill = "#AA4499",
                        size = 0.5) +
          geom_line(aes(x = get(paste0("Wavelength_F", cell_number)),
                        y = get(paste0("CD_F", cell_number, "_", "90C"))),
                        color = "#AA4499", linetype = "solid", size = 0.5)
          
      # Scale Star Shape and Labels:
      p <- p +
        scale_starshape_manual(values = c("10C" = 15, "20C" = 11, "30C" = 23,
                                          "40C" = 13, "50C" = 5, "60C" = 6,
                                          "70C" = 7, "80C" = 8, "90C" = 9), 
                               labels = c("10 \u00b0C", "20 \u00b0C",
                                          "30 \u00b0C", "40 \u00b0C",
                                          "50 \u00b0C", "60 \u00b0C",
                                          "70 \u00b0C", "80 \u00b0C",
                                          "90 \u00b0C"))
      
      # Plot Labels and Limits:
      p <- p +
        labs(
          x = "Wavelength (nm)",
          y = expression(bold(paste("MRE", " \u00d7 10"^-3, " (deg cm"^2,
                                    "dmol"^-1, ")")))
        ) +
        scale_x_continuous(expand = c(0.004, 0),
                           breaks = function(x) seq(cd_spectra_x_lower, 
                                                    cd_spectra_x_upper, by = 5),
                           minor_breaks = NULL,
                           labels = function(x) ifelse(x %% 10 == 0,
                                                       as.character(x), ""), 
                           # Keep Labels for Major Ticks Only
                           limits = c(cd_spectra_x_lower, cd_spectra_x_upper)) +
        # Set X-Axis Limits
        scale_y_continuous(expand = c(0.004, 0),
                           breaks = function(y) seq(cd_spectra_y_lower, 
                                                    cd_spectra_y_upper, by = 5),
                           minor_breaks = function(y) seq(cd_spectra_y_lower, 
                                                    cd_spectra_y_upper, by = 5),
                           # Set Minor Tick Positions for Y-Axis
                           labels = function(y) ifelse(y %% 10 == 0,
                                                       as.character(y), ""), 
                           # Keep Labels for Major Ticks Only
                           limits = c(cd_spectra_y_lower, cd_spectra_y_upper)) +
        # Set Y-Axis Limits
        guides(starshape = guide_legend(title = NULL, 
                                        override.aes = list(size = 2))) 
                                        # Adjust Legend Point Size
      
      # Define File Name:
      filename <- paste0(output_directory, "/Cell_", cell_number,
                         "_Spectra_Forward.svg")
      
      # Export Plot with Custom Width and Height:
      ggsave(filename = filename, plot = p, device = "svg", 
             width = 4.5, height = 4)
      
    } else {
      print(paste("Cell", cell_number, "CD file not found. 
                  Skipping processing and plotting."))
    }
  }
  
  # Process Cells:
  for (cell_number in 1:6) {
    process_cell_f(cell_number = cell_number)
  }
    
################################################################################

# Plotting Melting Curves - Increasing:
  process_melt_f <- function(cell_number, directory) {
    # Getting Ready for Exporting:
    melting_plots_f <- "Melting Plots - Increasing"
    
    # Initialize p:
    p <- NULL
    
    # Check if Folder Already Exists:
    if (!file.exists(file.path(directory, melting_plots_f))) {
      # If the Folder Does Not Exist, Create It and Set Output Directory:
      dir.create(file.path(directory, melting_plots_f), recursive = TRUE)
      output_directory <- paste(directory, melting_plots_f, sep = "/")
    } else {
      # If Folder Already Exists, Set Output Directory:
      output_directory <- paste(directory, melting_plots_f, sep = "/")
    }
    
    # Locating Melting Data File:
    melt_file <- list.files(".", paste0("1-Cell ", cell_number, ".csv"))
    
    # If File Exists:
    if (length(melt_file) > 0) {
      # Reading Melting Data File:
      melt_f_data <- read.csv(melt_file, skip = 15, header = FALSE)
      melt_f_data <- melt_f_data[1:81, 1:3]
      
      # Naming Columns:
      header_melt_f <- paste0("Temperature_M", cell_number)
      header_melt_f <- c(header_melt_f, paste0("CD_M", cell_number), 
                         paste0("HT_M", cell_number))
      
      # Ensure Consistency by Setting Column Names Directly:
      colnames(melt_f_data) <- header_melt_f
      
      # Converting Data Columns to Numeric:
      melt_f_data <- as.data.frame(lapply(melt_f_data, as.numeric))
      
      # Defining Blank Value at 222 nm:
      if (exists(paste0("C", cell_number, "B"))) {
          blank_222 <- get(paste0("C", cell_number, "B"))[29, 2]
      } else {
          blank_222 <- 0
      }
      
      # Subtracting Blank and Calculating MRE for the CD Data Columns:
      if (exists("melt_f_data")) {
        melt_f_data[, 2] <- (melt_f_data[, 2] - 
                               blank_222) / 
          (get(paste0("C", cell_number, "P")) * get(paste0("C", 
                                                           cell_number, "C"))
                                  * get(paste0("C", cell_number, "R"))) * 0.1
      }
      
      # Save Data Frame as CSV:
      write.csv(melt_f_data[, c(paste0("Temperature_M", cell_number), 
                                paste0("CD_M", cell_number))], 
                file = paste0(output_directory, "/Cell_", cell_number,
                              "_Melting_Forward.csv"), 
                row.names = FALSE)
      
      # Set Theme:
      theme_set(theme_bw())
      theme_update(
        text = element_text(family = "Gill Sans MT", size = 14, 
                            color = "black"), # Set Font
        axis.title.x = element_text(face = "bold"), # Bold X Axis Title
        axis.title.y = element_text(face = "bold", angle = 90), 
        # Bold Y Axis Title
        axis.text.x = element_text(face = "bold", size = 12, color = "black",
                                   margin = margin(t = 5)), #Bold X Axis Text
        axis.text.y = element_text(face = "bold", size = 12, color = "black",
                                   margin = margin(r = 5)), # Bold Y Axis Text
        panel.grid.major = element_blank(),  # Remove Major Gridlines
        panel.grid.minor = element_blank(),   # Remove Minor Gridlines
        panel.border = element_rect(color = "black", linewidth = 3), 
        # Set Border to Black
        plot.background = element_rect(fill = "transparent", color = NA), 
        # Set Plot Background to Transparent
        panel.background = element_rect(fill = "transparent", color = NA), 
        # Set Panel Background to Transparent
        axis.ticks.length = unit(c(0.1, 0.05),
                                 "inches"), 
        # Set Lengths for Major and Minor Tick Marks
        axis.ticks = element_line(linewidth = c(1,
                                                0.5), color = "black"), 
        # Set Thicknesses for Major and Minor Tick Marks
        aspect.ratio = 1, # Make Plot Square
        legend.background = element_blank(), # Remove Legend Background
        legend.position = c(0.85, 0.72), # Position Legend
        legend.box.background = element_rect(color = "black", size = 1), 
        # Introduce Legend Panel
        legend.key.size = unit(0.8, "lines"),  # Adjust Legend Size
        legend.text = element_text(size = 10) # Adjust Legend Text Size
      )
      
      # Plotting:
      p <- ggplot(melt_f_data, aes(x = get(paste0("Temperature_M",
                                                  cell_number)),
                                   y = get(paste0("CD_M", cell_number)))) +
                                   geom_point(size = 1.5)
      
      # Plot Labels and Limits:
      p <- p +
        labs(
          x = "Temperature (\u00b0C)",
          y = expression(bold(paste("MRE"[222], " \u00d7 10"^-3, " (deg cm"^2,
                                    "dmol"^-1, ")")))
        ) +
        scale_x_continuous(expand = c(0.0045, 0),
                           breaks = function(x) seq(melting_curve_x_lower, 
                                                melting_curve_x_upper, by = 5),
                           minor_breaks = NULL,
                           labels = function(x) ifelse(x %% 10 == 0,
                                                       as.character(x), ""), 
                           # Keep Labels for Major Ticks Only
                           limits = c(melting_curve_x_lower, 
                                      melting_curve_x_upper)) + 
        # Set X-Axis Limits
        scale_y_continuous(expand = c(0.0045, 0),
                           breaks = function(y) seq(melting_curve_y_lower, 
                                    melting_curve_y_upper, by = 2.5),
                           minor_breaks = function(y) seq(melting_curve_y_lower,
                                          melting_curve_y_upper, by = 0.5),
                           # Set Minor Tick Positions for Y-Axis
                           labels = function(y) ifelse(y %% 5 == 0,
                                                       as.character(y), ""), 
                           # Keep Labels for Major Ticks Only
                           limits = c(melting_curve_y_lower, 
                                      melting_curve_y_upper)) 
                                      # Set Y-Axis Limits
      
      # Define File Name:
      filename <- paste0(output_directory, "/Cell_", cell_number,
                         "_Melting_Forward.svg")
      
      # Export Plot with Custom Width and Height:
      ggsave(filename = filename, plot = p, device = "svg", 
             width = 4, height = 4)
      
    } else {
      print(paste("Cell", cell_number,
                  "melting file not found.Skipping processing and plotting."))
    }
  }
  
  # Process Cells:
  for (cell_number in 1:6) {
    process_melt_f(cell_number = cell_number, directory = ".")
  }

################################################################################

  # Combined Plot - Increasing:
  melting_combined <- function(directory) {
    
    # Define Initial List of Colors for Each Series:
    initial_colors <- c("#EE3377", "#CC3311", "#EE7733", "#009988", "#33BBEE",
                        "#0077BB")
    
    # Define Initial List of Shapes for Data Points:
    initial_shapes <- c(21, 22, 23, 24, 25, 8)
    
    # Search Directory for Melting Files:
    melting_files <- list.files(directory, 
                                pattern = "Cell_\\d+_Melting_Forward.csv", 
                                full.names = TRUE)
    
    # Initialize List to Store Data Frames:
    melting <- list()
    
    # Loop Through Each File, Extract Data, and Store in List:
    for (file in melting_files) {
      # Read Data From File:
      data <- read.csv(file, header = TRUE)  # Ensure Header is Read
      
      # Store in List:
      melting[[file]] <- data
    }
    
    # Combine Data Frames into a Single Data Frame using dplyr's bind_rows
    combined_melting <- bind_cols(melting, .id = "series")
    
    # Set Theme:
    theme_set(theme_bw())
    theme_update(
      text = element_text(family = "Gill Sans MT", size = 14, 
                          color = "black"), # Set Font
      axis.title.x = element_text(face = "bold"), # Bold X Axis Title
      axis.title.y = element_text(face = "bold", angle = 90), 
      # Bold Y Axis Title
      axis.text.x = element_text(face = "bold", size = 12, color = "black",
                                 margin = margin(t = 5)), #Bold X Axis Text
      axis.text.y = element_text(face = "bold", size = 12, color = "black",
                                 margin = margin(r = 5)), # Bold Y Axis Text
      panel.grid.major = element_blank(),  # Remove Major Gridlines
      panel.grid.minor = element_blank(),   # Remove Minor Gridlines
      panel.border = element_rect(color = "black", linewidth = 3), 
      # Set Border to Black
      plot.background = element_rect(fill = "transparent", color = NA), 
      # Set Plot Background to Transparent
      panel.background = element_rect(fill = "transparent", color = NA), 
      # Set Panel Background to Transparent
      axis.ticks.length = unit(c(0.1, 0.05),
                               "inches"), 
      # Set Lengths for Major and Minor Tick Marks
      axis.ticks = element_line(linewidth = c(1,
                                              0.5), color = "black"), 
      # Set Thicknesses for Major and Minor Tick Marks
      aspect.ratio = 1,  # Make Plot Square
      legend.background = element_blank(), # Remove Legend Background
      legend.box.background = element_rect(color = "black", size = 1), 
      # Introduce Legend Panel
      legend.key.size = unit(0.85, "lines"),  # Adjust Legend Size
      legend.text = element_text(size = 10), # Adjust Legend Text Size
    )
    
    # Plotting:
    p <- ggplot() +
      labs(
        x = "Temperature (\u00b0C)",
        y = expression(bold(paste("MRE"[222], " \u00d7 10"^-3,
                                  " (deg cm"^2, "dmol"^-1, ")")))
      ) +
      scale_x_continuous(expand = c(0.0045, 0),
                         breaks = seq(melting_curve_x_lower, 
                                      melting_curve_x_upper, by = 5),
                         minor_breaks = NULL,
                         labels = function(x) ifelse(x %% 10 == 0,
                                                     as.character(x), ""),  
                         limits = c(melting_curve_x_lower, 
                                    melting_curve_x_upper)) + 
      scale_y_continuous(expand = c(0.0045, 0),
                         breaks = seq(melting_curve_y_lower, 
                                      melting_curve_y_upper, by = 2.5),
                         minor_breaks = seq(melting_curve_y_lower, 
                                            melting_curve_y_upper, by = 0.5),  
                         labels = function(y) ifelse(y %% 5 == 0,
                                                     as.character(y), ""),  
                         limits = c(melting_curve_y_lower, 
                                    melting_curve_y_upper))
    
    # Initialize Lists to Store Used Colors, Labels, and Shapes:
    used_colors <- character(0)
    used_labels <- character(0)
    used_shapes <- numeric(0)
    
    # Initialize Index Variable:
    actual_index <- 1
    
    # Loop Through Each Set of Data Columns and Plot Each Series:
    for (i in 1:6) {
      temperature_col <- paste0("Temperature_M", i)
      cd_col <- paste0("CD_M", i)
      
      # Check if Both Temperature and CD Columns Exist for Current Series:
      if (temperature_col %in% colnames(combined_melting) && 
          cd_col %in% colnames(combined_melting)) {
        series_data <- combined_melting %>% 
          filter(!is.na(!!rlang::sym(temperature_col)) &
                   !is.na(!!rlang::sym(cd_col)))
        
        # Plot Data for Current Series:
        p <- p + geom_line(data = series_data, aes_string(x = temperature_col,
                            y = cd_col, color = as.factor(actual_index),
                            fill = as.factor(actual_index),
                            shape = as.factor(actual_index)), size = 0.5) + 
                 geom_point(data = series_data, aes_string(x = temperature_col,
                            y = cd_col, color = as.factor(actual_index),
                            fill = as.factor(actual_index),
                            shape = as.factor(actual_index)), size = 1.5)
        
        # Add Actual Color, Label, and Shape to Used Lists:
        used_colors <- c(used_colors, initial_colors[i])
        used_labels <- c(used_labels, get(paste0("C", i, "N")))
        used_shapes <- c(used_shapes, initial_shapes[i])
        
        # Increment Actual Index:
        actual_index <- actual_index + 1
      } else {
        print(paste("Columns not found for series", i))
      }
    }
    
    # Add Legend Using Saved Information from Iteration:
    p <- p + scale_color_manual(values = used_colors, labels = used_labels,
                                name = NULL) +
      scale_fill_manual(values = used_colors, labels = used_labels,
                        name = NULL) +
      scale_shape_manual(values = used_shapes, labels = used_labels,
                         name = NULL)
    
    # Define File Name:
    filename <- paste0(directory, "/Combined_Melting_Forward.svg")
    
    # Export Plot with Custom Width and Height:
    ggsave(filename = filename, plot = p, device = "svg", width = 5,
           height = 4)
  }
  
  # Calling Combined Melting Curve Function:
  melting_combined(directory = paste0(".", "/Melting Plots - Increasing"))

################################################################################
  
# Plotting Circular Dichroism Spectra - Decreasing:
  process_cell_r <- function(cell_number) {
    # Getting Ready for Exporting:
    cd_spectra <- "CD Spectra - Decreasing"
    
    # Initialize p:
    p <- NULL
    
    # Check if Folder Already Exists:
    if (!file.exists(file.path(directory, cd_spectra))) {
      # If the Folder Does Not Exist, Create It and Set Output Directory:
      dir.create(file.path(directory, cd_spectra), recursive = TRUE)
      output_directory <- paste(directory, cd_spectra, sep = "/")
    } else {
      # If Folder Already Exists, Set Output Directory:
      output_directory <- paste(directory, cd_spectra, sep = "/")
    }
    
    # Locating CD Data File:
    cd_file <- list.files(".", paste0(")-Cell ", cell_number, "-R.csv"))
    
    # If File Exists:
    if (length(cd_file) > 0) {
      # Reading CD Data File:
      cd_r_data <- read.csv(cd_file, skip = 20, header = FALSE)
      cd_r_data <- cd_r_data[1:final_row, 1:10]
      
      # Naming Columns:
      header_cd_r <- paste0("Wavelength_R", cell_number)
      header_cd_r <- c(header_cd_r, paste0("CD_R", cell_number, "_", c("90C",
                      "80C", "70C", "60C", "50C", "40C", "30C", "20C", "10C")))
      names(cd_r_data) <- header_cd_r
      
      # Combining Raw CD - Reverse Data Row-Wise:
      cd_r_combined <- cbind(cd_r_data)
      
      # Converting Data Columns to Numeric:
      cd_r_combined <- as.data.frame(lapply(cd_r_combined, as.numeric))
      
      # Subtracting Blank and Calculating MRE for the CD Data Columns:
      if (exists(paste0("C", cell_number, "B"))) {
        blank_data <- get(paste0("C", cell_number, "B"))
        cd_r_combined[, -1] <- (cd_r_combined[, -1] - 
                                  blank_data[[paste0("CD_B", cell_number)]]) / 
          (get(paste0("C", cell_number, "P")) * get(paste0("C", 
          cell_number, "C")) * get(paste0("C", cell_number, "R"))) * 0.1
      }
      
      # Read High Tension Voltage:
      ht_r_data <- read.csv(cd_file, skip = HT_start, header = FALSE)
      ht_r_data <- ht_r_data[1:HT_end, 1:10]
      
      # Naming Columns:
      header_ht_r <- paste0("Wavelength_R", cell_number)
      header_ht_r <- c(header_ht_r, paste0("HT_R", cell_number, "_", c("10C",
                     "20C", "30C", "40C", "50C", "60C", "70C", "80C", "90C")))
      names(ht_r_data) <- header_ht_r
      
      # Combining Raw HT - Reverse Data Row-Wise:
      ht_r_combined <- cbind(ht_r_data)
      
      # Converting Data Columns to Numeric:
      ht_r_combined <- as.data.frame(lapply(ht_r_combined, as.numeric))
      
      # Filter CD Data Based on High Tension Voltage Greater Than 600 V:
      cd_r_combined <- cd_r_combined[ht_r_combined[[paste0("HT_R", cell_number,
                       "_10C")]] < 600 | is.na(ht_r_combined[[paste0("HT_R",
                       cell_number, "_10C")]]), ]
      
      # Save Data Frame as CSV:
      write.csv(cd_r_combined[, ], 
                file = paste0(output_directory, "/Cell_", cell_number,
                              "_Spectra_Reverse.csv"), 
                row.names = FALSE)
      
      # Set Theme:
      theme_set(theme_bw())
      theme_update(
        text = element_text(family = "Gill Sans MT", size = 14, 
                            color = "black"), # Set Font
        axis.title.x = element_text(face = "bold"), # Bold X Axis Title
        axis.title.y = element_text(face = "bold", angle = 90), 
        # Bold Y Axis Title
        axis.text.x = element_text(face = "bold", size = 12, color = "black",
                                   margin = margin(t = 5)), #Bold X Axis Text
        axis.text.y = element_text(face = "bold", size = 12, color = "black",
                                   margin = margin(r = 5)), # Bold Y Axis Text
        panel.grid.major = element_blank(),  # Remove Major Gridlines
        panel.grid.minor = element_blank(),   # Remove Minor Gridlines
        panel.border = element_rect(color = "black", linewidth = 3), 
        # Set Border to Black
        plot.background = element_rect(fill = "transparent", color = NA), 
        # Set Plot Background to Transparent
        panel.background = element_rect(fill = "transparent", color = NA), 
        # Set Panel Background to Transparent
        axis.ticks.length = unit(c(0.1, 0.05),
                                 "inches"), 
        # Set Lengths for Major and Minor Tick Marks
        axis.ticks = element_line(linewidth = c(1,
                                                0.5), color = "black"), 
        # Set Thicknesses for Major and Minor Tick Marks
        aspect.ratio = 1, # Make Plot Square
        legend.background = element_blank(), # Remove Legend Background
        legend.position = c(0.85, 0.72), # Position Legend
        legend.box.background = element_rect(color = "black", size = 1), 
        # Introduce Legend Panel
        legend.key.size = unit(0.8, "lines"),  # Adjust Legend Size
        legend.text = element_text(size = 10) # Adjust Legend Text Size
      )
      
      # Plotting:
      p <- ggplot(cd_r_combined) +
        geom_hline(yintercept = 0, linetype = "solid", color = "black", 
                  size = 1) +  # Add Horizontal Line at y = 0
        geom_star(aes(x = get(paste0("Wavelength_R", cell_number)), 
                      y = get(paste0("CD_R", cell_number, "_", "90C")), 
                      starshape = "90C"), color = "#AA4499", fill = "#AA4499",
                  size = 0.5) +
        geom_line(aes(x = get(paste0("Wavelength_R", cell_number)),
                      y = get(paste0("CD_R", cell_number, "_", "90C"))),
                  color = "#AA4499", linetype = "solid", size = 0.5) +
        geom_star(aes(x = get(paste0("Wavelength_R", cell_number)), 
                      y = get(paste0("CD_R", cell_number, "_", "80C")), 
                      starshape = "80C"), color = "#882255", fill = "#882255",
                  size = 0.5) +
        geom_line(aes(x = get(paste0("Wavelength_R", cell_number)),
                      y = get(paste0("CD_R", cell_number, "_", "80C"))),
                  color = "#882255", linetype = "solid", size = 0.5) +
        geom_star(aes(x = get(paste0("Wavelength_R", cell_number)), 
                      y = get(paste0("CD_R", cell_number, "_", "70C")), 
                      starshape = "70C"), color = "#CC6677", fill = "#CC6677",
                  size = 0.5) +
        geom_line(aes(x = get(paste0("Wavelength_R", cell_number)),
                      y = get(paste0("CD_R", cell_number, "_", "70C"))),
                  color = "#CC6677", linetype = "solid", size = 0.5) +
        geom_star(aes(x = get(paste0("Wavelength_R", cell_number)), 
                      y = get(paste0("CD_R", cell_number, "_", "60C")), 
                      starshape = "60C"), color = "#DDCC77", fill = "#DDCC77",
                  size = 0.5) +
        geom_line(aes(x = get(paste0("Wavelength_R", cell_number)),
                      y = get(paste0("CD_R", cell_number, "_", "60C"))),
                  color = "#DDCC77", linetype = "solid", size = 0.5) +
        geom_star(aes(x = get(paste0("Wavelength_R", cell_number)), 
                      y = get(paste0("CD_R", cell_number, "_", "50C")), 
                      starshape = "50C"), color = "#999933", fill = "#999933",
                  size = 0.5) +
        geom_line(aes(x = get(paste0("Wavelength_R", cell_number)),
                      y = get(paste0("CD_R", cell_number, "_", "50C"))),
                  color = "#999933", linetype = "solid", size = 0.5) +
        geom_star(aes(x = get(paste0("Wavelength_R", cell_number)), 
                      y = get(paste0("CD_R", cell_number, "_", "40C")), 
                      starshape = "40C"), color = "#117733", fill = "#117733",
                  size = 0.5) +
        geom_line(aes(x = get(paste0("Wavelength_R", cell_number)),
                      y = get(paste0("CD_R", cell_number, "_", "40C"))),
                  color = "#117733", linetype = "solid", size = 0.5) +
        geom_star(aes(x = get(paste0("Wavelength_R", cell_number)), 
                      y = get(paste0("CD_R", cell_number, "_", "30C")), 
                      starshape = "30C"), color = "#44AA99", fill = "#44AA99",
                  size = 0.5) +
        geom_line(aes(x = get(paste0("Wavelength_R", cell_number)),
                      y = get(paste0("CD_R", cell_number, "_", "30C"))),
                  color = "#44AA99", linetype = "solid", size = 0.5) +
        geom_star(aes(x = get(paste0("Wavelength_R", cell_number)), 
                      y = get(paste0("CD_R", cell_number, "_", "20C")), 
                      starshape = "20C"), color = "#88CCEE", fill = "#88CCEE",
                  size = 0.5) +
        geom_line(aes(x = get(paste0("Wavelength_R", cell_number)),
                      y = get(paste0("CD_R", cell_number, "_", "20C"))),
                  color = "#88CCEE", linetype = "solid", size = 0.5) +
        geom_star(aes(x = get(paste0("Wavelength_R", cell_number)), 
                      y = get(paste0("CD_R", cell_number, "_", "10C")), 
                      starshape = "10C"), color = "#332288", fill = "#332288",
                  size = 0.5) +
        geom_line(aes(x = get(paste0("Wavelength_R", cell_number)),
                      y = get(paste0("CD_R", cell_number, "_", "10C"))),
                  color = "#332288", linetype = "solid", size = 0.5)
      
      # Scale Star Shape and Labels:
      p <- p +
        scale_starshape_manual(values = c("90C" = 9, "80C" = 8, "70C" = 7,
                                          "60C" = 6, "50C" = 5, "40C" = 13,
                                          "30C" = 23, "20C" = 11, "10C" = 15), 
                               labels = c("10 \u00b0C", "20 \u00b0C",
                                          "30 \u00b0C", "40 \u00b0C",
                                          "50 \u00b0C", "60 \u00b0C",
                                          "70 \u00b0C", "80 \u00b0C",
                                          "90 \u00b0C"))
      
      # Plot Labels and Limits:
      p <- p +
        labs(
          x = "Wavelength (nm)",
          y = expression(bold(paste("MRE", " \u00d7 10"^-3, " (deg cm"^2,
                                    "dmol"^-1, ")")))
        ) +
        scale_x_continuous(expand = c(0.004, 0),
                           breaks = function(x) seq(cd_spectra_x_lower, 
                                                    cd_spectra_x_upper, by = 5),
                           minor_breaks = NULL,
                           labels = function(x) ifelse(x %% 10 == 0,
                                                       as.character(x), ""), 
                           # Keep Labels for Major Ticks Only
                           limits = c(cd_spectra_x_lower, cd_spectra_x_upper)) +
        # Set X-Axis Limits
        scale_y_continuous(expand = c(0.004, 0),
                           breaks = function(y) seq(cd_spectra_y_lower, 
                                                    cd_spectra_y_upper, by = 5),
                           minor_breaks = function(y) seq(cd_spectra_y_lower, 
                                                    cd_spectra_y_upper, by = 5),
                           # Set Minor Tick Positions for Y-Axis
                           labels = function(y) ifelse(y %% 10 == 0,
                                                       as.character(y), ""), 
                           # Keep Labels for Major Ticks Only
                           limits = c(cd_spectra_y_lower, cd_spectra_y_upper)) +
        # Set Y-Axis Limits
        guides(starshape = guide_legend(title = NULL, 
                                        override.aes = list(size = 2))) 
      # Adjust Legend Point Size
      
      # Define File Name:
      filename <- paste0(output_directory, "/Cell_", cell_number,
                         "_Spectra_Reverse.svg")
      
      # Export Plot with Custom Width and Height:
      ggsave(filename = filename, plot = p, device = "svg", 
             width = 4.5, height = 4)
      
    } else {
      print(paste("Cell", cell_number, "CD file not found. 
                  Skipping processing and plotting."))
    }
  }
  
  # Process Cells:
  for (cell_number in 1:6) {
    process_cell_r(cell_number = cell_number)
  }
  
  ##############################################################################
  
  # Plotting Melting Curves - Decreasing:
  process_melt_r <- function(cell_number, directory) {
    # Getting Ready for Exporting:
    melting_plots_r <- "Melting Plots - Decreasing"
    
    # Initialize p:
    p <- NULL
    
    # Check if Folder Already Exists:
    if (!file.exists(file.path(directory, melting_plots_r))) {
      # If the Folder Does Not Exist, Create It and Set Output Directory:
      dir.create(file.path(directory, melting_plots_r), recursive = TRUE)
      output_directory <- paste(directory, melting_plots_r, sep = "/")
    } else {
      # If Folder Already Exists, Set Output Directory:
      output_directory <- paste(directory, melting_plots_r, sep = "/")
    }
    
    # Locating Melting Data File:
    melt_file <- list.files(".", paste0("1-Cell ", cell_number, "-R.csv"))
    
    # If File Exists:
    if (length(melt_file) > 0) {
      # Reading Melting Data File:
      melt_r_data <- read.csv(melt_file, skip = 15, header = FALSE)
      melt_r_data <- melt_r_data[1:81, 1:3]
      
      # Naming Columns:
      header_melt_r <- paste0("Temperature_M", cell_number)
      header_melt_r <- c(header_melt_r, paste0("CD_M", cell_number), 
                         paste0("HT_M", cell_number))
      
      # Ensure Consistency by Setting Column Names Directly:
      colnames(melt_r_data) <- header_melt_r
      
      # Converting Data Columns to Numeric:
      melt_r_data <- as.data.frame(lapply(melt_r_data, as.numeric))
      
      # Defining Blank Value at 222 nm:
      if (exists(paste0("C", cell_number, "B"))) {
        blank_222 <- get(paste0("C", cell_number, "B"))[29, 2]
      } else {
        blank_222 <- 0
      }
      
      # Subtracting Blank and Calculating MRE for the CD Data Columns:
      if (exists("melt_r_data")) {
        melt_r_data[, 2] <- (melt_r_data[, 2] - 
                               blank_222) / 
          (get(paste0("C", cell_number, "P")) * get(paste0("C", 
           cell_number, "C")) * get(paste0("C", cell_number, "R"))) * 0.1
      }
      
      # Save Data Frame as CSV:
      write.csv(melt_r_data[, c(paste0("Temperature_M", cell_number), 
                                paste0("CD_M", cell_number))], 
                file = paste0(output_directory, "/Cell_", cell_number,
                              "_Melting_Reverse.csv"), 
                row.names = FALSE)
      
      # Set Theme:
      theme_set(theme_bw())
      theme_update(
        text = element_text(family = "Gill Sans MT", size = 14, 
                            color = "black"), # Set Font
        axis.title.x = element_text(face = "bold"), # Bold X Axis Title
        axis.title.y = element_text(face = "bold", angle = 90), 
        # Bold Y Axis Title
        axis.text.x = element_text(face = "bold", size = 12, color = "black",
                                   margin = margin(t = 5)), #Bold X Axis Text
        axis.text.y = element_text(face = "bold", size = 12, color = "black",
                                   margin = margin(r = 5)), # Bold Y Axis Text
        panel.grid.major = element_blank(),  # Remove Major Gridlines
        panel.grid.minor = element_blank(),   # Remove Minor Gridlines
        panel.border = element_rect(color = "black", linewidth = 3), 
        # Set Border to Black
        plot.background = element_rect(fill = "transparent", color = NA), 
        # Set Plot Background to Transparent
        panel.background = element_rect(fill = "transparent", color = NA), 
        # Set Panel Background to Transparent
        axis.ticks.length = unit(c(0.1, 0.05),
                                 "inches"), 
        # Set Lengths for Major and Minor Tick Marks
        axis.ticks = element_line(linewidth = c(1,
                                                0.5), color = "black"), 
        # Set Thicknesses for Major and Minor Tick Marks
        aspect.ratio = 1, # Make Plot Square
        legend.background = element_blank(), # Remove Legend Background
        legend.position = c(0.85, 0.72), # Position Legend
        legend.box.background = element_rect(color = "black", size = 1), 
        # Introduce Legend Panel
        legend.key.size = unit(0.8, "lines"),  # Adjust Legend Size
        legend.text = element_text(size = 10) # Adjust Legend Text Size
      )
      
      # Plotting:
      p <- ggplot(melt_r_data, aes(x = get(paste0("Temperature_M",
                                                  cell_number)),
                                   y = get(paste0("CD_M", cell_number)))) +
        geom_point(size = 1.5)
      
      # Plot Labels and Limits:
      p <- p +
        labs(
          x = "Temperature (\u00b0C)",
          y = expression(bold(paste("MRE"[222], " \u00d7 10"^-3, " (deg cm"^2,
                                    "dmol"^-1, ")")))
        ) +
        scale_x_continuous(expand = c(0.0045, 0),
                           breaks = function(x) seq(melting_curve_x_lower, 
                                                melting_curve_x_upper, by = 5),
                           minor_breaks = NULL,
                           labels = function(x) ifelse(x %% 10 == 0,
                                                       as.character(x), ""), 
                           # Keep Labels for Major Ticks Only
                           limits = c(melting_curve_x_lower, 
                                      melting_curve_x_upper)) + 
        # Set X-Axis Limits
        scale_y_continuous(expand = c(0.0045, 0),
                           breaks = function(y) seq(melting_curve_y_lower, 
                                              melting_curve_y_upper, by = 2.5),
                           minor_breaks = function(y) seq(melting_curve_y_lower, 
                                          melting_curve_y_upper, by = 0.5),
                           # Set Minor Tick Positions for Y-Axis
                           labels = function(y) ifelse(y %% 5 == 0,
                                                       as.character(y), ""), 
                           # Keep Labels for Major Ticks Only
                           limits = c(melting_curve_y_lower, 
                                      melting_curve_y_upper)) 
      # Set Y-Axis Limits
      
      # Define File Name:
      filename <- paste0(output_directory, "/Cell_", cell_number,
                         "_Melting_Reverse.svg")
      
      # Export Plot with Custom Width and Height:
      ggsave(filename = filename, plot = p, device = "svg", 
             width = 4, height = 4)
      
    } else {
      print(paste("Cell", cell_number,
                  "melting file not found.Skipping processing and plotting."))
    }
  }
  
  # Process Cells:
  for (cell_number in 1:6) {
    process_melt_r(cell_number = cell_number, directory = ".")
  }
  
  ##############################################################################
  
  # Combined Plot - Decreasing:
  melting_combined <- function(directory) {
    
    # Define Initial List of Colors for Each Series:
    initial_colors <- c("#EE3377", "#CC3311", "#EE7733", "#009988", "#33BBEE",
                        "#0077BB")
    
    # Define Initial List of Shapes for Data Points:
    initial_shapes <- c(21, 22, 23, 24, 25, 8)
    
    # Search Directory for Melting Files:
    melting_files <- list.files(directory, 
                                pattern = "Cell_\\d+_Melting_Reverse.csv", 
                                full.names = TRUE)
    
    # Initialize List to Store Data Frames:
    melting <- list()
    
    # Loop Through Each File, Extract Data, and Store in List:
    for (file in melting_files) {
      # Read Data From File:
      data <- read.csv(file, header = TRUE)  # Ensure Header is Read
      
      # Store in List:
      melting[[file]] <- data
    }
    
    # Combine Data Frames into a Single Data Frame using dplyr's bind_rows
    combined_melting <- bind_cols(melting, .id = "series")
    
    # Set Theme:
    theme_set(theme_bw())
    theme_update(
      text = element_text(family = "Gill Sans MT", size = 14, 
                          color = "black"), # Set Font
      axis.title.x = element_text(face = "bold"), # Bold X Axis Title
      axis.title.y = element_text(face = "bold", angle = 90), 
      # Bold Y Axis Title
      axis.text.x = element_text(face = "bold", size = 12, color = "black",
                                 margin = margin(t = 5)), #Bold X Axis Text
      axis.text.y = element_text(face = "bold", size = 12, color = "black",
                                 margin = margin(r = 5)), # Bold Y Axis Text
      panel.grid.major = element_blank(),  # Remove Major Gridlines
      panel.grid.minor = element_blank(),   # Remove Minor Gridlines
      panel.border = element_rect(color = "black", linewidth = 3), 
      # Set Border to Black
      plot.background = element_rect(fill = "transparent", color = NA), 
      # Set Plot Background to Transparent
      panel.background = element_rect(fill = "transparent", color = NA), 
      # Set Panel Background to Transparent
      axis.ticks.length = unit(c(0.1, 0.05),
                               "inches"), 
      # Set Lengths for Major and Minor Tick Marks
      axis.ticks = element_line(linewidth = c(1,
                                              0.5), color = "black"), 
      # Set Thicknesses for Major and Minor Tick Marks
      aspect.ratio = 1,  # Make Plot Square
      legend.background = element_blank(), # Remove Legend Background
      legend.box.background = element_rect(color = "black", size = 1), 
      # Introduce Legend Panel
      legend.key.size = unit(0.85, "lines"),  # Adjust Legend Size
      legend.text = element_text(size = 10), # Adjust Legend Text Size
    )
    
    # Plotting:
    p <- ggplot() +
      labs(
        x = "Temperature (\u00b0C)",
        y = expression(bold(paste("MRE"[222], " \u00d7 10"^-3,
                                  " (deg cm"^2, "dmol"^-1, ")")))
      ) +
      scale_x_continuous(expand = c(0.0045, 0),
                         breaks = seq(melting_curve_x_lower, 
                                      melting_curve_x_upper, by = 5),
                         minor_breaks = NULL,
                         labels = function(x) ifelse(x %% 10 == 0,
                                                     as.character(x), ""),  
                         limits = c(melting_curve_x_lower, 
                                    melting_curve_x_upper)) + 
      scale_y_continuous(expand = c(0.0045, 0),
                         breaks = seq(melting_curve_y_lower, 
                                      melting_curve_y_upper, by = 2.5),
                         minor_breaks = seq(melting_curve_y_lower, 
                                            melting_curve_y_upper, by = 0.5),  
                         labels = function(y) ifelse(y %% 5 == 0,
                                                     as.character(y), ""),  
                         limits = c(melting_curve_y_lower, 
                                    melting_curve_y_upper))
    
    # Initialize Lists to Store Used Colors, Labels, and Shapes:
    used_colors <- character(0)
    used_labels <- character(0)
    used_shapes <- numeric(0)
    
    # Initialize Index Variable:
    actual_index <- 1
    
    # Loop Through Each Set of Data Columns and Plot Each Series:
    for (i in 1:6) {
      temperature_col <- paste0("Temperature_M", i)
      cd_col <- paste0("CD_M", i)
      
      # Check if Both Temperature and CD Columns Exist for Current Series:
      if (temperature_col %in% colnames(combined_melting) && 
          cd_col %in% colnames(combined_melting)) {
        series_data <- combined_melting %>% 
          filter(!is.na(!!rlang::sym(temperature_col)) &
                   !is.na(!!rlang::sym(cd_col)))
        
        # Plot Data for Current Series:
        p <- p + geom_line(data = series_data, aes_string(x = temperature_col,
                                y = cd_col, color = as.factor(actual_index),
                                fill = as.factor(actual_index),
                                shape = as.factor(actual_index)), size = 0.5) + 
          geom_point(data = series_data, aes_string(x = temperature_col,
                                y = cd_col, color = as.factor(actual_index),
                                fill = as.factor(actual_index),
                                shape = as.factor(actual_index)), size = 1.5)
        
        # Add Actual Color, Label, and Shape to Used Lists:
        used_colors <- c(used_colors, initial_colors[i])
        used_labels <- c(used_labels, get(paste0("C", i, "N")))
        used_shapes <- c(used_shapes, initial_shapes[i])
        
        # Increment Actual Index:
        actual_index <- actual_index + 1
      } else {
        print(paste("Columns not found for series", i))
      }
    }
    
    # Add Legend Using Saved Information from Iteration:
    p <- p + scale_color_manual(values = used_colors, labels = used_labels,
                                name = NULL) +
      scale_fill_manual(values = used_colors, labels = used_labels,
                        name = NULL) +
      scale_shape_manual(values = used_shapes, labels = used_labels,
                         name = NULL)
    
    # Define File Name:
    filename <- paste0(directory, "/Combined_Melting_Reverse.svg")
    
    # Export Plot with Custom Width and Height:
    ggsave(filename = filename, plot = p, device = "svg", width = 5,
           height = 4)
  }
  
  # Calling Combined Melting Curve Function:
  melting_combined(directory = paste0(".", "/Melting Plots - Decreasing"))
