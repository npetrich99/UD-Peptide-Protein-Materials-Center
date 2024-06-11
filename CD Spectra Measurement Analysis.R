
# Circular Dichroism (CD) Spectroscopy Spectra Measurement Data Analysis

# Code written by Nolan Petrich

# This script takes raw CD data from the Jasco J-1500 and plots a CD spectrum.
# This includes background subtraction of blank samples. Due to the naming 
# convention of files, the name of the blank sample CSV files must be saved of 
# the form Sample-Comment-No. The comment must be labeled as "Blank1", "Blank2",
# etc. with no spaces between the word and number. If this convention is not 
# used, then the name of the blank file must be edited, starting at line 196. 
# Additionally, a similar naming convention must be used for the actual sample
# CSV files. They should be saved of the form Sample-Comment-No. with the 
# comment labeled as "Sample1", "Sample2", etc. with no spaces between the word 
# and number. This convention must be used for this code. When running the 
# script, users will be prompted to select the folder that contains the CSV data
# files, as well as specify the sample conditions. When providing the wavelength 
# range, ensure that the values are those that were measured in the CD data. 
# There are also options to edit the lower and upper bounds for the x and y axis 
# for the CD spectra plots. Default values for all inputs are provided.

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
    
    # Save Starting and Ending Wavelengths to Variables in Global Environment:
    assign("starting_wavelength", wavelengths$starting, envir = .GlobalEnv)
    assign("ending_wavelength", wavelengths$ending, envir = .GlobalEnv)
    
    # Save CD Spectra Bounds to Variables in Global Environment:
    assign("cd_spectra_x_lower", cd_spectra_bounds$x_lower, envir = .GlobalEnv)
    assign("cd_spectra_x_upper", cd_spectra_bounds$x_upper, envir = .GlobalEnv)
    assign("cd_spectra_y_lower", cd_spectra_bounds$y_lower, envir = .GlobalEnv)
    assign("cd_spectra_y_upper", cd_spectra_bounds$y_upper, envir = .GlobalEnv)
    
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

# Plotting Circular Dichroism Spectrum:
process_cell <- function(cell_number) {
  # Getting Ready for Exporting:
  cd_spectrum <- "CD Spectra"
  
  # Initialize p:
  p <- NULL
  
  # Check if Folder Already Exists:
  if (!file.exists(file.path(directory, cd_spectrum))) {
    # If the Folder Does Not Exist, Create It and Set Output Directory:
    dir.create(file.path(directory, cd_spectrum), recursive = TRUE)
    output_directory <- paste(directory, cd_spectrum, sep = "/")
  } else {
    # If Folder Already Exists, Set Output Directory:
    output_directory <- paste(directory, cd_spectrum, sep = "/")
  }
  
  # Locating CD Data File:
  cd_file <- list.files(".", paste0("Sample", cell_number, "-1.csv"))
  
  # If File Exists:
  if (length(cd_file) > 0) {
    # Reading CD Data File:
    cd_data <- read.csv(cd_file, skip = 20, header = FALSE)
    cd_data <- cd_data[1:final_row, 1:3]
    
    # Naming Columns:
    header_cd <- paste0("Wavelength_", cell_number)
    header_cd <- c(header_cd, paste0("CD_", cell_number), 
                   paste0("HT_", cell_number))
    names(cd_data) <- header_cd
    
    # Combining Raw CD - Data Row-Wise:
    cd_combined <- cbind(cd_data)
    
    # Converting Data to Numeric:
    cd_combined <- as.data.frame(lapply(cd_combined, as.numeric))
    
    # Subtracting Blank and Calculating MRE for the CD Data:
    if (exists(paste0("C", cell_number, "B"))) {
      blank_data <- get(paste0("C", cell_number, "B"))
      cd_combined[, 2] <- (cd_combined[, 2] - 
                            blank_data[[paste0("CD_B", cell_number)]]) / 
        (get(paste0("C", cell_number, "P")) * get(paste0("C", 
        cell_number, "C")) * get(paste0("C", cell_number, "R"))) * 0.1
    }
    
    # Filter CD Data Based on High Tension Voltage Greater Than 600 V:
    cd_combined <- cd_combined[cd_combined$HT < 600 | is.na(cd_combined$HT), ]
    
    # Save Data Frame as CSV:
    write.csv(cd_combined[, ], 
              file = paste0(output_directory, "/Cell_", cell_number,
                            "_Spectrum.csv"), 
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
      legend.box.background = element_rect(color = "black", linewidth = 1), 
      # Introduce Legend Panel
      legend.key.size = unit(0.8, "lines"),  # Adjust Legend Size
      legend.text = element_text(size = 10) # Adjust Legend Text Size
    )
    
    # Plotting:
    p <- ggplot(cd_combined) +
      geom_hline(yintercept = 0, linetype = "solid", color = "black", 
                 size = 1) +  # Add Horizontal Line at y = 0
      geom_point(aes(x = get(paste0("Wavelength_", cell_number)), 
                    y = get(paste0("CD_", cell_number))), shape = 21, 
                    color = "#00539F", fill = "#00539F", size = 0.5) +
      geom_line(aes(x = get(paste0("Wavelength_", cell_number)),
                    y = get(paste0("CD_", cell_number))), color = "#00539F",
                    linetype = "solid", size = 0.5)
    
    # Label:
    p <- p +
      scale_starshape_manual(labels = c(C1N))
    
    # Plot Label and Limits:
    p <- p +
      labs(
        x = "Wavelength (nm)",
        y = expression(bold(paste("MRE", " \u00d7 10"^-3, " (deg cm"^2,
                                  "dmol"^-1, ")")))
      ) +
      scale_x_continuous(expand = c(0.0045, 0),
                         breaks = function(x) seq(cd_spectra_x_lower, 
                                                  cd_spectra_x_upper, by = 5),
                         minor_breaks = NULL,
                         labels = function(x) ifelse(x %% 10 == 0,
                                                     as.character(x), ""), 
                         # Keep Labels for Major Ticks Only
                         limits = c(cd_spectra_x_lower, cd_spectra_x_upper)) +
      # Set X-Axis Limits
      scale_y_continuous(expand = c(0.0045, 0),
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
    
    ########The line above may throw an error.
    
    # Define File Name:
    filename <- paste0(output_directory, "/Cell_", cell_number,
                       "_Spectrum.svg")
    
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
  process_cell(cell_number = cell_number)
}
