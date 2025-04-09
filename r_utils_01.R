Setup <- R6Class("Setup",
  public = list(
    
    # Stores the final merged configuration settings
    final_config = NULL,
    
    # Class constructor
    initialize = function(config_file) {
      # Step 1: Load default configuration file
      default_config <- tryCatch({
        yaml::read_yaml("default_settings.yml")  # Load default settings from YAML
      }, error = function(e) {
        message("Error loading default config: ", e)
        list()  # Return empty list if there's an error
      })

      # Step 2: Load custom configuration file
      custom_config <- tryCatch({
        yaml::read_yaml(config_file)  # Load custom settings from specified YAML file
      }, error = function(e) {
        message("Error loading custom config: ", e)
        list()  # Return empty list if there's an error
      })

      # Step 3: Merge default and custom configurations (custom will overwrite default)
      final_config <- modifyList(default_config, custom_config)

      # Step 4: Add date to the configuration (current date and time)
      final_config$date <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")

      # Step 5: Add ID if not already present (generates a new UUID)
      if (is.null(final_config$id) || final_config$id == "") {
        final_config$id <- toupper(substr(UUIDgenerate(), 1, 10))  # Generate a new ID
      }

      # Step 6: Print initial status about the analysis with ID and creation date
      message(sprintf("New analysis with id: %s; created: %s", 
                      final_config$id, final_config$date))

      # Step 7: Check that required settings are provided in the config
      required_settings <- c(
        "output_column", "comparison", "train_ancestry", "infer_ancestry", 
        "ancestry_column", "seed", "data_path", "data_type", "tech", "output_directory"
      )
      
      for (setting in required_settings) {
        stopifnot(!is.null(final_config[[setting]]))  # Ensure required settings exist
      }

      # Step 8: Ensure covariate, if present, is a list
      if (!is.null(final_config$covariate)) {
        stopifnot(is.list(final_config$covariate))
      }

      # Step 9: Handle classification setup based on comparison settings
      if (length(final_config$comparison) > 2) {
        final_config$multiclass <- TRUE
        message("Multiclass problem comparing: ", paste(final_config$comparison, collapse=", "))
      } else if (length(final_config$comparison) == 2) {
        final_config$multiclass <- FALSE
        message("Binaryclass problem comparing: ", paste(final_config$comparison, collapse=", "))
      } else {
        stop("Can't classify with only one class in comparison.")
      }

      # Step 10: Save final configuration and check settings
      self$final_config <- final_config
      self$check_settings()
    },
    
    # Step 11: Check validity of the configuration settings
    check_settings = function() {
      stopifnot(is.numeric(self$final_config$seed))  # Ensure seed is a number
      stopifnot(file.exists(self$final_config$data_path))  # Ensure data path exists
      stopifnot(grepl("\\.h5ad$", self$final_config$data_path))  # Ensure file has .h5ad extension
      stopifnot(self$final_config$tech %in% c("methylation", "transcriptomics", "proteomics"))  # Ensure tech is valid

      # Classification checks
      stopifnot(is.character(self$final_config$ancestry_column))
      stopifnot(is.character(self$final_config$train_ancestry))
      stopifnot(is.character(self$final_config$infer_ancestry))

      # Machine learning settings validation
      stopifnot(is.numeric(self$final_config$nfolds))
      stopifnot(self$final_config$nfolds >= 2)  # At least 2 folds for cross-validation
    },

    # Step 12: Generate output file paths based on the output directory
    out = function(x) {
      file.path(self$final_config$output_directory, x)  # Generate full file path
    },

    # Step 13: Create the output directory (if it doesn't exist) and initialize log
    make_output_directory = function() {
      dir.create(self$final_config$output_directory, recursive = TRUE, showWarnings = FALSE)  # Create the directory
      message("Output will be saved to ", normalizePath(self$final_config$output_directory))  # Print the save location
      
      # Create a log file with headers
      writeLines("Step\tMemory_MB\tTime", self$out("Log.tsv"))
    },

    # Step 14: Log process steps (with memory usage and timestamp)
    log = function(text) {
      log_line <- paste0(
        text, "\t",
        round(memory.size(), 2), "\t",  # Memory usage in MB
        format(Sys.time(), "%Y-%m-%d %H:%M:%S")  # Timestamp
      )
      write(log_line, file = self$out("Log.tsv"), append = TRUE)  # Append log entry
    },

    # Step 15: Return the current configuration settings
    return_settings = function() {
      return(self$final_config)  # Return the final config
    },

    # Step 16: Add new settings to the configuration (with optional overwrite)
    add_settings = function(new_settings, overwrite = FALSE) {
      for (key in names(new_settings)) {
        if (!overwrite && key %in% names(self$final_config)) {
          message(sprintf("Setting '%s' already exists. Use overwrite=TRUE to modify it.", key))
        } else {
          self$final_config[[key]] <- new_settings[[key]]  # Update the setting
        }
      }
    }
  )
)

