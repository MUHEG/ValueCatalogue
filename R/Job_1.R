source("Private/paths.R")
rmarkdown::render("RMD/Parent_RMDs/Initial_Report.Rmd",
                  params = list(
                                file_names_ls = file_names_ls, # Defined in "Private/paths.R"
                                path_to_input_1L_chr = path_to_input_1L_chr, # Defined in "Private/paths.R"
                                path_to_output_1L_chr = path_to_output_1L_chr, # Defined in "Private/paths.R"
                                sheet_names_ls = sheet_names_ls # Defined in "Private/paths.R"
                                ),
                  output_dir = path_to_output_1L_chr)