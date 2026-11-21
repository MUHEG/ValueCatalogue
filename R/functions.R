# Functions go here prior to export to libraries
get_variant_segment <- function(text_chr,
                                category_chr,
                                variants_ls,
                                delimiter_chr = ":",
                                end_patterns_chr = NULL) {
  
  if (is.na(text_chr) ||
      is.na(category_chr) ||
      !category_chr %in% names(variants_ls)) {
    
    return(
      list(
        segment_chr = text_chr,
        remainder_chr = NA_character_
      )
    )
    
  }
  
  start_patterns_chr <- variants_ls[[category_chr]]
  
  start_positions_tb <- purrr::map_dfr(
    start_patterns_chr,
    function(pattern_chr) {
      
      loc_int <- stringr::str_locate(
        string = stringr::str_to_lower(text_chr),
        pattern = stringr::fixed(
          stringr::str_to_lower(pattern_chr)
        )
      )[1]
      
      tibble::tibble(
        pattern_chr = pattern_chr,
        start_int = loc_int
      )
      
    }
  ) %>%
    dplyr::filter(
      !is.na(start_int)
    )
  
  if (nrow(start_positions_tb) == 0) {
    
    return(
      list(
        segment_chr = text_chr,
        remainder_chr = NA_character_
      )
    )
    
  }
  
  start_pattern_chr <- start_positions_tb %>%
    dplyr::arrange(start_int) %>%
    dplyr::slice(1) %>%
    dplyr::pull(pattern_chr)
  
  start_loc_mtx <- stringr::str_locate(
    string = stringr::str_to_lower(text_chr),
    pattern = stringr::fixed(
      stringr::str_to_lower(start_pattern_chr)
    )
  )
  
  segment_start_int <- start_loc_mtx[2]
  
  remainder_chr <- substr(
    text_chr,
    segment_start_int + 1,
    nchar(text_chr)
  )
  
  delimiter_loc_mtx <- stringr::str_locate(
    remainder_chr,
    stringr::fixed(delimiter_chr)
  )
  
  if (!is.na(delimiter_loc_mtx[1])) {
    
    remainder_chr <- substr(
      remainder_chr,
      delimiter_loc_mtx[2] + 1,
      nchar(remainder_chr)
    )
    
  }
  
  terminal_patterns_chr <- unlist(
    variants_ls[
      names(variants_ls) != category_chr
    ],
    use.names = FALSE
  )
  
  if (!is.null(end_patterns_chr)) {
    
    terminal_patterns_chr <- c(
      terminal_patterns_chr,
      end_patterns_chr
    )
    
  }
  
  terminal_positions_tb <- purrr::map_dfr(
    terminal_patterns_chr,
    function(pattern_chr) {
      
      loc_int <- stringr::str_locate(
        string = stringr::str_to_lower(remainder_chr),
        pattern = stringr::fixed(
          stringr::str_to_lower(pattern_chr)
        )
      )[1]
      
      tibble::tibble(
        pattern_chr = pattern_chr,
        start_int = loc_int
      )
      
    }
  ) %>%
    dplyr::filter(
      !is.na(start_int)
    )
  
  if (nrow(terminal_positions_tb) == 0) {
    
    return(
      list(
        segment_chr = stringr::str_squish(
          remainder_chr
        ),
        remainder_chr = NA_character_
      )
    )
    
  }
  
  first_terminal_int <- terminal_positions_tb %>%
    dplyr::arrange(start_int) %>%
    dplyr::slice(1)
  
  segment_chr <- stringr::str_squish(
    substr(
      remainder_chr,
      1,
      first_terminal_int$start_int - 1
    )
  )
  
  terminal_pattern_chr <- first_terminal_int$pattern_chr
  
  terminal_loc_mtx <- stringr::str_locate(
    string = remainder_chr,
    pattern = stringr::fixed(
      terminal_pattern_chr
    )
  )
  
  trailing_chr <- stringr::str_squish(
    substr(
      remainder_chr,
      terminal_loc_mtx[2] + 1,
      nchar(remainder_chr)
    )
  )
  
  list(
    segment_chr = segment_chr,
    remainder_chr = dplyr::if_else(
      trailing_chr == "",
      NA_character_,
      trailing_chr
    )
  )
  
}
make_split_cases <- function(x) {
  
  if (length(x) == 0L ||
      is.na(x) ||
      !nzchar(trimws(x))) {
    return(
      tibble::tibble(
        value_chr = character(0),
        split_pattern_chr = character(0)
      )
    )
  }
  
  marker_pattern <- "\\([0-9]+\\)"
  
  markers <- gregexpr(
    marker_pattern,
    x,
    perl = TRUE
  )[[1]]
  
  # No numbered items
  if (markers[1] == -1) {
    
    return(
      tibble::tibble(
        value_chr = x,
        split_pattern_chr = NA_character_
      )
    )
    
  }
  
  marker_lengths <- attr(
    markers,
    "match.length"
  )
  
  split_patterns_chr <- regmatches(
    x,
    gregexpr(
      marker_pattern,
      x,
      perl = TRUE
    )
  )[[1]]
  
  prefix <- trimws(
    substr(
      x,
      1,
      markers[1] - 1
    )
  )
  
  items_chr <- vapply(
    seq_along(markers),
    function(i) {
      
      start <- markers[i] + marker_lengths[i]
      
      end <- if (i < length(markers)) {
        markers[i + 1] - 1
      } else {
        nchar(x)
      }
      
      trimws(
        substr(
          x,
          start,
          end
        )
      )
      
    },
    FUN.VALUE = character(1)
  )
  
  tibble::tibble(
    value_chr = paste(
      prefix,
      items_chr
    ),
    split_pattern_chr = split_patterns_chr
  )
  
}
transform_to_expand_multiple_matches <- function(
    .data,
    source_col,
    multiple_matches_col,
    separator_chr = ";") {
  
  source_col <- rlang::ensym(source_col)
  multiple_matches_col <- rlang::ensym(
    multiple_matches_col
  )
  
  source_col_nm_chr <- rlang::as_name(
    source_col
  )
  
  multiple_matches_col_nm_chr <- rlang::as_name(
    multiple_matches_col
  )
  
  .data %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      
      .expanded_values_chr = list(
        
        if (is.na(
          .data[[multiple_matches_col_nm_chr]]
        )) {
          
          as.character(
            .data[[source_col_nm_chr]]
          )
          
        } else {
          
          trimws(
            unlist(
              strsplit(
                .data[[multiple_matches_col_nm_chr]],
                separator_chr,
                fixed = TRUE
              )
            )
          )
          
        }
        
      )
      
    ) %>%
    dplyr::ungroup() %>%
    tidyr::unnest_longer(
      .expanded_values_chr
    ) %>%
    dplyr::mutate(
      !!source_col := .expanded_values_chr
    ) %>%
    dplyr::select(
      -.expanded_values_chr
    )
  
}
transform_to_extract_split_segment <- function(
    .data,
    split_pattern_col,
    target_col,
    end_patterns_chr = NULL,
    remainder_col = NULL) {
  
  split_pattern_col <- rlang::ensym(
    split_pattern_col
  )
  
  target_col <- rlang::ensym(
    target_col
  )
  
  if (!is.null(remainder_col)) {
    remainder_col <- rlang::ensym(
      remainder_col
    )
  }
  
  target_col_nm_chr <- rlang::as_name(
    target_col
  )
  
  outp_tb <- .data %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      
      .segment_chr = {
        
        text_chr <- as.character(
          !!target_col
        )
        
        split_pattern_chr <- as.character(
          !!split_pattern_col
        )
        
        if (is.na(split_pattern_chr)) {
          
          text_chr
          
        } else {
          
          split_num_int <- as.integer(
            gsub(
              "[()]",
              "",
              split_pattern_chr
            )
          )
          
          next_pattern_chr <- paste0(
            "\\(",
            split_num_int + 1,
            "\\)"
          )
          
          start_loc_int <- stringr::str_locate(
            text_chr,
            stringr::fixed(
              split_pattern_chr
            )
          )[2]
          
          if (is.na(start_loc_int)) {
            
            text_chr
            
          } else {
            
            remaining_chr <- trimws(
              substr(
                text_chr,
                start_loc_int + 1,
                nchar(text_chr)
              )
            )
            
            next_loc_mtx <- stringr::str_locate(
              remaining_chr,
              next_pattern_chr
            )
            
            if (!is.na(next_loc_mtx[1])) {
              
              trimws(
                substr(
                  remaining_chr,
                  1,
                  next_loc_mtx[1] - 1
                )
              )
              
            } else {
              
              remaining_chr
              
            }
            
          }
          
        }
        
      },
      
      .remainder_chr = {
        
        if (is.null(end_patterns_chr)) {
          
          NA_character_
          
        } else {
          
          terminal_pattern_chr <- end_patterns_chr[
            purrr::detect_index(
              end_patterns_chr,
              ~ stringr::str_detect(
                .segment_chr,
                stringr::fixed(.x)
              )
            )
          ]
          
          if (length(terminal_pattern_chr) == 0 ||
              identical(
                terminal_pattern_chr,
                character(0)
              )) {
            
            NA_character_
            
          } else {
            
            split_mtx <- stringr::str_split_fixed(
              .segment_chr,
              stringr::fixed(
                terminal_pattern_chr
              ),
              n = 2
            )
            
            trimws(
              split_mtx[1, 2]
            )
            
          }
          
        }
        
      },
      
      .segment_chr = {
        
        if (is.null(end_patterns_chr)) {
          
          .segment_chr
          
        } else {
          
          terminal_pattern_chr <- end_patterns_chr[
            purrr::detect_index(
              end_patterns_chr,
              ~ stringr::str_detect(
                .segment_chr,
                stringr::fixed(.x)
              )
            )
          ]
          
          if (length(terminal_pattern_chr) == 0 ||
              identical(
                terminal_pattern_chr,
                character(0)
              )) {
            
            .segment_chr
            
          } else {
            
            trimws(
              stringr::str_split_fixed(
                .segment_chr,
                stringr::fixed(
                  terminal_pattern_chr
                ),
                n = 2
              )[1, 1]
            )
            
          }
          
        }
        
      }
      
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      !!target_col := .segment_chr
    ) %>%
    dplyr::select(
      -.segment_chr
    )
  
  if (!is.null(remainder_col)) {
    
    outp_tb <- outp_tb %>%
      dplyr::mutate(
        !!remainder_col := .remainder_chr
      )
    
  }
  
  outp_tb %>%
    dplyr::select(
      -dplyr::any_of(".remainder_chr")
    )
  
}
transform_to_extract_variant_segment <- function(
    .data,
    source_col,
    text_col,
    multiple_matches_col,
    variants_ls,
    delimiter_chr = ":",
    end_patterns_chr = NULL,
    remainder_col = NULL) {
  
  source_col <- rlang::ensym(source_col)
  text_col <- rlang::ensym(text_col)
  multiple_matches_col <- rlang::ensym(
    multiple_matches_col
  )
  
  if (!is.null(remainder_col)) {
    
    remainder_col <- rlang::ensym(
      remainder_col
    )
    
  }
  
  outp_tb <- .data %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      
      .result_ls = list(
        
        if (is.na(
          !!multiple_matches_col
        )) {
          
          list(
            segment_chr = as.character(
              !!text_col
            ),
            remainder_chr = NA_character_
          )
          
        } else {
          
          get_variant_segment(
            text_chr = as.character(
              !!text_col
            ),
            category_chr = as.character(
              !!source_col
            ),
            variants_ls = variants_ls,
            delimiter_chr = delimiter_chr,
            end_patterns_chr = end_patterns_chr
          )
          
        }
        
      )
      
    ) %>%
    dplyr::ungroup()
  
  outp_tb <- outp_tb %>%
    dplyr::mutate(
      
      !!text_col := purrr::map_chr(
        .result_ls,
        "segment_chr"
      )
      
    )
  
  if (!is.null(remainder_col)) {
    
    outp_tb <- outp_tb %>%
      dplyr::mutate(
        
        !!remainder_col := purrr::map_chr(
          .result_ls,
          "remainder_chr"
        )
        
      )
    
  }
  
  outp_tb %>%
    dplyr::select(
      -.result_ls
    )
  
}
transform_to_filtered_patterns <- function(
    .data,
    col,
    allowable_patterns_chr = NULL,
    banned_patterns_chr = NULL) {
  
  col <- rlang::ensym(col)
  
  .data %>%
    dplyr::filter(
      
      purrr::map_lgl(
        
        as.character(!!col),
        
        function(x) {
          
          if (is.na(x))
            return(FALSE)
          
          x_lower_chr <- stringr::str_to_lower(x)
          
          allowable_lgl <- if (is.null(
            allowable_patterns_chr
          )) {
            
            TRUE
            
          } else {
            
            any(
              stringr::str_detect(
                string = x_lower_chr,
                pattern = stringr::fixed(
                  stringr::str_to_lower(
                    allowable_patterns_chr
                  )
                )
              )
            )
            
          }
          
          banned_lgl <- if (is.null(
            banned_patterns_chr
          )) {
            
            FALSE
            
          } else {
            
            any(
              stringr::str_detect(
                string = x_lower_chr,
                pattern = stringr::fixed(
                  stringr::str_to_lower(
                    banned_patterns_chr
                  )
                )
              )
            )
            
          }
          
          allowable_lgl & !banned_lgl
          
        }
        
      )
      
    )
  
}
transform_to_permissible_values <- function(.data,
                                            source_col,
                                            permissible_values_chr,
                                            original_value_col,
                                            multiple_matches_col,
                                            residual_value_col = NULL,
                                            unmatched_action_chr = c(
                                              "error",
                                              "leave",
                                              "set_na"
                                            )) {
  
  unmatched_action_chr <- match.arg(
    unmatched_action_chr
  )
  
  source_col <- rlang::ensym(source_col)
  original_value_col <- rlang::ensym(original_value_col)
  multiple_matches_col <- rlang::ensym(multiple_matches_col)
  
  if (!is.null(residual_value_col)) {
    residual_value_col <- rlang::ensym(
      residual_value_col
    )
  }
  
  source_col_nm_chr <- rlang::as_name(
    source_col
  )
  
  outp_tb <- .data %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      
      .matches_chr = list({
        
        value_chr <- as.character(
          !!source_col
        )
        
        exact_idx <- match(
          tolower(value_chr),
          tolower(permissible_values_chr)
        )
        
        if (!is.na(exact_idx)) {
          
          permissible_values_chr[
            exact_idx
          ]
          
        } else {
          
          permissible_values_chr[
            stringr::str_detect(
              string = stringr::str_to_lower(
                value_chr
              ),
              pattern = stringr::fixed(
                stringr::str_to_lower(
                  permissible_values_chr
                )
              )
            )
          ]
          
        }
        
      }),
      
      .n_matches_int = length(
        .matches_chr
      ),
      
      .replacement_chr = {
        
        if (.n_matches_int == 1) {
          .matches_chr[[1]]
        } else {
          NA_character_
        }
        
      },
      
      !!original_value_col := dplyr::case_when(
        !is.na(.replacement_chr) &
          .replacement_chr != as.character(
            !!source_col
          ) ~ as.character(
            !!source_col
          ),
        TRUE ~ NA_character_
      ),
      
      !!multiple_matches_col := dplyr::case_when(
        .n_matches_int > 1 ~ paste(
          sort(.matches_chr),
          collapse = "; "
        ),
        TRUE ~ NA_character_
      ),
      
      !!source_col := dplyr::coalesce(
        .replacement_chr,
        as.character(
          !!source_col
        )
      )
      
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(
      -.matches_chr,
      -.n_matches_int,
      -.replacement_chr
    )
  
  original_col_nm_chr <- rlang::as_name(
    original_value_col
  )
  
  multiple_matches_col_nm_chr <- rlang::as_name(
    multiple_matches_col
  )
  
  unmatched_lgl <- is.na(
    outp_tb[[original_col_nm_chr]]
  ) &
    is.na(
      outp_tb[[multiple_matches_col_nm_chr]]
    ) &
    !(
      tolower(
        outp_tb[[source_col_nm_chr]]
      ) %in%
        tolower(
          permissible_values_chr
        )
    )
  
  if (unmatched_action_chr == "error") {
    
    unmatched_chr <- unique(
      outp_tb[[source_col_nm_chr]][
        unmatched_lgl
      ]
    )
    
    if (length(unmatched_chr) > 0) {
      
      stop(
        paste0(
          "No permissible value found for: ",
          paste(
            unmatched_chr,
            collapse = "; "
          )
        ),
        call. = FALSE
      )
      
    }
    
  }
  
  if (unmatched_action_chr == "set_na") {
    
    outp_tb <- outp_tb %>%
      dplyr::mutate(
        
        !!original_value_col := dplyr::if_else(
          unmatched_lgl &
            is.na(.data[[original_col_nm_chr]]),
          as.character(.data[[source_col_nm_chr]]),
          .data[[original_col_nm_chr]]
        ),
        
        !!source_col := dplyr::if_else(
          unmatched_lgl,
          NA_character_,
          as.character(.data[[source_col_nm_chr]])
        )
        
      )
    
  }
  
  if (!is.null(
    residual_value_col
  )) {
    
    outp_tb <- outp_tb %>%
      dplyr::mutate(
        
        !!residual_value_col := dplyr::if_else(
          
          !is.na(
            .data[[original_col_nm_chr]]
          ),
          
          {
            residual_chr <- stringr::str_remove_all(
              string = .data[[
                original_col_nm_chr
              ]],
              pattern = stringr::regex(
                .data[[
                  source_col_nm_chr
                ]],
                ignore_case = TRUE
              )
            )
            
            residual_chr <- stringr::str_squish(
              residual_chr
            )
            
            dplyr::if_else(
              residual_chr == "",
              NA_character_,
              residual_chr
            )
            
          },
          
          NA_character_
          
        )
        
      )
    
  }
  
  outp_tb
  
}
transform_to_split_cases <- function(.data,
                                     col,
                                     split_pattern_col = NULL) {
  
  col <- rlang::ensym(col)
  
  if (!is.null(split_pattern_col)) {
    split_pattern_col <- rlang::ensym(
      split_pattern_col
    )
  }
  
  outp_tb <- .data %>%
    dplyr::mutate(
      .split = purrr::map(
        !!col,
        make_split_cases
      )
    ) %>%
    tidyr::unnest(
      .split
    ) %>%
    dplyr::mutate(
      !!col := value_chr
    ) %>%
    dplyr::select(
      -value_chr
    )
  
  if (!is.null(split_pattern_col)) {
    
    outp_tb <- outp_tb %>%
      dplyr::rename(
        !!split_pattern_col :=
          split_pattern_chr
      )
    
  } else {
    
    outp_tb <- outp_tb %>%
      dplyr::select(
        -split_pattern_chr
      )
    
  }
  
  outp_tb
  
}
transform_to_split_cols <- function(.data,
                                    source_col,
                                    before_col,
                                    after_col,
                                    split_pattern_chr = ":") {
  
  source_col <- rlang::ensym(source_col)
  before_col <- rlang::ensym(before_col)
  after_col <- rlang::ensym(after_col)
  
  .data %>%
    dplyr::mutate(
      .split = stringr::str_split_fixed(
        string = as.character(!!source_col),
        pattern = stringr::fixed(split_pattern_chr),
        n = 2
      ),
      !!before_col := trimws(
        dplyr::if_else(
          .split[, 2] == "",
          as.character(!!source_col),
          .split[, 1]
        )
      ),
      !!after_col := trimws(
        dplyr::if_else(
          .split[, 2] == "",
          as.character(!!source_col),
          .split[, 2]
        )
      )
    ) %>%
    dplyr::select(-.split)
  
}
