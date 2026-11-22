# Functions go here prior to export to libraries
get_variant_segment <- function(
    text_chr,
    category_chr,
    variants_ls,
    delimiter_chr = ":",
    end_patterns_chr = NULL,
    continuation_patterns_chr = NULL,
    fallback_to_end_patterns_lgl = TRUE) {
  
  make_return_ls <- function(segment_chr,
                             remainder_chr = NA_character_,
                             matched_start_pattern_chr = NA_character_,
                             matched_terminal_pattern_chr = NA_character_,
                             used_fallback_lgl = FALSE) {
    
    list(
      segment_chr = segment_chr,
      remainder_chr = remainder_chr,
      matched_start_pattern_chr = matched_start_pattern_chr,
      matched_terminal_pattern_chr = matched_terminal_pattern_chr,
      used_fallback_lgl = used_fallback_lgl
    )
    
  }
  
  if (is.na(text_chr) ||
      is.na(category_chr) ||
      !category_chr %in% names(variants_ls)) {
    
    return(
      make_return_ls(
        segment_chr = text_chr
      )
    )
    
  }
  
  start_patterns_chr <- paste0(
    variants_ls[[category_chr]],
    delimiter_chr
  )
  
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
    
    if (isTRUE(
      fallback_to_end_patterns_lgl
    ) &&
    !is.null(
      end_patterns_chr
    )) {
      
      end_positions_tb <- purrr::map_dfr(
        
        end_patterns_chr,
        
        function(pattern_chr) {
          
          loc_int <- stringr::str_locate(
            string = stringr::str_to_lower(
              text_chr
            ),
            pattern = stringr::fixed(
              stringr::str_to_lower(
                pattern_chr
              )
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
      
      if (nrow(end_positions_tb) > 0) {
        
        first_end_tb <- end_positions_tb %>%
          dplyr::arrange(start_int) %>%
          dplyr::slice(1)
        
        segment_chr <- stringr::str_squish(
          substr(
            text_chr,
            1,
            first_end_tb$start_int - 1
          )
        )
        
        remainder_chr <- stringr::str_squish(
          substr(
            text_chr,
            first_end_tb$start_int,
            nchar(text_chr)
          )
        )
        
        return(
          
          make_return_ls(
            segment_chr = segment_chr,
            remainder_chr = dplyr::if_else(
              remainder_chr == "",
              NA_character_,
              remainder_chr
            ),
            matched_terminal_pattern_chr =
              first_end_tb$pattern_chr,
            used_fallback_lgl = TRUE
          )
          
        )
        
      }
      
    }
    
    return(
      make_return_ls(
        segment_chr = text_chr
      )
    )
    
  }
  
  first_start_tb <- start_positions_tb %>%
    dplyr::arrange(start_int) %>%
    dplyr::slice(1)
  
  working_chr <- substr(
    text_chr,
    first_start_tb$start_int +
      nchar(first_start_tb$pattern_chr),
    nchar(text_chr)
  )
  
  terminal_patterns_chr <- paste0(
    unlist(
      variants_ls[
        names(variants_ls) != category_chr
      ],
      use.names = FALSE
    ),
    delimiter_chr
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
      
      locs_mtx <- stringr::str_locate_all(
        string = stringr::str_to_lower(
          working_chr
        ),
        pattern = stringr::fixed(
          stringr::str_to_lower(
            pattern_chr
          )
        )
      )[[1]]
      
      if (nrow(locs_mtx) == 0) {
        
        return(
          tibble::tibble(
            pattern_chr = character(),
            start_int = integer()
          )
        )
        
      }
      
      valid_lgl <- rep(
        TRUE,
        nrow(locs_mtx)
      )
      
      if (!is.null(
        continuation_patterns_chr
      )) {
        
        for (i in seq_len(
          nrow(locs_mtx)
        )) {
          
          following_chr <- substr(
            working_chr,
            locs_mtx[i, 2] + 1,
            nchar(working_chr)
          )
          
          valid_lgl[i] <- !stringr::str_detect(
            following_chr,
            paste0(
              "^\\s*(",
              paste(
                continuation_patterns_chr,
                collapse = "|"
              ),
              ")"
            )
          )
          
        }
        
      }
      
      tibble::tibble(
        pattern_chr = pattern_chr,
        start_int = locs_mtx[
          valid_lgl,
          1
        ]
      )
      
    }
    
  ) %>%
    dplyr::filter(
      !is.na(start_int)
    )
  
  if (nrow(
    terminal_positions_tb
  ) == 0) {
    
    return(
      
      make_return_ls(
        segment_chr = stringr::str_squish(
          working_chr
        ),
        matched_start_pattern_chr =
          first_start_tb$pattern_chr
      )
      
    )
    
  }
  
  first_terminal_tb <- terminal_positions_tb %>%
    dplyr::arrange(start_int) %>%
    dplyr::slice(1)
  
  segment_chr <- stringr::str_squish(
    substr(
      working_chr,
      1,
      first_terminal_tb$start_int - 1
    )
  )
  
  remainder_chr <- stringr::str_squish(
    substr(
      working_chr,
      first_terminal_tb$start_int,
      nchar(working_chr)
    )
  )
  
  make_return_ls(
    segment_chr = segment_chr,
    remainder_chr = dplyr::if_else(
      remainder_chr == "",
      NA_character_,
      remainder_chr
    ),
    matched_start_pattern_chr =
      first_start_tb$pattern_chr,
    matched_terminal_pattern_chr =
      first_terminal_tb$pattern_chr
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
transform_to_clean_text <- function(
    .data,
    source_col,
    original_value_col = NULL,
    keep_after_patterns_chr = NULL,
    except_if_contains_patterns_chr = NULL,
    remove_trailing_patterns_chr = NULL,
    remove_keep_after_pattern_lgl = TRUE) {
  
  source_col <- rlang::ensym(
    source_col
  )
  
  if (!is.null(
    original_value_col
  )) {
    
    original_value_col <- rlang::ensym(
      original_value_col
    )
    
  }
  
  outp_tb <- .data %>%
    
    dplyr::mutate(
      
      .cleaned_chr =
        as.character(
          !!source_col
        )
      
    )
  
  if (!is.null(
    original_value_col
  )) {
    
    outp_tb <- outp_tb %>%
      
      dplyr::mutate(
        
        !!original_value_col :=
          .cleaned_chr
        
      )
    
  }
  
  if (!is.null(
    keep_after_patterns_chr
  )) {
    
    outp_tb <- outp_tb %>%
      
      dplyr::rowwise() %>%
      
      dplyr::mutate(
        
        .cleaned_chr = {
          
          x_chr <- .cleaned_chr
          
          exception_found_lgl <- FALSE
          
          if (!is.null(
            except_if_contains_patterns_chr
          )) {
            
            exception_found_lgl <- any(
              
              purrr::map_lgl(
                
                except_if_contains_patterns_chr,
                
                function(pattern_chr) {
                  
                  stringr::str_detect(
                    
                    string =
                      stringr::str_to_lower(
                        x_chr
                      ),
                    
                    pattern =
                      stringr::fixed(
                        stringr::str_to_lower(
                          pattern_chr
                        )
                      )
                    
                  )
                  
                }
                
              )
              
            )
            
          }
          
          if (!exception_found_lgl) {
            
            positions_tb <- purrr::map_dfr(
              
              keep_after_patterns_chr,
              
              function(pattern_chr) {
                
                loc_mtx <- stringr::str_locate(
                  
                  string =
                    stringr::str_to_lower(
                      x_chr
                    ),
                  
                  pattern =
                    stringr::fixed(
                      stringr::str_to_lower(
                        pattern_chr
                      )
                    )
                  
                )
                
                tibble::tibble(
                  
                  pattern_chr =
                    pattern_chr,
                  
                  start_int =
                    loc_mtx[1],
                  
                  end_int =
                    loc_mtx[2]
                  
                )
                
              }
              
            ) %>%
              
              dplyr::filter(
                !is.na(start_int)
              )
            
            if (nrow(
              positions_tb
            ) > 0) {
              
              first_match_tb <-
                positions_tb %>%
                dplyr::arrange(
                  start_int
                ) %>%
                dplyr::slice(1)
              
              if (isTRUE(
                remove_keep_after_pattern_lgl
              )) {
                
                x_chr <- substr(
                  
                  x_chr,
                  
                  first_match_tb$end_int + 1,
                  
                  nchar(x_chr)
                  
                )
                
              } else {
                
                x_chr <- substr(
                  
                  x_chr,
                  
                  first_match_tb$start_int,
                  
                  nchar(x_chr)
                  
                )
                
              }
              
            }
            
          }
          
          stringr::str_squish(
            x_chr
          )
          
        }
        
      ) %>%
      
      dplyr::ungroup()
    
  }
  
  if (!is.null(
    remove_trailing_patterns_chr
  )) {
    
    outp_tb <- outp_tb %>%
      
      dplyr::mutate(
        
        .cleaned_chr = purrr::map_chr(
          
          .cleaned_chr,
          
          function(x_chr) {
            
            if (is.na(x_chr))
              return(
                NA_character_
              )
            
            for (pattern_chr in
                 remove_trailing_patterns_chr) {
              
              if (stringr::str_ends(
                
                string = x_chr,
                
                pattern =
                stringr::fixed(
                  pattern_chr
                )
                
              )) {
                
                x_chr <- substr(
                  
                  x_chr,
                  
                  1,
                  
                  nchar(x_chr) -
                    nchar(
                      pattern_chr
                    )
                  
                )
                
              }
              
            }
            
            stringr::str_squish(
              x_chr
            )
            
          }
          
        )
        
      )
    
  }
  
  outp_tb %>%
    
    dplyr::mutate(
      
      !!source_col :=
        .cleaned_chr
      
    ) %>%
    
    dplyr::select(
      -.cleaned_chr
    )
  
}
# transform_to_expand_ce_perspectives <- function(
#     .data,
#     source_col,
#     variants_ls,
#     perspective_col = "Perspective",
#     perspective_text_col = "CE_Perspective_Text",
#     original_text_col = NULL,
#     matched_pattern_col = NULL) {
#   
#   source_col <- rlang::ensym(
#     source_col
#   )
#   
#   perspective_col <- rlang::sym(
#     perspective_col
#   )
#   
#   perspective_text_col <- rlang::sym(
#     perspective_text_col
#   )
#   
#   if (!is.null(
#     original_text_col
#   )) {
#     
#     original_text_col <- rlang::sym(
#       original_text_col
#     )
#     
#   }
#   
#   if (!is.null(
#     matched_pattern_col
#   )) {
#     
#     matched_pattern_col <- rlang::sym(
#       matched_pattern_col
#     )
#     
#   }
#   
#   get_perspective_blocks <- function(
#     text_chr,
#     variants_ls) {
#     
#     if (is.na(text_chr)) {
#       
#       return(
#         
#         tibble::tibble(
#           
#           perspective_chr =
#             NA_character_,
#           
#           matched_pattern_chr =
#             NA_character_,
#           
#           perspective_text_chr =
#             NA_character_
#           
#         )
#         
#       )
#       
#     }
#     
#     matches_tb <- purrr::imap_dfr(
#       
#       variants_ls,
#       
#       function(
#     patterns_chr,
#     perspective_chr
#       ) {
#         
#         purrr::map_dfr(
#           
#           patterns_chr,
#           
#           function(pattern_chr) {
#             
#             locs_mtx <-
#               stringr::str_locate_all(
#                 
#                 string =
#                   stringr::str_to_lower(
#                     text_chr
#                   ),
#                 
#                 pattern =
#                   stringr::fixed(
#                     stringr::str_to_lower(
#                       pattern_chr
#                     )
#                   )
#                 
#               )[[1]]
#             
#             if (
#               nrow(
#                 locs_mtx
#               ) == 0
#             ) {
#               
#               return(
#                 
#                 tibble::tibble(
#                   
#                   perspective_chr =
#                     character(),
#                   
#                   matched_pattern_chr =
#                     character(),
#                   
#                   start_int =
#                     integer(),
#                   
#                   end_int =
#                     integer()
#                   
#                 )
#                 
#               )
#               
#             }
#             
#             tibble::tibble(
#               
#               perspective_chr =
#                 perspective_chr,
#               
#               matched_pattern_chr =
#                 pattern_chr,
#               
#               start_int =
#                 locs_mtx[, 1],
#               
#               end_int =
#                 locs_mtx[, 2]
#               
#             )
#             
#           }
#           
#         )
#         
#       }
#     
#     )
#     
#     if (
#       nrow(
#         matches_tb
#       ) == 0
#     ) {
#       
#       return(
#         
#         tibble::tibble(
#           
#           perspective_chr =
#             NA_character_,
#           
#           matched_pattern_chr =
#             NA_character_,
#           
#           perspective_text_chr =
#             stringr::str_squish(
#               text_chr
#             )
#           
#         )
#         
#       )
#       
#     }
#     
#     matches_tb <- matches_tb %>%
#       
#       dplyr::mutate(
#         
#         match_length_int =
#           end_int -
#           start_int +
#           1
#         
#       ) %>%
#       
#       dplyr::group_by(
#         start_int
#       ) %>%
#       
#       dplyr::slice_max(
#         
#         order_by =
#           match_length_int,
#         
#         n = 1,
#         
#         with_ties = FALSE
#         
#       ) %>%
#       
#       dplyr::ungroup() %>%
#       
#       dplyr::select(
#         -match_length_int
#       ) %>%
#       
#       dplyr::arrange(
#         start_int
#       )
#     
#     outp_ls <- vector(
#       
#       mode = "list",
#       
#       length =
#         nrow(
#           matches_tb
#         )
#       
#     )
#     
#     for (i in seq_len(
#       nrow(matches_tb)
#     )) {
#       
#       block_start_int <-
#         matches_tb$end_int[
#           i
#         ] + 1
#       
#       if (
#         
#         i <
#         nrow(
#           matches_tb
#         )
#         
#       ) {
#         
#         block_end_int <-
#           
#           matches_tb$start_int[
#             i + 1
#           ] - 1
#         
#       } else {
#         
#         block_end_int <-
#           
#           nchar(
#             text_chr
#           )
#         
#       }
#       
#       outp_ls[[i]] <-
#         
#         tibble::tibble(
#           
#           perspective_chr =
#             
#             matches_tb$
#             perspective_chr[
#               i
#             ],
#           
#           matched_pattern_chr =
#             
#             matches_tb$
#             matched_pattern_chr[
#               i
#             ],
#           
#           perspective_text_chr =
#             
#             stringr::str_squish(
#               
#               substr(
#                 
#                 text_chr,
#                 
#                 block_start_int,
#                 
#                 block_end_int
#                 
#               )
#               
#             )
#           
#         )
#       
#     }
#     
#     dplyr::bind_rows(
#       outp_ls
#     )
#     
#   }
#   
#   outp_tb <- .data %>%
#     
#     dplyr::mutate(
#       
#       .original_text_chr =
#         
#         as.character(
#           !!source_col
#         )
#       
#     ) %>%
#     
#     dplyr::rowwise() %>%
#     
#     dplyr::mutate(
#       
#       .perspectives_ls = list(
#         
#         get_perspective_blocks(
#           
#           text_chr =
#             
#             as.character(
#               !!source_col
#             ),
#           
#           variants_ls =
#             variants_ls
#           
#         )
#         
#       )
#       
#     ) %>%
#     
#     dplyr::ungroup() %>%
#     
#     tidyr::unnest(
#       .perspectives_ls
#     ) %>%
#     
#     dplyr::mutate(
#       
#       !!perspective_col :=
#         perspective_chr,
#       
#       !!perspective_text_col :=
#         perspective_text_chr
#       
#     )
#   
#   if (!is.null(
#     matched_pattern_col
#   )) {
#     
#     outp_tb <- outp_tb %>%
#       
#       dplyr::mutate(
#         
#         !!matched_pattern_col :=
#           matched_pattern_chr
#         
#       )
#     
#   }
#   
#   if (!is.null(
#     original_text_col
#   )) {
#     
#     outp_tb <- outp_tb %>%
#       
#       dplyr::mutate(
#         
#         !!original_text_col :=
#           .original_text_chr
#         
#       )
#     
#   }
#   
#   outp_tb %>%
#     
#     dplyr::select(
#       
#       -tidyselect::any_of(
#         c(
#           "perspective_chr",
#           "matched_pattern_chr",
#           "perspective_text_chr",
#           ".original_text_chr"
#         )
#       )
#       
#     )
#   
# }
transform_to_expand_ce_statistics <- function(
    .data,
    source_col,
    statistic_patterns_ls,
    statistic_col = "Statistic_Type",
    statistic_text_col = "Statistic_Text",
    original_text_col = NULL,
    keep_unmatched_rows_lgl = TRUE) {
  
  source_col <- rlang::ensym(
    source_col
  )
  
  statistic_col <- rlang::sym(
    statistic_col
  )
  
  statistic_text_col <- rlang::sym(
    statistic_text_col
  )
  
  if (!is.null(
    original_text_col
  )) {
    
    original_text_col <- rlang::sym(
      original_text_col
    )
    
  }
  
  extract_statistics <- function(
    text_chr) {
    
    if (
      is.na(text_chr)
    ) {
      
      return(
        
        tibble::tibble(
          
          statistic_type_chr =
            NA_character_,
          
          statistic_text_chr =
            NA_character_
          
        )
        
      )
      
    }
    
    matches_tb <- purrr::imap_dfr(
      
      statistic_patterns_ls,
      
      function(
    patterns_chr,
    statistic_chr
      ) {
        
        purrr::map_dfr(
          
          patterns_chr,
          
          function(pattern_chr) {
            
            locs_mtx <-
              
              stringr::str_locate_all(
                
                string =
                  stringr::str_to_lower(
                    text_chr
                  ),
                
                pattern =
                  stringr::fixed(
                    stringr::str_to_lower(
                      pattern_chr
                    )
                  )
                
              )[[1]]
            
            if (
              nrow(
                locs_mtx
              ) == 0
            ) {
              
              return(
                tibble::tibble()
              )
              
            }
            
            tibble::tibble(
              
              statistic_type_chr =
                statistic_chr,
              
              pattern_chr =
                pattern_chr,
              
              start_int =
                locs_mtx[,1],
              
              end_int =
                locs_mtx[,2]
              
            )
            
          }
          
        )
        
      }
    
    )
    
    if (
      nrow(matches_tb) == 0
    ) {
      
      if (
        keep_unmatched_rows_lgl
      ) {
        
        return(
          
          tibble::tibble(
            
            statistic_type_chr =
              NA_character_,
            
            statistic_text_chr =
              stringr::str_squish(
                text_chr
              )
            
          )
          
        )
        
      }
      
      return(
        tibble::tibble()
      )
      
    }
    
    matches_tb <- matches_tb %>%
      
      dplyr::mutate(
        
        match_width =
          end_int -
          start_int
        
      ) %>%
      
      dplyr::group_by(
        start_int
      ) %>%
      
      dplyr::slice_max(
        match_width,
        n = 1,
        with_ties = FALSE
      ) %>%
      
      dplyr::ungroup() %>%
      
      dplyr::arrange(
        start_int
      )
    
    out_ls <- vector(
      "list",
      nrow(matches_tb)
    )
    
    for (
      i in seq_len(
        nrow(matches_tb)
      )
    ) {
      
      start_pos <-
        matches_tb$end_int[i] + 1
      
      end_pos <-
        
        if (
          i < nrow(matches_tb)
        ) {
          
          matches_tb$start_int[
            i + 1
          ] - 1
          
        } else {
          
          nchar(
            text_chr
          )
          
        }
      
      out_ls[[i]] <-
        
        tibble::tibble(
          
          statistic_type_chr =
            
            matches_tb$
            statistic_type_chr[i],
          
          statistic_text_chr =
            
            stringr::str_squish(
              
              substr(
                
                text_chr,
                
                start_pos,
                
                end_pos
                
              )
              
            )
          
        )
      
    }
    
    dplyr::bind_rows(
      out_ls
    )
    
  }
  
  outp_tb <- .data %>%
    
    dplyr::mutate(
      
      .original_text_chr =
        
        as.character(
          !!source_col
        )
      
    ) %>%
    
    dplyr::rowwise() %>%
    
    dplyr::mutate(
      
      .stats_ls = list(
        
        extract_statistics(
          
          as.character(
            !!source_col
          )
          
        )
        
      )
      
    ) %>%
    
    dplyr::ungroup() %>%
    
    tidyr::unnest(
      .stats_ls
    ) %>%
    
    dplyr::mutate(
      
      !!statistic_col :=
        statistic_type_chr,
      
      !!statistic_text_col :=
        statistic_text_chr
      
    )
  
  if (!is.null(
    original_text_col
  )) {
    
    outp_tb <- outp_tb %>%
      
      dplyr::mutate(
        
        !!original_text_col :=
          
          .original_text_chr
        
      )
    
  }
  
  outp_tb %>%
    
    dplyr::select(
      
      -tidyselect::any_of(
        
        c(
          ".original_text_chr",
          "statistic_type_chr",
          "statistic_text_chr"
        )
        
      )
      
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
transform_to_expand_numbered_items <- function(
    .data,
    source_col,
    item_nbr_col = "Item_Nbr",
    item_seq_col = "Item_Seq",
    item_text_col = "Item_Text",
    pre_item_text_col = "Pre_Item_Text",
    post_item_text_col = "Post_Item_Text",
    original_text_col = NULL,
    keep_unmatched_rows_lgl = TRUE) {
  
  source_col <- rlang::ensym(
    source_col
  )
  
  item_nbr_col <- rlang::sym(
    item_nbr_col
  )
  
  item_seq_col <- rlang::sym(
    item_seq_col
  )
  
  item_text_col <- rlang::sym(
    item_text_col
  )
  
  pre_item_text_col <- rlang::sym(
    pre_item_text_col
  )
  
  post_item_text_col <- rlang::sym(
    post_item_text_col
  )
  
  if (!is.null(
    original_text_col
  )) {
    
    original_text_col <- rlang::sym(
      original_text_col
    )
    
  }
  
  split_numbered_items <- function(
    text_chr) {
    
    if (is.na(
      text_chr
    )) {
      
      return(
        
        tibble::tibble(
          
          item_nbr_chr =
            NA_character_,
          
          item_seq_int =
            NA_integer_,
          
          item_text_chr =
            NA_character_,
          
          pre_item_text_chr =
            NA_character_,
          
          post_item_text_chr =
            NA_character_
          
        )
        
      )
      
    }
    
    locs_mtx <-
      
      stringr::str_locate_all(
        
        text_chr,
        
        "\\([0-9]+\\)"
        
      )[[1]]
    
    if (nrow(
      locs_mtx
    ) == 0) {
      
      if (
        keep_unmatched_rows_lgl
      ) {
        
        return(
          
          tibble::tibble(
            
            item_nbr_chr =
              NA_character_,
            
            item_seq_int =
              1L,
            
            item_text_chr =
              stringr::str_squish(
                text_chr
              ),
            
            pre_item_text_chr =
              stringr::str_squish(
                text_chr
              ),
            
            post_item_text_chr =
              NA_character_
            
          )
          
        )
        
      }
      
      return(
        
        tibble::tibble(
          
          item_nbr_chr =
            character(),
          
          item_seq_int =
            integer(),
          
          item_text_chr =
            character(),
          
          pre_item_text_chr =
            character(),
          
          post_item_text_chr =
            character()
          
        )
        
      )
      
    }
    
    pre_item_text_chr <-
      
      stringr::str_squish(
        
        substr(
          
          text_chr,
          
          1,
          
          locs_mtx[
            1,
            1
          ] - 1
          
        )
        
      )
    
    post_item_text_chr <-
      
      NA_character_
    
    outp_ls <- vector(
      
      mode = "list",
      
      length =
        nrow(
          locs_mtx
        )
      
    )
    
    for (
      
      i in seq_len(
        nrow(
          locs_mtx
        )
      )
      
    ) {
      
      item_nbr_chr <-
        
        substr(
          
          text_chr,
          
          locs_mtx[
            i,
            1
          ],
          
          locs_mtx[
            i,
            2
          ]
          
        )
      
      if (
        
        i <
        nrow(
          locs_mtx
        )
        
      ) {
        
        item_end_int <-
          
          locs_mtx[
            i + 1,
            1
          ] - 1
        
      } else {
        
        item_end_int <-
          
          nchar(
            text_chr
          )
        
      }
      
      item_text_chr <-
        
        substr(
          
          text_chr,
          
          locs_mtx[
            i,
            2
          ] + 1,
          
          item_end_int
          
        )
      
      if (
        
        i ==
        nrow(
          locs_mtx
        )
        
      ) {
        
        post_item_text_chr <-
          
          stringr::str_squish(
            
            item_text_chr
            
          )
        
      }
      
      outp_ls[[i]] <-
        
        tibble::tibble(
          
          item_nbr_chr =
            item_nbr_chr,
          
          item_seq_int =
            i,
          
          item_text_chr =
            stringr::str_squish(
              item_text_chr
            ),
          
          pre_item_text_chr =
            pre_item_text_chr,
          
          post_item_text_chr =
            ifelse(
              i ==
                nrow(
                  locs_mtx
                ),
              post_item_text_chr,
              NA_character_
            )
          
        )
      
    }
    
    dplyr::bind_rows(
      outp_ls
    )
    
  }
  
  outp_tb <- .data %>%
    
    dplyr::mutate(
      
      .original_text_chr =
        
        as.character(
          !!source_col
        )
      
    ) %>%
    
    dplyr::rowwise() %>%
    
    dplyr::mutate(
      
      .items_ls = list(
        
        split_numbered_items(
          
          as.character(
            !!source_col
          )
          
        )
        
      )
      
    ) %>%
    
    dplyr::ungroup() %>%
    
    tidyr::unnest(
      .items_ls
    ) %>%
    
    dplyr::mutate(
      
      !!item_nbr_col :=
        item_nbr_chr,
      
      !!item_seq_col :=
        item_seq_int,
      
      !!item_text_col :=
        item_text_chr,
      
      !!pre_item_text_col :=
        pre_item_text_chr,
      
      !!post_item_text_col :=
        post_item_text_chr
      
    )
  
  if (!is.null(
    original_text_col
  )) {
    
    outp_tb <-
      
      outp_tb %>%
      
      dplyr::mutate(
        
        !!original_text_col :=
          
          .original_text_chr
        
      )
    
  }
  
  outp_tb %>%
    
    dplyr::select(
      
      -tidyselect::any_of(
        c(
          "item_nbr_chr",
          "item_seq_int",
          "item_text_chr",
          "pre_item_text_chr",
          "post_item_text_chr",
          ".original_text_chr"
        )
      )
      
    )
  
}
transform_to_expand_unnumbered_perspectives <- function( # Could this be more general (e.g., pass a perspectives fn if req, etc)
    .data,
    source_col,
    variants_ls,
    perspective_col = "Perspective",
    perspective_text_col = NULL) {
  
  source_col <- rlang::ensym(
    source_col
  )
  
  perspective_col <- rlang::sym(
    perspective_col
  )
  
  if (!is.null(
    perspective_text_col
  )) {
    
    perspective_text_col <- rlang::sym(
      perspective_text_col
    )
    
  }
  
  extract_blocks <- function(
    text_chr,
    variants_ls) {
    
    if (
      is.na(text_chr)
    ) {
      
      return(
        
        tibble::tibble(
          
          perspective_chr =
            NA_character_,
          
          perspective_text_chr =
            NA_character_
          
        )
        
      )
      
    }
    
    matches_tb <- purrr::imap_dfr(
      
      variants_ls,
      
      function(
    patterns_chr,
    perspective_chr
      ) {
        
        purrr::map_dfr(
          
          patterns_chr,
          
          function(pattern_chr) {
            
            locs_mtx <-
              
              stringr::str_locate_all(
                
                string =
                  stringr::str_to_lower(
                    text_chr
                  ),
                
                pattern =
                  stringr::fixed(
                    stringr::str_to_lower(
                      pattern_chr
                    )
                  )
                
              )[[1]]
            
            if (
              nrow(
                locs_mtx
              ) == 0
            ) {
              
              return(
                tibble::tibble()
              )
              
            }
            
            tibble::tibble(
              
              perspective_chr =
                perspective_chr,
              
              pattern_chr =
                pattern_chr,
              
              start_int =
                locs_mtx[,1],
              
              end_int =
                locs_mtx[,2]
              
            )
            
          }
          
        )
        
      }
    
    )
    
    if (
      nrow(
        matches_tb
      ) == 0
    ) {
      
      return(
        
        tibble::tibble(
          
          perspective_chr =
            NA_character_,
          
          perspective_text_chr =
            text_chr
          
        )
        
      )
      
    }
    
    matches_tb <- matches_tb %>%
      
      dplyr::mutate(
        
        width_int =
          end_int -
          start_int
        
      ) %>%
      
      dplyr::group_by(
        start_int
      ) %>%
      
      dplyr::slice_max(
        width_int,
        n = 1,
        with_ties = FALSE
      ) %>%
      
      dplyr::ungroup() %>%
      
      dplyr::arrange(
        start_int
      )
    
    out_ls <- vector(
      "list",
      nrow(matches_tb)
    )
    
    for (
      
      i in seq_len(
        nrow(matches_tb)
      )
      
    ) {
      
      start_pos <-
        matches_tb$end_int[i] + 1
      
      end_pos <-
        
        if (
          i < nrow(matches_tb)
        ) {
          
          matches_tb$
            start_int[
              i + 1
            ] - 1
          
        } else {
          
          nchar(
            text_chr
          )
          
        }
      
      out_ls[[i]] <-
        
        tibble::tibble(
          
          perspective_chr =
            
            matches_tb$
            perspective_chr[
              i
            ],
          
          perspective_text_chr =
            
            stringr::str_squish(
              
              substr(
                
                text_chr,
                
                start_pos,
                
                end_pos
                
              )
              
            )
          
        )
      
    }
    
    dplyr::bind_rows(
      out_ls
    )
    
  }
  
  numbered_tb <- .data %>%
    
    dplyr::filter(
      !is.na(
        CE_Item_Nbr
      )
    )
  
  unnumbered_tb <- .data %>%
    
    dplyr::filter(
      is.na(
        CE_Item_Nbr
      )
    ) %>%
    
    dplyr::rowwise() %>%
    
    dplyr::mutate(
      
      .perspectives_ls = list(
        
        extract_blocks(
          
          as.character(
            !!source_col
          ),
          
          variants_ls
          
        )
        
      )
      
    ) %>%
    
    dplyr::ungroup() %>%
    
    tidyr::unnest(
      .perspectives_ls
    ) %>%
    
    dplyr::mutate(
      
      !!perspective_col :=
        
        dplyr::coalesce(
          
          .data[[rlang::as_string(
            perspective_col
          )]],
          
          perspective_chr
          
        )
      
    )
  
  if (
    
    !is.null(
      perspective_text_col
    )
    
  ) {
    
    unnumbered_tb <-
      
      unnumbered_tb %>%
      
      dplyr::mutate(
        
        !!perspective_text_col :=
          
          perspective_text_chr
        
      )
    
  } else {
    
    unnumbered_tb <-
      
      unnumbered_tb %>%
      
      dplyr::mutate(
        
        !!source_col :=
          
          perspective_text_chr
        
      )
    
  }
  
  dplyr::bind_rows(
    
    numbered_tb,
    
    unnumbered_tb
    
  ) %>%
    
    dplyr::select(
      
      -tidyselect::any_of(
        c(
          "perspective_chr",
          "perspective_text_chr",
          ".perspectives_ls"
        )
      )
      
    )
  
}
transform_to_extract_ce_perspective <- function(
    .data,
    source_col,
    variants_ls,
    perspective_col = "Perspective",
    original_text_col = NULL,
    matched_pattern_col = NULL) {
  
  source_col <- rlang::ensym(
    source_col
  )
  
  perspective_col <- rlang::sym(
    perspective_col
  )
  
  if (!is.null(
    original_text_col
  )) {
    
    original_text_col <- rlang::sym(
      original_text_col
    )
    
  }
  
  if (!is.null(
    matched_pattern_col
  )) {
    
    matched_pattern_col <- rlang::sym(
      matched_pattern_col
    )
    
  }
  
  outp_tb <- .data %>%
    
    dplyr::mutate(
      
      .working_chr =
        
        as.character(
          !!source_col
        ),
      
      .perspective_chr =
        NA_character_,
      
      .matched_pattern_chr =
        NA_character_
      
    )
  
  for (
    perspective_chr in
    names(
      variants_ls
    )
  ) {
    
    for (
      pattern_chr in
      variants_ls[[perspective_chr]]
    ) {
      
      matched_lgl <-
        
        is.na(
          outp_tb$.perspective_chr
        ) &
        
        stringr::str_detect(
          
          string =
            stringr::str_to_lower(
              outp_tb$.working_chr
            ),
          
          pattern =
            stringr::fixed(
              stringr::str_to_lower(
                pattern_chr
              )
            )
          
        )
      
      outp_tb <-
        
        outp_tb %>%
        
        dplyr::mutate(
          
          .perspective_chr =
            
            dplyr::if_else(
              
              matched_lgl,
              
              perspective_chr,
              
              .perspective_chr
              
            ),
          
          .matched_pattern_chr =
            
            dplyr::if_else(
              
              matched_lgl,
              
              pattern_chr,
              
              .matched_pattern_chr
              
            )
          
        )
      
    }
    
  }
  
  if (!is.null(
    original_text_col
  )) {
    
    outp_tb <- outp_tb %>%
      
      dplyr::mutate(
        
        !!original_text_col :=
          
          .working_chr
        
      )
    
  }
  
  outp_tb <- outp_tb %>%
    
    dplyr::mutate(
      
      !!source_col :=
        
        purrr::map2_chr(
          
          .working_chr,
          
          .matched_pattern_chr,
          
          function(
    text_chr,
    pattern_chr
          ) {
            
            if (
              !is.na(pattern_chr)
            ) {
              
              text_chr <-
                
                stringr::str_remove(
                  
                  string =
                    stringr::str_trim(
                      text_chr
                    ),
                  
                  pattern =
                    stringr::fixed(
                      pattern_chr,
                      ignore_case = TRUE
                    )
                  
                )
              
            }
            
            stringr::str_squish(
              text_chr
            )
            
          }
    
        ),
    
    !!perspective_col :=
      
      .perspective_chr
    
    )
  
  if (!is.null(
    matched_pattern_col
  )) {
    
    outp_tb <- outp_tb %>%
      
      dplyr::mutate(
        
        !!matched_pattern_col :=
          
          .matched_pattern_chr
        
      )
    
  }
  
  outp_tb %>%
    
    dplyr::select(
      
      -.working_chr,
      -.perspective_chr,
      -.matched_pattern_chr
      
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
transform_to_extract_trailing_text <- function(
    .data,
    source_col,
    trailing_text_col = "CE_Trailing_Text",
    section_patterns_chr =
      c(
        "Incremental NMB:",
        "ICER - ",
        "ICER:",
        "ROI:",
        "INMB"
      )) {
  
  source_col <- rlang::ensym(
    source_col
  )
  
  trailing_text_col <- rlang::sym(
    trailing_text_col
  )
  
  outp_tb <- .data
  
  outp_tb[[rlang::as_string(
    trailing_text_col
  )]] <- NA_character_
  
  for (
    
    i in seq_len(
      nrow(outp_tb)
    )
    
  ) {
    
    txt_chr <-
      
      outp_tb[[rlang::as_string(
        source_col
      )]][i]
    
    if (
      is.na(txt_chr)
    ) {
      
      next
    }
    
    locs_int <- integer()
    
    for (
      pattern_chr in
      section_patterns_chr
    ) {
      
      loc_int <-
        
        stringr::str_locate(
          txt_chr,
          stringr::fixed(
            pattern_chr
          )
        )[1]
      
      if (
        !is.na(loc_int)
      ) {
        
        locs_int <-
          c(
            locs_int,
            loc_int
          )
        
      }
      
    }
    
    if (
      length(
        locs_int
      ) == 0
    ) {
      
      next
    }
    
    split_int <-
      min(
        locs_int
      )
    
    if (
      split_int < 25
    ) {
      
      next
    }
    
    outp_tb[[rlang::as_string(
      trailing_text_col
    )]][i] <-
      
      stringr::str_squish(
        
        substr(
          
          txt_chr,
          
          split_int,
          
          nchar(
            txt_chr
          )
          
        )
        
      )
    
    outp_tb[[rlang::as_string(
      source_col
    )]][i] <-
      
      stringr::str_squish(
        
        substr(
          
          txt_chr,
          
          1,
          
          split_int - 1
          
        )
        
      )
    
  }
  
  outp_tb
  
}
transform_to_extract_variant_segment <- function(
    .data,
    source_col,
    text_col,
    multiple_matches_col,
    variants_ls,
    delimiter_chr = ":",
    end_patterns_chr = NULL,
    continuation_patterns_chr = NULL,
    fallback_to_end_patterns_lgl = TRUE,
    remainder_col = NULL,
    original_text_col = NULL,
    debug_lgl = FALSE) {
  
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
  
  if (!is.null(original_text_col)) {
    
    original_text_col <- rlang::ensym(
      original_text_col
    )
    
  }
  
  outp_tb <- .data %>%
    
    dplyr::rowwise() %>%
    
    dplyr::mutate(
      
      .original_text_chr =
        as.character(
          !!text_col
        ),
      
      .result_ls = list(
        
        if (is.na(
          !!multiple_matches_col
        )) {
          
          list(
            segment_chr = as.character(
              !!text_col
            ),
            remainder_chr = NA_character_,
            matched_start_pattern_chr = NA_character_,
            matched_terminal_pattern_chr = NA_character_,
            used_fallback_lgl = FALSE
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
            
            end_patterns_chr =
              end_patterns_chr,
            
            continuation_patterns_chr =
              continuation_patterns_chr,
            
            fallback_to_end_patterns_lgl =
              fallback_to_end_patterns_lgl
            
          )
          
        }
        
      )
      
    ) %>%
    
    dplyr::ungroup() %>%
    
    dplyr::mutate(
      
      !!text_col := purrr::map_chr(
        .result_ls,
        "segment_chr"
      )
      
    )
  
  if (!is.null(remainder_col)) {
    
    outp_tb <- outp_tb %>%
      
      dplyr::mutate(
        
        !!remainder_col :=
          purrr::map_chr(
            .result_ls,
            "remainder_chr"
          )
        
      )
    
  }
  
  if (!is.null(original_text_col)) {
    
    outp_tb <- outp_tb %>%
      
      dplyr::mutate(
        
        !!original_text_col :=
          .original_text_chr
        
      )
    
  }
  
  if (isTRUE(debug_lgl)) {
    
    outp_tb <- outp_tb %>%
      
      dplyr::mutate(
        
        matched_start_pattern_chr =
          purrr::map_chr(
            .result_ls,
            "matched_start_pattern_chr"
          ),
        
        matched_terminal_pattern_chr =
          purrr::map_chr(
            .result_ls,
            "matched_terminal_pattern_chr"
          ),
        
        used_fallback_lgl =
          purrr::map_lgl(
            .result_ls,
            "used_fallback_lgl"
          )
        
      )
    
  }
  
  outp_tb %>%
    
    dplyr::select(
      -.result_ls,
      -.original_text_chr
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
