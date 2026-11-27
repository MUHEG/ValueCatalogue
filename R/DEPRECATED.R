# get_temporal_segments <- function(
    #     text_chr,
#     patterns_chr) {
#   
#   if (
#     is.na(text_chr)
#   ) {
#     
#     return(
#       tibble::tibble(
#         Time_Horizon = NA_character_,
#         Temporal_Result_Text = NA_character_
#       )
#     )
#     
#   }
#   
#   matches_tb <-
#     
#     purrr::map_dfr(
#       
#       patterns_chr,
#       
#       function(pattern_chr) {
#         
#         loc_mat <-
#           
#           stringr::str_locate_all(
#             
#             stringr::str_to_lower(
#               text_chr
#             ),
#             
#             stringr::fixed(
#               stringr::str_to_lower(
#                 pattern_chr
#               )
#             )
#             
#           )[[1]]
#         
#         if (
#           nrow(loc_mat) == 0
#         ) {
#           
#           return(NULL)
#           
#         }
#         
#         tibble::tibble(
#           
#           Horizon_Label =
#             pattern_chr,
#           
#           Start =
#             loc_mat[, "start"]
#           
#         )
#         
#       }
#       
#     )
#   
#   if (
#     nrow(matches_tb) <= 1
#   ) {
#     
#     return(
#       
#       tibble::tibble(
#         
#         Time_Horizon =
#           
#           if (
#             nrow(matches_tb) == 1
#           ) {
#             
#             matches_tb$Horizon_Label
#             
#           } else {
#             
#             NA_character_
#             
#           },
#         
#         Temporal_Result_Text =
#           text_chr
#         
#       )
#       
#     )
#     
#   }
#   
#   matches_tb <-
#     
#     matches_tb %>%
#     
#     dplyr::arrange(
#       Start
#     )
#   
#   purrr::map_dfr(
#     
#     seq_len(
#       nrow(matches_tb)
#     ),
#     
#     function(i) {
#       
#       start_pos <-
#         
#         matches_tb$Start[i]
#       
#       end_pos <-
#         
#         if (
#           
#           i <
#           nrow(matches_tb)
#           
#         ) {
#           
#           matches_tb$Start[i + 1] - 1
#           
#         } else {
#           
#           nchar(
#             text_chr
#           )
#           
#         }
#       
#       tibble::tibble(
#         
#         Time_Horizon =
#           
#           matches_tb$Horizon_Label[i],
#         
#         Temporal_Result_Text =
#           
#           stringr::str_squish(
#             
#             substr(
#               text_chr,
#               start_pos,
#               end_pos
#             )
#             
#           )
#         
#       )
#       
#     }
#     
#   )
#   
# }

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
# transform_to_extract_temporal_segments <- function(
    #     .data,
#     source_col,
#     horizon_patterns_ls) {
#   
#   source_col <- rlang::ensym(
#     source_col
#   )
#   
#   all_patterns_chr <-
#     
#     unique(
#       
#       unlist(
#         horizon_patterns_ls,
#         use.names = FALSE
#       )
#       
#     )
#   
# 
#   
#   .data %>%
#     
#     dplyr::mutate(
#       
#       Temporal_Segments =
#         
#         purrr::map(
#           
#           !!source_col,
#           
#           ~ get_temporal_segments(
#             .x,
#             all_patterns_chr
#           )
#           
#         )
#       
#     )
#   
# }