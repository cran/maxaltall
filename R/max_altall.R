# You can learn more about package authoring with RStudio at:
#
#   https://r-pkgs.org
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

#roxygen2 code to import external packages
#' @importFrom data.table as.data.table
#' @importFrom data.table :=
#' @importFrom dplyr case_when
#' @importFrom dplyr group_by
#' @importFrom dplyr mutate
#' @importFrom dplyr rename
#' @importFrom dplyr row_number
#' @importFrom dplyr select
#' @importFrom magrittr %>%
#' @importFrom tidyr pivot_wider

max_altall <- function(file, type, thres_1 = 0.8, thres_2 = 0.2, export = "FALSE", export_dir = tempdir()) {

  #Checks to make sure the argument values are valid

  if (file.exists(file) == FALSE) {
    stop("Error: the specified 'file' does not exist")
  }

  if (type != "aa" && type != "nuc") {
    stop("Error: 'type' must be either 'aa' or 'nuc'")
  }

  if (type == "aa" && (thres_1 > 1 | thres_1 < 0.05)) {
    stop("Error: 'thres_1' must be in the range of 0.05 to 1 for amino acid states")
  }

  if (type == "nuc" && (thres_1 > 1 | thres_1 < 0.25)) {
    stop("Error: 'thres_1' must be in the range of 0.25 to 1 for nucleotide states")
  }

  if (!is.character(export) || export != "FALSE" && export != "TRUE") {
    stop("Error: 'export' must be either 'FALSE' or 'TRUE'")
  }

  if (file.exists(export_dir) == FALSE) {
    stop("Error: the specified 'export_dir' does not exist")
  }

  #Introduce a tolerance value to avoid floating-point comparison issues

  tolerance <- sqrt(.Machine$double.eps)

  #Sanity checks on the values specified for the arguments thres_1 and thres_2

  if (type == "aa" && (thres_2 - tolerance > thres_1 | thres_2 - tolerance > (1 - thres_1) | thres_2 + tolerance < (1 - thres_1 )/19)) {
    stop("Error: 'thres_2' not in permissible range for amino acid states - see documentation")
  }

  if (type == "nuc" && (thres_2 - tolerance > thres_1 | thres_2 - tolerance > (1 - thres_1) | thres_2 + tolerance < (1 - thres_1)/3)) {
    stop("Error: 'thres_2' not in permissible range for nucleotide states - see documentation")
  }

  #Code to execute if the .state file contains amino acid states

  if (type == "aa") {

    #Read table and assign to variable "nodes"
    nodes <- read.table(file = file, header=TRUE)

    #Check to ensure that the correct type (aa) was selected by checking the number of columns in the inputted file
    if (ncol(nodes) != 23) {
      stop("Error: input file does not match 'aa' format.")
    }

    #rename column headers as one letter amino acid abbreviations
    nodes <- nodes %>% rename("A" = "p_A", "R" = "p_R", "N" = "p_N", "D" = "p_D", "C" = "p_C", "Q" = "p_Q","E" = "p_E","G" = "p_G","H" = "p_H","I" = "p_I","L" = "p_L","K" = "p_K","M" = "p_M","F" = "p_F","P" = "p_P","S" = "p_S","T" = "p_T","W" = "p_W","Y" = "p_Y","V" = "p_V")


    #Create a new column named "max" which holds the state (amino acid or nucleotide residue) with the highest probability at each site for each node
    nodes$max <- colnames(nodes[,(4:23)])[apply(nodes[,(4:23)],1,which.max)]

    #Use vectorized operations and convert the nodes table into a data table to improve code efficiency
    x <- c()
    y <- c()
    nodes_dt <- as.data.table(nodes)

    #Create two new columns, one with the probability of the max (highest probability) state and the other with the probability of the alt (second highest probability) state
    #In the event that there are ties for the highest or second highest state, the selected state defaults to whichever one appears earliest in the state list.
    nodes_dt[, c("max_prob", "max2_prob") := {vals <- as.matrix(.SD); list(apply(vals, 1, function(row) if (all(is.na(row))) NA else sort(row, decreasing = TRUE, na.last = TRUE)[1]), apply(vals, 1, function(row) if (length(na.omit(row)) < 2) NA else sort(row, decreasing = TRUE, na.last = TRUE)[2]))}, .SDcols = 4:23]

    #Create a new column named "max2" which holds the actual state (not probability as in the previous command) that has the second highest probability
    #A tie in the second highest probabilities will result in the state that appears earlier in the list to be selected as max2
    value_columns <- names(nodes_dt)[4:23]; nodes_dt[, max2 := {values <- as.matrix(.SD); second_max_col <- function(values_row) {if (sum(!is.na(values_row)) < 2) return(NA); sorted_indices <- order(-values_row); value_columns[sorted_indices[2]]}; apply(values, 1, second_max_col)}, .SDcols = value_columns]
  }


  #Code to execute if the .state file contains nucleotide states
  if (type == "nuc") {

    #Read table and assign to variable "nodes"
    nodes <- read.table(file = file, header=TRUE)

    #Check to ensure that the correct type (nuc) was selected by checking the number of columns in the input file
    if (ncol(nodes) != 7) {
      stop("Error: input file does not match 'nuc' format.")
    }

    #rename column headers as one letter nucleotide abbreviations
    nodes <- nodes %>% rename("A" = "p_A", "C" = "p_C", "G" = "p_G","T" = "p_T")

    #Create a new column named "max" which holds the state (specific nucleotide) with the highest probability at each site for each node
    nodes$max <- colnames(nodes[,(4:7)])[apply(nodes[,(4:7)],1,which.max)]

    #Use vectorized operations and convert the nodes table into a data table to improve code efficiency
    x <- c()
    y <- c()
    nodes_dt <- as.data.table(nodes)

    #Create two new columns, one with the probability of the max (highest probability) state and the other with the probability of the alt (second highest probability) state

    #In the event that there are ties for the highest or second highest state probability, the selected state defaults to whichever one appears earlier in the state list

    nodes_dt[, c("max_prob", "max2_prob") := {vals <- as.matrix(.SD); list(apply(vals, 1, function(row) if (all(is.na(row))) NA else sort(row, decreasing = TRUE, na.last = TRUE)[1]), apply(vals, 1, function(row) if (length(na.omit(row)) < 2) NA else sort(row, decreasing = TRUE, na.last = TRUE)[2]))}, .SDcols = 4:7]

    #Create a new column named "max2" which holds the actual state with the second highest probability
    #A tie in the second highest probabilities will result in the state that appears earlier in the list to be selected as max2
    value_columns <- names(nodes_dt)[4:7]; nodes_dt[, max2 := {values <- as.matrix(.SD); second_max_col <- function(values_row) {if (sum(!is.na(values_row)) < 2) return(NA); sorted_indices <- order(-values_row); value_columns[sorted_indices[2]]}; apply(values, 1, second_max_col)}, .SDcols = value_columns]
  }


  #Scan max_prob and max2_prob column to determine if substitution of states will occur
  #If the probability of the highest probability state (max) is less than or equal to thres_1 (0.8 by default) and the probability of the next highest probability state (max2) is greater than or equal to thres_2 (0.2 by default), then place max2 into new column "alt".
  #If probability of max state > thres_1 or the probability of max2 state < thres_2, then place the highest probability state into new column "alt" so that alt = max
  nodes <- nodes_dt %>% mutate(alt = case_when(max_prob <= thres_1 & max2_prob >= thres_2 ~ max2, TRUE ~ max))


  #Create modified data tables with only two columns from the original table: the node number and max or alt state
  nodes_mod1 <- subset(nodes, select = c(Node, max))
  nodes_mod2 <- subset(nodes, select = c(Node, alt))

  #Reshape the previously created modified data tables so that each node has its own column of max or alt states arranged in order of position (essentially the complete sequence but written vertically)
  nodes_mod3 <- nodes_mod1 %>% group_by(Node) %>% mutate(row = row_number()) %>% pivot_wider(names_from = Node, values_from = max) %>% select(-row)
  nodes_mod4 <- nodes_mod2 %>% group_by(Node) %>% mutate(row = row_number()) %>% pivot_wider(names_from = Node, values_from = alt) %>% select(-row)

  if (export == FALSE) {

  #Deposit complete maximum likelihood sequence for each node
  for (col_name in colnames(nodes_mod3)) {

    #Extract column values of the max states
    column_values <- nodes_mod3[[col_name]]

    #Create new variable “nodename” to identify each string to be created
    nodename <- paste(">",col_name,"_ML", "\n", sep = "")

    #Create a string of all the max column values without spaces
    values_string <- paste(column_values, collapse = "")

    #Add each new string to the text file “node_sequences_all.txt”
    cat(nodename, sep = "")
    cat(values_string, "\n\n", sep = "")
  }

  #Deposit complete alt sequence of each node
  for (col_name in colnames(nodes_mod4)) {

    #Extract column values of the alt states
    column_values <- nodes_mod4[[col_name]]

    #Create new variable “nodename” to identify each string to be created
    nodename <- paste(">",col_name,"_altall","\n", sep = "")

    #Create a string of all the alt column values without spaces
    values_string <- paste(column_values, collapse = "")

    #Adds each new string to the same text file “node_sequences_all.txt”
    cat(nodename, sep = "")
    cat(values_string, "\n\n", sep = "")
  }
  }
  if (export == TRUE) {

  #Create a new text file for depositing complete sequences of each node, which will contain both the max and alt sequences
  txt_file_path <- paste0(export_dir, "/node_sequences_all.txt")
  con <- file(txt_file_path, "w")

  #Deposit complete maximum likelihood sequence for each node
  for (col_name in colnames(nodes_mod3)) {

    #Extract column values of the max states
    column_values <- nodes_mod3[[col_name]]

    #Create new variable “nodename” to identify each string to be created
    nodename <- paste(">",col_name,"_ML", "\n", sep = "")

    #Create a string of all the max column values without spaces
    values_string <- paste(column_values, collapse = "")

    #Add each new string to the text file “node_sequences_all.txt”
    cat(nodename, file = con, sep = "")
    cat(values_string, "\n\n", file = con, sep = "")
  }


  #Deposit complete alt sequence of each node
  for (col_name in colnames(nodes_mod4)) {

    #Extract column values of the alt states
    column_values <- nodes_mod4[[col_name]]

    #Create new variable “nodename” to identify each string to be created
    nodename <- paste(">",col_name,"_altall","\n", sep = "")

    #Create a string of all the alt column values without spaces
    values_string <- paste(column_values, collapse = "")

    #Adds each new string to the same text file “node_sequences_all.txt”
    cat(nodename, file = con, sep = "")
    cat(values_string, "\n\n", file = con, sep = "")
  }

  close(con)
  }
}
