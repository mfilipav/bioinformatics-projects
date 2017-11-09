# Infer the phylogeny tree from multiple DNA sequences
# Clustering by UPGMA algorithm, and manhattan distances

# Inspired by Computational Biology course at ETH Zürich
# By: Modestas Filipavicius
# Date: 2017-11-13

# Ape handles phylogenetic trees
library(ape)

# Given edges, vector of lengths and node info dataframe, make a tree Ape object
# the dataframe  containing node descriptions.
convert_to_phylo = function(sequences, connected_edges, edge_lengths, node_info) {
    
    edges = connected_edges
    
    for (name in rownames(node_info)) {
        index = which(rownames(node_info) == name)
        edges[which(edges == name)] = as.numeric(index)
    }
    edges = matrix(as.numeric(edges), ncol = 2)
    
    edges[which(edges > length(sequences))] = - edges[which(edges > length(sequences))]
    root = setdiff(edges[,1], edges[,2])
    edges[which(edges==root)] = length(sequences) + 1
    
    k = length(sequences) + 2
    for(x in unique(edges[which(edges < 0)])) {
        edges[which(edges==x)] = k
        k = k + 1
    }
    
    # Construct object
    tree = list()
    class(tree) = "phylo"
    tree$edge = edges
    tree$edge.length = edge_lengths
    tree$Nnode = as.integer(length(sequences) - 1)
    tree$tip.label = names(sequences)
    
    # Return the tree as abject in Ape's phylo class
    return(tree)
}

# Plot the phylogenetic tree with node labels and edge lengths.
tree_plot = function(tree) {
    plot(tree)
    edgelabels(format(tree$edge.length, digits = 2))
}

# Initialize the dataframe which holds the node info
initialize_node_info = function(sequences) {

    N = length(sequences)
    node_names = names(sequences)
    node_sizes = rep(1, times = N)
    node_heights = rep(0, times = N)
    node_info = data.frame(node_sizes, node_heights)
    rownames(node_info) = node_names

    # Return a data frame that contains information on currently existing tip nodes.
    # node_info: a dataframe containing information on the currently existing nodes.
    #                   The row names are the names of the currently existing tips, i.e.
    #                   are the same as the names in the sequence list, node_height is
    #                   0 and node_size is 1 as the all the currently existing nodes are tips.
    return(node_info)
}


#### add_new_node ####
add_new_node = function(node_info, merging_nodes) {
    # Add new merged node to the node description data frame.
    # The new node is a combination of the nodes supplied in the merging_nodes,
    # e.g. if one needs to merge nodes "bird" and "fish", the new node in the
    # dataframe will be called "bird.fish".
    #    node_info: the dataframe created by initialize_node_info, containing
    #                      current node sizes and node heights
    #    merging_nodes: a vector of two names of the nodes being merged

    new_node_name = paste(merging_nodes, collapse = ".")
    new_node_row = data.frame(node_sizes = 0, node_heights = 0)
    rownames(new_node_row) = new_node_name
    new_node_info = rbind(node_info, new_node_row)
    
    # Return the node_info dataframe with a row for the new node added, and
    # the new node name.
    #    node_info: the dataframe where the rows are labelled by current nodes and columns
    #                      contain the node heights and sizes.
    #    new_node_name: the name of the newly added node, created from names in merging_nodes.
    return(list(node_info = new_node_info, new_node_name = new_node_name))
}

#### Distance functions ####
# Returns hamming distance given two sequences
get_hamming_distance = function(sequence1, sequence2) {

    seq1 = unlist(strsplit(sequence1, ""))
    seq2 = unlist(strsplit(sequence2, ""))
    
    distance = 0
    for (i in 1: length(seq1)) {

        if (seq1[i] != seq2[i]) {
            distance = distance + 1
        }
    }
    
  # Return the numerical value of the distance
  return(distance)
}


get_JC69_distance = function(sequence1, sequence2) {

    hamm_dist = get_hamming_distance(sequence1, sequence2)
    distance = -0.75*log(1 - 1.3333333333* hamm_dist/nchar(sequence1)) 
    return(distance)
}


get_K80_distance = function(sequence1, sequence2) {

    s1 = unlist(strsplit(sequence1, ""))
    s2 = unlist(strsplit(sequence2, ""))
    
    S = 0
    V = 0
    for (i in 1:length(s1)) {
        
        if ((s1[i]=="T" && s2[i]=="C") || (s1[i]=="C" && s2[i]=="T") ||
            (s1[i]=="A" && s2[i]=="G") || (s1[i]=="G" && s2[i]=="A")) {
            S = S +1
        }
        
        if ((s1[i]=="A" && s2[i]=="C") || (s1[i]=="C" && s2[i]=="A") ||
            (s1[i]=="A" && s2[i]=="T") || (s1[i]=="T" && s2[i]=="A") ||
            (s1[i]=="G" && s2[i]=="C") || (s1[i]=="C" && s2[i]=="G") ||
            (s1[i]=="G" && s2[i]=="T") || (s1[i]=="T" && s2[i]=="G")) {
            V = V +1
        }
    }
    
    V = V/length(s1)
    S = S/length(s1)
    distance = -0.5 * log(1-2*S-V) - 0.25*log(1-2*V)
    
    return(distance)
}

#### precompute_dist_matrix ####
# Compute the initial distance matrix using one of the distance measures, for
# nxn matrix, where n = # of seqs
precompute_dist_matrix = function(sequences, distance_measure) {
    
        N = length(sequences)
        distance_matrix = matrix(nrow = N, ncol = N)
        diag(distance_matrix) = Inf
        rownames(distance_matrix) = names(sequences)
        colnames(distance_matrix) = names(sequences)
          
  
      # Assign one of three dist functions to be used, depending on which
      # string representation of that function was called in the 
      # 'precompute_dist_matrix' function head... 
      if (distance_measure == "hamming") {
          dist_fun = get_hamming_distance
      } else if (distance_measure == "JC69") {
          dist_fun = get_JC69_distance
      } else if (distance_measure == "K80") {
          dist_fun = get_K80_distance
      }
      
      # Traverse all matrix without diagonal:
      for (row in 1:N) {
          for (col in 1:N) {
              if (row != col) {
                  distance_matrix[row,col] = dist_fun(sequences[[row]], sequences[[col]])
              }
          }
      }
      
      # Return the NxN matrix of inter-sequence distances with Inf on the diagonal
      return(distance_matrix)
}


#### get_merge_node_distance ####
get_merge_node_distance = function(merging_nodes, existing_node, node_info, distance_matrix) {

    node_sizes_ls = list()
    #node_sizes_ls[1] = df$node_sizes[which(row.names(df) == existing_node)]
    
    node_sizes_ls[1] = node_info$node_sizes[which(row.names(node_info) == merging_nodes[1])]
    node_sizes_ls[2] = node_info$node_sizes[which(row.names(node_info) == merging_nodes[2])]
    
    # Collect distances between merging nodes and existing nodes
    node_distances_ls = list()
    existing_node_idx = which(colnames(distance_matrix) == existing_node)
    for (i in 1:length(merging_nodes)) {
        
        node_distances_ls[i] = distance_matrix[existing_node_idx, which(colnames(distance_matrix) == merging_nodes[i])]
    }
    
    # Calculate d(existing_node, mergining_nodes)
    new_distance = (node_sizes_ls[[1]]*node_distances_ls[[1]] + 
                        node_sizes_ls[[2]]*node_distances_ls[[2]]) /
                        (sum(unlist(node_sizes_ls)))

    return(new_distance)
}

#### update_distance_matrix ####
update_distance_matrix = function(distance_matrix, merging_nodes, new_node_name, node_info) {
  
    
    # Create two lists with indeces for sequences that are merging and 
    # those that are not merging
    not_merging_node_idx_ls = which(!rownames(distance_matrix) %in% merging_nodes)
    merging_node_idx_ls = which(rownames(distance_matrix) %in% merging_nodes)
    
    # For every sequence in dist matrix, calculate its distance to the newly
    # merged sequence. Collect Infs as well. Return a list.
    new_merged_dist_col_ls = list()
    for (i in 1:length(distance_matrix[1, ])) {
        existing_node = colnames(distance_matrix)[i] # i=1 is "a", i=2 is "b", ..
        new_dist = get_merge_node_distance(merging_nodes, existing_node, 
                                           node_info, distance_matrix)
        #print(new_dist)
        new_merged_dist_col_ls = append(new_merged_dist_col_ls, new_dist)
    }
    
    # Final cleanup for list of new distances: unlist and add terminal "Inf"
    new_merged_dist_col_ls = unlist(new_merged_dist_col_ls)
    new_merged_dist_col_ls = append(new_merged_dist_col_ls, Inf)
    
    # Update the distance matrix by renaming rows and cols,
    # and retaining only the unmerged rows and cols
    #updated_distance_matrix = distance_matrix[not_merging_node_idx_ls, not_merging_node_idx_ls]
    updated_distance_matrix = cbind(distance_matrix, new_merged_dist_col_ls)
    updated_distance_matrix = rbind(updated_distance_matrix, new_merged_dist_col_ls)
    
    rownames(updated_distance_matrix)[length(rownames(updated_distance_matrix))] = new_node_name
    colnames(updated_distance_matrix)[length(colnames(updated_distance_matrix))] = new_node_name
    
    remove = which(rownames(updated_distance_matrix) %in% merging_nodes)
    
    # must add "as.matrix", otherwise when matrix shrinks to a single "Inf",
    # it is converted to a class "numeric"
    updated_distance_matrix = as.matrix(updated_distance_matrix[-remove, -remove])
    # print(updated_distance_matrix)
    
    return(updated_distance_matrix)
}


#### upgma_single_iter ####
# Single iteration of UPGMA algorithm
upgma_single_iter = function(distance_matrix, edges, edge_lengths, node_info) {
 
    # Find coordinates for lowest distance between two nodes
    # in both upper and lower matrix triangles
    coords = which(distance_matrix == min(unlist(distance_matrix), na.rm = TRUE), 
                   arr.ind = TRUE)
    
    # Pick the coordinates in upper matrix triangle (for consistency)
    # Ex. [1] 3 4
    min_coords = as.vector(coords[which(coords[, 1] < coords[, 2]), ])
    
    # Distance between to merging nodes:
    merging_nodes_dist = unlist(distance_matrix[min_coords[1], min_coords[2]])
    
    # Make a list with row (1st element in min_coords), col(2nd) names of 
    # corresponding lowest distance. Ex. [1] "c" "d"
    merging_nodes = c(rownames(distance_matrix)[min_coords[1]],
                      colnames(distance_matrix)[min_coords[2]])

    # Make a new merged node name, "a" + "b" = "a.b"
    new_node_name = add_new_node(node_info, merging_nodes)[[2]]
    
    node_des_ixd1 = which(rownames(node_info) == merging_nodes[1])
    node_des_ixd2 = which(rownames(node_info) == merging_nodes[2])
    branch1 = merging_nodes_dist/2 - node_info$node_heights[node_des_ixd1]
    branch2 = merging_nodes_dist/2 - node_info$node_heights[node_des_ixd2]
    
    # originally initiated in build_upgma_tree
    #edges = matrix(nrow = 0, ncol = 2)
    #edge_lengths = c()
    
    new_edge1 = c(new_node_name, merging_nodes[[1]])
    new_edge2 = c(new_node_name, merging_nodes[[2]])
    edges = rbind(edges, new_edge1)
    edges = rbind(edges, new_edge2)
    rownames(edges) = NULL
    
    # update edge lengths
    edge_lengths = append(edge_lengths, branch1)
    edge_lengths = append(edge_lengths, branch2)
    
    
    # NEW NODE DESCRIPTION
    # Adds a new line to node_info, makes a new merged node name
    # return(list(node_info = new_node_info, new_node_name = new_node_name))
    new_node_info = add_new_node(node_info, merging_nodes)[[1]]
    
    # Update node description df: add up merging nodes sizes, recalculate heights
    merging_node_index = which(rownames(new_node_info) == new_node_name)
    
    
    node_info_idx1 = which(rownames(node_info) == merging_nodes[1])
    node_info_idx2 = which(rownames(node_info) == merging_nodes[2])
    
    new_node_size = node_info$node_sizes[node_info_idx1] +
        node_info$node_sizes[node_info_idx2]
    
    new_node_height = merging_nodes_dist/2
    
    new_node_info[merging_node_index, ] = c(new_node_size, new_node_height)
    
    
    # Remove newly merged nodes from node description

    # new_node_info = new_node_info[-c(node_des_ixd1, node_des_ixd2), ]
    
    
    # Now we can recalculate distance matrix
    distance_matrix = update_distance_matrix(distance_matrix, merging_nodes, 
                                                     new_node_name, node_info)
    
    # Re-assign node_info with new_node_info, once the dist
    # matrix is updated
    node_info = new_node_info
    
    return(list(distance_matrix = distance_matrix, edges = edges, edge_lengths = edge_lengths, node_info = node_info))
}

#### build_upgma_tree ####
build_upgma_tree = function(sequences, distance_measure) {
    
    N = length(sequences)
    node_info = initialize_node_info(sequences)
    edges = matrix(nrow = 0, ncol = 2)
    edge_lengths = vector(mode = "numeric", length = 0)
  
    # Initialize distance matrix:
    distance_matrix = precompute_dist_matrix(sequences, distance_measure)
  
    while (dim(distance_matrix)[1] > 1) {
        result = upgma_single_iter(distance_matrix, edges, 
                                 edge_lengths, node_info)
        
        distance_matrix = result$distance_matrix
        edges = result$edges
        edge_lengths = result$edge_lengths
        node_info = result$node_info
    }
    # Return the UPGMA tree of sequences
    tree = convert_to_phylo(sequences, edges, edge_lengths, node_info)
    return(tree)
}

test_tree_building = function() {
  sequences = list(orangutan = "TCACACCTAAGTTATATCTATATATAATCCGGGCCGGG",
                   chimp =     "ACAGACTTAAAATATACATATATATAATCCGGGCCGGG",
                   human =     "AAAGACTTAAAATATATATATATATAATCCGGGCCGGG",
                   gorilla =   "ACACACCCAAAATATATATATATATAATCCGGGCCGGG",
                   unicorn =   "CCACACCCAAAATACCCCCGCGTATAATCCGGGCCCCC")
  
  distance_measure = 'hamming'
  #distance_measure = 'JC69'
  #distance_measure = 'K80'
  
  tree = build_upgma_tree(sequences, distance_measure)
  tree_plot(tree)
}

test_tree_building()
