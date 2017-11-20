# Simulates evolution on input DNA sequences. Nucleotide mutations are encoded in
# the Substitution model (matrix Q), which reflects the relative rate of 
# transitions (A<-->G, C<-->T) to transversions (A<-->C, G<-->T)

library(RUnit)
library(Matrix)
library(ape)

nucleotides = c("T", "C", "A", "G")

convert_to_nucleotides = function(sequence) {
    # Convert from numbers to letters, 1234 --> ACGT
    nucl_sequence = paste0(lapply(sequence, get_nucleotide_from_number), collapse = "")
    return(nucl_sequence)
}

get_nucleotide_from_number = function(i) {
    # Transform a nucleotide number into the appropriate letter
    return(nucleotides[i])
}

get_first_nucleotide = function(pi) {
    # Randomly generate a nucleotide given initial nucleotide probabilities
    # Make a cumulative sum vector from pi vector
    pi_cum = cumsum(pi)
    
    # Generate one number between 0 and 1
    rand = runif(1, min = 0, max = 1)
    
    # Take the first occurance of pi_cum element which is bigger or equal 
    # to generated random number rand
    nucleotide = which(pi_cum >= rand)[1]    
    
    # Return the sampled nucleotide
    # nucleotide: integer nucleotide value
    return(nucleotide)
}


initialize_substitution_matrix = function(pi, alpha1, alpha2, beta) {
    # Create transition rate matrix Q
    
    t = pi * c(0, alpha1, beta, beta)
    
    ab_matrix = rbind(c(0, alpha1, beta, beta), c(alpha1, 0, beta, beta), 
                      c(beta, beta, 0, alpha2), c(beta, beta, alpha2, 0))
    
    # (1X4) .x (4x4) = (4x4)
    Q = t(pi * ab_matrix)
    
    # Fill up the diagonals so that rows sum to zero
    for (i in 1:4) {
            Q[i, i] = -sum(Q[i, ])
    }
    
    # Return the transition rate matrix
    # Q: 4 by 4 matrix of rates
    return(Q)
}

get_starting_sequence = function(pi, N) {
    # Sample a starting sequence of length N
    
    starting_sequence = vector()
    
    # Fill up vector with values from {1, 2, 3, 4}
    for (i in 1:N) {
        starting_sequence[i] = get_first_nucleotide(pi)
    }
    
    # Return the sampled sequence
    # starting_sequence: vector of integer nucleotide values
    return(starting_sequence)
}

evolve_single_nucleotide = function(nucleotide, branch_length, Q) {
    # Evolve a given nucleotide along a branch of specified length.
    # nucleotide: nucleotide at the beginning of the branch
    # branch_length: the length of the branch along which evolution happens
    # Q: the transition rate matrix
    
    # debugging
    # Q = matrix
    # nucleotide = 3
    # branch_length = 1
    # 
    # Eq 1.16 and 1.17 explain spectral decomposition of Q, where Lambda is
    # L = diag{exp(lambda1*t), ... , exp(lamda4*t)}, t = branch_length
    # P(t) = exp(Qt) = U * L * U_inv
    # Sauce:
    # http://www.stats.ox.ac.uk/__data/assets/pdf_file/0014/5063/Iceland_Yang.pdf
    # Lambda = diag(exp(eigen(Q)$values*branch_length))
    # 
    # U = eigen(Q)$vectors
    # U_inv = solve(U)
    # # U  %*% U_inv
    # 
    # P = U %*% Lambda %*% U_inv
    # # sum(P[4, ]) == 1 # better return True!!!
    
    # Easier way to get the P(Qt)
    P = expm(Q*branch_length)
    
    # Make a cumulative sum vector from P vector
    P_cum = cumsum(P[nucleotide, ])
    
    # Generate one number between 0 and 1
    rand = runif(1, min = 0, max = 1)
    
    # Take the first occurance of P_cum element which is bigger or equal 
    # to generated random number rand
    evolved_nucleotide = which(P_cum >= rand)[1]
    
    # Return the nucleotide resulting from the evolutionary process along the branch
    # evolved_nucleotide: the new integer nucleotide value at the end of the branch
    return(evolved_nucleotide)
}


get_evolved_sequence = function(sequence, branch_length, Q) {
    # Evolve a given nucleotide sequence along a branch of specified length (bya)
    # sequence: nucleotide sequence at the beginning of the branch
    # branch_length: the length of the branch along which evolution happens
    # Q: the transition rate matrix

    evolved_sequence = vector()
    for (i in 1:length(sequence)) {
        evolved_sequence[i] = evolve_single_nucleotide(sequence[i], branch_length, Q)
    }
    
    # Return the nucleotide sequence after all positions have evolved along the given branch.
    # evolved_sequence: the vector of new integer nucleotide values at the end of the branch
    return(evolved_sequence)
}


simulate_evolution = function(newick_tree, pi, alpha1, alpha2, beta, N) {
    # Simulate evolution along the given tree.
    # newick_tree: the tree in newick text format

    # Transfrom the tree from text format to an object of the phylo class 
    # which represents the tree in R
    tree = read.tree(text = newick_tree)
    
    # Reorder the tree for easier traversing
    tree = reorder(tree, order = "cladewise")
    # > tree
    # Phylogenetic tree with 4 tips and 3 internal nodes.
    # Tip labels:
    #     [1] "orangutan" "gorilla"   "human"     "chimp"    
    # Rooted; includes branch lengths.
    
    # Set up the Q matrix
    Q = initialize_substitution_matrix(pi, alpha1, alpha2, beta)

    # Set up the starting sequence @ the root of the tree
    starting_sequence = get_starting_sequence(pi, N)
    
    # Prepare a list to store evolved sequences at each node
    sequence_per_node = list()
    sequence_per_node[[tree$edge[1,1]]] = starting_sequence
    
    # Walk the tree while evolving sequences along appropriate branches
    for (i in 1:length(tree$edge.length)) {
        # 6 tree edges: 13.00  2.75 10.25  4.75  5.50  5.50
        # > tree$edge
        # [,1] paren child length
        # [1,]    5    1    13
        # [2,]    5    6    2.75
        # [3,]    6    2    10.25
        # [4,]    6    7    4.75
        # [5,]    7    3    5.5
        # [6,]    7    4    5.5

        node_parent = tree$edge[i, 1]
        node_child = tree$edge[i, 2]
        branch_length = tree$edge.length[i]
        parent_sequence = sequence_per_node[[node_parent]]
        
        child_sequence = get_evolved_sequence(parent_sequence, branch_length, Q)

        sequence_per_node[[node_child]] = child_sequence
    }

    # Transform the alignment from nucleotide indices to nucleotide characters
    # and filter out the sequences at the tips
    alignment = list()
    
    # loop through 1 to # surviving descendants
    for (i in 1:length(tree$tip.label)) {
        alignment[[tree$tip.label[i]]] = convert_to_nucleotides(sequence_per_node[[i]])
    }
    
    # Return the simulated alignment.
    # The alignment should be in the form of a list where the tip label corresponds to the
    # appropriate simulated sequence, e.g. alignment$human = ACTG
    return(alignment)
}


test_simulation = function() {
    newick_tree = "(orangutan:13,(gorilla:10.25,(human:5.5,chimp:5.5):4.75):2.75);"
    N = 20
    beta = 1.35
    alpha1 = 44.229
    alpha2 = 21.781
    pi = c(0.22, 0.26, 0.33, 0.19)
    result = simulate_evolution(newick_tree, pi, alpha1, alpha2, beta, N)
    #return(result)
    # Check alignment
    print(result$orangutan)
    print(result$gorilla)
    print(result$human)
    print(result$chimp)
}

test_simulation()

# # measure aligned distances
# aligned = test_simulation() # modify to return sequences
# library(stringdist)
# bla = stringdistmatrix(parent_sequence, child_sequence)
# heatmap(bla) # look for pi ratios
# length(which(parent_sequence != child_sequence)) # HAMMING DISTANCE
# stringdistmatrix(aligned[1], aligned[2])
