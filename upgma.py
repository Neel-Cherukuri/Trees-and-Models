import random

def findSmallest(dM):
    # Find the coordinates of the shortest distance between organisms
    min_value = float('inf')
    min_positions = []

    for i in range(len(dM)):
        for j in range(len(dM[i])):
            if i != j and dM[i][j] < min_value:
                min_value = dM[i][j]
                min_positions = [(i, j)]

    return random.choice(min_positions) if min_positions else None, min_value

def updateMatrix(dM, row, col):
    # Determine which index to keep (the smaller one) and which to remove
    keep, remove = min(row, col), max(row, col)

    # Update the matrix with new distances
    for i in range(len(dM)):
        if i != keep and i != remove:
            dM[i][keep] = dM[keep][i] = (dM[i][keep] + dM[i][remove]) / 2

    # Remove the column and row of the merged cluster
    for i in range(len(dM)):
        del dM[i][remove]
    del dM[remove]

    return dM

def updateSpecies(sp, r, c, dist):
    # Updates the species list with distance
    merged_species = "({}, {}: {})".format(sp[r], sp[c], dist)
    sp[r] = merged_species
    del sp[c]
    return sp

def UPGMA(dM, sp):
    while len(dM) > 1:
        (leastRow, leastCol), min_value = findSmallest(dM)  # Finds the smallest non-zero distance and its value
        dM = updateMatrix(dM, leastRow, leastCol)  # Updates the distance matrix
        sp = updateSpecies(sp, leastRow, leastCol, min_value)  # Updates the species list with distance
    return sp[0]  # Return the final phylogenetic tree

# Example Data
distanceMatrix = [
    [0, 12, 12, 13, 15, 15],
    [12, 0, 2, 6, 8, 8],
    [12, 2, 0, 6, 9, 9],
    [13, 6, 6, 0, 8, 8],
    [15, 8, 9, 8, 0, 4],
    [15, 8, 9, 8, 4, 0]
]

speciesList = ["M_Spacii", "T_Pain", "G_Unit", "Q_Doba", "R_Mani", "A_Finch"]
