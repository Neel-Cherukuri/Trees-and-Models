import math

# List of codons
modelCodons = ['TTT', 'TTC', 'TTA', 'TTG', 'CTT', 'CTC', 'CTA',
               'CTG', 'ATT', 'ATC', 'ATA', 'ATG', 'GTT', 'GTC',
               'GTA', 'GTG', 'TCT', 'TCC', 'TCA', 'TCG', 'AGT',
               'AGC', 'CCT', 'CCC', 'CCA', 'CCG', 'ACT', 'ACC',
               'ACA', 'ACG', 'GCT', 'GCC', 'GCA', 'GCG', 'TAT',
               'TAC', 'CAT', 'CAC', 'CAA', 'CAG', 'AAT', 'AAC',
               'AAA', 'AAG', 'GAT', 'GAC', 'GAA', 'GAG', 'TGT',
               'TGC', 'CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG',
               'GGT', 'GGC', 'GGA', 'GGG', 'TGG', 'TAA', 'TAG',
               'TGA']

def scoreModels():
    codingMatrix = getProbs("./codingModel.tab")  # 64x64 matrix of probabilities for coding triplet mutations
    noncodingMatrix = getProbs("./noncodingModel.tab")  # 64x64 matrix of probabilities for non-coding triplet mutations
    id2ancestorSeq = getSeq("./Ancestor.fa")  # reads the ancestor sequences
    id2spaciiSeq = getSeq("./Spacii.fa")  # reads the spacii sequences
    allID = list(id2ancestorSeq.keys())  # the two sequences above share dictionary indexes

    for ID in allID:
        cScore = 0  # variable to contain the score with the coding model
        nScore = 0  # variable to contain the score with the non-coding model

        ancestorSeq = id2ancestorSeq[ID]
        spaciiSeq = id2spaciiSeq[ID]

        # Iterate over codons in the sequences
        for i in range(0, len(ancestorSeq), 3):
            ancestorCodon = ancestorSeq[i:i + 3]
            spaciiCodon = spaciiSeq[i:i + 3]

            if len(ancestorCodon) == 3 and len(spaciiCodon) == 3:
                ancestorIndex = modelCodons.index(ancestorCodon)
                spaciiIndex = modelCodons.index(spaciiCodon)

                # Update scores using log probabilities
                cScore += math.log(codingMatrix[ancestorIndex][spaciiIndex])
                nScore += math.log(noncodingMatrix[ancestorIndex][spaciiIndex])

        # Print the result for each ID in the specified format
        if cScore > nScore:
            print(
                f"Sequence '{ID}'is classified as likely coding with a coding score of {cScore:.2f} and a non-coding score of {nScore:.2f}.")
        else:
            print(
                f"Sequence '{ID}'is classified as likely non-coding with a coding score of {cScore:.2f} and a non-coding score of {nScore:.2f}.")

def getProbs(f1):
    f = open(f1)
    pMatrix = []
    for line in f:
        tmp = line.rstrip().split("\t")
        tmp = [float(i) for i in tmp]
        pMatrix.append(tmp)
    f.close()
    return pMatrix

def getSeq(filename):
    f = open(filename)
    id2seq = {}
    currkey = ""
    for line in f:
        if line.find(">") == 0:
            currkey = (line[1:].split("|")[0])
            id2seq[currkey] = ""
        else:
            id2seq[currkey] += line.rstrip()
    f.close()
    return id2seq

# Call the function to perform the analysis
scoreModels()
