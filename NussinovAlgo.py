# -*- coding: utf-8 -*-
"""
Created on Tue Apr 19 02:57:31 2022

@author: Ahmed Nasser
"""
import numpy as np


def checkMatch(x, y):
    matchGC = [('G', 'C'), ('C', 'G')]
    matchAU = [('A', 'U'), ('U', 'A')]
    matchUG = [('U', 'G'), ('G', 'U')]

    if (x, y) in matchGC: return scoreGC
    elif (x, y) in matchAU: return scoreAU
    elif wobblePair[0] == 'Y' and (x, y) in matchUG: return 1
    else: return 0
    return 0


def maxBifurcate(i, j):
    localMax = idx = 0
    for k in range(i+1, j):
        nav = Mat[i][k] + Mat[k+1][j]
        if nav >= localMax:
            localMax = nav
            idx = k
    return localMax, idx


def fillMatrix():
    for j in range(1, N):   #col
        for i in range(N-j):   #row
            left = Mat[i][(i+j)-1]
            down = Mat[i+1][i+j]
            diag = Mat[i+1][(i+j)-1] + checkMatch(rnaDict[i+1], rnaDict[i+j+1])
            bifurc = maxBifurcate(i, i+j)
            Mat[i][i+j] = max(down, left, diag, bifurc[0])
    return


def traceBack(i, j):
    if i < j:
        if Mat[i][j] == Mat[i][j-1]: traceBack(i, j-1)
        elif Mat[i][j] == Mat[i+1][j]: traceBack(i+1, j)
        elif Mat[i][j] == Mat[i+1][j-1] + checkMatch(rnaSeq[i], rnaSeq[j]):
            dotBracket[min(i, j)] = '('
            dotBracket[max(i, j)] = ')'
            traceBack(i+1, j-1)
        else:
            kIdx = maxBifurcate(i, j)[1]
            traceBack(i, kIdx)
            traceBack(kIdx+1, j)
    return



''' Test Cases:-
                        (1)                     (2)                 (3)
                ===================     ===================     =============
                | CGGACCCAGACUUUC |     | UAACGUACUGGAGUA |     | GGGAAAUCC |
    Wobble:     |=================|     |=================|     |===========|
        (Y)=>   | .((.)).(((.))). |     | ()(()).(()).(). |     | .((..())) |
        (N)=>   | ()((...((.))).) |     | ()(())(((..))). |     | .((..())) |
                ===================     ===================     =============
'''
if __name__ == "__main__":
    #Input-> ..Step
    rnaSeq = input("Please Enter the RNA Sequence.. ").upper()
    wobblePair = input("Wobble Pair allowed? (Y/N).. ").upper()

    #Scoring-> ..
    scoreGC = 1
    scoreAU = 1

    #Initialization-> .. Step
    N = len(rnaSeq)
    Mat = np.zeros((N, N))   #Initialize Matrix with ZERO
    dotBracket = ['.'] * N   #Set FoldList as dotBracket[N]
    rnaDict = {}   #Dictionary of Indexes

    for idx in range(1, N+1):   #Enumerate rnaSeq into Dictionary[0:N-1]
        rnaDict[idx] = rnaSeq[idx-1]

    #Processing-> ..Step
    fillMatrix()   #Build the Scoring Matrix
    traceBack(0, N-1)   #Traceback from the Top-Right cell

    print('\n', Mat, '\n\n\t', rnaSeq, '\n\t', ''.join(dotBracket))
