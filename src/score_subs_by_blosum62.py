#! /usr/bin/python3

import sys
import argparse
import numpy as np
from parse_fasta import parse_fasta


BLOSUM62 = """4 -1 -2 -2 0 -1 -1 0 -2 -1 -1 -1 -1 -2 -1 1 0 -3 -2 0
-1 5 0 -2 -3 1 0 -2 0 -3 -2 2 -1 -3 -2 -1 -1 -3 -2 -3
-2 0 6 1 -3 0 0 0 1 -3 -3 0 -2 -3 -2 1 0 -4 -2 -3
-2 -2 1 6 -3 0 2 -1 -1 -3 -4 -1 -3 -3 -1 0 -1 -4 -3 -3
0 -3 -3 -3 9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1
-1 1 0 0 -3 5 2 -2 0 -3 -2 1 0 -3 -1 0 -1 -2 -1 -2
-1 0 0 2 -4 2 5 -2 0 -3 -3 1 -2 -3 -1 0 -1 -3 -2 -2
0 -2 0 -1 -3 -2 -2 6 -2 -4 -4 -2 -3 -3 -2 0 -2 -2 -3 -3
-2 0 1 -1 -3 0 0 -2 8 -3 -3 -1 -2 -1 -2 -1 -2 -2 2 -3
-1 -3 -3 -3 -1 -3 -3 -4 -3 4 2 -3 1 0 -3 -2 -1 -3 -1 3
-1 -2 -3 -4 -1 -2 -3 -4 -3 2 4 -2 2 0 -3 -2 -1 -2 -1 1
-1 2 0 -1 -3 1 1 -2 -1 -3 -2 5 -1 -3 -1 0 -1 -3 -2 -2
-1 -1 -2 -3 -1 0 -2 -3 -2 1 2 -1 5 0 -2 -1 -1 -1 -1 1
-2 -3 -3 -3 -2 -3 -3 -3 -1 0 0 -3 0 6 -4 -2 -2 1 3 -1
-1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4 7 -1 -1 -4 -3 -2
1 -1 1 0 -1 0 0 0 -1 -2 -2 0 -1 -2 -1 4 1 -3 -2 -2
0 -1 0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1 1 5 -2 -2 0
-3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1 1 -4 -3 -2 11 2 -3
-2 -2 -2 -3 -2 -1 -2 -3 2 -1 -1 -2 -1 3 -3 -2 -2 2 7 -1
0 -3 -3 -3 -1 -2 -2 -3 -3 3 1 -2 1 -1 -2 -2 0 -3 -1 4"""


class ScoringMatrix():
    def __init__(self, scores=BLOSUM62):
        self.scoreMatrix = np.array([[int(i) for i in row.split(" ")] for
                                     row in scores.splitlines()])
        self.states = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L",
                       "K", "M", "F", "P", "S", "T", "W", "Y", "V"]

    def getScore(self, state1, state2):
        ind1 = self.states.index(state1)
        ind2 = self.states.index(state2)
        score = self.scoreMatrix[ind1, ind2]
        return score

    def getScoreList(self, subList):
        scores = {}
        for i in subList:
            pos = int(i[1:-1])
            scores[pos] = self.getScore(i[0], i[-1])
        return scores


if __name__ == "__main__":
    if len(sys.argv[1:]) == 0:
        sys.argv.append("-h")

    parser = argparse.ArgumentParser()
    parser.add_argument("subs", help="space-separated list of substitutions, \
                                      formatted A263V",
                        nargs="+")
    args = parser.parse_args()

    b62 = ScoringMatrix()
    possScores = [11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0, -1, -2, -3, -4]
    hexPal = ["#0000FF", "#1717FF", "#2E2EFF", "#4545FF", "#5C5CFF", "#7373FF",
              "#8B8BFF", "#A2A2FF", "#B9B9FF", "#D0D0FF" "#E7E7FF" "#FFFFFF",
              "#FFBFBF", "#FF7F7F", "#FF3F3F", "#FF0000"]
    rgbPal = [(0, 0, 255), (23, 23, 255), (46, 46, 255), (69, 69, 255),
              (92, 92, 255), (115, 115, 255), (139, 139, 255), (162, 162, 255),
              (185, 185, 255), (208, 208, 255), (231, 231, 255),
              (255, 255, 255), (255, 191, 191), (255, 127, 127), (255, 63, 63),
              (255, 0, 0)]
    scores = b62.getScoreList(args.subs)
    print("pos\tscore\trgb")
    for k, v in scores.items():
        print(str(k) + "\t" + str(v) + "\t" + str(rgbPal[possScores.index(v)]))
