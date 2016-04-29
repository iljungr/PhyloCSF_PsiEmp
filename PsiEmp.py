#!/usr/bin/env python
# Copyright 2016 Irwin Jungreis
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at http://www.apache.org/licenses/LICENSE-2.0
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""
PsiEmp.py
Evaluator for PhyloCSF Psi_EMP.
"""

from __future__ import division
from __future__ import print_function
import os
from SomeUtilities import myopen
from Distr import Distribution, Density
from math import log, exp, sqrt
import cPickle

MaxLenEmpDistr = 10  # Maximum number of codons for using empirical distribution
MaxLenEmpMean  = 60  # Maximum number of codons for using empirical mean and stdev

decibanFactor = 10 / log(10)

class PsiEmpEvaluator(object) :
    "Class for computing PhyloCSF-PsiEMP. Usage: instance(rawScore, numCodons)"
    """
    PhyloCSF-PsiEMP is like PhyloCSF-Psi except that instead of length-scaled
        normal approximation to distribution it uses:
            if NumCodons <= MaxLenEmpDistr :
                empirical distribution
            elif NumCodons <= MaxLenEmpMean :
                empirical distribution of MaxLenEmpDistr codons scaled
                using empirical mean and stdev for length NumCodons
            else :
                empirical distribution of MaxLenEmpDistr codons scaled
                using best linear fit to empirical means or log-scale stdevs of regions
                with MaxLenEmpMean // 2 <= NumCodons <= MaxLenEmpMean
    """
    def __init__(self, scoreFile) :
        """
        scoreFile is a cPickle file containing raw PhyloCSF scores of coding and
            non-coding regions of regions of various lengths.
        Format: (codingScores, noncodingScores) where each is a dictionary 
            {NumCodons : [score1, score2,...], ...}. The keys must be
            range(1, MaxLenEmpMean + 1) and there must be a sufficient number of scores
            for each value of NumCodons to compute an empirical distribution.
        """
        scores = cPickle.load(myopen(scoreFile, 'rb')) # (codingScores, nonCodingScores)
        self.means          = [None, None]#0: coding, 1: noncoding (typ)
        self.stdevs         = [None, None]
        self.densities      = [None, None]
        self.meanRegress    = [None, None]#[(slope,intercept) for coding & noncoding mean]
        self.stdevRegress   = [None, None]
        self.densityToScale = [None, None]#Normalized density for NumCodons=MaxLenEmpDistr
        for ind in range(2) : # 0: coding, 1: noncoding
            assert(sorted(scores[ind].keys()) == range(1, MaxLenEmpMean + 1))
            self.means[ind]  = [0] + [mean(scores[ind][ii])
                                      for ii in range(1, MaxLenEmpMean + 1)]
            self.stdevs[ind] = [0] + [stdev(scores[ind][ii])
                                      for ii in range(1, MaxLenEmpMean + 1)]
            self.densities[ind] = [0] + [Density(Distribution(scores[ind][ii],
                                                              absTol = -1))
                                         for ii in range(1, MaxLenEmpDistr + 1)]
            self.meanRegress[ind] = linear_regression(
                [(ii, self.means[ind][ii])
                 for ii in range(MaxLenEmpMean // 2, MaxLenEmpMean + 1)])
            self.stdevRegress[ind] = linear_regression(
                [(log(ii), log(self.stdevs[ind][ii]))
                 for ii in range(MaxLenEmpMean // 2, MaxLenEmpMean + 1)])
                
            meanForScaling = self.means[ind][MaxLenEmpDistr]
            stdevForScaling = self.stdevs[ind][MaxLenEmpDistr]
            normalizedScoresForScaling = [(score - meanForScaling) / stdevForScaling
                                          for score in scores[ind][MaxLenEmpDistr]]
            self.densityToScale[ind] = Density(Distribution(normalizedScoresForScaling,
                                                            absTol = -1))

    def eval(self, rawScore, numCodons) :
        "Return PhyloCSF-PsiEMP for given PhyloCSF raw score and region length in codons."
        densities = [self._eval_density(rawScore, numCodons, ind) for ind in [0, 1]]
        
        # Handle 0 density, and evaluation beyond reliable range:
        if densities[0] >= densities[1] * 1e5 :
            return 50
        if densities[0] <= densities[1] * 1e-5 :
            return -50
        
        """
        The following should be adjusted to avoid pathalogical results when densities[1]
        is too low because we are so far from the center of the distribution. For
        example, the 4-codon region of agam4.2 chr2R:2095242-2095253+ has extremely low
        score -217.2695, yet it gets a positive psi-emp of 9.526356958337086. 
        To prevent this, we should limit the psi-emp score when the raw score is less 
        than the non-coding mean...
        """
        
        return decibanFactor * log(densities[0] / densities[1])

    __call__ = eval

    def _get_mean_stdev(self, numCodons, ind) :
        if numCodons <= MaxLenEmpMean :
            m = self.means[ind][numCodons]
            s = self.stdevs[ind][numCodons]
        else :
            m = self.meanRegress[ind][0] * numCodons + self.meanRegress[ind][1]
            s = exp(log(numCodons) * self.stdevRegress[ind][0] +
                    self.stdevRegress[ind][1])
        return m, s

    def _eval_density(self, rawScore, numCodons, ind) : # ind = 0 for coding, 1 noncoding
        if numCodons <= MaxLenEmpDistr :
            return self.densities[ind][numCodons].eval(rawScore)
        else : # Scale the distribution for MaxLenEmpDistr
            m, s = self._get_mean_stdev(numCodons, ind)
            return self.densityToScale[ind].eval((rawScore - m) / s) / s

    def dump(self, outFileName) :
        "Save this PsiEmpEvaluator as a cPickle file."
        cPickle.dump(self, myopen(outFileName, 'wb'), protocol = 2)

    @staticmethod
    def load(inFileName) :
        "Return a cPickled PsiEmpEvaluator read from the specified file name."
        psiEmpEval = cPickle.load(myopen(inFileName, 'rb'))
        return psiEmpEval
        
def mean(sampleValues) :
    return sum(sampleValues) / len(sampleValues)

def stdev(sampleValues) :
    ave = mean(sampleValues)
    return sqrt(sum((x - ave) * (x - ave) for x in sampleValues) / (len(sampleValues)-1))

def linear_regression(pairs) :
    # Given a set of (x, y) pairs, return m and b (slope and intercept)
    # that give the best fit for y = m x + b
    nn = len(pairs)
    xsum = sum([x for (x, y) in pairs])
    ysum = sum([y for (x, y) in pairs])
    xysum = sum([x * y for (x, y) in pairs])
    xsqrsum = sum([x * x for (x, y) in pairs])
    m = (nn * xysum - xsum * ysum) / (nn * xsqrsum - xsum * xsum)
    b = (ysum - m * xsum) / nn
    return m, b 
