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
Distr.py
Classes Distribution and Density.

Note: Creating a Density instance requires rpy (not rpy2), which can be obtained
    here: https://sourceforge.net/projects/rpy/files/rpy/
"""

from __future__ import division
from __future__ import print_function
import itertools
import bisect

class Distribution(object) :
    """
    Class for representing the cumulative distribution of a set of numbers.
    """
    def __init__(self, numbersIt, weightsIt = None, relTol = None, absTol = None) :
        # Create a Distribution from the iterable of numbers and weights.
        # If absTol or relTol is not None, numbers within absTol (if not None) or
        #   relTol * (max - min), may be combined.
        #   relTol is ignored if absTol is not None.
        #   Reasonable positive value for relTol might be 0.0001
        #   If absTol < 0 (or relTol < 0 and absTol == None) no numbers will be combined.
        _set_Distribution_fields(self, numbersIt, weightsIt, relTol, absTol)

    def add(self, numbersIt, weightsIt = None) :
        # Add more numbers and weights. Weights must be commensurate with original.
        # Adding one number at a time will be slow.
        _add(self, numbersIt, weightsIt)
    
    def cum_prob(self, upTo, strictlyLess = False) :
        # Return the total probability < upTo if strictlyLess, otherwise <= upTo
        # Numbers within 1e-9 are considered equal to account for floating point errors.
        # If upTo is between two input numbers within absTol, results are fuzzy.
        return _cum_prob(self, upTo, strictlyLess)

    def dump(self, fileName) :
        # Save the distribution in a file with the specified name
        _dump(self, fileName)
    
    @staticmethod
    def load(fileName) :
        # Return a Distribution read from the file with specified name
        return _load(fileName)

    @staticmethod
    def from_text_file(fileName, **keywordArgs) :
        # Return a Distribution from the numbers in a text file (white-space delimited).
        # Keyword args can include relTol and absTol for Distribution constructor.
        textFile = open(abspath(expanduser(fileName)))
        return Distribution(map(float, textFile.read().split()), **keywordArgs)

class Density(object) :
    """
    Class for representing the PDF (probability density function) for a distribution.
    Warning: For performance reasons, Distribution approximates a set of close points
        as a single point with a weight. When this is passed to the R density function
        it estimates a different bandwidth from what it would get if it were given
        the original points. Difference is not due to the slight shifting of the points;
        rather it is due to representing several equal points as one point with a weight.
        The weights representation seem to get a bandwidth roughly twice the original
        in cases I tried. So it might be appropriate to divide "adjust" by 2.
    """
    def __init__(self, distr, bandwidth = None, adjust = None) :
        """
        Use the R density function to interpolate the PDF from the distribution.
        The bandwidth and adjust args are the bw and adjust args of the R function.
        """
        _set_density_fields(self, distr, bandwidth, adjust)
        
    def eval(self, x) :
        # Return the density at x. Evaluate by linear interpolation.
        return _density_eval(self, x)
    
    __call__ = eval
    
    @staticmethod
    def from_text_file(fileName, **keywordArgs) :
        # Return a Density from the numbers in a text file (white-space delimited).
        # keywordArgs are for Density and Distribution constructors.
        densKeywords = ['bandwidth', 'adjust']
        densKeywordArgs = dict(item for item in keywordArgs.items()
                                    if item[0] in densKeywords)
        distrKeywordArgs = dict(item for item in keywordArgs.items()
                                     if item[0] not in densKeywords)
        distr = Distribution.from_text_file(fileName, **distrKeywordArgs)
        return Density(distr, **densKeywordArgs)
    

#-------------------------------------------#
# Implementation of Distribution:

"""
Fields of Distribution:
    absTol
    numbers
    cumProbs
    totalWeight # Remember for use in add()
"""

import cPickle, sys
from os.path import abspath, expanduser

def _set_Distribution_fields(distr, numbersIt, weightsIt, relTol, absTol) :
    if weightsIt == None :
        weightsIt = itertools.repeat(1)
    numbers = list(numbersIt)
    weights = list(itertools.islice(weightsIt, len(numbers)))

    # for unknown reason, it crashes in sorted sometimes but not others
    # so if it fails try, try, again
    for attempt in range(1, 11) :
        try :
            sortedInds = sorted(range(len(numbers)), key = lambda ind : numbers[ind])
        except SystemError :
            print('*** SystemError in Distr.py. Trying again. ***', file = sys.stderr)
        else :
            break
    else :
        print('*** SystemError in Distr.py. Giving up. ***', file = sys.stderr)

    if absTol == None :
        if relTol != None :
            absTol = relTol * (numbers[sortedInds[-1]] - numbers[sortedInds[0]])
        else :
            absTol = -1
    distr.absTol = absTol
    distr.numbers = []
    distr.cumProbs = []
    distr.totalWeight = sum(weights) # Remember for use in add()
    if len(numbers) == 0 :
        return
    newWts = []
    for ind in sortedInds :
        num = numbers[ind]
        wt = weights[ind]
        if len(distr.numbers) == 0 or \
           num >= distr.numbers[-1] + absTol : # Always true if absTol < 0
            distr.numbers.append(num)
            newWts.append(wt)
        else :
            # Absorb new number and weight into previous cluster without changing mean
            newWt = newWts[-1] + wt
            distr.numbers[-1] = (distr.numbers[-1] * newWts[-1] + num * wt) / newWt
            newWts[-1] = newWt
    weightSoFar = 0
    for wt in newWts :
        weightSoFar += wt
        distr.cumProbs.append(weightSoFar / distr.totalWeight)
        
def _add(distr, numbersIt, weightsIt) :
    _set_Distribution_fields(distr,
        itertools.chain(distr.numbers, numbersIt),
        itertools.chain([distr.cumProbs[0] * distr.totalWeight], 
                        [(distr.cumProbs[ii + 1] - distr.cumProbs[ii]) * distr.totalWeight
                         for ii in range(len(distr.cumProbs) - 1)],
                        weightsIt),
        relTol = None,
        absTol = distr.absTol)

def _cum_prob(distr, upTo, strictlyLess = False) :
    # Return the total probability < upTo if strictlyLess else <= upTo
    # Numbers within 1e-9 are considered equal to account for floating point errors.
    # If upTo is between two input numbers within absTol, larger one is included.
    index = bisect.bisect(distr.numbers, upTo - 1e-9 if strictlyLess else upTo + 1e-9)
    return distr.cumProbs[index - 1] if index > 0 else 0

def _dump(distr, fileName) :
    cPickle.dump(distr, open(abspath(expanduser(fileName)), 'wb'), protocol = 2)

def _load(fileName) :
    return cPickle.load(open(abspath(expanduser(fileName)), 'rb'))
                
#-------------------------------------------#
# Implementation of Density:

"""
Fields of Density:
    xs
    ys
"""

def _set_density_fields(dens, distr, bandwidth, adjust) :
    import rpy
    weights = [distr.cumProbs[0] if ii == 0 else
               distr.cumProbs[ii] - distr.cumProbs[ii - 1]
               for ii in range(len(distr.cumProbs))]
    argList = [('', distr.numbers), ('weights', weights)]
    if bandwidth != None :
        argList.append(('bw', bandwidth))
    if adjust != None :
        argList.append(('adjust', adjust))
    rpyDensity = rpy.r.density.lcall(argList)
    dens.xs = rpyDensity['x']
    dens.ys = rpyDensity['y']

def _density_eval(dens, x) :
    # Evaluate by linear interpolation
    xs, ys = dens.xs, dens.ys
    index = bisect.bisect(xs, x)
    if index < 1 or index >= len(xs) :
        return 0
    return ys[index - 1] + \
           (ys[index] - ys[index - 1]) * (x - xs[index - 1]) / (xs[index] - xs[index - 1])

        