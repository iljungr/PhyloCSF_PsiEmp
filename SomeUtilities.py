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
SomeUtilities.py
"""
from __future__ import division
import sys, os, shutil

class MyException(Exception) :
    "Use for any error that shouldn't be handled, but can print OS-level string."
    pass

def err_msg(message) :
    print >>sys.stderr, message
    sys.stderr.flush()

def get_absolute_path(path) :
    "Convert ., .., and ~"
    return os.path.abspath(os.path.expanduser(path))

pjoin = os.path.join

def myopen(path, *args) :
    path = get_absolute_path(path)
    if not file_exists(path) and file_exists(path + '.gz') :
        path = path + '.gz'
    if path[-3:] == '.gz' :
        import gzip
        return gzip.open(path, *args)
    elif path[-4:] == '.bgz' :
        import bgzf # Make sure sys.path has appropriate dir
        return bgzf.open(path, *args)
    else :
        return open(path, *args)

def file_exists(fileName) :
    return os.path.exists(get_absolute_path(fileName))

def is_directory(fileName) :
    return os.path.isdir(get_absolute_path(fileName))

isdir = is_directory

def mv(source, destination, dump = False, onlyDump = False) :
    if dump or onlyDump :
        err_msg('mv %s %s' % (source, destination))
    if onlyDump :
        return
    os.rename(get_absolute_path(source), get_absolute_path(destination))

def cp(source, destination, dump = False, onlyDump = False) :
    if dump or onlyDump :
        err_msg('cp %s %s' % (source, destination))
    if onlyDump :
        return
    shutil.copy(get_absolute_path(source), get_absolute_path(destination))

def rmdir(path, dump = False, onlyDump = False) :
    if dump or onlyDump :
        err_msg('rmdir %s' % path)
    if onlyDump :
        return
    os.rmdir(get_absolute_path(path))
    
def rm_minus_r(path, dump = False, onlyDump = False) :
    "Like rm -rf path"
    if dump or onlyDump :
        err_msg('rm -rf %s' % path)
    if onlyDump :
        return
    shutil.rmtree(path)
    
def rm(path, dump = False, onlyDump = False) :
    if dump or onlyDump :
        err_msg('rm %s' % path)
    if onlyDump :
        return
    os.remove(get_absolute_path(path))

def mkdir(path, dump = False, onlyDump = False) :
    if dump or onlyDump :
        err_msg('mkdir ' + path)
    if onlyDump :
        return
    return os.mkdir(get_absolute_path(path))

def assure_dir(path, dump = False) :
    if not is_directory(path) :
        mkdirs(path, dump)

def mkdirs(path, dump = False):
    "Create a directory and all parent directories."
    path = get_absolute_path(path)
    if os.path.isdir(path):
        pass
    elif os.path.isfile(path):
        raise OSError("cannot create directory, file already exists: '%s'" % path)
    else:
        head, tail = os.path.split(path)
        if head and not os.path.isdir(head):
            mkdirs(head, dump)
        if tail:
            mkdir(path, dump)

def same_file(path1, path2) :
    "Return True if the two paths point to the same file."
    return os.path.samefile(get_absolute_path(path1), get_absolute_path(path2))

def parse_arg(argName, shortArgName = None) :
    # If --argName is in sys.argv, return True and remove it.
    # If shortArgName (which should be a single character) is in an argument of the 
    #   form -string, where string doesn't start with -, then return True and remove it.
    # Otherwise return False
    if '--' + argName in sys.argv[1:] :
        del sys.argv[sys.argv.index('--' + argName, 1)]
        return True
    elif shortArgName != None :
        for pos, arg in enumerate(sys.argv) :
            if pos == 0 or arg[0] != '-' or (len(arg) > 1 and arg[1] == '-') :
                continue
            if shortArgName in arg[1:] :
                sys.argv[pos] = sys.argv[pos].replace(shortArgName, '')
                if sys.argv[pos] == '-' :
                    del sys.argv[pos]
                return True
    return False
        
def parse_arg_value(argName, argType = str, default = None) :
    # If sys.argv includes a string of the form --argName=value, 
    #   return the value, casted to type argType, and remove it from sys.argv.
    # Otherwise return default.
    for pos, arg in enumerate(sys.argv) :
        if pos == 0 or not arg.startswith('--' + argName + '=') :
            continue
        del sys.argv[pos]
        return argType(arg[len(argName) + 3 :])
    return default
                    
def without_trailing_digits(word) :
    while word != '' and word[-1].isdigit() :
        word = word[: -1]
    return word
remove_trailing_digits = without_trailing_digits
    
def without_trailing_char(string, char) :
    if string != '' and string[-1] == char :
        return string[:-1]
    else :
        return string
        
def without_trailing_nl(line) :
    return without_trailing_char(line, '\n')
    
remove_trailing_nl = without_trailing_nl
strip_nl = without_trailing_nl

def neighbors(iter) :
    # Return a list of the consecutive pairs of elements in the iterator.
    return list(neighbors_iter(iter))

def neighbors_iter(iter) :
    # Generate all consecutive pairs of elements in the iterator.
    firstTime = True
    for item in iter :
        if firstTime :
            firstTime = False
        else :
            yield (prev, item)
        prev = item

def antistrand(strand) :
    if strand == '+' :
        return '-'
    elif strand == '-' :
        return '+'
    else :
        return strand

dnaComplementDict = {'G':'C', 'C':'G', 'A':'T', 'T':'A',
                     'g':'c', 'c':'g', 'a':'t', 't':'a'}
rnaComplementDict = {'G':'C', 'C':'G', 'A':'U', 'U':'A',
                     'g':'c', 'c':'g', 'a':'u', 'u':'a'}

def reverse_complement(string, rna = False) :
    """
    Reverse complement the string, complementing upper or lower case a, c, g, and t (or u 
        if rna is True), preserving case, and leaving all other characters the same, e.g., 
        'N' or '-'.
    Fail if string contains t or T if rna or if it contains u or U and not rna.
    """
    if rna :
        assert('t' not in string and 'T' not in string), string
        complementDict = rnaComplementDict
    else :
        assert('u' not in string and 'U' not in string), string
        complementDict = dnaComplementDict
    return ''.join(complementDict.get(char, char) for char in string[::-1])
    
def stop(status = 1):
    sys.exit(status)
    
def regionString_to_triples(regionStr) :
    # Parse a string of form chrom:START-END[+chrom:START-END]*, 
    #    where the numbers can include commas.
    # Return [(chrom, START, END)]
    triples = []
    for regionStr in regionStr.split('+') :
        split1 = regionStr.split(':')
        split2 = split1[1].replace(',', '').split('-')
        triples.append((split1[0], int(split2[0]), int(split2[1])))
    return triples

def get_intervals_length(intervals, includesChrom = True) :
    # Return number of bases in a list of intervals: [(CHROM, START, END)...] 
    # (or [(START, END)...] if not includesChrom)
    if includesChrom :
        startIndex = 1
    else :
        startIndex = 0
    return sum(interval[startIndex + 1] - interval[startIndex] + 1
               for interval in intervals)
    
def intervals_prefix(intervals, strand, numBases, includesChrom = True) :
    # Return subset of intervals numBases long starting at beginning relative to strand.
    # Intervals is a list: [(CHROM, START, END)]
    # (or [(START, END)...] if not includesChrom)
    if numBases < 0 or numBases > get_intervals_length(intervals, includesChrom) :
        raise ValueError, ('intervals_prefix: invalid numBases (%s) for length %s' %
                            (numBases, get_intervals_length(intervals, includesChrom)))
    result = []
    reverse = strand == '-'
    startInd = [0, 1][includesChrom]
    for interval in intervals[::[1, -1][reverse]] :
        intervalLen = interval[startInd + 1] - interval[startInd] + 1
        if numBases >= intervalLen :
            result.insert([len(result), 0][reverse], interval)
            numBases -= intervalLen
            continue
        if numBases > 0 :
            if reverse :
                result.insert(0, interval[:startInd] +
                                 (interval[startInd + 1] - numBases + 1,
                                  interval[startInd + 1]))
            else :
                result.insert(len(result), interval[:startInd] +
                                           (interval[startInd],
                                            interval[startInd] + numBases - 1))
        return result
    return result

def intervals_suffix(intervals, strand, numBases, includesChrom = True) :
    # Return subset of intervals numBases long ending at end relative to strand.
    # Intervals is a list: [(CHROM, START, END)]
    # (or [(START, END)...] if not includesChrom)
    return intervals_prefix(intervals, ['-', '+'][strand == '-'], numBases, includesChrom)

def bed_line_to_intervals(line) :
    """
    Input a line in a BED format file.
    Return (name, chrom, intervals, strand),
       where intervals = [[start1, end1], [start2, end2], ...] in the input order, which
       should satisfy start1 <= end1 < start2 <= end2 < ...
    Note that the resulting intervals treat chromosomes as starting at position 1,
       and include both endpoints in the interval (as is done in gtf/gff files,
       the genome browser, and PhyloCSF), rather than treating chromosomes as starting at
       position 0 and ending beyond the end of the interval, as is done in input BED file.
       In other words, the starts will differ by 1 from the numbers in the BED file. 
       For example, an interval containing just the first base of a chromosome would be 
       0, 1 in the BED file, and 1, 1 here. (End is ignored for BED12, but not for BED6.)
    For Bed6 files, return one interval for the whole thing.
    """
    words = map(str.strip, line.split('\t')) # Remove trailing \n, as well as spaces
    chrom = words[0]
    firstStart = int(words[1]) + 1
    name = words[3]
    strand = words[5]
    if len(words) > 11 : # BED12
        numExons = int(words[9])
        exonSizes = map(int, without_trailing_char(words[10], ',').split(','))
        exonRelStarts = map(int, without_trailing_char(words[11], ',').split(','))
        intervals = [[firstStart + relStart, firstStart + relStart + size - 1]
                     for size, relStart in zip(exonSizes, exonRelStarts)]
        # Checks:
        if len(exonSizes) != numExons or len(exonRelStarts) != numExons :
            raise AssertionError, 'Mismatched counts: ' + line
        if int(words[2]) != intervals[-1][1] :
            raise AssertionError, 'Mismatched ends: ' + line
    elif len(words) <= 9 : # BED6
        firstEnd = int(words[2])
        intervals = [[firstStart, firstEnd]]
    else :
        raise AssertionError, 'Line with invalid number of fields: ' + line   
    return name, chrom, intervals, strand

def is_bed_comment(line) :
    # Return true if line from bed file is a comment or browser metadata 
    return line == '' or line[0] == '#' or line.split()[0] in ['browser', 'track']

def csv_to_columns(file) :
    # Read a comma-separated-value file.
    # Return a dictionary and a list: {columnName:[values]}, [columnNames in order].
    # columnNames and values will be stripped of leading and trailing spaces.
    valuesDict = None
    for line in file :
        if line == '' or line == '\n' or line[0] == '#' :
            continue # Skip comments and blank lines
        words = [word.strip() for word in line.split(',')] # strip removes trailing \n
        if valuesDict == None : # Header row
            valuesDict = dict((word, []) for word in words)
            columnNames = words
            numCols = len(words)
        else :
            numWords = len(words)
            if numWords < numCols :
                words += '' * (numCols - numWords)
            elif numWords > numCols :
                err_msg('read_csv: ignoring extra words.')
                del words[numCols :]
            for word, columnName in zip(words, columnNames) :
                valuesDict[columnName].append(word)
    return valuesDict, columnNames

def csv_to_table(file, indexColumnName) :
    # Read a comma-separated-value file.
    # Return a Table, where row names are taken from the column with name indexColumnName.
    # columnNames and values will be stripped of leading and trailing spaces.
    table = None
    for line in file :
        if line == '' or line == '\n' or line[0] == '#' :
            continue # Skip comments and blank lines
        words = [word.strip() for word in line.split(',')] # strip removes trailing \n
        if table == None : # Header row
            indexColumn = words.index(indexColumnName)
            table = Table(words[: indexColumn] + words[indexColumn + 1 :])
            numCols = len(words)
        else :
            numWords = len(words)
            if numWords < numCols :
                words += '' * (numCols - numWords)
            elif numWords > numCols :
                raise ValueError, 'csv_to_table: too many words (%d vs %d).' % (numWords,
                                                                                numCols)
            table.add_row(words[indexColumn],
                          words[: indexColumn] + words[indexColumn + 1 :])
    return table
    
def table_to_csv(file, table, indexColumnName = 'name') :
    # Print the Table as a comma-separated file (with index column first)
    columnNames = table.column_names()
    print >>file, ','.join([indexColumnName] + list(columnNames))
    for rowName in table.row_names() :
        print >>file, ','.join([rowName] + [table[rowName][columnName]
                                            for columnName in columnNames])

class Table(object) :
    """
    Class for storing a table with (ordered) row and column names so that a cell can be
        accessed as atable[rowName][columnName].
    Only access by methods -- implementation is subject to change.
    """
    def __init__(self, columnNames) :
        self.colNames = tuple(columnNames)
        self.colDict = dict((colName, colInd)
                            for colInd, colName in enumerate(self.colNames))
        self.rowNames = []
        self.rowDict = {} # {rowName : orderedValues, ...}
    def add_row(self, rowName, orderedValues) :
        # Values must be in order of column names
        self.rowNames.append(rowName)
        self.rowDict[rowName] = tuple(orderedValues)
    def delete_row(self, rowName) :
        self.rowNames.remove(rowName)
        del self.rowDict[rowName]
    def row_names(self) :
        return self.rowNames
    def column_names(self) :
        return self.colNames
    def __getitem__(self, rowName) :
        return _TableRow(self, self.rowDict[rowName])
        
class _TableRow(object) :
    # For internal use by Table. Implementation subject to change.
    def __init__(self, table, row) :
        self.table = table
        self.row = row
    def __getitem__(self, colName) :
        return self.row[self.table.colDict[colName]]

def split(strToSplit, sep = None) :
    """Equivalent to strToSplit.split(sep) except that split('') -> [] instead of [''].
       This makes it the inverse of join: split(sep.join(listOfStrs), sep) -> listOfStrs.
    """
    return strToSplit.split(sep) if strToSplit != '' else []

def split_str(strToSplit, maxLength = 70) :
    """
    Split string into lines of length at most maxLength. Useful for writing fasta files.
    If maxLength is None, just return the string without splitting it.
    """
    if maxLength == None :
        yield strToSplit
        return
    start = 0
    while start < len(strToSplit) - maxLength :
        yield strToSplit[start : start + maxLength]
        start += maxLength
    yield strToSplit[start : ]
    
def write_to_fasta(outFile, header, sequence, maxLength = 70) :
    print >>outFile, '>' + header
    for line in split_str(sequence, maxLength) :
        print >>outFile, line

def iter_fasta(faFile) :
    # Iterate through a fasta file returning pairs: sequence name, sequence
    seqName = None
    seq = ''
    for line in faFile :
        line = strip_nl(line)
        if len(line) == 0 :
            continue
        if line[0] == '>' :
            if seqName != None :
                yield seqName, seq
            seqName = line[1:]
            seq = ''
        else :
            seq += line
    if seqName != None :
        yield seqName, seq
