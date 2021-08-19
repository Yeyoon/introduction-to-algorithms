#!/usr/bin/env python2.7

import unittest
from dnaseqlib import *

### Utility classes ###

# Maps integer keys to a set of arbitrary values.
class Multidict:
    # Initializes a new multi-value dictionary, and adds any key-value
    # 2-tuples in the iterable sequence pairs to the data structure.
    def __init__(self, pairs=[]):
        self.data = {}
        for p in pairs:
            k,v = p[0],p[1]
            self.put(k,v)
        
    # Associates the value v with the key k.
    def put(self, k, v):
        if k in self.data:
            self.data[k].append(v)
        else:
            self.data[k] = [v]
            
    # Gets any values that have been associated with the key k; or, if
    # none have been, returns an empty sequence.
    def get(self, k):
        if k in self.data:
            return self.data[k]
        return []

def hash_func1(seq):
    v = 0
    for c in seq:
        v = v + ord(c)

    return v % len(seq)

# Given a sequence of nucleotides, return all k-length subsequences
# and their hashes.  (What else do you need to know about each
# subsequence?)
def subsequenceHashes(seq, k):
    base_seq = ''
    tmp = []
    print("seq,k -> ",seq,k)
    for i in range(k):
        try:
            tmp.append(next(seq))
        except StopIteration:
            return []

    base_seq = ''.join(tmp)

    h = RollingHash(base_seq)
    print("first yield :",base_seq)
    yield (base_seq,h.current_hash(),0)

    i = 0
    for c in seq:
        previtm = tmp[0]
        nextitm = c
        tmp = tmp[1:]
        tmp.append(c)
        base_seq = ''.join(tmp)
        i += 1
        print("i is ",i,base_seq)
        yield (base_seq,h.slide(previtm,nextitm),i)



# Similar to subsequenceHashes(), but returns one k-length subsequence
# every m nucleotides.  (This will be useful when you try to use two
# whole data files.)
def intervalSubsequenceHashes(seq, k, m):
    raise Exception("Not implemented!")

# Searches for commonalities between sequences a and b by comparing
# subsequences of length k.  The sequences a and b should be iterators
# that return nucleotides.  The table is built by computing one hash
# every m nucleotides (for m >= k).
def getExactSubmatches(a, b, k, m):
    mh = Multidict()
    for s in subsequenceHashes(a,k):
        mh.put(s[1],s)

    for s in subsequenceHashes(b,k):
        values = mh.get(s[1])
        for v in values:
            if v[0] == s[0]:
                yield (v[2],s[2])
                

if __name__ == '__main__':
    if len(sys.argv) != 4:
        print('Usage: {0} [file_a.fa] [file_b.fa] [output.png]'.format(sys.argv[0]))
        sys.exit(1)

    # The arguments are, in order: 1) Your getExactSubmatches
    # function, 2) the filename to which the image should be written,
    # 3) a tuple giving the width and height of the image, 4) the
    # filename of sequence A, 5) the filename of sequence B, 6) k, the
    # subsequence size, and 7) m, the sampling interval for sequence
    # A.
    compareSequences(getExactSubmatches, sys.argv[3], (500,500), sys.argv[1], sys.argv[2], 8, 100)
