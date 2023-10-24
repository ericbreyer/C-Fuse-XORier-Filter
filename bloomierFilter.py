from sklearn.utils import murmurhash3_32
from bitarray import bitarray
import math
from random import randint, random, seed
import csv
import pandas as pd
import sys
import matplotlib.pyplot as plt
from collections import defaultdict
import functools
import array
import time
from joblib import Parallel, delayed

seed("Baker Comes First!")

def get_size(obj, seen=None):
    """Recursively finds size of objects"""
    size = sys.getsizeof(obj)
    if seen is None:
        seen = set()
    obj_id = id(obj)
    if obj_id in seen:
        return 0
    # Important mark as seen *before* entering recursion to gracefully handle
    # self-referential objects
    seen.add(obj_id)
    if isinstance(obj, dict):
        size += sum([get_size(v, seen) for v in obj.values()])
        size += sum([get_size(k, seen) for k in obj.keys()])
    elif hasattr(obj, '__dict__'):
        size += get_size(obj.__dict__, seen)
    elif hasattr(obj, '__iter__') and not isinstance(obj, (str, bytes, bytearray)):
        size += sum([get_size(i, seen) for i in obj])
    return size

def hash_function_factory(m, seed):
    return lambda x : int(murmurhash3_32(x, seed) % m)

class bloom_filter:
    def init_cn(self, c, n):
        r_over_n = math.log(c, .618)
        r = n * r_over_n
        k = max(1, round(r_over_n * math.log(2)))
        self.hashes = []
        for _ in range(k):
            self.hashes.append(hash_function_factory(r, random()))
        self.bitarr = bitarray((0 for _ in range(2**math.ceil(math.log(r,2)))))
        return self
    
    def init_nr(self, n, r):
        r_over_n = r / n
        k = max(1, round(r_over_n * math.log(2)))
        self.hashes = []
        for _ in range(k):
            self.hashes.append(hash_function_factory(r, random()))
        self.bitarr = bitarray((0 for _ in range(r)))
        return self

    def insert ( self , key ):
        for hash in self.hashes:
            self.bitarr[hash(key)] = 1
    def test ( self , key ):
        for hash in self.hashes:
            if self.bitarr[hash(key)] == 0:
                return False
        return True
    
    def __sizeof__(self) -> int:
        return sys.getsizeof(self.bitarr) + sys.getsizeof(self.hashes)


# e - fp rate
# k - num hash funcs
# q - bits per entry
# m - size of table
# c - size of table / num keys 
# S - support
# D - domain
# R - Range

class bloomier_filter:

    def __init__(self, elems, c, k, q):
        start = time.time()
        self.m = int(c * len(elems))
        self.k = k
        self.q = q

        S = {elem for elem in elems}

        self.table1 = array.array("B", [0 for _ in range(self.m)])
        self.table2 = [0 for _ in range(self.m)]
        
        res = None
        tries = 0
        while res == None:
            print("finding matching")
            hashes = []
            for _ in range(k):
                hashes.append(hash_function_factory(self.m, random() * (2 ** 30) + tries))
            hashes.append(hash_function_factory(2**q, random() * (2 ** 30)  + tries))

            self.hashes = tuple(hashes)

            res = self.findMatch(S)
            tries += 1
            if tries > 10:
                raise Exception("too many tries with k=", k)


        PI, matching = res
        for t in PI:
            v = elems[t]
            hashes = self.hashAll(t, self.hashes)
            neighborhood = hashes[:len(hashes)-1]
            M = hashes[len(hashes)-1]
            l = matching[t]
            L = neighborhood[l]
            

            self.table1[L] = l ^ M ^ functools.reduce(lambda a,b : a ^ b, [self.table1[neighborhood[i]] for i in range(self.k)], self.table1[L])
            self.table2[L] = v
        # print("made table")
        self.buildTime = time.time() - start

    def findMatch(self, S):
        print("find match", len(S))
        E = set()
        PI = []
        matching = {}
        singletons = self.singletons(S)

        for t in S:
            j = self.tweak(t, singletons)
            if j == None:
                continue
            E.add(t)
            matching[t] = j

        if len(E) == 0:
            return None
        
        PIprime, matchingPrime = [], {}
        H = S.difference(E)
        if len(H) > 0:
            res = self.findMatch(H)
            if res == None:
                return None
            PIprime, matchingPrime = res

        PI = PIprime
        for t in E:
            PI.append(t)
        return PI, matching | matchingPrime

    def tweak(self, t, singletons):

        neighborhood = self.hashAll(t, self.hashes)
        neighborhood = neighborhood[:len(neighborhood)-1]

        for i in range(len(neighborhood)):
            if neighborhood[i] in singletons:
                return i
        return None

    def singletons(self, S):
        locFreqs = defaultdict(lambda : 0)
        for t in S:
            neighborhood = self.hashAll(t, self.hashes)
            neighborhood = neighborhood[:len(neighborhood)-1]
            for n in neighborhood:
                locFreqs[n] += 1
        singles = set()
        for k, v in locFreqs.items():
            if v == 1:
                singles.add(k)
        return singles

    @functools.cache
    def hashAll(self, t, hashes):
        return [hashes[i](t) for i in range(len(hashes))]
    
    @functools.cache
    def findPlace(self, t):
        hashes = self.hashAll(t, self.hashes)
        neighborhood = hashes[:len(hashes)-1]
        M = hashes[len(hashes)-1]
        
        l = functools.reduce(lambda a,b : a ^ b, [self.table1[neighborhood[i]] for i in range(self.k)], M)
        if l < self.k:
            L = neighborhood[l]
            return L
        return None

    def lookup(self, t):
        L = self.findPlace(t)
        if L == None:
            return None
        return self.table2[L]
    
    def set_value(self, t, v):
        L = self.findPlace(t)
        if L == None:
            return False
        self.table2[L] = v
        return True
    
    # def __sizeof__(self) -> int:
    #     return sys.getsizeof(self.table1) + sys.getsizeof(self.table2) + sys.getsizeof(self.hashes)

# fps = []
# bftimes = []
# hashtimes = []
# bfstorages = []
# dictstorages = []
# topofrange = 1000000
# steps = range(15, 2, -1)

# for step in steps:
#     print(step)
#     trials = 1

#     start = time.time()
#     function = {k : k for k in range(0, topofrange, step)}
#     filledFunction = dict(function)
#     # for k in range(0, topofrange):
#     #     if not (k in filledFunction):
#     #         filledFunction[k] = 0
#     hashbuildtime = time.time() - start

#     bfs = [(lambda : bloomier_filter(function,int(1.25 * topofrange), 4, 8))() for _ in range(trials)]

#     fp = 0
#     neg = 0
#     for i in range(topofrange):
#         for bf in bfs:
#             if bf.lookup(i) != None and not i in filledFunction:
#                 fp += 1
#         if  not i in filledFunction:
#             neg += trials

#     fps.append(fp/neg)
#     bftimes.append(bf.buildTime)
#     hashtimes.append(hashbuildtime)
#     bfstorages.append(get_size(bf))
#     dictstorages.append(get_size(filledFunction))

# for i in range(len(steps)):
#     print("1/" + str(steps[i]), fps[i], bftimes[i], bfstorages[i])

# fig, axs = plt.subplots(2,2)

# axs[0][0].plot(steps, fps, label="fprate")
# axs[0][0].set_title("fp rate")

# axs[0][1].plot(steps, [b/h for b, h in zip(bfstorages, dictstorages)], label="storage")
# axs[0][1].set_title("bf storage/hash storage")

# axs[1][0].plot(steps, [b/h for b, h in zip(bftimes, hashtimes)], label="build times")
# axs[1][0].set_title("bf time/hash time")

# fig.show()
# plt.show()
# print(get_size(bf) / get_size(function))


data = pd.read_csv('user-ct-test-collection-01.txt', sep="\t")
querylist = data.Query.dropna()

function = defaultdict(int)
# testSet = defaultdict(int)

for i, q in enumerate(querylist):
    function[q] += 1

total_keys_size = sum([get_size(key) for key in function])


bf = bloomier_filter(function, 1.25, 3, 8)

print(get_size(bf)/total_keys_size, get_size(function)/total_keys_size)


# for k in range(querylist):
#     if bf.lookup(k) != None and not k in function:
#         fp += 1
#     if  not i in filledFunction:
#         neg += trials

# test_set = [str(random()) for _ in range(TEST_NUM)]