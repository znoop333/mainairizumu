# -*- coding: utf-8 -*-
"""
Created on Tue Mar  7 09:35:46 2017

@author: znoop333
"""

"""
Pavel Curtis

@julzerator

Look up - Knuth's sudoku "Dancing 9" algorithm.


Mainarizumu (John Lopez @godeke;  @julzerator from Australia, Nick @foton, @DoctaStooge, @jcooper0) 
- usually 7x7 grid. (17x17 is "Olymic" version)
- like Sudoku (inspired by crossword puzzles, which don't work in Japanese)

Row constraint: every row must have the integers 1-N, exactly once

Column constraint: every column must have the integers 1-N, exactly once

1. Greater than constraint: left cell must be greater than right cell. 
(e.g., so left cell can't be 1)

2. Less than constraint: ibid.

3. Difference constraint: absolute difference between adjacent cells is fixed.

Known value constraint: a cell has a known value (e.g., given for children's
version of puzzle.)


"Satisfying" puzzle: brute forcing and guessing should not be required. 

Unique puzzle: there should be only 1 solution.

Satisfiability solver - mathematical model for constraints and outputs.

"""


"""
Evaluation of state for a solution:
    if every cell has exactly one possible value:
        the puzzle is uniquely solved.
    else: 
        assign a cell then recursively Evaluate

When assigning a value to a cell:
    record history of when it was assigned (depth of search) and why (which rule caused it)
    
"""


"""
Input - text format:
         col 01234

row 0:       _>_>_
boundary 1:  _____
row 2:       __2__
boundary 3:  1____
row 4:       _____

>    greater than constraints: cell[0,0] > cell[0,1] and cell[0,1] > cell[0,2] 
2    Known value in cell[1,1]
1    Difference constraint (vertical) in [2,0]: cell[1,0] - cell[2,0] = +/- 1 

separate arrays: S - PossibleValues (the current state of the board)
[2,3,4]   [1,3,4]   [1,2,3,4]
[1,3,4,5] [2]       [1,3,4,5]
[1,3,4,5] [1,3,4,5] [1,2,3,4,5]

H - ConstraintsArrayHorizontal: horizontal constraints (always has empty left and right columns to make it N+1xN+1)
[]    [>]   [>]  []
[]    []    []   []
[]    []    []   []
[]    []    []   []

V - ConstraintsArrayVertical: vertical constraints (always has an empty top row to make it N+1xN+1)
[]    []    []   []
[]    []    []   []
[1]   []    []   []
[]    []    []   []


A puzzle is solved when S has exactly one value in every cell and all constraints are satisfied.

A puzzle is inconsistent when either: any cell of S is empty (no possible solutions) 
or a constraint cannot be satisfied 

A puzzle has multiple solutions if there exists S1 and S2 with distinct values that have exactly 
one value in each cell, and all constraints are satisfied

A solution generator takes S1 and creates S2 by assigning exactly one value to S1[r,c] = v
(out of the possible values for this cell). Constraints are applied to S2 to eliminate other values
in cells in this order:

    Row constraint (r): the other cells in this row cannot have value v.
         Eliminate values in S1[r,C] = v for all C != c. Mark them as modified.

    Column constraint (c): the other cells in this column cannot have value v
         Eliminate values in S1[R,c] = v for all R != r. Mark them as modified.

    H constraints: h1 = H[r,c], h2 = H[r,c+1]
    V constraints: v1 = V[r,c], v2 = H[r+1,c]
    
    For each constraint K in [h1, h2, v1, v2]:
        if K = '>' (greater than constraint): 
           Eliminate values of S1 that violate this constraint and mark those cells as modified.

    let S2 = modified S1

    for each modified cell:
        if no values remain - the puzzle is inconsistent. exit with error.
        else  -- if only one value remains - check if the puzzle is solved
            recurse on S2 to apply constraints
    

    Greater/Less than constraint - H: 
        if H[r,c] = '>', S1[r,c-1] must be > S1[r,c]. 
           Eliminate values of S1[r,c-1] that violate this constraint and mark S1[r,c-1] as modified.

        if H[r,c+1] = '>', S1[r,c] < S1[r,c+1]. 
           Eliminate values of S1[r,c+1] that violate this constraint and mark S1[r,c+1] as modified.
        
        Repeat for V[r,c] = '>' and V[r+1,c] = '>'
        
    Greater/Less than constraint - V:
        same as above but use r +/- 1 instead of c +/- 1
    
    Difference constraint - H:
        if H[]

"""

import numpy as np

#a = np.array([[1, 2], [3,4]], dtype='byte')
#a = np.array([['a', 'b'], ['c','c']], dtype='S1')
#print(a[0,:])
#print(str(a[0,:]))

#np.array()

"""
logical size: N = 3 above
storage grid size: M = 2N - 1

grid layout: N x N x N+1 so that last dimension can be accessed as 1-N not 0-(N-1)

dimension 3: 
    
"""

N = 5
M = N + 1

S = np.ones([N,N,N+1], 'byte')
H = np.zeros([M,M], 'S1')
V = np.zeros([M,M], 'S1')

# put in known value constraints
S[0,1,:] = 0
S[0,1,3] = 1
 
# put in less/greater than constraints
H[0,1] = '<'
H[1,3] = '<'
H[2,2] = '<'
H[3,3] = '>'

# use 'v' for lower cell is less than upper cell. (no good character available)
V[2,2] = 'v'
V[3,2] = 'v'

# use '^' for lower cell is greater than upper cell

# a cell is solved when there is exactly one possible integer 1:N
# which is marked as a possible solution (1). 
# also, element 0 in the third dimension should be changed from 1 (unsolved)
# to 0 when the cell is solved.
def isCellSolved(S_, N_, ii, jj):
    solved = S_[ii,jj,1:N_+1].nonzero()
    if solved[0].size == 1:
        return 1
    return 0

#print( isCellSolved(S, N, 0,0) )
#print( isCellSolved(S, N, 0,1) )

# a cell cannot be solved if all of the possible integers 1-N are marked 0
def canCellBeSolved(S_, N_, ii, jj):
    return S_[ii,jj,1:N_+1].any()

#print( canCellBeSolved(S, N, 0,0) )
#print( canCellBeSolved(S, N, 0,1) )

# if a cell is solved, return the solution.
# if a cell is not solved but it has more than 1 possible value, return the first value
# else return 0
def firstCellPossibleSolution(S_, N_, ii, jj):
    solved = S_[ii,jj,1:N_+1].nonzero()
    if solved[0].size >= 1:
        return solved[0][0]+1
    return 0

#print( firstCellPossibleSolution(S, N, 0,0) )
#print( firstCellPossibleSolution(S, N, 0,1) )

# a 2D text grid of whether each cell can be solved (first possible solution)
# or 0 for cells that can't be solved
def PuzzleStatusString(S_, N_):
    t = ''
    for ii in range(N_):
        for jj in range(N_):
            t += str(firstCellPossibleSolution(S_, N_, ii, jj)) + ' '
        t += '\n'
    return t

print('PuzzleStatusString:')

print(PuzzleStatusString(S, N))

print('PuzzleCountString:')
print(PuzzleCountString(S2, N))

S2 = S.copy()
S2[1,1,:] = 0

print('PuzzleStatusString (no solution):')
print(PuzzleStatusString(S2, N))

def countPossibleSolutions(S_, N_, ii, jj):
    solved = S_[ii,jj,1:N_+1].nonzero()
    return solved[0].size

def PuzzleCountString(S_, N_):
    t = ''
    for ii in range(N_):
        for jj in range(N_):
            t += str(countPossibleSolutions(S_, N_, ii, jj)) + ' '
        t += '\n'
    return t

print('PuzzleCountString (no solution):')
print(PuzzleCountString(S2, N))

    
# if any cell in the puzzle can't be solved, the whole puzzle is now inconsistent
def canPuzzleBeSolved(S_, N_):
    for ii in arange(N_):
        for jj in arange(N_):
            if not canCellBeSolved(S_, N_, ii, jj):
                return 0
    return 1
            
print(canPuzzleBeSolved(S, N))
print(canPuzzleBeSolved(S2, N))

#    for ii in range(N_):
#        if np.max()
#    pass

# arguments are mutable! remember to copy() if you don't want to change them.
#def mutate(S_):
#    S_ = S_.copy()
#    S_[0,0,:]=0
#
#mutate(S2)
#print(PuzzleStatusString(S2, N))

def AssignValue(S_, N_, rr, cc, vv):
    S_[rr, cc, :] = 0
    S_[rr, cc, vv] = 1

# return a new puzzle where the value vv has been removed from all columns in row rr.
# if no changes have been made, (S_ == S2).all() is true.
# this puzzle may be inconsistent.
def ApplyRowConstraint(S_, N_, rr, cc, vv):
    if cc > 0 and S_[rr, 0:cc, vv].any() \
                    or cc < N_-1 and S_[rr, cc+1:N_, vv].any():
        S2 = S_.copy() 
        if cc > 0:
            S2[rr, 0:cc, vv] = 0
            for c in range(cc):
                S2[rr, c, 0] = S2[rr, c, :].any()
        if cc < N_-1:
            S2[rr, cc+1:N_, vv] = 0
            for c in range(cc+1, N_):
                S2[rr, c, 0] = S2[rr, c, :].any()
        return S2
    return S_

# same but along rows in a single column
def ApplyColConstraint(S_, N_, rr, cc, vv):
    if rr > 0 and S_[0:rr, cc, vv].any() \
                    or rr < N_-1 and S_[rr+1:N_, cc, vv].any():
        S2 = S_.copy() 
        if rr > 0:
            S2[0:rr, cc, vv] = 0
            for r in range(rr):
                S2[r, cc, 0] = S2[r, cc, :].any()
        if rr < N_-1:
            S2[rr+1:N_, cc, vv] = 0
            for r in range(rr+1, N_):
                S2[r, cc, 0] = S2[r, cc, :].any()
        return S2
    return S_
      
S2 = S.copy()
AssignValue(S2, N, 3, 0, 2)
S2[:,0,1] = 0
#S2[3,0,:] = 0
#S2[3,0,4] = 1

print('ApplyColConstraint')
print(PuzzleStatusString(S2, N))
S3 = ApplyColConstraint(S2, N, 3, 0, 2)  
  
print(PuzzleStatusString(S3, N))
print((S3 == S2).all())

# ApplyColConstraint() with the same arguments should be idempotent
#S4 = ApplyColConstraint(S2, N, 3, 0, 2)  

S4 = ApplyRowConstraint(S3, N, 0, 0, 3)
print(PuzzleStatusString(S4, N))
print((S3 == S4).all())
  
print(canPuzzleBeSolved(S4, N))

# if a new puzzle is created by ApplyColConstraint(), ApplyRowConstraint()
# must be run for each row that changed in the column to propagate the
# new constraints, possibly leading to a unique solution or no solution.

# if there is a Horizontal constraint that applies to cell[rr,cc], 
# apply it to the adjacent cell and generate a new state S2 such that
# illegal values have been eliminated. if a new state is generated, the
# row and column constraints must be reapplied to the puzzle.
#
# there can be up to 2 H constraints applying to cell[rr,cc]: 
# one where cell[rr,cc] OP cell[rr,cc+1] (so that the reference cell is LEFT)
# another where cell[rr,cc-1] OP cell[rr,cc] (so that the reference cell is RIGHT)
# where OP is "<", ">", "i" for i < N (difference).
# the constraints are applied as LEFT then RIGHT.
def ApplyHConstraint(S_, N_, rr, cc):
    rightOP = H[rr, cc]
    leftOP = H[rr, cc+1]

    if leftOP == b'' and rightOP == b'':
        return S_ # no changes
    
    S2 = S_.copy() 
    
    # set(L) < set(R) implies:
    # l < max(R) for all l in L 
    # and
    # r > min(L) for all r in R

    # e.g.: L[3,4,5] < R[2,3,4] implies
    # max(R) = 4 
    # L < 4 so L becomes [3] (exclude 4,5)
    # min(L) = 3
    # R > 3 so R becomes [4] (exclude 2,3)
    # now max(R) = 4
    # L < 4 so L becomes [3] (exclude 4)
    # min(L) = 3, no further changes. L[3] < R[4]
    if rightOP == b'<' and cc > 0:
        minL = S2[rr, cc-1, 1:].min()
        maxR = S2[rr, cc, 1:].max()

        S2[rr, cc-1, maxR:] = 0
        S2[rr, cc, :minL+1] = 0
          
        # restore the invariants
        S2[rr, cc-1, 0] = S2[rr, cc-1, :].any()
        S2[rr, cc, 0] = S2[rr, cc, :].any()
        
    # set(L) > set(R) implies
    # same as above but reversed L/R
    if rightOP == b'>' and cc > 0:
        maxL = S2[rr, cc-1, 1:].max()
        minR = S2[rr, cc, 1:].min()

        S2[rr, cc, maxL:] = 0
        S2[rr, cc-1, :minR+1] = 0
          
        # restore the invariants
        S2[rr, cc-1, 0] = S2[rr, cc-1, :].any()
        S2[rr, cc, 0] = S2[rr, cc, :].any()

    # [rr, cc]  leftOP  [rr, cc+1]
    if leftOP == b'<' and cc < N:
        minL = S2[rr, cc, 1:].min()
        maxR = S2[rr, cc+1, 1:].max()

        S2[rr, cc, maxR:] = 0
        S2[rr, cc+1, :minL+1] = 0
          
        # restore the invariants
        S2[rr, cc, 0] = S2[rr, cc, :].any()
        S2[rr, cc+1, 0] = S2[rr, cc+1, :].any()
    
    if leftOP == b'>' and cc < N:
        maxL = S2[rr, cc, 1:].max()
        minR = S2[rr, cc+1, 1:].min()

        S2[rr, cc, maxL:] = 0
        S2[rr, cc+1, :minR+1] = 0
          
        # restore the invariants
        S2[rr, cc, 0] = S2[rr, cc, :].any()
        S2[rr, cc+1, 0] = S2[rr, cc+1, :].any()

            
    pass

# np.nonzero
# r_


aa = np.array([0,1,0,0,1])
bb = np.array([0,0,1,1,0])
a_ = aa[:7]
a1_ = a_.argmax()


#b_ = aa[:2]
#print(a_)
#print(b_)
#print('ApplyLTConstraint')
#ApplyLTConstraint(a_, b_)
#print(a_)
#print(b_)
#print(aa)

# modifies L_, R_, both, or neither to enforce L_ < R_ for all
# elements of L_ and R_, excluding 0th element (isNotSolved invariant)
def ApplyLTConstraint(L_, R_):
    N_, = L_.shape
    # argmax - Only the first occurrence of the max is returned, but I want
    # the last occurance. So reverse the array, call argmax, and fix the offset
#    minL = N_ - L_[::-1].argmin()
#    maxR = N_ - R_[::-1].argmax()

    # better solution: since it's a bitfield, it's sorted.
    minL = L_.nonzero()[0][0]
    maxR = R_.nonzero()[0][-1]
    
#    print('minL {0}, maxR {1}'.format(minL, maxR))
    if L_[maxR:].any() or R_[1:minL+1].any():
        L_[maxR:] = 0
        R_[1:minL+1] = 0
          
        # restore the invariants
        L_[0] = L_[1:].any()
        R_[0] = R_[1:].any()
        
        return True # change(s) were made. rerunning suggested until this function returns False
    return False # no changes

print('0 aa {0}, bb {1}'.format(aa, bb))
changed = ApplyLTConstraint(aa, bb)
print('1 aa {0}, bb {1}, changed {2}'.format(aa, bb, changed))
changed = ApplyLTConstraint(aa, bb)
print('2 aa {0}, bb {1}, changed {2}'.format(aa, bb, changed))


# difference constraint k - 
# S[rr, cc] - S[rr, cc+1] = +/- k
# where k is [1,N-1]

# set(L) - set(R) = +/- k 
# for all l in L, if exists (l+k) or (l-k) in R, then keep l. else discard.
# for all r in R, if exists (r+k) or (r-k) in L, then keep r. else discard.

# e.g.   L[1,2,3,4] - R[1,2,4,6] = +/- 3                                
# L + 3 = [4,5,6,7] which matches at R[4,6]
# only keep L[1,3]
# L - 3 = [-2, -2, 0, 1] which matches at R[1].
# only keep L[4]
# therefore L[1,3,4] (union)

# R + 3 = [4,5,7,9] which matches at L[4]
# keep R[1]
# R - 3 = [-2, -1, 1, 3] which matches at L[1,3]
# keep R[4,6]
# therefore R[1,4,6] (union)
def ApplyDifferenceConstraint(L_, R_, k_):
    N_, = L_.shape
    modified = False
    for ii in range(N_):
        if L_[ii] == 0:
            continue

        lPlus = L_[ii] + k_
        lMinus = L_[ii] - k_
        if lPlus < N_ and R_[lPlus] == 1 \
                            or lMinus > 0 and R_[lMinus] == 1:
            pass
        else:
            L_[ii] = 0
            modified = True

    for ii in range(N_):
        if R_[ii] == 0:
            continue

        lPlus = R_[ii] + k_
        lMinus = R_[ii] - k_
        if lPlus < N_ and L_[lPlus] == 1 \
                            or lMinus > 0 and L_[lMinus] == 1:
            pass
        else:
            R_[ii] = 0
            modified = True
                     
    return modified 
    
#    l = L_[1:].nonzero()[0]
#    lPlus = l + k_
#    lMinus = l - k_
#    
#    print('lPlus: {0}'.format(lPlus))
#    print('lMinus: {0}'.format(lMinus))
#
#    lPlus = lPlus[lPlus > 0 and lPlus < N_]
#    lMinus = lMinus[lMinus > 0 and lMinus < N_]
#    
#    print('lPlus: {0}'.format(lPlus))
#    print('lMinus: {0}'.format(lMinus))
#    
#    lSearch = np.concatenate(lPlus, lMinus)
#    print('lSearch: {0}'.format(lSearch))
#    
#    lSearch = np.unique( lSearch )

    pass

aa = np.array([0,1,1,1,1,0,0]) # L[1,2,3,4]
bb = np.array([0,1,1,0,1,0,1]) # R[1,2,4,6]
changed = False
print('0 aa {0}, bb {1}, changed {2}'.format(aa, bb, changed))
changed = ApplyDifferenceConstraint(aa, bb, 2)
print('1 aa {0}, bb {1}, changed {2}'.format(aa, bb, changed))
changed = ApplyDifferenceConstraint(aa, bb, 2)
print('2 aa {0}, bb {1}, changed {2}'.format(aa, bb, changed))






