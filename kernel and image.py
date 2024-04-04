import numpy as np
import math
import sympy
from sympy import symbols, diff, Poly
from sympy.polys.monomials import itermonomials
from sympy.polys.orderings import monomial_key
from sympy import *
import itertools

x1, x2, x3, x4 = symbols(' x1 x2 x3 x4', real=True)
init_printing(use_unicode = True)

b = 4 #beta parameter

def picker(vals, val, k): #assigns the value to the k'th element of vals
    temp = np.empty(len(vals), sympy.core.power.Pow)
    for i in range(len(vals)):
        if k == i:
            temp[i]  = val
        else:
            temp[i] = 0
    return temp


def dp0peter(vals):
    f = vals[0]
    p1 = 0
    p2 = -1*x1*diff(f, x4)
    p3 = -1*x2*diff(f, x4)
    p4 = (x1*diff(f, x2) + x2*diff(f, x3))
    temp = np.array([p1, p2, p3, p4])
    polytemp = np.array([])
    for i in range(4):
        temp[i] = expand(temp[i])
        polytemp = np.append(polytemp, Poly(temp[i], x1, x2, x3, x4))
    return polytemp 

def dp1peter(vals):
    X_1 = vals[0]
    X_2 = vals[1]
    X_3 = vals[2]
    X_4 = vals[3]
    pp12= -1*x1*diff(X_1,x4)
    pp13= -1*x2*diff(X_1,x4)
    pp14 = (x1*diff(X_1,x2)+x2*diff(X_1,x3))
    pp23 = -1*x2*diff(X_2,x4)
    pp24 = (-2*X_1+x1*diff(X_2,x2)+x1*diff(X_4,x4)+x2*diff(X_2,x3))
    pp34 = (x1*diff(X_3,x2)-X_2+x2*diff(X_3,x3)+x2*diff(X_4,x2))
    temp = np.array([pp12, pp13, pp14, pp23, pp24, pp34])
    polytemp = np.array([])
    for i in range(6):
        temp[i] = expand(temp[i])
        polytemp = np.append(polytemp, Poly(temp[i], x1, x2, x3, x4))
    return polytemp

def dpi0(vals):
    f = vals[0]
    p1 = (-1*b*x1*diff(f, x4))
    p2 = (-1*x2*diff(f, x4))
    p3 = (-(x2 + x3)*diff(f, x4))
    p4 = (b*x1*diff(f, x1) + x2*diff(f, x2) + (x2 + x3)*diff(f, x3))
    temp = np.array([p1, p2, p3, p4])
    polytemp = np.array([])
    for i in range(4):
        #temp[i] = expand(temp[i])
        polytemp = np.append(polytemp, Poly(temp[i], x1, x2, x3, x4))
    return polytemp

def dpi1(vals): #input polynomial array, output polynomial array
    X_1 = vals[0]
    X_2 = vals[1]
    X_3 = vals[2]
    X_4 = vals[3]
    pp12 =(b*x1*diff(X_2, x4) - x2*diff(X_1, x4))
    pp13 = (b*x1*diff(X_3, x4)) - (x2 + x3)*diff(X_1, x4)
    pp14 = (b*x1*(diff(X_1, x1) + diff(X_4, x4)) - b*X_1 + x2*diff(X_1, x2) + (x2 + x3)*diff(X_1, x3))
    pp23 = (x2*diff(X_3, x4) - (x2 + x3)*diff(X_2, x4))
    pp24 = (b*x1*diff(X_2, x1) + x2*(diff(X_2, x2) + diff(X_4, x4)) - X_2 + (x2 + x3)*diff(X_2, x3))
    pp34 = (b*x1*diff(X_3, x1) + x2*diff(X_3, x2) + (x2 + x3)*(diff(X_3, x3) + diff(X_4, x4)) - (X_2 + X_3))
    temp = np.array([pp12, pp13, pp14, pp23, pp24, pp34])
    polytemp = np.array([])
    for i in range(6):
        #temp[i] = expand(temp[i])
        polytemp = np.append(polytemp, Poly(temp[i], x1, x2, x3, x4))
    return polytemp 

def dpi2(vals): #input polynomial array, output polynomial array
    X_12 = vals[0]
    X_13 = vals[1]
    X_14 = vals[2]
    X_23 = vals[3]
    X_24 = vals[4]
    X_34 = vals[5]
    ppp123 = (b*x1 * diff(X_23, x4) - x2 * diff(X_13, x4) + (x2 + x3)*diff(X_12, x4))
    ppp124 = (-1*b*x1*diff(X_12, x1) + b*X_12 + b*x1*diff(X_24, x4) - x2*diff(X_12, x2) + X_12 - x2*diff(X_14, x4) - (x2 + x3)*diff(X_12, x3))
    ppp134 = (-1*b * x1 * diff(X_13, x1) + b*X_13 + b*x1*diff(X_34, x4) - x2*diff(X_13, x2) + X_12 + X_13 - (x2 + x3)*(diff(X_13, x3) + diff(X_14, x4)))
    ppp234 = (-b*x1*diff(X_23, x1) - x2*diff(X_23, x2) + X_23 + x2*diff(X_34, x4) + X_23 - (x2 + x3)*(diff(X_23, x3) + diff(X_24, x4)))
    temp = np.array([ppp123, ppp124, ppp134, ppp234])
    polytemp = np.array([])
    for i in range(4):
        #temp[i] = expand(temp[i])
        polytemp = np.append(polytemp, Poly(temp[i], x1, x2, x3, x4))
    return polytemp

def dpi3(vals):
    X_123 = vals[0]
    X_124 = vals[1]
    X_134 = vals[2]
    X_234 = vals[3]
    pppp1234 = b*x1*diff(X_123, x1) + b*x1*diff(X_234, x4) - b*X_123 + x2*(diff(X_123, x2) - diff(X_134, x4))  - 2*X_123 + (x2 + x3)*(diff(X_123, x3) + diff(X_124, x4))
    temp = np.array([pppp1234])
    polytemp = np.array([])
    for i in range(1):
        #temp[i] = expand(temp[i])
        polytemp = np.append(polytemp, Poly(temp[i], x1, x2, x3, x4))
    return polytemp

def kerdn(n, degree, map): #dpi_n and map need to agree (will get outofbounds though)
    numInp = math.comb(4, n) #number of input partials
    inputFunc = np.empty(numInp)

    terms = math.comb(degree + 4 - 1, 4 - 1) #number of monomials of degree "degree"
    basis = sorted(itermonomials([x1, x2, x3, x4], degree, degree), key=monomial_key('grlex', [x4, x3, x2, x1])) #ordered basis

    outputVec = np.array([]) 
    for k in range(numInp): #for each partial
        for poly in basis:
            chosen = picker(inputFunc, poly, k) #assign it to one of them
            applied = map(chosen) #apply whatever map
            numOut = len(applied) #num of partials in the output

            for i in range(numOut): #big list of all of the coefficients then by partials
                for polz in basis:
                    coeff = applied[i].coeff_monomial(polz) #extract the relevant coefficient
                    outputVec = np.append(outputVec, coeff)
    N = returnNullspace(outputVec, terms, numInp, numOut) 
    processNullspace(N, degree, basis, numInp, map)
    return True

def returnNullspace(bigVec, terms, inp, out): #reshapes big list into proper matrix and returns nullspace
    temp = bigVec.reshape(inp * terms, out * terms)
    temp = temp.transpose()
    M = Matrix(temp)
    return M.nullspace()


def processNullspace(N, degree, basis, numInp, map): #does the printing output
    terms = len(basis)
    for j in range(len(N)): #for each vector in the nullspace
        for i in range(numInp):
            k = N[j][terms*i:terms*(i+1)] #print it out per partial
            print(k)

        testTerms = np.empty(numInp, sympy.core.power.Pow) #make sure that it truly is in the nullspace
        for h in range(0, numInp):
            asPoly = 0
            for l in range(terms):
                asPoly += N[j][l + terms*h]*basis[l]
            print(asPoly)
            testTerms[h] = asPoly
        checking = map(testTerms)
        if np.all(checking == 0):
            print("\r\nverified")
        else:
            print("uh oh")
        print("_______")
    print("number of terms: " + str(len(N)))
    print("all done")
    return True

#just copied dmitry code
def imdn(n, degree, map): #dpi_n and map need to agree (will get outofbounds though) 
    numInp = math.comb(4, n) #number of input partials
    inputFunc = np.empty(numInp)

    terms = math.comb(degree + 4 - 1, 4 - 1) #number of monomials of degree "degree"
    basis = sorted(itermonomials([x1, x2, x3, x4], degree, degree), key=monomial_key('grlex', [x4, x3, x2, x1])) #ordered basis
    
    outputVec = np.array([]) 
    for k in range(numInp): #for each partial
        for poly in basis:
            chosen = picker(inputFunc, poly, k) #assign it to one of them
            applied = map(chosen) #apply whatever map
            numOut = len(applied) #num of partials in the output

            for i in range(numOut): #big list of all of the coefficients then by partials
                for polz in basis:
                    coeff = applied[i].coeff_monomial(polz) #extract the relevant coefficient
                    outputVec = np.append(outputVec, coeff)
    N = returnImage(outputVec, terms, numInp, numOut) 
    processImage(N, degree, basis, numOut, map)
    return True

def returnImage(bigVec, terms, inp, out): #reshapes big list into proper matrix and returns image
    temp = bigVec.reshape(inp * terms, out * terms)
    temp = temp.transpose()
    M = Matrix(temp)
    return M.columnspace()

def processImage(N, degree, basis, numInp, map): #does the printing output
    terms = len(basis)
    for j in range(len(N)): #for each vector in the image
        for i in range(numInp):
            k = N[j][terms*i:terms*(i+1)] #print it out per partial
            print(k)
        temp = np.empty(numInp, sympy.core.power.Pow)
        for h in range(0, numInp):
            asPoly = 0
            for l in range(terms):
                asPoly += N[j][l + terms*h]*basis[l]
            print(asPoly)
        print("_______")
    print("number of terms: " + str(len(N)))
    print("all done")
    return True

def checkComposite(first, second, degree, f):
    partials = math.comb(4, f)
    inputFunc = np.empty(partials)
    basis = sorted(itermonomials([x1, x2, x3, x4], degree, degree), key=monomial_key('grlex', [x4, x3, x2, x1]))
    allGood = True
    for k in range(partials):
        for poly in basis:
            chosen = picker(inputFunc, poly, k)
            singleApp = first(chosen)
            doubleApp = second(singleApp)
            if np.all(doubleApp != 0):
                allGood = False
    if allGood:
        print("checks composition")
    else:
        print("BAD")

def bracketAndMaybeCheck(SC, degree):
    n = degree
    wumbo = list(itertools.combinations([1, 2, 3, 4], 2))
    print(wumbo)
    C = [[[0 for k in range(n)] for j in range(n)] for i in range(n)]
    pi = [[0 for j in range(n)] for i in range(n)]
    p = 0
    for pair in itertools.combinations([1, 2, 3, 4], 2):
        i = pair[0] - 1
        j = pair[1] - 1
        for k in range(degree):
            C[i][j][k] = SC[p][k]
            C[j][i][k] = -1*SC[p][k]
            pi[i][j] += SC[p][k]
        p += 1
    print(pi)
    return True

def generateBivector(SC, dim):
    pairs = list(itertools.combinations([1, 2, 3, 4], 2)) #4 depends on the dim
    pi = [0 for i in range(len(pairs))]
    for l in range(len(pairs)):
        i = pairs[l][0] - 1
        j = pairs[l][1] - 1
        for k in range(dim):
            pi[l] += SC[l][k]
    print(pi)

    return pi

def checkPoisson(bivred, dim):
    bivec = [[0 for j in range(dim)] for i in range(dim)]
    pairs = list(itertools.combinations([1, 2, 3, 4], 2))
    for l in range(len(pairs)):
        i = pairs[l][0] - 1
        j = pairs[l][1] - 1
        bivec[i][j] = bivred[l]
        bivec[j][i] = -1*bivred[l]
    print(bivec)
    triples = list(itertools.combinations([1, 2, 3, 4], 3)) #4 depends on the dim
    v = [x1, x2, x3, x4]
    for triple in triples:
        i = triple[0] - 1
        j = triple[1] - 1
        k = triple[2] - 1
        sum = 0
        for l in range(dim):
            xl = v[l]
            sum += bivec[i][l]*diff(bivec[j][k], xl) + bivec[j][l]*diff(bivec[k][i], xl) + bivec[k][l]*diff(bivec[i][j], xl)
        if sum == 0:
            print("hooray")
        else:
            print("bad")

#checkComposite(dpi2, dpi3, 4, 2)
#kerdn(1, 4, dpi1)
#imdn(0, 4, dpi0)

SC1 = [[0 ,0, 0, 0],
      [0, 0, 0, 0],
      [b*x1, 0, 0, 0],
      [0, 0, 0, 0],
      [0, 1*x2, 0, 0],
      [0, 1*x2, 1*x3, 0]]
SC2 = [[0 ,0, 0, 0], #random guy online (g, 4.4)
       [0 ,0, 0, 0],
       [x1 ,0, 0, 0],
       [0 ,0, 0, 0],
       [0 ,x1 + x2, 0, 0],
       [0 ,0, x2 + x3, 0],]

pi = generateBivector(SC2, 4)
checkPoisson(pi, 4)
