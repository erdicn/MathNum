import numpy as np

#elimination de gauss 
# M matrice (n+1)*n car le +1 cest les solutions 
def ResolutionGauss(M):
    """returns a matrix with all the solutions of each variable [x0, x1, x3,...,xn]"""
    # print(M)
    sol_M = [[] for _ in range(len(M))]
    for i in range(len(M)):
        zero_counter = True # becomes false after the first non zero encounter
        sol_M[i].append(0)
        for j in range(len(M[i])-1):
            if zero_counter and M[i][j] == 0:
                sol_M[i][0] += 1
            elif zero_counter: # if the first condition is not satisfied and the zero counter is still true 
                zero_counter = False
            sol_M[i].append(M[i][j])
        sol_M[i].append(M[i][-1])
    
    # We need to sort our solution matrix acording to the pivots 
    sol_M.sort()
    
    # print(sol_M)
    
    def G_Listadd(L1, L2):
        """returns the L1 + L2 but it doesnt add the 1st elements(since they are where the 0s start
           precondition: len(L1) == len(L2)"""
        new_L = [0]
        until_pivot = True
        for i in range(1, len(L1)-1):
            added = L1[i]+L2[i]
            if until_pivot and added == 0:
                new_L[0] += 1
            elif until_pivot:
                until_pivot = False
            new_L.append(added)
        new_L.append(L1[-1]+ L2[-1])
        return new_L
    assert(G_Listadd([1, 0, 1, 2], [0, 1,2,3]) == [0, 1, 3, 5])
    
    def G_Listsubstract(L1,L2):
        """returns the L1 - L2 but it doesnt substract the 1st elements(since they are where the 0s start
           precondition: len(L1) == len(L2)"""
        new_L = [0]
        until_pivot = True
        for i in range(1, len(L1)-1):
            substracted = L1[i]-L2[i]
            if until_pivot and substracted == 0:
                new_L[0] += 1
            elif until_pivot:
                until_pivot = False
            new_L.append(substracted)
        new_L.append(L1[-1]- L2[-1])
        return new_L
    assert(G_Listsubstract([1, 0, 1, 2], [0, 1,2,3]) == [0, -1, -1, -1])
    
    
    def G_ListxScalar(L, s):# Dilatation ?
        """returns the L multiplied by s however doenst multiply the 0th element since it is the until pivot"""
        if len(L) == 0: return []
        return [L[0]] + [L[i]*s for i in range(1, len(L))]
    assert(G_ListxScalar([1, 2, 3], 2) == [1, 4, 6])
    assert(G_ListxScalar([0, 2, 3], 2) == [0, 4, 6])
    assert(G_ListxScalar([], 2) == [])
    
    # TODO here i am not sure if i have all the test cases 
    for i in range(len(sol_M)):
        for j in range(len(sol_M)):
            if sol_M[i][0] != sol_M[j][0] or j == i: continue # we break when we cant get the pivots further or we continue since we might hace the same pivots after  
            pivot = sol_M[i][0] + 1 
            multiplier = sol_M[j][pivot] / sol_M[i][pivot]
            # print(multiplier)
            sol_M[j] = G_Listsubstract(sol_M[j], G_ListxScalar(sol_M[i], multiplier))
    # print(sol_M)
    
    # now we go back up 
    sol_M.sort()
    
    before = -1
    for i in sol_M:
        if i[0] != before+1:
            print("Error the gauss algorithm does not work properly or the matrix does not have a solurion probaply need a more thourough look at above")
        before += 1
    if sol_M[-1][0] != len(sol_M) - 1: 
        print("solution does not exist")
        return []# solutions does not exist
        
    for i in range(1, len(sol_M)+1):
        for j in range(i+1, len(sol_M)+1):
            multiplier = sol_M[-j][-i-1] / sol_M[-i][-i-1]
            sol_M[-j] = G_Listsubstract(sol_M[-j] , G_ListxScalar(sol_M[-i], multiplier))
            
    solutions = [0 for _ in range(len(M))]
    for i in range(1, len(sol_M)+1):
        solutions[-i] = sol_M[-i][-1] / sol_M[-i][-i-1]
    #     print(solutions[-i], sol_M[-i][-1] / sol_M[-i][-i-1], sol_M[-i][-1], sol_M[-i][-i-1])
    # print(solutions)
    # print(len(sol_M), len(sol_M[0]))
    return solutions
            
    
    
test_M = [[1,3,5],
          [0,1,1]]
test_M1= [[1,3,5],
          [2,7,11]]
test_M2= [[0,1,1],
          [1,3,5]]
test_M3 = [[3.0, 2.0, -4.0, 3.0], [2.0, 3.0, 3.0, 15.0], [5.0, -3, 1.0, 14.0]]
assert(ResolutionGauss(test_M3) == [3, 1, 2])
print("solutions of ", test_M,"=" , ResolutionGauss(test_M))
print("solutions of ", test_M2,"=" , ResolutionGauss(test_M2))
print("solutions of ", test_M3,"=" , ResolutionGauss(test_M3))
print(ResolutionGauss([[4,8, 56, 21], [8, 56, 272, 56], [56, 272, 1568, 320]]))
print("real sols ",29/5, 71/60, -5/24)
#assert(ResolutionGauss([[4,8, 56, 21], [8, 56, 272, 56], [56, 272, 1568, 320]])== [29/5, 71/60, -5/24])

# Symbolic polynomial calculations using lists and indices as implicit constants as ax^i where i is the indice and a is the constant before which is the element in the ithe place
def PolynomialTimesConstant(P, c): 
    """recomende de utiliser PolynomialTimesPolynomial avec seulement P, [s] a la place de ceci si on nest pas 100% sur qu c va rester une constante"""
    return [i*c for i in P]

assert(PolynomialTimesConstant([1,2,3], 3) == [3,6,9])
assert(PolynomialTimesConstant([1,2,3], 0) == [0,0,0])
assert(PolynomialTimesConstant([0,0,0], 0) == [0,0,0])
assert(PolynomialTimesConstant([], 0) == [])
assert(PolynomialTimesConstant([], 1) == [])

def PolynomialDivisionConstant(P, c): 
    """Precondition: c != 0"""
    return PolynomialTimesConstant(P, 1/c)

assert(PolynomialDivisionConstant([3,6,9], 3) == [1,2,3])

def PolynomialTimesPolynomial(P1, P2):
    new_P = [0 for i in range(len(P1)+len(P2))] # multiplication de deux polynomale peux resulter maximum a un polynome de degree de leur somme 
    for i in range(len(P1)):
        for j in range(len(P2)):
            new_P[i+j] += P1[i] * P2[j]
    while len(new_P) and new_P[-1] == 0:
        new_P.pop()
    return new_P 

assert(PolynomialTimesPolynomial([1, 2, 3, 4], [1,1]) == [1, 3, 5, 7, 4])
assert(PolynomialTimesPolynomial([1, 2, 3, 4], [1,1, 1, 1, 1, 1]) == [1,3,6,10, 10, 10, 9, 7, 4])
assert(PolynomialTimesPolynomial([1], [2,1]) == [2,1])
assert(PolynomialTimesPolynomial([2], [1,1]) == [2,2])
assert(PolynomialTimesPolynomial([2], [1,1, 0 , 0]) == [2,2])
assert(PolynomialTimesPolynomial([0, 2], [1,1, 0 , 0]) == [0, 2,2])

# TODO faire long division https://en.wikipedia.org/wiki/Polynomial_long_division pour au futur si on a besoin de diviser des polynomes 

def PolynomialAdition(P1, P2):
    # TODO ici je suis pas sur si cest la meilleur moyen de faire pour determiner le plus gand et petit degree polynomiale
    # Peut etre les parcourir en meme temps et les append a la liste ?
    big_P = []
    small_P = []
    if len(P1) > len(P2): 
        big_P = P1
        small_P = P2
    else: 
        big_P = P2 
        small_P = P1
        
    new_P = [i for i in big_P]
    for i in range(len(small_P)):
        new_P[i] += small_P[i]
        
    while len(new_P) != 0 and new_P[-1] == 0:
        new_P.pop()
    
    return new_P 

assert(PolynomialAdition([1,2,3], [1,2]) == [2,4, 3])
assert(PolynomialAdition([1,2,3], [0,0,0,0,0,0,0]) == [1,2, 3])
assert(PolynomialAdition([1,2,0], [1,2]) == [2,4])
assert(PolynomialAdition([1,2,0], [-1,-2]) == [])


def PolynomialSubstraction(P1, P2):
    """returns P1 - P2
       multiplies P2 by -1 then adds them together""" # TODO not sure if it is optimal or not
    return PolynomialAdition(P1, PolynomialTimesConstant(P2, -1))

assert(PolynomialSubstraction([1,2,3], [1,2]) == [0,0, 3])
assert(PolynomialSubstraction([1,2,3], [0,0,0,0,0,0,0]) == [1,2, 3])
assert(PolynomialSubstraction([1,2,0], [1,2]) == [])
assert(PolynomialSubstraction([1,2,0], [-1,-2]) == [2,4])

def ListToFunction(L, x):
    if len(L) == 0:
        return 0
    tot = L[0]
    for i in range(1, len(L)):
        tot += L[i] * (x)**i
    return tot

def ListToStr(L):
    if len(L) == 0:
        return 0
    s = "" + str(L[0])
    for i in range(1, len(L)):
        if L[i] < 0 : s +=  "{:.2f}x^{}".format(L[i], i)
        else :        s += "+{:.2f}x^{}".format(L[i], i) #si cest plus grand que 0 le + ca sajoute pas automatiquement donc il faut l'ajouter
    return s

def LangragePolynomial(data, degree = -1):
    """:param data :is an array of tuples with (x,f(x))
       :param degree: the degree which we want to approach (n)
       TODO returns a list of constants before the x^n ex: [1, 3 , 2, 0, 5] => 1x^0 + 3x^1 + 2x^2 + 0*x^3 + 5*x^4 """
    if degree == -1: degree = len(data) - 1
    
    #  TODO x f(x) here can only be constants maybe change them with polynomial multiplication ???
    #  Comme pour linstant on a pas la fonctionalite de division on les stock apart et a la fin on les mets ensemble
    L_up   = [[1]]* (len(data)) # on fait degree + 1 car dans les listes on prends on compte la constante aussi
    L_down = [1]  * (len(data)) 
    L      = [[]] * (len(data)) 
    P      = [0]  * (len(data)) 
    for i in range(len(data)):
        for j in range(len(data)):
            if i == j:
                continue
            L_up[i]    = PolynomialTimesPolynomial(L_up[i], [-data[j][0], 1])
            L_down[i] *= (data[i][0]-data[j][0])                                    # here all data x are constants atention  
        L[i] = PolynomialDivisionConstant(L_up[i], L_down[i])                       # Atention here we assume that all x_i - x_j  are constants
        P    = PolynomialAdition(P, PolynomialTimesConstant( L[i] , data[i][1]))      
    while len(P) != 0 and P[-1] == 0:
        P.pop()
    return P    
assert(LangragePolynomial([(0, 1), (1, 0), (2,-1)]) == [1, -1])
assert(LangragePolynomial([(-2, 3), (0, 5), (4,8), (6, 5)]) == [5.0, 1.25, 0.04166666666666671, -0.04166666666666667])
assert(LangragePolynomial([(1,1), (2,4), (3,9)]) == [0,0,1])
assert(LangragePolynomial([(1,1), (2,8), (3,27), (4,64)]) == [0,0,0,1])
## Dans ceci on a des floating point calculation errors at the end we have [0.0, -1.4210854715202004e-14, 1.0, -1.7763568394002505e-15] which is really close to y = x^2
#assert(LangragePolynomial([(1,1), (2,4), (3,9), (4, 16)]) == [0,0,1])
#assert(LangragePolynomial([(1,1), (2,4), (3,9), (4, 16), (5, 25)]) == [0,0,1])# 


def MoindresCarees(data, degree):
    """returns the polynome as a list [a0, a1, a2, ..., am]"""
    if len(data) < 2: 
        print("To aproximate you need at least 2 points nothing has been done")
        return []
    if degree > len(data) - 2: 
        degree = len(data) - 2 # car le degree du polynome de degreee max nb points -2
        print("The degree was given higher than number of points - 2")
    M_to_solve = [[0 for i in range(degree+2)] for j in range(degree+1)] # +1 to include degree and +2 to include degree plus the result
    m = degree
    n = len(data)-1
    
    list_with_calculations = [0 for i in range(0, 2*m + 1)] # +1 since 2m is included 
    for i in range(0, 2*m + 1):
        for xj, fj in data:
            list_with_calculations[i] += xj**i 
    print(list_with_calculations)
    for row in range(len(M_to_solve)):
        for col in range(len(M_to_solve)):
            M_to_solve[row][col] = list_with_calculations[row+col]
        fj_xj_multiplication_sum = 0
        for xj, fj in data:
            fj_xj_multiplication_sum += fj* (xj**row) 
        M_to_solve[row][-1] = fj_xj_multiplication_sum
    solution = ResolutionGauss(M_to_solve)
    return solution

test_data = [(-2, 3), (0, 5), (4,8), (6, 5)]
assert(MoindresCarees(test_data, 1)==[4.55, 0.35])

def S_degree_moindre_carre_calcul(data, sol_moindres_carees):
    error = 0
    for xj, fj in data:
        error += ((ListToFunction(sol_moindres_carees,xj))-fj)**2
    return error

assert(S_degree_moindre_carre_calcul(test_data, MoindresCarees(test_data, 1)) == 7.850000000000001)