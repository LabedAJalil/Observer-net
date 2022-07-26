# -*- coding: utf-8 -*-
"""
@author: JALIL
"""


import time
import numpy as np 
import networkx as nx


# n transitions
n = 8
# m places
m = 6

Post = np.array([
[0,0,0,0,0,0,0,2],
[1,0,0,0,0,0,0,0],
[0,1,0,0,0,0,0,0],
[0,0,1,0,0,0,1,0],
[0,0,0,1,0,0,0,0],
[0,0,0,0,1,1,0,0]])

Pre = np.array([
[1,1,0,0,0,0,0,0],
[0,0,1,0,0,0,0,0],
[0,0,0,1,0,0,0,0],
[0,0,0,0,1,0,0,0],
[0,0,0,0,0,1,0,0],
[0,0,0,0,0,0,1,2]])

# initial marking
M = np.array([120,0,0,0,0,0])

# label for observable transition  \u03B5
transition_labeling = {'t(1)':'a','t(2)':'0','t(3)':'0','t(4)':'0',
                    't(5)':'b','t(6)':'a','t(7)':'0','t(8)':'b'}

transitions_with_label = {'a':[0,5],'b':[4,7],'None':[1,2,3,6],'obs':[0,4,5,7]}

N = {'Places':m,'Transitions':n,'Pre':Pre,'Post':Post}
LPN = {'N':N, 'M0': M, 'E': transition_labeling, 'transition': transitions_with_label}


def name(mphi):
    SUB = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")
    for m in mphi:
        s = ""
        for i,p in enumerate(m):
            if p >0:
                t = "P"+str(i+1).translate(SUB)
                if p ==1:
                    s += t+"+"
                else: 
                    s += str(p)+t+"+"
        print(s[:-1])
               
# labeling function
def l(t):
    return transition_labeling.get('t('+str(t+1)+')')

def T(e):
    return transitions_with_label.get(e)

# verifiyer if transition t is enabled at M 
def tras_enabled(M, t):
    return np.greater_equal(M, Pre[:,t]).all()

# get the set of all transistion enabled at M 
def setof_enabled(M):
    set_tran = []
    for t in range(n):
        if tras_enabled(M,t): set_tran.append(True)
        else: set_tran.append(False)
    return set_tran

# fired t if transition t is enabled at M 
def tras_fired(M, t):
    C =N['Post']-N['Pre']
    return (M +  C[:,t])

def setof_unt_enabled(M):
    set_tran = []
    for t in T('None'):
        if tras_enabled(M,t): set_tran.append(t)
    return set_tran

# the set of markings reachable from M by firing only unobservable transitions.
def unobservable_reach(M):
    # set of unobservable_reach 
    U=[]
    U.append(M)
    for m in U:
        for t in setof_unt_enabled(m):
            # nM new marking
            nM = tras_fired(m,t)   
            if not any((nM == m).all() for m in U): 
                U.append(nM)
    return U

# enabled  transitions  at MiΦ with  label e
#(Sigma.append(t) for t in observable_transitions if l(t)==e and tras_enabled(M, t))
def enabled_transitions(M,e):  
    Sigma = []
    for t in T(e):
        if tras_enabled(M, t):
            Sigma.append(t)
    return Sigma

# feasible transitions at MΦ with label e
def feasible_transitions(M_phi,e):
    SIGMA = []
    for m in M_phi:
        if enabled_transitions(m,e):
            return True
    return False
    #     if enabled_transitions(m,e): SIGMA.extend(enabled_transitions(m,e))
    # return SIGMA
    
# set of all reachable markings, when all feasible transitions in Σ(MΦ,e) are fired from every MiΦ in MΦ
def reachable_markings(M_phi,e):
    MU = []
    for m in M_phi:
        for t in enabled_transitions(m,e):
            nM=tras_fired(m,t)
            if not any((nM == m_).all() for m_ in MU):
                MU.append(nM)
    return MU

# The checker state transition function of Given the checker markingMΦ and the event e
def state_transition_function(M_phi,e):
    STATE = []
    for m in reachable_markings(M_phi,e):
        STATE.extend(unobservable_reach(m))
    if len(STATE)>1:
        STATE = np.unique(STATE, axis=0) 
    return STATE

def check_word(M,word):
    M_phi_0 = unobservable_reach(M)
    M_phi = M_phi_0
    for e in word:
        M_phi = state_transition_function(M_phi,e)

def exist(CM,M):
    for mphi in CM:
        if mphi.shape == M.shape:
            if np.equal(mphi,M).all(): 
                return True
    return False

def checker_states(M):
    CM= []
    event = ['a','b']
    M_phi_0 = np.asarray(unobservable_reach(M))
    CM.append(M_phi_0)
    for M_phi in CM:
        for e in event:
            if feasible_transitions(M_phi, e):
                nM_phi=state_transition_function(M_phi,e)
                if not exist(CM, nM_phi):
                    CM.append(nM_phi)
                
    print("STATE",len(CM))

def rechability(M):
    R=[]
    R.append(M)
    RG = nx.DiGraph()
    RGdict = {}
    path = {}
    for m in R:
        for t in range(n):
            if tras_enabled(m, t):
                Mn=tras_fired(m, t)
                if not any((Mn == ma).all() for ma in R):
                    R.append(Mn)
                RG.add_edge(str(m),str(Mn),label=l(t))
                if l(t) in path:
                    path[l(t)].append(str(Mn))
                else:
                    path[l(t)] = [str(Mn)]
              
        RGdict[str(m)] = path
        path = {}
    # dic1 = nx.to_dict_of_lists(RG,nodelist=None)
    # for k,v in RGdict.items():
    #     print(k,v)
    return RGdict

def existNFA(S,q_):
    for q in S:
        if q == q_:
            return True
    return False

def closure(nfa,q):
    S = []
    S.append(q)
    for s in S: 
        if "0" in nfa[s]:
            l = (nfa[s]).get("0")
            for q_ in l: 
                if not existNFA(S,q_):
                    S.append(q_)
    return S
    
def existObs(STATES,Q):
    for S in STATES:            
        if len(S) == len(Q):
            i = 0
            for q in S:
                for q_ in Q:
                    if q == q_:
                        i += 1
            if(i==len(S)):
                return True
            else:
                i = 0 
    return False

def NFA_TO_DFA(nfa):
    STATES = []
    event = ['a','b',]
    q0 = list(nfa.keys())[0] 
    STATES.append(closure(nfa,q0))
    for Q in STATES:
        for e in event:
            Q_ = move_event(nfa, Q, e)
            if not existObs(STATES,Q_) and Q_:
                STATES.append(Q_)
    print("STATE",len(STATES))
                            
def move_event(nfa, STATE, e):
    S = []
    for q in STATE:
        if e in nfa[q]:
            states = nfa[q].get(e)
            for q_ in states:
                for cq in closure(nfa,q_):
                    if not existNFA(S, cq):
                        S.append(cq)
    return S   
    
def Ymin(M,N,t):
    A = np.array([(M-N['Pre'][:,t]).transpose()])
    B = np.array([np.zeros(len(T('None')),dtype=int).transpose()])
    C = N['Post']-N['Pre']
    Cu = (C[:,T('None')]).transpose()
    Id = np.identity(len(T('None')))
    
    while A[A<0].size != 0 :
        i_s = np.where(A<0)[0][0]
        j_s = np.where(A<0)[1][0]
        I = [ indx for indx,row in enumerate(Cu[:,j_s])if row>0]
        for i in I:
            A = np.vstack([A, A[i_s]+Cu[i]])
            B = np.vstack([B, B[i_s]+(Id[i]).astype(int)])
        A = np.delete(A, i_s, 0)
        B = np.delete(B, i_s, 0)
    return B

def basis_rechability(M,N):
    BR=[] # Basis Reachability set 
    C = N['Post']-N['Pre'] # incidence matrix
    Cu = C[:,T('None')]
    BR.append(M) # insialization of basisi reachability set add the first basis marking to the set 
    # BRG = nx.DiGraph() # basis reachability graph 
    BRGdict = {}
    path = {}

    for Mb in BR:
        for t in T("obs"):
            if Ymin(Mb,N,t).size != 0:
                for yu in Ymin(Mb,N,t):
                    Mn = Mb + Cu.dot(yu) + C[:,t]
                    if not any((Mn == ma).all() for ma in BR):
                        BR.append(Mn)                    
                    if l(t) in path:
                        path[l(t)].append(str(Mn))
                    else:
                        path[l(t)] = [str(Mn)]
                  
                BRGdict[str(Mb)] = path
        path = {}
                
    # for k,v in BRGdict.items():
    #     print("Mb",k)
    #     print(v)
    print("States in BR : ",len(BR))
    return BRGdict


K = [np.array([2,0,0,0,0,0]),  np.array([4,0,0,0,0,0]), 
      np.array([8,0,0,0,0,0]),  np.array([10,0,0,0,0,0]),
      np.array([20,0,0,0,0,0]), np.array([40,0,0,0,0,0]),  np.array([60,0,0,0,0,0]),]



############################################################ FIRST TEST ############################################################ 
for M in K:
    print("---------------------------------------------")
    print(M)
    print("---------------------------------------------")
    
    print('Observer net')
    start = time.time()
    checker_states(M)
    check_word(M,"aba")
    end = time.time()
    print(end - start)
    
    print('Reachability')
    start = time.time()  
    nfa=rechability(M)
    NFA_TO_DFA(nfa)
    end = time.time()           
    print(end - start)
    
    print('Basis Reachability')
    start = time.time()  
    nfa=basis_rechability(M, N)
    end = time.time()   
    
    print("basis time",end-start)               
    NFA_TO_DFA(nfa)
    end = time.time()           
    print(end - start)
    
############################################################ SECOND TEST  ############################################################
    
start = time.time()
print(len('a'))
check_word(M,"a")
end = time.time()
print(end - start)


start = time.time()
print(len('abababababababababababababababababababab'))
check_word(M,"abababababababababababababababababababab")
end = time.time()
print(end - start)
    
start = time.time()
print(len('abababababababababababababababababababababababababababababab'))
check_word(M,"abababababababababababababababababababababababababababababab")
end = time.time()
print(end - start)
    
    
start = time.time()
print(len('abababababababababababababababababababababababababababababababababababababababab'))
check_word(M,"abababababababababababababababababababababababababababababababababababababababab")
end = time.time()
print(end - start)
    
start = time.time()
print(len('abababababababababababababababababababababababababababababababababababababababababababababababababab'))
check_word(M,"abababababababababababababababababababababababababababababababababababababababababababababababababab")
end = time.time()
print(end - start)
    
start = time.time()
print(len('abababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababab'))
check_word(M,"abababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababababab")
end = time.time()
print(end - start)
