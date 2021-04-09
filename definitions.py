#Dependencies.

import copy
import math
import networkx as nx
import itertools as it
from sympy.solvers import solve
from sympy.solvers import solveset
from sympy import Symbol
import multiprocessing
from joblib import Parallel, delayed
from tqdm import tqdm

#The class for the grahps.
#
#The distinguished vertex should be 0 and its neighbours 1,2,3,4 in counter-clockwise order.
#
#As in plantri, G[i] contains the vertices andyacent to vertex i in counter-clockwise order. If i is 1,2,3 or 4, G[i] should start with i-1,0,i+1.
#
#The faces should begin with the outer faces 012, 023, 034, 041, followed by the corner faces containing 12, 23, 34, 41 in that order. Then the faces containing sides 1,2,3 or 4.
#
#The corner faces should begin with 21, 32, 43 or 14.
#
#faces_containing_vertex[i] gives the faces containing vertex i in counter-clockwisse order. If i is 1,2,3,4, then it should begin with the 2 outer faces.

class Tiling_graph:
    def __init__(self, G):
        self.G = G
        self.faces = []
        self.facesContainingVertex = [[] for g in self.G]
    
    def ady(self,i,j): #Decides if tile i and j are adyacent
        if i in self.G[j]:
            return True
        return False
    
    def reorder(self,P): #Reorders labels according to the permutation P, P[i] gets relabeled to i
        GG = []
        for i in P:
            VV = [P.index(j) for j in self.G[i]]
            GG.append(VV)
        self.G = copy.deepcopy(GG)
        for i in range(len(self.G)):
            for k in range(5):
                if k in self.G[i]:
                    ind = self.G[i].index(k)
                    self.G[i] = self.G[i][ind:] + self.G[i][:ind]
                    break
        for i in range(1,5):
            ind = self.G[i].index(0)
            self.G[i] = self.G[i][ind-1:] + self.G[i][:ind-1]
        self.build_faces()
        
    def is_pyramid(self,v): #Determines if the vertex v is the apex of a square pyramid
        if len(self.G[v]) != 4:
            return False
        for i in range(4):
            j = (i+1)%4
            if not self.ady(self.G[v][i],self.G[v][j]):
                return False
        return True
        
    def pyramid_apexes(self): #Lists all vertices that serve as pyramid apex
        gV = []
        for i,d in enumerate([len(g) for g in self.G]):
            if d != 4:
                continue
            if self.is_pyramid(i):
                gV.append(i)
        return gV

    def nx_graph(self): #returns a networkx graph
        H=nx.Graph()
        H.add_nodes_from(range(len(self.G)))
        for i,g in enumerate(self.G):
            for j in range(len(self.G)):
                if i==j:
                    continue
                if j in g:
                    H.add_edge(i,j)
        H.nodes[0]['0']='0'
        return H
    
    def is_0iso(self,L): #Decides if there is an isomorphism which fixes 0 between self and a graph in L
        for H in L:
            G = self.nx_graph()
            if nx.algorithms.isomorphism.is_isomorphic(G, H.nx_graph(),
                                                       node_match = nx.algorithms.isomorphism.categorical_node_match('0',None)):
                return True
        return False
    
    def build_faces(self): #Builds the faces of G
        Marks = [[False]*len(g) for g in self.G]
        OuterFaces = [None]*4
        CornerFaces = [None]*4
        OtherFaces = []
        for i,g in enumerate(self.G):
            for j in g: #ij is a directed edge
                if Marks[i][self.G[i].index(j)]:
                    continue
                face = [i]
                old = i
                new = j
                Marks[old][self.G[old].index(new)] = True #puts a mark on ij
                while new != i:
                    face.append(new)
                    ind = self.G[new].index(old)
                    old = new
                    new = self.G[old][ind-1]
                    Marks[old][self.G[old].index(new)] = True
                AddToOtherFace = True
                if 0 in face:
                    for k in range(4):
                        if 1+k in face and 1+(1+k)%4 in face:
                            OuterFaces[k] = face
                            AddToOtherFace = False
                            break
                else:
                    for k in range(4):
                        if 1+k in face and 1+(1+k)%4 in face:
                            CornerFaces[k] = face
                            AddToOtherFace = False
                            break
                if AddToOtherFace:
                    OtherFaces.append(face)
        #Rearrange CornerFaces
        for i in range(4):
            ind = CornerFaces[i].index(i+1)
            CornerFaces[i] = CornerFaces[i][ind-1:] + CornerFaces[i][:ind-1]
        self.faces = OuterFaces + CornerFaces + OtherFaces
        #Update facesContainingVertex
        self.facesContainingVertex = [[None]*len(g) for g in self.G]
        for i,F in enumerate(self.faces):
            for j,v in enumerate(F):
                prev = F[j-1]
                ind = self.G[v].index(prev)
                self.facesContainingVertex[v][ind] = i
        for i in range(1,5):
            ind = self.facesContainingVertex[i].index((i-1)%4)
            self.facesContainingVertex[i] = self.facesContainingVertex[i][ind:] + self.facesContainingVertex[i][:ind]
    
    def draw2(self): #Draws using networkx.draw_planar
        G=self.nx_graph()
        G.remove_node(0)
        nx.draw_planar(G, with_labels = True)
        
    def draw(self): #Draws using networkx's spring:layout starting with a fixed square.
        G=self.nx_graph()
        G.remove_node(0)
        pos = {}
        pos[1]=(0,0); pos[2]=(0,1); pos[3]=(1,1); pos[4]=(1,0)
        for i in range(5,len(self.G)):
            pos[i]=(0.5,0.5)
        nx.draw(G,pos=nx.spring_layout(G,pos=pos,fixed=[1,2,3,4],iterations=7),with_labels=True)

#Functions to load and save data.

def load_plantri(filename): #Returns a list of all good matrices witn N+5 vertices
    Data = []
    string = ''
    with open(filename,'r', newline='') as file:
        for l in file:
            string += l
    assert string[:15] == '>>planar_code<<'
    lines = string[15:].split(chr(0))[:-1]
    first_line_of_graph = True
    G = []
    num_verts = 0
    for l in lines:
        if first_line_of_graph:
            G = []
            num_verts = ord(l[0])
            first_line_of_graph = False
            l=l[1:]            
        G.append([ord(c)-1 for c in l])
        if len(G) == num_verts:
            Data.append(Tiling_graph(G))
            first_line_of_graph = True
    return Data

def save_data(Data,filename):
    with open(filename,'w') as file:
        for G in Data:
            file.write(str(len(G.G))+'\n')
            for v in G.G:
                file.write(str(v)[1:-1]+'\n')

def load_data(filename):
    Data = []
    with open(filename,'r') as file:
        First_line = True
        N = 0
        G = []
        for l in file:
            if First_line:
                N=int(l)
                First_line = False
                G = []
                continue
            G.append([int(i) for i in l.split(',')])
            if len(G) == N:
                PG = Tiling_graph(G)
                PG.build_faces()
                Data.append(PG)
                First_line = True
    return Data

#A function that identifies the possible distinguished vertices of a graph and returns a list of the graphs with this distinguished vertex.

def add_distinguished(G): #Constructs the grahps with an apex labeled as 0
    graphs = []
    for v in G.pyramid_apexes():
        P = [] #Construct permutation with apex and neighbours first
        P.append(v)
        for w in G.G[v]:
            P.append(w)
        [P.append(ii) for ii in range(len(G.G)) if ii!=v and ii not in G.G[v]]
        G0 = copy.deepcopy(G)
        G0.reorder(P)
        if not G0.is_0iso(graphs):
            graphs.append(G0)
    return graphs

#Functions for filtering the list of graphs.

def first_filter(G):
    sideTiles = [None]*4
    for i in range(4):
        sideTiles[i] = {t for j in G.facesContainingVertex[1+i] for t in G.faces[j]}.difference({0,1,2,3,4})
    for i in range(4):
        problems = sideTiles[i].intersection(sideTiles[(i+1)%4]).difference(G.faces[4+i][2:])
        if len(problems) >= 1:
            return False
    return True

def second_filter(G):
    for i in range(2):
        side1 = G.G[1+i]
        side3 = G.G[3+i]
        problems = set(side1).intersection(side3).difference({0,1,2,3,4})
        if len(problems) >= 1: 
            #print(list(problems))
            return False
    return True

def has_min_deg_5(G):
    mindeg = min([len(g) for g in G.G[5:]])
    if mindeg >=5:
        return True
    return False

#Auxiliary functions for the angle search

def is_corner(f,G):
    F = G.faces[f]
    if F[1]<=4:
        return True
    return False

def is_side(f,G):
    F = G.faces[f]
    if F[0]<=4 and not is_corner(f,G):
        return True
    return False

#Returns faces f1,f2 such that (v,f1), (v,f), (v,f2) are consecutive angles
def adj_angles(v,f,G):
    Fs = G.facesContainingVertex[v]
    ind = Fs.index(f)
    return ((Fs[ind-1],Fs[(ind+1)%len(Fs)]))

#Tile containing opposite side
def adj_side(v,f1,f2,G):
    F1 = G.faces[f1]
    F2 = G.faces[f2]
    ind1 = F1.index(v)
    ind2 = F2.index(v)
    if F1[ind1-1] == F2[(ind2+1)%len(F2)]:
        return F1[ind1-1]
    elif F2[ind2-1] == F1[(ind1+1)%len(F1)]:
        return F2[ind2-1]

def clean_solution(Solution, PermSize):
    newSol = str(Solution).split(',')
    newSol[0] = newSol[0][1:]
    newSol[-1] = newSol[-1][:-1]
    R = ['Not Solved']*PermSize
    for S in newSol:
        Var, Res = S.split(':')
        if(Var[0] == ' '):
            Var = Var[1:]
        Res = Res[1:]
        if not ('x' in Res):
            i = int(Var[1:])
            R[i] = float(Res)
    return R
    
def construct_angles(G):
    visitedTiles = set()
    Angles = []
    #Corner Angles
    for f,F in enumerate(G.faces[4:8],4):
        for v in F:
            if v >= 5:
                Angles.append((v,f))
                #visitedTiles.add(v)
    #Side Angles (adjacent, not just touching)
    for side in range(1,5):
        for i,v in enumerate(G.G[side][3:],3):
            if i != 3:
                Angles.append((v,G.facesContainingVertex[side][i-1]))
            if i != len(G.G[side])-1:
                Angles.append((v,G.facesContainingVertex[side][i]))
            visitedTiles.add(v)
    #Inner Angles
    for v in visitedTiles:
        Angles += [(v,f) for f in G.facesContainingVertex[v] if (v,f) not in Angles]
    for v in set(range(5,len(G.G))).difference(visitedTiles):
        Angles += [(v,f) for f in G.facesContainingVertex[v]]
    #See when each side and face was filled
    IndexFilledSide = [0]*4
    IndexFilledFace = [0]*len(G.faces)
    for i,(v,f) in enumerate(Angles):
        for side in range(1,5):
            if side in G.G[v] and side in G.faces[f]:
                IndexFilledSide[side-1] = i
            IndexFilledFace[f] = i
    return Angles, set(IndexFilledSide), set(IndexFilledFace[4:])

#Functions to construct, save and load the angle permutations

def cyclic_perms(L):
    Res = []
    for i in range(len(L)):
        Res.append(L[i:]+L[:i])
    L.reverse()
    for i in range(len(L)):
        Res.append(L[i:]+L[:i])
    return Res

def have_common_data(L1, L2): 
    for i in L1:
        for j in L2:
            if i == j:
                return True
    return False

def has_consec_ones(L):
    for i in range(len(L)):
        if L[i-1] == L[i] == 1:
            return True
    return False

def angle_perms(numSides):
    MinAngle = [1,90,91]
    MaxAngle = [89,90,179]
    AngPerms = []
    for I in it.product(range(3),repeat=numSides):
        if numSides == 4 and has_consec_ones(I):
            continue
        P=cyclic_perms(list(I))
        if have_common_data(AngPerms,P):
            continue
        if sum([MinAngle[i] for i in I]) > (numSides - 2) *180: #Aquí teníamos desigualdades no estrictas, ¿por qué funcionaba?
            continue
        if sum([MaxAngle[i] for i in I]) < (numSides - 2) *180:
            continue
        AngPerms.append(list(I))
    return AngPerms

def save_ang_perms(numSides,filename):
    AngPerms = angle_perms(numSides)
    with open(filename,'w') as file:
        for p in AngPerms:
            S = ''
            for i in p:
                S += f'{i},'
            file.write(S[:-1] + '\n')

def load_ang_perms(filename):
    AngPerms = []
    with open(filename,'r') as file:
        for l in file:
            AngPerms.append([int(i) for i in l.split(',')])
    return AngPerms

#The class for the nodes
#angle-types: a -> 0, r -> 1, o -> 2, p->3

class NodeState:
    #G: graph G
    #P: the angle-type permutation
    def __init__(self, ig, ip, G, P, Angles,PrintProof = False):
        self.iP = ip #index of angle-type permutation
        self.iG = ig #index of G
        
        self.PrintProof = PrintProof
        
        self.AngleEqs = ['x0+x1+x2+x3-360.0']
        self.SideEqs = []
        
        self.LastAngSol = []
        self.LastSideSol = []
        
        self.indexAngle = 0
        
        #Assignation of the angles for each tile
        self.Assigned = {a:-1 for a in Angles}
        self.AssignedAT = {a:-1 for a in Angles}
        
        #ToAssign[v] contains the set of possible assignations. len(P) means that it need a plan angle.
        self.ToAssign = [{i for i in range(len(P))} for T in G.G]
        for i,T in enumerate(G.G):
            if len(T) > len(P):
                self.ToAssign[i].add(len(P))
        
        #ToAssignAT counts how many angles of each angle type have not been assigned
        at0 = P.count(0); at1 = P.count(1); at2 = P.count(2)
        self.ToAssignAT = [[at0,at1,at2,len(V)-at0-at1-at2] for V in G.G]
        
        #The difference between the total face angle and the maximum and minimum possible angle sum in that face
        self.MaxFaceSums = [0] * len(G.faces)
        self.MinFaceSums = [0] * len(G.faces)
        for f,F in enumerate(G.faces):
            if is_corner(f,G):
                angsum = 90
            elif is_side(f,G):
                angsum = 180
            else:
                angsum = 360
            self.MaxFaceSums[f] -= angsum
            self.MinFaceSums[f] -= angsum
            maxsum = 0
            minsum = 0
            for v1 in F:
                if v1 < 5:
                    continue
                self.MaxFaceSums[f] += 180
                self.MinFaceSums[f] += 1
    
    #-----------------
    def next(self, Angles, P, G, IndexFilledSide, IndexFilledFace):
        Res = [] #list of branches
        
        v,f = Angles[self.indexAngle]
        
        if self.PrintProof:
            print(f'Assigning angle {v}, {G.faces[f]}.')
        
        PossibleAssignments = self.poss_angs(Angles,G,P)
        
        if self.PrintProof:
            print(f'Possible assignments: {PossibleAssignments}')
        
        for assign in PossibleAssignments:
            NewBranch = copy.deepcopy(self)
            NewBranch.indexAngle += 1
            NewBranch.Assigned[(v,f)] = assign
            maxangle = {0:89,1:90,2:179,3:180,-1:180}
            minangle = {0:1,1:90,2:91,3:180,-1:1}
            if assign != len(P):
                NewBranch.ToAssignAT[v][P[assign]] -= 1
                NewBranch.ToAssign[v].remove(assign)
                self.MaxFaceSums[f] = self.MaxFaceSums[f] - maxangle[-1] + maxangle[P[assign]]
                self.MinFaceSums[f] = self.MinFaceSums[f] - minangle[-1] + minangle[P[assign]]
            else:
                NewBranch.ToAssignAT[v][3] -= 1
                if NewBranch.ToAssignAT[v][3]==0:
                    NewBranch.ToAssign[v].remove(len(P))
                self.MaxFaceSums[f] = self.MaxFaceSums[f] - maxangle[-1] + maxangle[3]
                self.MinFaceSums[f] = self.MinFaceSums[f] - minangle[-1] + minangle[3]
            if self.indexAngle in IndexFilledFace:
                if not NewBranch.add_face_eq(f,G,P):
                    if self.PrintProof:
                        print(f'Assignmment {assign} failed face equations.')
                    continue
                NewBranch.FaceFilled = False
            FilledSides = self.filled_sides(Angles,IndexFilledSide,G,P)
            if len(FilledSides) > 0:
                if not NewBranch.add_side_eq(FilledSides,G,P):
                    if self.PrintProof:
                        print(f'Assignmment {assign} failed side equations.')
                    continue
            Res.append(NewBranch)
        if self.PrintProof:
            print('---')
        return Res

    #-----------------
    def poss_angs(self,Angles,G,P):
        Res = []
        v,f = Angles[self.indexAngle]
        for assign in self.ToAssign[v]:
            if assign < len(P):
                at = P[assign]
            else:
                at = 3
            if self.PrintProof:
                print(f'Trying to assign {assign}, angle-type {at}')
            if self.check_tile(v,f,assign,at,G,P) and self.check_face(v,f,at,G,P):
                Res.append(assign)
        return Res    

    #-----------------
    def check_tile(self,v,f,assign,at,G,P):
        if at == 3:
            return True
        if len(P) == sum(self.ToAssignAT[v][:3]):
            return True
        indf = G.facesContainingVertex[v].index(f)
        i = (indf+1)%len(G.facesContainingVertex[v])
        while True:
            f1 = G.facesContainingVertex[v][i]
            a = self.Assigned[v,f1]
            if 0 <= a < len(P):
                if (assign - a)%len(P) not in {1,len(P)-1}:
                    if self.PrintProof:
                        print(f'Failed tile test.')
                    return False
                break
            if a == -1:
                break
            i = (i+1)%len(G.facesContainingVertex[v])
        i = (indf-1)%len(G.facesContainingVertex[v])
        while True:
            f1 = G.facesContainingVertex[v][i]
            a = self.Assigned[v,f1]
            if 0 <= a < len(P):
                if (assign - a)%len(P) not in {1,len(P)-1}:
                    if self.PrintProof:
                        print(f'Failed tile test.')
                    return False
                break
            if a == -1:
                break
            i = (i-1)%len(G.facesContainingVertex[v])
        if self.PrintProof:
            print(f'Passed tile test.')
        return True
    
    #-----------------
    def check_face(self,v,f,at,G,P):
        #This works because we have few faces and sides
        maxangle = {0:89,1:90,2:179,3:180,-1:180}
        minangle = {0:1,1:90,2:91,3:180,-1:1}
        MaxSum = self.MaxFaceSums[f] - maxangle[-1] + maxangle[at]
        MinSum = self.MinFaceSums[f] - minangle[-1] + minangle[at]
        if MinSum > 0:
            if self.PrintProof:
                print(f'Angle-type sum too large by at least {MinSum}')
            return False
        if MaxSum < 0:
            if self.PrintProof:
                print(f'Angle-type sum too small by at least {-MaxSum}')
            return False
        if self.PrintProof:
            print('Passed face test.')
        return True
    
    #-----------------
    def filled_sides(self,Angles,IndexFilledSide,G,P):
        Res = []
        v,f = Angles[self.indexAngle]
        f1,f2 = adj_angles(v,f,G)
        v1 = adj_side(v,f,f1,G)
        v2 = adj_side(v,f,f2,G)
        for vv,ff in [(v1,f1),(v2,f2)]:
            if 0 <= self.Assigned[v,ff] < len(P):
                if vv <= 4:
                    if self.indexAngle in IndexFilledSide:
                        Res.append((vv,))
                    continue
                if 0 <= self.Assigned[vv,f] < len(P) and 0 <= self.Assigned[vv,ff] < len(P):
                    Res.append((v,vv,f,ff))
        return Res
    
    #-----------------
    def add_face_eq(self, f, G, P):
        NewEq = ''
        straightCounter = 0
        for v in G.faces[f]:
            if v < 5:
                continue
            if self.Assigned[v,f] < len(P):
                NewEq += f'+x{self.Assigned[v,f]}'
            else:
                straightCounter += 1
        if is_corner(f,G): #Corner
            NewEq += '-90'
        elif is_side(f,G): #Side
            NewEq += '-180'
        else: #Inside
            NewEq += f'-{360 - straightCounter*180}'
        self.AngleEqs.append(NewEq)
        if self.PrintProof:
            print(f'New face equation: {NewEq}')
        #Checks if the system can still be solved. Or if it's solved with a non acceptable answer. 
        eqSolution = solve(self.AngleEqs, rational=None)
        if self.PrintProof:
            print(f'Reduced equations: {eqSolution}')
        if len(eqSolution) == 0:
            return False
        #if it has a possible solution yet. Checks if it's acceptable.
        possiblyRes = clean_solution(eqSolution, len(P))
        for i in range(len(possiblyRes)): #Checks for each variable (-1:= no exact answer)
            if possiblyRes[i] == 'Not Solved':
                continue
            if P[i] == 0:
                if possiblyRes[i] >= 90 or possiblyRes[i] <= 0:
                    return False
            if P[i] == 1:
                if not possiblyRes[i] == 90:
                    return False
            if P[i] == 2:
                if possiblyRes[i] <= 90 or possiblyRes[i] >=180:
                    return False
        self.LastAngSol = possiblyRes
        return True
    
    #-----------------
    #Side 0 is between angles 0 and 1
    def side_var(self,v,f1,f2,P):
        i = min(self.Assigned[v,f1],self.Assigned[v,f2])
        if i == 0 and max(self.Assigned[v,f1],self.Assigned[v,f2]) == len(P)-1:
             i = len(P)-1
        return f'x{i}'
    
    #-----------------
    def add_side_eq(self, sides, G, P):
        for s in sides:
            #Eqs for a side of the square
            if len(s)==1:
                side = s[0]
                NewEq = ''
                for i,v in enumerate(G.G[side][3:],3):
                    f1 = G.facesContainingVertex[side][i-1]
                    f2 = G.facesContainingVertex[side][i]
                    NewEq += f'+{self.side_var(v,f1,f2,P)}'
                NewEq += '-1.0'
                self.SideEqs.append(NewEq)
            else:
                v1,v2,f1,f2 = s
                NewEq = f'{self.side_var(v1,f1,f2,P)}-{self.side_var(v2,f1,f2,P)}'
                self.SideEqs.append(NewEq)
        if self.PrintProof:
            print(f'New side equation: {NewEq}')
        #Checks if the system can still be solved. Or if it's solved with a non acceptable answer. 
        eqSolution = solve(self.SideEqs, rational=None)
        if self.PrintProof:
            print(f'Reduced equations: {eqSolution}')
        if len(eqSolution) == 0:
            return False
        #if it has a possible solution yet. Checks if it's acceptable.
        possiblyRes = clean_solution(eqSolution, len(P))
        sidesSolved = 0
        for i in range(len(possiblyRes)): #Checks for each variable.
            if possiblyRes[i] == 'Not Solved':
                continue
            if possiblyRes[i] > 1 or possiblyRes[i] <= 0:
                    return False
            sidesSolved += 1
        self.LastSideSol = possiblyRes
        
        if len(P) == 4: #Only works for 4-side tiles.
            if sidesSolved == len(P):
                DiagonalsToCheck = []
                if self.LastAngSol[0] != 'Not Solved' and self.LastAngSol[2] != 'Not Solved':
                    DiagonalsToCheck.append([0,2])
                if self.LastAngSol[1] != 'Not Solved' and self.LastAngSol[3] != 'Not Solved':
                    DiagonalsToCheck.append([1,3])
                for AngPair in DiagonalsToCheck:
                    a1 = self.LastSideSol[AngPair[0]-1]
                    b1 = self.LastSideSol[AngPair[0]]
                    alpha1 = self.LastAngSol[AngPair[0]]
                    d1 = a1*a1 + b1*b1 - 2*a1*b1*math.cos(math.radians(alpha1))
                    
                    a2 = self.LastSideSol[AngPair[1]-1]
                    b2 = self.LastSideSol[AngPair[1]]
                    alpha2 = self.LastAngSol[AngPair[1]]
                    d2 = a2*a2 + b2*b2 - 2*a2*b2*math.cos(math.radians(alpha2))
                    
                    if self.PrintProof:
                        print(f'Sides: {self.LastSideSol}')
                        print(f'Angles: {self.LastAngSol}')
                    
                    if not math.isclose(d1,d2,rel_tol = 1e-3):
                        if self.PrintProof:
                            print(f'Diagonals suqared did not match: {d1,d2}')
                        return False
                    
                    if self.PrintProof:
                        print(f'Passed diagonal test with: {d1,d2}')
                    
                    A1 = a1*b1*math.sin(math.radians(alpha1))/2
                    A2 = a2*b2*math.sin(math.radians(alpha2))/2
                    A = A1 + A2
                    if not math.isclose(A, 1.0/(len(G.G) - 5), rel_tol = 1e-3):
                        if self.PrintProof:
                            print(f'Area is not 1/{len(G.G) - 5}: {A}')
                        return False
                    
                    if self.PrintProof:
                        print(f'Passed area test: {A}')
        return True

#The function which explores the tree

def search(ig, G, numSides,PrintProof = False): #index of graph, graph, number of sides
    Res = [] #Contains final result
    Perms = load_ang_perms(str(numSides)+'_perms.txt') #List of possible angle-types a tile may have
    
    #An angle is a pair (v,F) where v is a tile and F is a face of G containing v
    Angles, IndexFilledSide, IndexFilledFace = construct_angles(G)
    
    for ip, P in enumerate(Perms):
        if PrintProof:
            print(f'\nPermutation: {P}\n')
        Stack = [NodeState(ig, ip, G, P, Angles,PrintProof)]
        while len(Stack) > 0:
            Branch = Stack.pop()
            if(Branch.indexAngle == len(Angles)):
                Res.append(Branch)
                if PrintProof:
                    print(f'\nFOUND A POSSIBLE SOLUTION!\n')
                continue
            NewBranches = Branch.next(Angles, P, G, IndexFilledSide, IndexFilledFace)
            Stack.extend(NewBranches)
        if PrintProof:
            print(f'\nFound {len(Res)} possible solutions.\n')
    return Res
