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
import matplotlib.pyplot as plt

#The class for the graphs.
#
#The distinguished vertex should be 0 and its neighbours 1,2,3,4 in counter-clockwise order.
#
#As in plantri, G[i] contains the vertices andyacent to vertex i in counter-clockwise order. If i is 1,2,3 or 4, G[i] should start with i-1,0,i+1.
#
#The faces should begin with the outer faces 012, 023, 034, 041, followed by the corner faces containing 12, 23, 34, 41 in that order. Then the faces containing sides 1,2,3 or 4.
#
#The corner faces should begin with 21, 32, 43 or 14.
#
#facesContainingVertex[i] gives the faces containing vertex i in counter-clockwisse order. If i is 1,2,3,4, then it should begin with the 2 outer faces.

class Tiling_graph:
    def __init__(self, G):
        self.G = G
        self.faces = []
        self.facesContainingVertex = [[] for g in self.G]
    
    #Decides if tile i and j are adyacent
    def ady(self,i,j): #Decides if tile i and j are adyacent
        if i in self.G[j]:
            return True
        return False
    
    #Reorders labels according to the permutation P, P[i] gets relabeled to i
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
        
    #Determines if the vertex v is the center of the wheel graph W_5
    def is_pyramid(self,v): #Determines if the vertex v is the apex of a square pyramid
        if len(self.G[v]) != 4:
            return False
        for i in range(4):
            if len(self.G[self.G[v][i]]) <= 3:
                return False
            j = (i+1)%4
            if not self.ady(self.G[v][i],self.G[v][j]):
                return False
        for i in range(2):
            j = (i+2)%4
            if self.ady(self.G[v][i],self.G[v][j]):
                return False
        return True
        
    #Lists all vertices that are the center of a wheel graph W_5
    def pyramid_apexes(self): #Lists all vertices that serve as pyramid apex
        gV = []
        for i,d in enumerate([len(g) for g in self.G]):
            if self.is_pyramid(i):
                gV.append(i)
        return gV

    #returns a networkx graph
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
    
    #Decides if there is an isomorphism which fixes 0 between self and a graph in L
    def is_0iso(self,L): #Decides if there is an isomorphism which fixes 0 between self and a graph in L
        for H in L:
            G = self.nx_graph()
            if nx.algorithms.isomorphism.is_isomorphic(G, H.nx_graph(),
                                                       node_match = nx.algorithms.isomorphism.categorical_node_match('0',None)):
                return True
        return False
    
    #Builds the faces of G
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
    
    #Draws using networkx.draw_planar
    def draw2(self): #Draws using networkx.draw_planar
        G=self.nx_graph()
        G.remove_node(0)
        nx.draw_planar(G, with_labels = True)
        plt.show()
        
    #Draws using networkx's spring layout starting with a fixed square. May not be planar.
    def draw(self,iterations=7): #Draws using networkx's spring:layout starting with a fixed square.
        G=self.nx_graph()
        G.remove_node(0)
        pos = {}
        pos[1]=(0,0); pos[2]=(0,1); pos[3]=(1,1); pos[4]=(1,0)
        for i in range(5,len(self.G)):
            pos[i]=(0.5,0.5)
        nx.draw(G,pos=nx.spring_layout(G,pos=pos,fixed=[1,2,3,4],iterations=7),with_labels=True)
        plt.show()
#Class ends here

#Takes output from plantri and returns a list of graphs
def load_plantri(filename): #Returns a list of all good matrices witn N+5 vertices
    Data = []
    with open(filename,'r', newline='') as file:
        string = ''.join(file.readlines())
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

#Saves a list of graphs
def save_data(Data,filename):
    with open(filename,'w') as file:
        for G in Data:
            file.write(str(len(G.G))+'\n')
            for v in G.G:
                file.write(str(v)[1:-1]+'\n')

                
#Loads a list of graphs
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

#Identifies the possible distinguished vertices of a graph and returns a list of the graphs with the distinguished vertex labeled as 0.
def add_distinguished(G): #Constructs the graphs with an apex labeled as 0
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

#Returns False iff G contains a vertex corresponding to a tile which intersects consecutive sides but doesn't contain the vertex inbetween
def first_filter(G):
    sideTiles = [None]*4
    for i in range(4):
        sideTiles[i] = {t for j in G.facesContainingVertex[1+i] for t in G.faces[j]}.difference({0,1,2,3,4})
    for i in range(4):
        problems = sideTiles[i].intersection(sideTiles[(i+1)%4]).difference(G.faces[4+i][2:])
        if len(problems) >= 1:
            return False
    return True

#Returns False iff G contains a tile touching opposite sides
def second_filter(G):
    for i in range(2):
        side1 = G.G[1+i]
        side3 = G.G[3+i]
        problems = set(side1).intersection(side3).difference({0,1,2,3,4})
        if len(problems) >= 1: 
            #print(list(problems))
            return False
    return True

#Checks if G has min degree at least d
def has_min_deg(G, d):
    mindeg = min([len(g) for g in G.G[5:]])
    if mindeg >=d:
        return True
    return False

#Checks if face f corresponds to a vertex of the square
def is_corner(f,G):
    F = G.faces[f]
    if F[1]<=4:
        return True
    return False


#Checks if face f corresponds to a vertex of the tiling on a side of the square which is not a vertex of the square
def is_side(f,G):
    F = G.faces[f]
    if F[0]<=4 and not is_corner(f,G):
        return True
    return False

#A pair (v,f), where v is a vertex and f is a face, corresponds to the angle of a tile

#Returns faces f1,f2 such that (v,f1), (v,f), (v,f2) are consecutive angles on the tile corresponding to v
def adj_angles(v,f,G):
    Fs = G.facesContainingVertex[v]
    ind = Fs.index(f)
    return ((Fs[ind-1],Fs[(ind+1)%len(Fs)]))

#Assuming f1 and f2 determine a side of the tile corresponding to v
#Returns the tile adjacent to v which also contains f1 and f2
def adj_side(v,f1,f2,G):
    F1 = G.faces[f1]
    F2 = G.faces[f2]
    ind1 = F1.index(v)
    ind2 = F2.index(v)
    if F1[ind1-1] == F2[(ind2+1)%len(F2)]:
        return F1[ind1-1]
    elif F2[ind2-1] == F1[(ind1+1)%len(F1)]:
        return F2[ind2-1]

    
#Takes the solution of the linear solver and turns it into a list with the solutions (perhaps symbolic) to each variable
def clean_solution(Solution, PermSize):
    newSol = str(Solution).split(',')
    newSol[0] = newSol[0][1:]
    newSol[-1] = newSol[-1][:-1]
    R = ['Not Solved']*PermSize
    RectSide = 'Not Solved'
    for S in newSol:
        Var, Res = S.split(':')
        if(Var[0] == ' '):
            Var = Var[1:]
        Res = Res[1:]
        if not ('x' in Res):
            i = int(Var[1:])
            if not i == PermSize:
                R[i] = float(Res)
            else:
                RectSide = float(Res)
    return R

#Gives a list of angles (vertex-face pairs) in the order in which they will be assigned
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
    return Angles, set(IndexFilledFace[4:])

#Checks if a list has two consecutive ones, this is used to check if a tile has two consecutive right angles
def has_consec_ones(L):
    for i in range(len(L)):
        if L[i-1] == L[i] == 1:
            return True
    return False


#Generates the possible angle-type permutations
def angle_perms(numSides, use5AT = False):
    AngPerms = []
    MinAngle = [1,90,91]
    MaxAngle = [89,90,179]
    numFirstAT = 2
    if use5AT:
        MinAngle = [1,45,46,90,91]
        MaxAngle = [44,45,89,90,179]
        numFirstAT = 4
    for first in range(numFirstAT):
        for I in it.product(range(len(MinAngle)),repeat=numSides-1):
            I = [first]+list(I)
            if I in AngPerms:
                continue
            if [I[0]] + list(reversed(I[1:])) in AngPerms:
                continue
            if sum([MinAngle[i] for i in I]) > (numSides - 2) *180:
                continue
            if sum([MaxAngle[i] for i in I]) < (numSides - 2) *180:
                continue
            if I == [1]*numSides:
                continue
            AngPerms.append(list(I))
    return AngPerms

#Generates and saves in a file the possible angle-type permutations
def save_ang_perms(numSides,filename, use5AT = False):
    AngPerms = angle_perms(numSides, use5AT)
    with open(filename,'w') as file:
        for p in AngPerms:
            S = ''
            for i in p:
                S += f'{i},'
            file.write(S[:-1] + '\n')

#Loads from a file the possible angle-type permutations
def load_ang_perms(filename):
    AngPerms = []
    with open(filename,'r') as file:
        for l in file:
            AngPerms.append([int(i) for i in l.split(',')])
    return AngPerms

#Dictionary for angle types
at_dict = []
#The class for the nodes
class NodeState:
        #ig: the index of the graph
    #ip: index of angle-type permutation
    #G: graph G
    #P: the angle-type permutation
    #Angles: list of angles (vertex-face pairs) in the order in which they will be assigned
    #use5AT: If True, uses labels for angle types cutting acute in 3 different cases: (0,45), [45,45], (45,90)
    #PrintProof: If True, prints all arguments used to analyse this graph
    def __init__(self, ig, ip, G, P, Angles, use5AT = False, PrintProof = False):
        
        self.iP = ip #index of angle-type permutation
        self.iG = ig #index of G
        
        self.PrintProof = PrintProof
        self.Case = ''
        
        #Angle equations
        eq=''
        for i in range(len(P)):
            eq += f'x{i}+'
        eq = eq[:-1]+'-'+str(180*(len(P)-2))+'.0'
        self.AngleEqs = [eq]
        
        #The values of the angles of the tile T
        self.LastAngSol = ['Not Solved']*len(P)
        
        #The index of the next angle to be assigned
        self.indexAngle = 0
       
        self.numAT = 3
        if use5AT:
            self.numAT = 5
        
        #Special conditions
        self.nonRight = []
        self.miniAcute = []
        
        self.non45 = False
        #Equiangular: for triangles only
        if len(P) == 3:
            for o in range(4):
                for v in G.G[1 + o]:
                    if v < 5:
                        continue
                    cR = set(G.facesContainingVertex[v]).intersection(G.facesContainingVertex[1 + (o+2)%4])
                    isInCorner = False
                    if len(cR) > 0:
                        for cF in cR:
                            if len({1,2,3,4}.intersection(G.faces[cF])) > 1:
                                isInCorner = True
                                break
                        cR = cR.pop()
                        self.nonRight.append((v, cR))
                        if isInCorner and len(G.G[1+o]) > 4:
                            self.miniAcute.append((v,cR))
        
        #Equiangular: for quadrilateral
        if len(P) == 4:
            for neighbours in G.G[5:]:
                if len({1,2,3,4}.intersection(neighbours)) >= 3:
                    self.non45 = True
                    break
                
        
        #Assignation of the angles for each tile
        self.Assigned = {a:-1 for a in Angles}
        
        #ToAssign[v] contains the list of possible assignations. len(P) represents a plain angle.
        self.ToAssign = [{i for i in range(len(P))} for T in G.G]
        for i,T in enumerate(G.G):
            if len(T) > len(P):
                self.ToAssign[i].add(len(P))
        
        #ToAssignAT counts how many angles of each angle type have not been assigned
        ats = [P.count(i) for i in range(self.numAT)]
        self.ToAssignAT = [ats + [len(V)-sum(ats)] for V in G.G]
        
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
    
    #Creates and returns the valid successor nodes
    def next(self, Angles, P, G, IndexFilledFace):
        Res = [] #list of branches
        
        v,f = Angles[self.indexAngle]
        
        if self.PrintProof:
            if self.indexAngle > 1:
                print(f'\nCase {self.Case} has the following possibilities:')
            elif self.indexAngle == 1:
                print(f'\nThis case has the following possibilities:')
            print(f'  Assigning the angle in tile {v} incident to tiles {G.faces[f]}')
        
        PossibleAssignments = self.poss_angs(Angles,G,P)
        
        if self.PrintProof:
            print(f'  Possible assignments: {PossibleAssignments}'.replace(str(len(P)),'p'))
        
        for case,assign in enumerate(PossibleAssignments):
            
            NewBranch = copy.deepcopy(self)
            if self.PrintProof:
                if self.indexAngle > 0:
                    NewBranch.Case += f'{case+1}.'
                    print(f'\nCase {NewBranch.Case}')
                else:
                    print('')
                print(f'  Assign {assign}'.replace(str(len(P)),'p') + f' to ({v}, {G.faces[f]})')
            
            #for triangles only
            if len(P) == 3 and assign != len(P):
                if P[assign] == 3 and ((v,f) in self.nonRight):
                    if self.PrintProof:
                        print(f'      Impossible: ({v},{f}) can not have a right angle, violates Lemma 9')
                    continue
                if P[assign] != 0 and ((v,f) in self.miniAcute):
                    if self.PrintProof:
                        print(f'      Impossible: ({v},{f}) can not have an angle <= 45, violates Lemma 9')
                    continue
            
            if NewBranch.step(assign, Angles, P, G, IndexFilledFace):
                Res.append(NewBranch)
        return Res
    
    #Checks if a new node is valid
    def step(self, assign, Angles, P, G, IndexFilledFace):
        v,f = Angles[self.indexAngle]
        self.Assigned[(v,f)] = assign
        #Update MaxFaceSums and MinFaceSums
        maxangle = {0:89,1:90,2:179,3:180,-1:180}
        minangle = {0:1,1:90,2:91,3:180,-1:1}
        if self.numAT == 5:
            maxangle = {0:44, 1:45, 2:89, 3:90, 4:179, 5:180, -1:180}
            minangle = {0:1, 1:45, 2:46, 3:90, 4:91, 5:180, -1:1}
        if assign != len(P):
            self.ToAssignAT[v][P[assign]] -= 1
            self.ToAssign[v].remove(assign)
            self.MaxFaceSums[f] = self.MaxFaceSums[f] - maxangle[-1] + maxangle[P[assign]]
            self.MinFaceSums[f] = self.MinFaceSums[f] - minangle[-1] + minangle[P[assign]]
        else:
            self.ToAssignAT[v][self.numAT] -= 1
            if self.ToAssignAT[v][self.numAT]==0:
                self.ToAssign[v].remove(len(P))
            self.MaxFaceSums[f] = self.MaxFaceSums[f] - maxangle[-1] + maxangle[self.numAT]
            self.MinFaceSums[f] = self.MinFaceSums[f] - minangle[-1] + minangle[self.numAT]
        OldLastAngSol = copy.copy(self.LastAngSol)
        #Check face equations
        if self.indexAngle in IndexFilledFace:
            if not self.add_face_eq(f,G,P):
#                if self.PrintProof:
#                    print(f'Assignmment {assign} failed face equations.')
                return False
            self.FaceFilled = False
        elif self.PrintProof:
            print(f'    No new angle equation')
        self.indexAngle += 1
        return True
    
    #Assigns an angle-type to the given tiling-vertex.
    def poss_angs(self,Angles,G,P):
        Res = []
        v,f = Angles[self.indexAngle]
        for assign in self.ToAssign[v]:
            if self.indexAngle == 0:
                assign = 0
            if assign < len(P):
                at = P[assign]
                if self.PrintProof:
                    print(f'    - to vertex {assign} of T which has angle-type \'{at_dict[at]}\'')
                    if self.indexAngle == 0:
                        print(f'      Only this vertex of T is considered, cyclic permutations of this labeling are considered separately')
            else:
                at = self.numAT
                if self.PrintProof:
                    print(f'    Angle-type \'p\'')
            if self.check_tile(v,f,assign,at,G,P) and self.check_face(v,f,at,G,P):
                Res.append(assign)
            if self.indexAngle == 0:
                break
        return Res    

    #Check if the given angle type at can be successful assigned to the tiling-vertex (v,f) without breaking tile conditions. 
    def check_tile(self,v,f,assign,at,G,P):
        if at == self.numAT:
            return True
        if len(P) == sum(self.ToAssignAT[v][:self.numAT]):
            return True
        indf = G.facesContainingVertex[v].index(f)
        i = (indf+1)%len(G.facesContainingVertex[v])
        while True:
            f1 = G.facesContainingVertex[v][i]
            a = self.Assigned[v,f1]
            if 0 <= a < len(P):
                if (assign - a)%len(P) not in {1,len(P)-1}:
                    if self.PrintProof:
                        print(f'      Impossible: The angles of this tile would be in the wrong order')
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
                        print(f'      The angles of this tile would be in the wrong order')
                    return False
                break
            if a == -1:
                break
            i = (i-1)%len(G.facesContainingVertex[v])
        return True
    
    #Check if the given angle type at can be successful assigned to the tiling-vertex (v,f) without breaking face conditions.
    def check_face(self,v,f,at,G,P):
        #This works because we have few faces and sides
        #Update MaxFaceSums and MinFaceSums
        maxangle = {0:89,1:90,2:179,3:180,-1:180}
        minangle = {0:1,1:90,2:91,3:180,-1:1}
        if self.numAT == 5:
            maxangle = {0:44, 1:45, 2:89, 3:90, 4:179, 5:180, -1:180}
            minangle = {0:1, 1:45, 2:46, 3:90, 4:91, 5:180, -1:1}
        MaxSum = self.MaxFaceSums[f] - maxangle[-1] + maxangle[at]
        MinSum = self.MinFaceSums[f] - minangle[-1] + minangle[at]
        if MinSum > 0:
            if self.PrintProof:
                print(f'      Impossible: Sum of angles at this tiling-vertex would be too large')
            return False
        if MaxSum < 0:
            if self.PrintProof:
                print(f'      Impossible: Sum of angles at this tiling-vertex would be too small')
            return False
        return True
    
    #Adds the new angle equations given by the current assignation. Returns false if the new system has no solution. True if it is still open.
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
            print(f'    New angle equation: {NewEq}=0')
        #Checks if the system can still be solved. Or if it's solved with a non acceptable answer. 
        eqSolution = solve(self.AngleEqs, rational=None)
        if self.PrintProof:
            print(f'    Reduced angle equations:\n' + f'      {eqSolution}'.replace(':',' =').replace('.0000000000000','.0').replace('.000000000000','.0'))
        if len(eqSolution) == 0:
            return False
        #if it has a possible solution yet. Checks if it's acceptable.
        possiblyRes = clean_solution(eqSolution, len(P))
        for i in range(len(possiblyRes)): #Checks for each variable
            #Special condition: quadrilateral
            if self.non45 and possiblyRes[i] == 45:
                if self.PrintProof:
                    print(f'      Impossible: ({v},{f}) can not have an angle of 45, violates Lemma 9')
                return False
            
            if possiblyRes[i] == 'Not Solved':
                continue
            isGood = True
            if self.numAT == 5:
                if (P[i] == 1 and possiblyRes[i] != 45) or (P[i] == 3 and possiblyRes[i] != 90):
                    isGood = False
                if P[i] == 0 and not (0 < possiblyRes[i] < 45):
                    isGood = False
                if P[i] == 2 and not (45 < possiblyRes[i] < 90):
                    isGood = False
                if P[i] == 4 and not (90 < possiblyRes[i] < 180):
                    isGood = False
            else:
                if P[i] == 1 and possiblyRes[i] != 90:
                    isGood = False
                if P[i] == 0 and not (0 < possiblyRes[i] < 90):
                    isGood = False
                if P[i] == 2 and not (90 < possiblyRes[i] < 180):
                    isGood = False
            if not isGood:
                if self.PrintProof:
                    print(f'Angle {i} out of bounds: {possiblyRes[i]}')
                return False
        self.LastAngSol = possiblyRes
        return True
#Class ends here

#The function which explores the tree
def search(ig, G, PermsFileName, use5AT = False, PrintProof = False): #index of graph, graph, number of sides
    Res = [] #Contains final result
    Perms = load_ang_perms(PermsFileName) #List of possible angle-types a tile may have
    
    #An angle is a pair (v,F) where v is a tile and F is a face of G containing v
    Angles, IndexFilledFace = construct_angles(G)
    
    #if using 5AT changes at_dict
    global at_dict
    if use5AT:
        at_dict = {0:'(sa)',1:'(ma)',2:'(la)',3:'r',4:'(so)',5:'(mo)',6:'(lo)',7:'p'}
    else:
        at_dict = {0:'a',1:'r',2:'o',3:'p'}
        
    for ip, P in enumerate(Perms):
        if PrintProof:
            print(f'\nAssuming T has labeling: {"".join([at_dict[p] for p in P])}\n')
        Stack = []
        N = NodeState(ig, ip, G, P, Angles, use5AT, PrintProof)
        NList = N.next(Angles, P, G, IndexFilledFace)
        #This assures you assign the first angletype to the first angle
        for N0 in NList:
            Stack.append(NList[0])
        while len(Stack) > 0:
            Branch = Stack.pop()
            if(Branch.indexAngle == len(Angles)):
                Res.append(Branch)
                if PrintProof:
                    print(f'\nFOUND A POSSIBLE SOLUTION!\n')
                continue
            NewBranches = Branch.next(Angles, P, G, IndexFilledFace)
            Stack.extend(NewBranches)
        if PrintProof:
            print(f'\nFound {len(Res)} possible solutions.\n\n-------------------------------------------------')
    return Res
