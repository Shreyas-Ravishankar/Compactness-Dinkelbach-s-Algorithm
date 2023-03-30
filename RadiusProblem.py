import sys
import numpy as np
import matplotlib.pyplot as plt
import math
import cplex
from cplex.exceptions import CplexError
from cplex.callbacks import LazyConstraintCallback as LCC
from cplex.callbacks import IncumbentCallback as ICC
import load_forest_data as lfd
import load_forest_data1 as lfd1
import collections
from pqdict import PQDict


sys.setrecursionlimit(10**8)

no_circle = 1
rows = 5
m = 1000
columns = 5
width = 10
nu = 8
ID = 1
g=30
U = np.array(list(range(1, nu + 1)))
area_p=0.2 #0.2,0.5
are_p = 0.07
min_comp = {}
for i in range(no_circle):
    min_comp[i+1] = 0
bud_p = 0.5 #0.5,1
limit = 1
dtol=0.0001 #distance cut tolerance
mip_gap=0
# min_area=7500
solve_new = True
use_lazy = True
use_random_data = False
use_inc = True
use_species = False
use_obj_cuts=True
use_compact=False
use_small_area=True
use_hole_callback = True
use_warmstart = False
use_heu = False
use_fractoint =False
use_connectivity = True
option = 1  ### 1: Best Bound, 0: DFS (default is 1)
global incr
incr = 0
global counter
do_Knapsack_count = True
counter = False
upperBound = 0
epsi = 1e-6
onal= {}
hole = []
sel = {}
w ={}
# Directories
data_dir = "C:/Users/shrey/OneDrive - Arizona State University/NRES2020/Forest-Conservation-Modified/Datasets/Carvajal Datasets/"
results_dir = "C:/Users/shrey/IdeaProjects/untitled3/Results/"



circ = []

for i in range(no_circle):
    circ.append(i+1)
#print(circ)

for i in circ:
    w[i] = {}

class Graph:

    # init function to declare class variables
    def __init__(self, V,nodes):
        self.V = V
        self.nodes=nodes
        #self.adj = [[] for i in nodes]
        self.adj=collections.defaultdict(list)

    def DFSUtil(self, temp, v, visited):

        # Mark the current vertex as visited
        #visited[v - 1] = True
        visited[v] = True

        # Store the vertex to list
        temp.append(v)

        # Repeat for all vertices adjacent
        # to this vertex v
        # for i in self.adj[v - 1]:
        #     if not visited[i - 1]:
        #         # Update the list
        #         temp = self.DFSUtil(temp, i, visited)
        for i in self.adj[v]:
            if not visited[i]:
                # Update the list
                temp = self.DFSUtil(temp, i, visited)
        return temp

        # method to add an undirected edge

    def addEdge(self, v, w):
        self.adj[v].append(w)
        self.adj[w].append(v)
        #self.adj[v-1].append(w)
        #self.adj[w-1].append(v)

        # Method to retrieve connected components

    # in an undirected graph
    def connectedComponents(self):
        visited = collections.defaultdict(bool)
        cc = []
        # for i in range(self.V):
        #     visited.append(False)
        #for v in range(1,self.V+1):
        for v in self.nodes:
            # #print(v)
            #if not visited[v - 1]:
            if not visited[v]:
                temp = []
                cc.append(self.DFSUtil(temp, v, visited))
            ##print(v,visited)
        return cc

# g=Graph(4,[1,5,6,7,9,10,20])
# g.addEdge(1,5)
# g.addEdge(6,7)
# g.addEdge(9,20)
# #print('test')

def generate_unit_vectors(nu, U):
    delta = 2 * math.pi / nu
    angle = 0
    u1 = {}
    u2 = {}
    for i in U:
        u1[i] = round(math.cos(angle), 6)
        u2[i] = round(math.sin(angle), 6)
        angle += delta
    return u1, u2


u1, u2 = generate_unit_vectors(nu, U)


def plot_solution(centers, width, selected, cost):
    for i in centers.keys():
        if i in selected:
            c = "r"
            a = 0.8
        else:
            c = "b"
            a = 0.08
        x = centers[i][0]
        y = centers[i][1]
        plt.plot([x - width / 2, x - width / 2], [y - width / 2, y + width / 2], alpha=a, color=c)
        plt.plot([x - width / 2, x + width / 2], [y + width / 2, y + width / 2], alpha=a, color=c)
        plt.plot([x + width / 2, x + width / 2], [y + width / 2, y - width / 2], alpha=a, color=c)
        plt.plot([x + width / 2, x - width / 2], [y - width / 2, y - width / 2], alpha=a, color=c)
        plt.text(x - width / 2, y, "({0})".format(round(cost[i], 2)), fontsize=10)
    plt.show()


def plot_solution_new(centers, width, selected, rad, centre, cost) -> object:
    fig, ax = plt.subplots()
    circ = {}
    for i in circ:
        cir[i] = plt.Circle(centre[i], radius=rad[i], color="b", alpha=0.2)
        ax.add_patch(cir[i])
        ax.plot()
    for i in centers.keys():
        if i in selected:
            c = "r"
            a = 0.8
        else:
            c = "b"
            a = 0.08
        x = centers[i][0]
        y = centers[i][1]
        plt.plot([x - width / 2, x - width / 2], [y - width / 2, y + width / 2], alpha=a, color=c)
        plt.plot([x - width / 2, x + width / 2], [y + width / 2, y + width / 2], alpha=a, color=c)
        plt.plot([x + width / 2, x + width / 2], [y + width / 2, y - width / 2], alpha=a, color=c)
        plt.plot([x + width / 2, x - width / 2], [y - width / 2, y - width / 2], alpha=a, color=c)
        plt.text(centers[i][0], centers[i][1], i)
    plt.show()


def generate_random_data(rows, columns, ID, width):
    def generate_polygon_centers(rows, columns, width):
        centers = {}
        x = width / 2
        y = width / 2
        ctr = 1
        for i in range(rows):
            for j in range(columns):
                centers[ctr] = (x + width * j, y + width * i)
                ctr += 1
        return centers

    def plot_grid(centers, width, cost):
        vertices = {}
        for i in centers.keys():
            vertices[i] = []
            x = centers[i][0]
            y = centers[i][1]
            plt.plot([x - width / 2, x - width / 2], [y - width / 2, y + width / 2], "b")
            vertices[i].append((x - width / 2, y - width / 2))
            plt.plot([x - width / 2, x + width / 2], [y + width / 2, y + width / 2], "b")
            vertices[i].append((x - width / 2, y + width / 2))
            plt.plot([x + width / 2, x + width / 2], [y + width / 2, y - width / 2], "b")
            vertices[i].append((x + width / 2, y + width / 2))
            plt.plot([x + width / 2, x - width / 2], [y - width / 2, y - width / 2], "b")
            vertices[i].append((x + width / 2, y - width / 2))
            vertices[i].append((x - width / 2, y - width / 2))
            plt.text(x - width / 2, y, "({0})".format(round(cost[i], 2)), fontsize=10)
        plt.show()
        return vertices

    centers = generate_polygon_centers(rows, columns, width)

    distance = {}
    for i in centers.keys():
        for j in centers.keys():
            distance[i, j] = ((centers[i][0] - centers[j][0]) ** 2 + (centers[i][1] - centers[j][1]) ** 2) ** 0.5

    neighbours = {}
    for i in centers.keys():
        neighbours[i] = []
        for j in centers.keys():
            if i != j and distance[i, j] <= width:
                neighbours[i].append(j)

    area = {}
    for i in centers.keys():
        area[i] = width ** 2

    filename = "grid_" + str(rows) + "x" + str(columns) + "_" + str(ID) + ".txt"

    def generate_costs(rows, columns, lb, ub):
        cost = {}
        for i in range(1, rows * columns + 1):
            cost[i] = np.random.uniform(lb, ub)
        return cost

    def write_data(filename, r, c, ub, lb):
        cost = generate_costs(r, c, lb, ub)
        f1 = open(filename, "w")
        for i in cost.keys():
            f1.write(str(i) + "\t" + str(cost[i]) + "\n")
        f1.close()
        return cost

    def read_data(filename):
        f1 = open(filename, "r")
        cost = {}
        for i in f1.readlines():
            j = i.strip().split("\t")
            cost[int(j[0])] = float(j[1])
        f1.close()
        return cost

    # cost=write_data(data_dir+filename,rows,columns,10,15)
    cost = read_data(data_dir + filename)
    vertices = plot_grid(centers, width, cost)

    circ = {}
    for i in centers.keys():
        circ[i] = round((area[i] / (2 * math.pi)) ** 0.5, 4)

    return centers, vertices, distance, neighbours, area, cost, circ


def use_subset_for_landscape(nodes,verticess,adj,area,cost,centers,ax,bx,ay,by):
    tbr=[]
    for i in nodes:
        trueval=True
        for j in verticess[i]:
            #print("hi")
            if (j[0]<ax or j[0]>bx) or (j[1]<ay or j[1]>by):
                trueval=False
                #print("Excluding ",i)
                break
        if not trueval:
            tbr.append(i)
    nodes=list(i for i in nodes if i not in tbr)
    for i in tbr:
        adj.pop(i)
        cost.pop(i)
        area.pop(i)
        verticess.pop(i)
        centers.pop(i)


    for i in adj.keys():
        adj[i]=list(i for i in adj[i] if i not in tbr)

    return area,cost,nodes,adj,verticess,centers

if use_random_data:
    nodes = list(range(1, rows * columns + 1))
    min_area = area_p * len(nodes) * width * width
    centers, verticess, distance, neighbours, area, cost, circ = generate_random_data(rows, columns, ID, width)
    xlim = columns * width
    ylim = rows * width
else:
    flgid = 5

    # A_max: Maximum area of the instance. Choose from one of these options
    # 1) 48.6 - for FLG9A instances.
    # 2) 40   - for Hardwicke and Shulkell and Random
    # 3) 80   - for NBCL5A instance.
    # 4) 120  - for El Dorado and Buttercreek.
    A_max = 40

    # T: Time periods: The value is 3 for the given model.
    T = 3

    # Yt: Number of years in each time period. Set this to 10 for El Dorado instance, 5 otherwise.
    Yt = 5

    #Data = lfd.FLG9A(A_max, T, flgid)  ## Import Data for FLG9A (5 instances)
    #Data=lfd.ButterCreek(A_max,T) ## Import Data for Buttercreek
    # Data=ld.NBCL5A(A_max,T)      ## Import Data for NBCL5A
    #Data=lfd.ElDorado(A_max,T)    ## Import Data for El Dorado
    #Data1 = lfd1.ElDorado(A_max, T)  ## Import Data for El Dorado
    # Data=lfd.ShulKellA(A_max,T)   ## Import Data for Shulkell
    Data = lfd.Hardwicke(A_max, T)  ## Import Data for Hardwicke
    Data1 = lfd1.Hardwicke(A_max, T)  ## Import Data for Hardwicke
    # Data=ld.Random50(A_max,T)     ## Import Data for random instance. Change the numbers among
    ## a) 50 b) 100 c) 200 d) 400 e) 900

    min_area = Data.A_min  ### - in same format
    nodes = Data.node_set  ### - in same format
    area = Data.area_set  ### - in same format

    min_area = area_p * sum(area.values())
    #print("min ", min_area)
    #print("total ", sum(area.values()))
    #print("ratio ", min_area / sum(area.values()))

    cost = Data.get_costs(Data.name, ID)
    for node in nodes:
        if node not in cost:
            cost[node] = 0
    neighbours = Data.adj  ### - in same format
    for i in neighbours.keys():
        for j in neighbours[i]:
            if i not in neighbours[j]:
                neighbours[j].append(i)

    verticess = Data.point_set  ### - in same format
    vertices = {}
    for i in verticess.keys():
        vertices[i] = []
        for vertex in verticess[i]:
            vertices[i].append((round(vertex[0], 2), round(vertex[1], 2)))

    centers = {}
    for i in vertices.keys():
        x = 0
        y = 0
        for j in vertices[i][:len(vertices[i]) - 1]:
            x += j[0]
            y += j[1]
        x = x / (len(vertices[i]) - 1)
        y = y / (len(vertices[i]) - 1)
        centers[i] = (x, y)

    distance = {}
    for i in centers.keys():
        for j in centers.keys():
            distance[i, j] = ((centers[i][0] - centers[j][0]) ** 2 + (centers[i][1] - centers[j][1]) ** 2) ** 0.5
    xlim = 175  ## 100 each for flg9a, (175,75) for hardwicke, (430,80) for el dorado
    ylim = 75
    #Data.plot_data(nodes)
    # input()

budget = bud_p * sum(list(cost.values()))

if use_small_area:
    ax = 00
    bx = 80
    ay = 00
    by = 80
    area, cost, nodes, neighbours, vertices,centers=use_subset_for_landscape(nodes,vertices,neighbours,area,cost,centers,ax,bx,ay,by)
    min_area = area_p * sum(area.values())
    mi_area = are_p * sum(area.values())
    budget = bud_p * sum(list(cost.values()))
    xlim = bx  ## 100 each for flg9a, (175,75) for hardwicke, (430,80) for el dorado
    ylim = by

## create variables names
x_bin = {}
for i in nodes:
    x_bin[i] = {}
for j in circ:
    for i in nodes:
        x_bin[i][j] = "x" + str(i) + str(j)

Radius = "r"
X = "x"
Y = "y"


#___________________________________________________________________Ratio Model_________________________________________________________________________________________________________
print("Radius Model===========================================================================================================================")

class IncumbentAnalysis(ICC):
    def __call__(self):
        global incr
        incr = 0
        global counter
        global upperBound
        counter = False
        radius = {}
        x={}
        y={}
        #print("Incumbent Callback")
        radius[1] = self.get_values(Radius)
        x[1] = self.get_values(X)
        y[1] = self.get_values(Y)
        print(f"X:{x[1]}, Y:{y[1]}, Radius:{radius[1]}")
        w = {}
        for j in circ:
            for i in nodes:
                w[(i,j)] = round(self.get_values(x_bin[i][j]))
        selected = {}
        semi = []
        compact = {}
        for j in circ:
            selected[j] = []
            for i in nodes:
                if w[(i,j)] == 1:
                    selected[j].append(i)
                if w[(i,j)] >0 and w[(i,j)] < 1:
                    semi.append(i)
        print("Selected Patches", selected)
        if use_random_data:
            plot_solution_new(centers, width, selected, radius[j], (x[j], y[j]), cost)
        else:
            Data.plot_soln_data(selected, [], "Radius model incumbent",nodes,(x,y),radius)

def dijkstra(G, start, end=None):
    start=str(start)

    inf= float('inf')
    D = {start: 0}  # mapping of nodes to their dist from start
    Q = PQDict(D)  # priority queue for tracking min shortest path
    P = {}  # mapping of nodes to their direct predecessors
    U = set(G.keys())  # unexplored nodes

    while U:  # nodes yet to explore
        (v, d) = Q.popitem()  # node w/ min dist d on frontier
        D[v] = d  # est dijkstra greedy score
        U.remove(v)  # remove from unexplored
        if v == end: break

        # now consider the edges from v with an unexplored head -
        # we may need to update the dist of unexplored successors
        for w in G[v]:  # successors to v
            if w in U:  # then w is a frontier node
                d = D[v] + G[v][w]  # dgs: dist of start -> v -> w
                if d < Q.get(w, inf):
                    Q[w] = d  # set/update dgs
                    P[w] = v  # set/update predecessor

    return D, P


def shortest_path(G, start, end):
    dist, pred = dijkstra(G, start, end)
    v = end
    path = [v]
    while v != start:
        v = pred[v]
        path.append(v)
    path.reverse()
    return path


def make_graph(filename):
    G = {}

    with open(filename) as file:
        for row in file:
            r = row.strip().split('\t')
            label = r.pop(0)
            neighbors = {v: int(length) for v, length in [e.split(',') for e in r]}
            G[label] = neighbors

    return G

def Radius_subProblem():
    class DistanceCuts(LCC):
        compactness = 0

        def __init__(self, env):
            LCC.__init__(self, env)

        def neighbourhood(self, z, adj): #neighbourhood of the connected reserves
            S = set()
            for i in z:
                for j in adj[i]:
                    S.add(j)
            return list(S - set(z))

        def get_patches(self, z, adj): #connected components (multi reserve)
            patch = []
            in_list = {}

            for i in z:
                in_list[i] = False
            for i in z:
                s = []
                l = []
                if in_list[i] == True:
                    continue
                else:
                    s.append(i)
                    in_list[i] = True
                    while len(s) > 0:
                        for j in s:
                            for k in adj[j]:
                                if (k in z) and (in_list[k] == False):
                                    in_list[k] = True
                                    s.append(k)
                            l.append(j)
                            s.remove(j)
                patch.append(l)
            #print("Patches mate!!!! Watch out!!!!!", patch)
            return patch

        def patch_areas(self, patch, area_set): #area of the reserves
            area = {}
            for i, j in enumerate(patch):
                a = 0
                for k in j:
                    a += area_set[k]

                area[i] = a
            return area

        def __call__(self):
            print("Lazy Callback")
            selected = {}
            x = {}
            y = {}
            RR = {}
            x[1] = self.get_values(X)
            y[1] = self.get_values(Y)
            RR[1] = self.get_values(Radius)
            print(f"X:{x[1]}, Y:{y[1]}, Radius:{RR[1]}")
            for k in circ:
                ask = {}
                for i in nodes:
                    for ind, j in enumerate(vertices[i][:len(vertices[i]) - 1]):
                        ask[i, ind] = ((j[0] - x[1]) ** 2 + (j[1] - y[1]) ** 2) ** 0.5

                w = {}
                for i in nodes:
                    w[i] = {}
                    w[i][k] = round(self.get_values(x_bin[i][k]))

                lets = True
                M = xlim/2
                for i in nodes:
                    if w[i][k] == 1:
                        #print(x_bin[i][k])
                        for ind, j in enumerate(vertices[i][:len(vertices[i]) - 1]):
                            if ask[i, ind] > RR[1] + dtol:
                                #print("2 This constraint is violated!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
                                lets = False
                                xt = x[1] - j[0]
                                yt = y[1] - j[1]
                                u1t = xt / ask[i, ind]
                                u2t = yt / ask[i, ind]
                                var = []
                                coeff1 = []
                                coeff2 = []
                                var.append(X)
                                var.append(Y)
                                coeff1.append(-u1t)
                                coeff1.append(-u2t)
                                coeff2.append(u1t)
                                coeff2.append(u2t)
                                var.append(x_bin[i][k])
                                coeff1.append(M)
                                coeff2.append(M)
                                var.append(Radius)
                                coeff1.append(-1)
                                coeff2.append(-1)
                                rhs1 = M - j[0] * u1t - j[1] * u2t
                                rhs2 = M + j[0] * u1t + j[1] * u2t
                                self.add(constraint=cplex.SparsePair(var, coeff1), sense="L", rhs=rhs1)
                                self.add(constraint=cplex.SparsePair(var, coeff2), sense="L", rhs=rhs2)

                if use_connectivity:
                    comp = False
                    selected[k] = []
                    for i in nodes:
                        if w[i][k] == 1:
                            selected[k].append(i)
                    patches = self.get_patches(selected[k], neighbours)
                    print("The length of patches is for circle:", k, "=", len(patches))
                    # portion of code to be corrected!! Something going on witj the IF statement, lets value. The ring inequalaties are being implemented only a few times.
                    if True:
                        # print("No cut added, this should be feasible")
                        if len(patches) > 1:
                            print("Adding ring inequalities", len(patches))
                            s = []
                            area_set = self.patch_areas(patches, area)
                            for i in range(len(patches)):
                                if area_set[i] < mi_area or 1:
                                    neigh = self.neighbourhood(patches[i], neighbours)
                                    if neigh != []:
                                        comp = True
                                        var = []
                                        coeff = []
                                        for j in patches[i]:
                                            var.append(x_bin[j][k])
                                            coeff.append(1)
                                        for j in neigh:
                                            var.append(x_bin[j][k])
                                            coeff.append(-1)
                                        rhs = len(patches[i]) - 1
                                        print(f"I am adding for {i}th patch of circ {k}")
                                        ccc = 0
                                        self.add(constraint=cplex.SparsePair(var, coeff), sense='L', rhs=rhs)

            #print("Here")

    #print("im here 1")
    try:
        prob = cplex.Cplex()
        prob.parameters.mip.display.set(4)
        prob.objective.set_sense(prob.objective.sense.minimize)
        prob.parameters.clocktype.set(2)
        prob.parameters.timelimit.set(7200)
        prob.parameters.mip.tolerances.mipgap.set(mip_gap)
        prob.parameters.barrier.display.set(2)
        prob.parameters.threads.set(1)
        prob.parameters.workmem.set(1024)
        prob.parameters.emphasis.memory.set(1)
        prob.parameters.workdir.set("C:/Users/shrey/OneDrive - Arizona State University/NRES2020/Forest-Conservation-Modified/Datasets/Carvajal Datasets/")
        prob.parameters.mip.strategy.file.set(2)
        prob.parameters.mip.strategy.nodeselect.set(option)
        print("Solving Radius Model")

        def add_compactness_by_rad_const1(inst):
            M=xlim/2
            inst.variables.add(obj=[1], lb=[(min_area/math.pi)**0.5], ub=[xlim/2], types="C", names=[Radius])
            inst.variables.add(obj=[0], lb=[0.0], ub=[xlim], types="C", names=[X])
            inst.variables.add(obj=[0], lb=[0.0], ub=[ylim], types="C", names=[Y])
            for j in circ:
                for i in nodes:
                    inst.variables.add(obj=[0], lb=[0], ub=[1], types="B", names=[x_bin[i][j]])
            #print("im here 1 in sub problem")
            for k in circ:
                for i in nodes:
                    for j in vertices[i][:len(vertices[i]) - 1]:
                        for u in U:
                            var = []
                            coeff = []
                            var.append(X)
                            var.append(Y)
                            coeff.append(-u1[u])
                            coeff.append(-u2[u])
                            var.append(x_bin[i][k])
                            coeff.append(M)
                            var.append(Radius)
                            coeff.append(-1)
                            rhs = M - j[0]*u1[u] - j[1]*u2[u]
                            inst.linear_constraints.add(lin_expr=[cplex.SparsePair(var, coeff)], senses="L", rhs=[rhs],
                                                        names=["radius constraint " + str(i) + "," + str(j) + "," + str(
                                                            u) + "," + str(k)])

        if solve_new:
            add_compactness_by_rad_const1(prob)
        else:
            add_variables_original(prob)

        ## constraint 3.2
        var = []
        coeff = []
        for j in circ:
            for i in nodes:
                var.append(x_bin[i][j])
                coeff.append(cost[i])
        prob.linear_constraints.add(lin_expr=[cplex.SparsePair(var, coeff)], senses="L", rhs=[budget], names=["budget"])

        var = []
        coeff = []
        for k in circ:
            for i in nodes:
                var.append(x_bin[i][k])
                coeff.append(area[i])
        prob.linear_constraints.add(lin_expr=[cplex.SparsePair(var, coeff)], senses="G", rhs=[min_area],
                                        names=["minarea"])


        if use_lazy:
            prob.register_callback(DistanceCuts)
        if use_inc:
            prob.register_callback(IncumbentAnalysis)
        # if use_heu:
        #     prob.register_callback(HeuristicAnalysis)
        # if use_fractoint:
        #     prob.register_callback(Fractoint)

        #add_compactness_by_rad_const1(prob)
        #### ---- feasibility checker of another solution ------ ###
        prob.solve()
        x={}
        y={}
        vall = {}
        radius ={}
        total_area = {}
        selected = {}
        sol = []
        obj = 0
        if prob.solution.get_status() == 101:
            selected = {}
            xvals={}
            for i in x_bin.keys():
                for j in x_bin[i].keys():
                    xvals[i]=prob.solution.get_values(x_bin[i][j])
            #xvals={j:prob.solution.get_values(x_bin[i][j]) for i in x_bin.keys() for j in x_bin[i].keys()}
            for k in circ:
                selected[k] = []
                semi = []
                dam = []
                for i in nodes:
                    count = round(prob.solution.get_values(x_bin[i][k]))
                    if count == 1:
                        selected[k].append(i)
                    if count > 0 and count < 1:
                        semi.append(i)

                total_area[k] = 0
                total_cost = 0
                for node in selected[k]:
                    total_area[k]+=area[node]
                objective = prob.solution.get_objective_value()
                radius[k] = prob.solution.get_values(Radius)
                x[k] = prob.solution.get_values(X)
                y[k] = prob.solution.get_values(Y)
                obj = objective
            circArea = math.pi*radius[1]*radius[1]
            print("The radius compactness is,", (total_area[1]/circArea))
            Data.plot_soln_data(selected, [], "Radius Model Solution (warmstart)", nodes, (x, y), radius)
            print("10101010 Solution completed!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!", objective)
            varr = prob.variables.get_names()
            for i in varr:
                if prob.variables.get_types(i) != "C":
                    vall[i] = round(prob.solution.get_values(i))
                else:
                    vall[i] = prob.solution.get_values(i)
            for node in selected[k]:
                sol.append(node)
            return sol, radius[1]
        else:
            return [], 100
    except CplexError as e:
        print(e)
        return [], 100

#val, sol = Radius_subProblem()