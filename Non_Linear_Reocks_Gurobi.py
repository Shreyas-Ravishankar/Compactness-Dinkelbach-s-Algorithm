### The aim of this code is to implement Dinkelbach’s algorithm

import numpy as np
import matplotlib.pyplot as plt
import math
import load_forest_data1 as lfd
import time
import sys
import pandas as pd
import collections
import pickle
from collections import defaultdict
import random
sys.setrecursionlimit(10**8)

import gurobipy as gp
from gurobipy import GRB


no_circle=1
nu = 8
ID = 1
U = np.array(list(range(1, nu + 1)))
area_p=0.1
min_comp = 0
bud_p = 0.15
dtol=0 #distance cut tolerance
mip_gap=0.00
costPerUnitArea = 12
use_connectivity = True
#### new

radiusDetails = {}
prev_sol = {}
sel = {}

# radiusMax = 100
# radiusMin = 0
#bigM = 100

#maxRadius = 0
for i in range(no_circle):
    prev_sol[i + 1] = []
# min_area=7500

### new until here

# Directories
data_dir = "C:/Users/shrey/OneDrive - Arizona State University/NRES2020/Forest-Conservation-Modified/Datasets/Carvajal Datasets/"
results_dir = "C:/Users/shrey/IdeaProjects/untitled3/Results/"
#data_dir = "./Datasets/"
#results_dir = "./Results/"
#
# _________________________________________________________________Smallest Enclosing circle______________________________________________________

# Python3 program to find the minimum enclosing
# circle for N integer points in a 2-D plane
from math import sqrt
from random import randint, shuffle

# Defining infinity
INF = 1e18


# Structure to represent a 2D point
class Point:
    def __init__(self, X=0, Y=0) -> None:
        self.X = X
        self.Y = Y


# Structure to represent a 2D circle
class Circle:
    def __init__(self, c=Point(), r=0.0, A=Point(0, 0), B=Point(0, 0), C=Point(0, 0)) -> None:
        self.C = c
        self.R = r
        self.P1 = A
        self.P2 = B
        self.P3 = C


# Function to return the euclidean distance
# between two points
def dist(a, b):
    return sqrt(pow(a.X - b.X, 2)
                + pow(a.Y - b.Y, 2))


# Function to check whether a point lies inside
# or on the boundaries of the circle
def is_inside(c, p):
    return dist(c.C, p) <= c.R


# The following two functions are used
# To find the equation of the circle when
# three points are given.

# Helper method to get a circle defined by 3 points
def get_circle_center(bx, by,
                      cx, cy):
    B = bx * bx + by * by
    C = cx * cx + cy * cy
    D = bx * cy - by * cx
    return Point((cy * B - by * C) / ((2 * D) + 0.000000000001),
                 (bx * C - cx * B) / ((2 * D) + 0.000000000001))


# Function to return the smallest circle
# that intersects 2 points
def circle_from1(A, B):
    # Set the center to be the midpoint of A and B
    C = Point((A.X + B.X) / 2.0, (A.Y + B.Y) / 2.0)

    # Set the radius to be half the distance AB
    return Circle(C, dist(A, B) / 2.0, Point(A.X, A.Y), Point(B.X, B.Y))


# Function to return a unique circle that
# intersects three points
def circle_from2(A, B, C):
    I = get_circle_center(B.X - A.X, B.Y - A.Y,
                          C.X - A.X, C.Y - A.Y)

    I.X += A.X
    I.Y += A.Y
    return Circle(I, dist(I, A), Point(A.X, A.Y), Point(B.X, B.Y), Point(C.X, C.Y))


# Function to check whether a circle
# encloses the given points
def is_valid_circle(c, P):
    # Iterating through all the points
    # to check whether the points
    # lie inside the circle or not
    for p in P:
        if (not is_inside(c, p)):
            return False
    return True


# Function to return the minimum enclosing
# circle for N <= 3
def min_circle_trivial(P):
    assert (len(P) <= 3)
    if not P:
        return Circle()

    elif (len(P) == 1):
        return Circle(P[0], 0)

    elif (len(P) == 2):
        return circle_from1(P[0], P[1])

    # To check if MEC can be determined
    # by 2 points only
    for i in range(3):
        for j in range(i + 1, 3):

            c = circle_from1(P[i], P[j])
            if (is_valid_circle(c, P)):
                return c

    return circle_from2(P[0], P[1], P[2])


# Returns the MEC using Welzl's algorithm
# Takes a set of input points P and a set R
# points on the circle boundary.
# n represents the number of points in P
# that are not yet processed.
def welzl_helper(P, R, n):
    # Base case when all points processed or |R| = 3
    if (n == 0 or len(R) == 3):
        return min_circle_trivial(R)

    # Pick a random point randomly
    idx = randint(0, n - 1)
    p = P[idx]
    boundry = []
    # Put the picked point at the end of P
    # since it's more efficient than
    # deleting from the middle of the vector
    P[idx], P[n - 1] = P[n - 1], P[idx]

    # Get the MEC circle d from the
    # set of points P - :p
    d = welzl_helper(P, R.copy(), n - 1)

    # If d contains p, return d
    if (is_inside(d, p)):
        return d

    # Otherwise, must be on the boundary of the MEC
    else:
        R.append(p)

    # Return the MEC for P - :p and R U :p
    return welzl_helper(P, R.copy(), n - 1)


def welzl(P):
    P_copy = P.copy()
    shuffle(P_copy)
    return welzl_helper(P_copy, [], len(P_copy))


# ------------------------------------------------------------------------------------------------------------------


class global_param:
    maxRadius = 0
    radiusMin=0
    bigM=100
    sol_queue={}
    Q_queue={}

    def __init__(self, area_p, bud_p):
        self.area_p = area_p
        self.bud_p = bud_p

    def store_solution(self,selected,iter,q):
        self.sol_queue[iter]=set(selected)
        self.Q_queue[iter]=q

    def check_repeat(self,selected):
        selected=set(selected)
        for i in self.sol_queue.keys():
            if selected==self.sol_queue[i]:
                return True
        return False



class getData:
    A_max={"FLG9A":48.6,"Hardwicke":40,
           "Shulkell":40,"NBCL5A":80,
           "ElDorado":120,"Buttercreek":120}

    XL = {"FLG9A":100,"Hardwicke":175,
           "ElDorado":430}  ## 100 each for flg9a, (175,75) for hardwicke, (430,80) for el dorado

    YL = {"FLG9A":100,"Hardwicke":75,
           "ElDorado":80}


    def __init__(self,type,small,ax=0,bx=0,ay=0,by=0):
        self.type=type[0]
        self.small=small
        if self.type=="Random":
            self.rows=type[1]
            self.columns=type[2]
            self.width=type[3]
            self.nodes = list(range(1, self.rows * self.columns + 1))
            self.min_area = area_p * len(self.nodes) * self.width * self.width
            self.generate_random_data(ID)
            self.xlim = self.columns * self.width
            self.ylim = self.rows * self.width
            self.budget=bud_p * sum(list(self.cost.values()))
            self.circ=[i for i in range(1,no_circle+1)]

        else:
            self.ax=ax
            self.ay=ay
            self.bx=bx
            self.by=by
            print(f"Getting datan for {self.type}")
            self.get_actual_data(ax,bx,ay,by)
            self.circ = [i for i in range(1, no_circle + 1)]


    def generate_random_data(self,ID):
        rows=self.rows
        columns=self.columns
        width=self.width
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

        self.centers=centers.copy()
        self.vertices=vertices.copy()
        self.distance=distance.copy()
        self.neighbours=neighbours.copy()
        self.area=area.copy()
        self.circ=circ.copy()
        self.cost=cost.copy()

    def get_actual_data(self,ax,bx,ay,by):
        flgid = 5
        A_max = self.A_max[self.type]
        T = 3 # T: Time periods: The value is 3 for the given model.
        Yt = 10 # Yt: Number of years in each time period. Set this to 10 for El Dorado instance, 5 otherwise.
        if self.type=="FLG9A":
            Data = lfd.FLG9A(A_max, T, flgid)
        elif self.type=="Buttercreek":
            Data = lfd.ButterCreek(A_max, T)
        elif self.type=="NBCL5A":
            Data = lfd.NBCL5A(A_max, T)
        elif self.type=="ElDorado":
            Data = lfd.ElDorado(A_max, T)
        elif self.type=="ShulKell":
            Data = lfd.ShulKellA(A_max, T)
        elif self.type=="Hardwicke":
            Data = lfd.Hardwicke(A_max, T)

        self.DATA=Data

        min_area = Data.A_min  ### - in same format
        nodes = Data.node_set  ### - in same format
        area = Data.area_set  ### - in same format

        min_area = area_p * sum(area.values())


        cost = Data.get_costs(Data.name, ID)
        for node in nodes:
            if node not in cost.keys():
                cost[node] = 0
            elif node in cost.keys():
                noise = random.normalvariate(5, 0.5)
                if noise > 0:
                    cost[node] = costPerUnitArea * area[node] + random.normalvariate(5, 0.5)
                else:
                    cost[node] = costPerUnitArea * area[node] + 0.001

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

        maxx = 0
        distance = {}
        distt = {}
        for i in centers.keys():
            for j in centers.keys():
                distance[i, j] = ((centers[i][0] - centers[j][0]) ** 2 + (centers[i][1] - centers[j][1]) ** 2) ** 0.5
        for i in nodes:
            distt[i] = {}
            for j in nodes:
                dis = []
                for k in vertices[j]:
                    dis.append(((centers[i][0] - k[0]) ** 2 + (centers[i][1] - k[1]) ** 2) ** 0.5)
                distt[i][j] = max(dis)
                if maxx < distt[i][j]:
                    maxx = distt[i][j]
        xlim = self.XL[self.type]  ## 100 each for flg9a, (175,75) for hardwicke, (430,80) for el dorado
        ylim = self.XL[self.type]
        print(f"Number of nodes originally {len(nodes)}")
        self.budget = bud_p * sum(list(cost.values()))

        if self.small:
            print("Getting subset datan")
            area, cost, nodes, neighbours, vertices, centers = self.use_subset_for_landscape(nodes, vertices, neighbours, area, cost,
                                                                                        centers, ax, bx, ay, by)
            min_area = area_p * sum(area.values())
            # mi_area = are_p * sum(area.values())
            self.budget = bud_p * sum(list(cost.values()))
            xlim = bx  ## 100 each for flg9a, (175,75) for hardwicke, (430,80) for el dorado
            ylim = by
            print(f"Number of nodes in small {len(nodes)}")


        self.centers = centers.copy()
        self.vertices = vertices.copy()
        self.distance = distance.copy()
        self.neighbours = neighbours.copy()
        self.cost=cost.copy()
        self.area = area.copy()
        self.nodes=nodes[:]
        self.xlim=xlim
        self.ylim=ylim
        self.min_area=min_area

        print("Minimum Area ", min_area)
        print("Total ", sum(area.values()))
        print("Ratio ", min_area / sum(area.values()))

    def use_subset_for_landscape(self,nodes, vertices, adj, area, cost, centers, ax, bx, ay, by):
        tbr = []
        for i in nodes:
            trueval = True
            for j in vertices[i]:
                #print("hi")
                if (j[0] < ax or j[0] > bx) or (j[1] < ay or j[1] > by):
                    trueval = False
                    #print("Excluding ", i)
                    break
            if not trueval:
                tbr.append(i)
        nodes = list(i for i in nodes if i not in tbr)
        for i in tbr:
            adj.pop(i)
            cost.pop(i)
            area.pop(i)
            vertices.pop(i)
            centers.pop(i)

        for i in adj.keys():
            adj[i] = list(i for i in adj[i] if i not in tbr)

        return area, cost, nodes, adj, vertices, centers


### graph class used
class Graph:

    # init function to declare class variables
    def __init__(self, V):
        self.V = V
        self.adj = [[] for i in range(V)]

    def DFSUtil(self, temp, v, visited):

        # Mark the current vertex as visited
        visited[v - 1] = True

        # Store the vertex to list
        temp.append(v)

        # Repeat for all vertices adjacent
        # to this vertex v
        for i in self.adj[v - 1]:
            if not visited[i - 1]:
                # Update the list
                temp = self.DFSUtil(temp, i, visited)
        return temp

        # method to add an undirected edge

    def addEdge(self, v, w):
        self.adj[v - 1].append(w)
        self.adj[w - 1].append(v)

        # Method to retrieve connected components

    # in an undirected graph
    def connectedComponents(self):
        visited = []
        cc = []
        for i in range(self.V):
            visited.append(False)
        for v in range(1,self.V+1):
            # print(v)
            if not visited[v - 1]:
                temp = []
                cc.append(self.DFSUtil(temp, v, visited))
            #print(v,visited)
        return cc

### generates unit vectors
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


## TODO: one of these two is redundant, need to find
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

def plot_solution_new(centers, width, selected, rad, centre, cost):
    fig, ax = plt.subplots()
    cir = plt.Circle(centre, radius=rad[1], color="b", alpha=0.2)
    ax.add_patch(cir)
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
        # plt.text(x - width / 2, y, "({0})".format(round(cost[i], 2)), fontsize=10)
    plt.show()

def plot_random_solution(centers,width,selected,rad,centre,cost):
    print("hi")


def Alternative_CCP_Knapsack_optimization(datan, circ_data, pointsInside, criticalPoints, u_upper_bound):
    probPoints = []
    for node in criticalPoints:
        if node not in pointsInside:
            pointsInside.append(node)
    for node in pointsInside:
        probPoints.append(node)
    for node in criticalPoints:
        if node not in probPoints:
            probPoints.append(node)

    cost = datan.cost
    area = datan.area
    budget = datan.budget
    circ = datan.circ
    min_area = datan.min_area
    nodes = datan.nodes
    x_bin={}
    for i in datan.nodes:
        x_bin[i] = {}
    for j in circ:
        for i in datan.nodes:
            x_bin[i][j] = "x" + str(i) + str(j)

    #print("im here 1", pointsInside, probPoints)
    try:
        prob = cplex.Cplex()
        prob.parameters.mip.display.set(4)
        prob.objective.set_sense(prob.objective.sense.maximize)
        prob.parameters.clocktype.set(2)
        prob.parameters.timelimit.set(7200)
        prob.parameters.mip.tolerances.mipgap.set(mip_gap)
        prob.parameters.barrier.display.set(2)
        prob.parameters.threads.set(1)
        prob.parameters.workmem.set(1024)
        prob.parameters.emphasis.memory.set(1)
        #prob.parameters.workdir.set("C:/Users/shrey/OneDrive - Arizona State University/NRES2020/Forest-Conservation-Modified/Datasets/Carvajal Datasets/")
        prob.parameters.mip.strategy.file.set(2)
        prob.parameters.mip.strategy.nodeselect.set(option)

        print("Inside Kanpsack alternative formulation", datan.min_area)
        def add_compactness_by_rad_const1(inst):
            for j in circ:
                for i in pointsInside:
                    #print(x_bin[i][j])
                    inst.variables.add(obj=[area[i]], lb=[0], ub=[1], types="B", names=[x_bin[i][j]])

        add_compactness_by_rad_const1(prob)
        ## constraint 3.2
        var = []
        coeff = []
        for j in circ:
            for i in pointsInside:
                var.append(x_bin[i][j])
                coeff.append(cost[i])
        prob.linear_constraints.add(lin_expr=[cplex.SparsePair(var, coeff)], senses="L", rhs=[budget], names=["budget"])

        var = []
        coeff = []
        for k in circ:
            for i in pointsInside:
                var.append(x_bin[i][k])
                coeff.append(area[i])
        prob.linear_constraints.add(lin_expr=[cplex.SparsePair(var, coeff)], senses="G", rhs=[min_area],
                                     names=["minarea"])


        var = []
        coeff = []
        for i in criticalPoints:
            var.append(x_bin[i][1])
            coeff.append(1)
        prob.linear_constraints.add(lin_expr=[cplex.SparsePair(var, coeff)], senses="G", rhs=[len(criticalPoints)],
                                    names=["CriticalPointsCon"])


        prob.solve()
        vall = {}
        total_area = {}
        selected = {}
        obj = 0
        if "optimal" in prob.solution.get_status_string():
            print(f"Problem status : {prob.solution.get_status_string()}")
            for k in circ:
                selected[k] = []
                semi = []
                dam = []
                for i in pointsInside:
                    count = prob.solution.get_values(x_bin[i][k])
                    if round(count) > 0:
                        selected[k].append(i)
                    if count > 0 and count < 1:
                        semi.append(i)

                total_area[k] = 0
                total_cost = 0
                for i in selected[k]:
                    total_area[k] += area[i]
                    total_cost += cost[i]
                objective = prob.solution.get_objective_value()
            obj = (objective/(math.pi*(circ_data[2]**2)))
            print("Objective:", objective)
            #datan.DATA.plot_soln_data(selected[1], "Knapsack SubProblem NEW", nodes, ([0], [0]), 0)
            #print("1010101010101010 Solution completed Knapsack!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
            print("Knapsack subprob solved")
            varr = prob.variables.get_names()
            for i in varr:
                if prob.variables.get_types(i) != "C":
                    vall[i] = round(prob.solution.get_values(i))
                else:
                    vall[i] = prob.solution.get_values(i)
            #print(vall)
            return obj, vall, objective
        else:
            return 0, [False], 0

    except CplexError as e:
        print(e)
        return 0, [False], 0

def solve_newtons_method(data,gpars):
    centers=data.centers
    cost=data.cost
    area=data.area
    budget=data.budget
    circ=data.circ
    vertices=data.vertices
    xlim=data.xlim
    ylim=data.ylim
    min_area=data.min_area
    distance=data.distance
    nodes=data.nodes
    neighbours=data.neighbours

    #### new
    filename=f"{data.type}_{data.ax}_{data.bx}_{data.ay}_{data.by}_a{int(100*gpars.area_p)}_b{int(100*gpars.bud_p)}.csv"
    #filename = "holdeData.csv"
    try:
        val = pd.read_csv(filename)
        flag = True
    except:
        print("File not present")
        flag = False
        ## create the datan to write
        data_to_write = {"NID": [], "R": []}
        for i in nodes:
            data_to_write["NID"].append(i)
            points = []
            for vertex in vertices[i]:
                points.append(Point(vertex[0], vertex[1]))
            circleDetails = welzl(points)
            data_to_write["R"].append(circleDetails.R)
            radiusDetails[i] = circleDetails.R
        data_f = pd.DataFrame(data_to_write)
        data_f.to_csv(filename)

    if flag:
        val = val.drop(['Unnamed: 0'], axis=1)
        for i in nodes:
            radiusDetails[i] = val.loc[val.NID == i].R.values[0]

    ### new until here
    def neighbourhood(z, adj):  # neighbourhood of the connected reserves
        S = set()
        for i in z:
            for j in adj[i]:
                S.add(j)
        return list(S - set(z))

    def get_patches(z, adj):  # connected components (multi reserve)
        patch = []
        in_list = {}
        p = {}
        for i in nodes:
            p[i] = False
        for i in z:
            p[i] = True
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
                            if p[k] and (in_list[k] == False):
                                in_list[k] = True
                                s.append(k)
                        l.append(j)
                        s.remove(j)
            patch.append(l)
        # print("Patches mate!!!! Watch out!!!!!", patch)
        return patch

    def patch_areas(patch, area_set):  # area of the reserves
        area = {}
        for i, j in enumerate(patch):
            a = 0
            for k in j:
                a += area_set[k]

            area[i] = a
        return area

    ## change from here

    m = gp.Model("qcp")

    ## create variables names
    x_bin = {}

    R ={}
    X = {}
    Y = {}

    ## new
    radi = []
    for i in nodes:
        radi.append(radiusDetails[i])
    print("Smallest radius", min(radi))
    print("Minimum area is", min_area)
    u_upper_bound = 1 / min_area


    for i in nodes:
        x_bin[i] = {}
    for j in circ:
        for i in nodes:
            x_bin[i][j] = m.addVar(lb=0,ub=1,vtype=GRB.BINARY,name=f"x{i},{j}")
        R[j] = m.addVar(lb=0,ub=1000,vtype=GRB.CONTINUOUS,name=f"r{j}")
        X[j] = m.addVar(lb=0,ub=1000,vtype=GRB.CONTINUOUS,name=f"X{j}")
        Y[j] = m.addVar(lb=0,ub=1000,vtype=GRB.CONTINUOUS,name=f"Y{j}")


    def mycallback(m, where):
        if where == GRB.Callback.MIPSOL:
            print("Lazy Callback")
            global maxRad
            selected = {}
            x = {}
            y = {}
            RR = {}
            m._vars = m.getVars()
            print(m._vars)
            x[1] = m.cbGetSolution(X[1])
            y[1] = m.cbGetSolution(Y[1])
            RR[1] = m.cbGetSolution(R[1])
            print(f"X:{x[1]}, Y:{y[1]}, Radius:{RR[1]}")
            w ={}
            for i in nodes:
                w[i,1] = m.cbGetSolution(x_bin[i][1])
            for k in circ:
                ask = {}
                for i in nodes:
                    for ind, j in enumerate(vertices[i][:len(vertices[i]) - 1]):
                        ask[i, ind] = ((j[0] - x[1]) ** 2 + (j[1] - y[1]) ** 2) ** 0.5


                if use_connectivity:
                    comp = False
                    selected[1] = []
                    for i in nodes:
                        if w[i, 1] == 1:
                            selected[1].append(i)
                    patches = get_patches(selected[1], neighbours)
                    print("The length of patches is for circle:", 1, "=", len(patches))
                    # portion of code to be corrected!! Something going on witj the IF statement, lets value. The ring inequalaties are being implemented only a few times.
                    if True:
                        # print("No cut added, this should be feasible")
                        if len(patches) > 1:
                            print("Adding ring inequalities", len(patches))
                            s = []
                            area_set = patch_areas(patches, area)
                            for i in range(len(patches)):
                                if area_set[i] < min_area or 1:
                                    neigh = neighbourhood(patches[i], neighbours)
                                    if neigh != []:
                                        comp = True
                                        rhs = len(patches[i]) - 1
                                        m.cbLazy(sum(x_bin[j][1] for j in patches[i]) <= sum(x_bin[k][1] for k in neigh) + rhs)



    try:
        m.Params.NonConvex=2
        ## using this function - distance constraints
        def add_compactness_by_rad_const(m):
            M = 100000
            for k in circ:
                for i in data.nodes:
                    for j in vertices[i][:len(vertices[i]) - 1]:
                        m.addConstr((X[k]-j[0])**2+(Y[k]-j[1])**2<=R[k]**2+M*(1-x_bin[i][k]),f"eu_dist_{i},{j},{k}")



        ## adding VI to strengthen model
        def add_patch_distance_const(m):
            for k in circ:
                for i in data.nodes:
                    for j in data.nodes:
                        if j>i:
                            m.addConstr(np.floor(distance[i,j])*(x_bin[i][k]+x_bin[j][k])<=4*R[k],f"VI_{i},{j},{k}")





        print("Adding distance constraints")
        add_compactness_by_rad_const(m)

        print("Adding budget constraint")
        m.addConstr(sum(x_bin[i][j]*cost[i] for j in circ for i in data.nodes)<=budget,"budget")

        print("Adding minimum area constraint")
        m.addConstr(sum(x_bin[i][j] * area[i] for j in circ for i in data.nodes) >= min_area, "min area")



        print("Adding patch dist VI")
        add_patch_distance_const(m)

        m.setObjective(R[1], GRB.MINIMIZE)

        ### iterative method starts here
        t1=time.time()
        m._vars = vars
        m.Params.LazyConstraints = 1
        m.optimize(mycallback)
        t2=time.time()
        print("Time Elapsed : ", round(t2 - t1, 4))

        selected={}
        for j in circ:
            selected[j] = []
            for i in data.nodes:
                if x_bin[i][j].X>0:
                    selected[j].append(i)

        data.DATA.plot_soln_data(selected,[],"Gurobi Sol Non-linear model",nodes,(X[1].X,Y[1].X), R[1].X)



        #################--------------------- final solution ------------------##########################


    except TypeError as e:
        print(e)



if __name__=="__main__":
    print("Starting code")
    personal=False
    if personal:
        use_enhance=True
        use_lazy = True                         ## switch to implement lazy callbacks
        use_warm = False                        ## switch to use warm start
        use_context = not use_lazy              ## switch to implement context callback (implements all constraints in lazy)
        use_random_data = False                 ## switch to use random datan
        use_inc = True and use_lazy            ## switch to use incumbent callback
        use_obj_cuts = True                     ## switch to use objective cutting planes
        use_small_area = False                   ## switch to use small area
        use_ch=True                             ## switch to use new convex hull based constraints
        option = 1  ### 1: Best Bound, 0: DFS (default is 1)
        no_circle = 1
        try_traditional_newtons=True           ## switch to use tradional newtons method. If false, the problem uses
                                                ## an alternative iterative method to get a feasible solution in each iteration
                                                ## instead of finding the optimal solution

        rows = 6
        columns = 7
        width = 10
        task=["Random",rows,columns,width]
        task = ["ElDorado"]
        data=getData(type=task,small=use_small_area,ax=00,bx=430,ay=00,by=80)
        gpars=global_param(area_p,bud_p)
        solve_newtons_method(data,gpars)
    else:
        #task_id = int(sys.argv[1])
        task_id=14
        #task_pykl = sys.argv[2]
        task_pykl="iise_instances.pykl"
        with open(task_pykl, 'rb') as fh:
            tsk = pickle.load(fh)[task_id]
        print(f"Processing task: {task_id} {tsk}")

        print("Hi")
        use_enhance = False
        use_lazy = True  ## switch to implement lazy callbacks
        use_context = False  ## switch to implement context callback (implements all constraints in lazy)
        use_random_data = False  ## switch to use random datan
        use_inc = False  ## switch to use incumbent callback
        use_obj_cuts = True  ## switch to use objective cutting planes

        if tsk[0]=="enhance":
            use_lazy=True
            use_inc=True
            use_enhance=True

        use_small_area = True  ## switch to use small area
        use_ch = False  ## switch to use new convex hull based constraints
        option = 1  ### 1: Best Bound, 0: DFS (default is 1)
        no_circle = 1
        try_traditional_newtons = True  ## switch to use tradional newtons method. If false, the problem uses
        ## an alternative iterative method to get a feasible solution in each iteration
        ## instead of finding the optimal solution

        rows = 6
        columns = 7
        width = 10
        #task = ["Random", rows, columns, width]
        task = ["FLG9A"]
        data = getData(type=task, small=use_small_area, ax=20, bx=60, ay=20, by=60)
        gpars = global_param(tsk[1][0], tsk[1][1])
        #solve_newtons_method(data, gpars)




