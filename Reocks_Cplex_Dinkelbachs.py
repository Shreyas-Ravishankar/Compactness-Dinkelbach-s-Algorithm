### The aim of this code is to implement Dinkelbach’s algorithm

import numpy as np
import matplotlib.pyplot as plt
import math
import cplex
from cplex.exceptions import CplexError
from cplex.callbacks import LazyConstraintCallback as LCC
from cplex.callbacks import IncumbentCallback as ICC
from cplex.callbacks import Context as CB
import time
from Alternate_CCP_veri import CCP_onal_veri
from CCP_Knapsack import Alternative_CCP_Knapsack_optimization
from RadiusProblem import Radius_subProblem
import pandas as pd
import load_forest_data as lfd

use_warmstart = 1
no_circle = 1
nu = 8
ID = 1
U = np.array(list(range(1, nu + 1)))
area_p = 0.2
min_comp = 0
bud_p = 0.5
dtol = 0  # distance cut tolerance
mip_gap = 0.05
incumbents = {}
incumbentCounter = 0
incumbentSol = {}
incumbentCenterRadius = {}
incumbentCriticalPoints = {}
ccpIntegerOptimizeProblemsol = {}
ccpIntegerOptimizeProblemobj = 0
criticalPointsSolution = {}
heuristic_counter = 0
heuristicCriticalPoints = {}
lazyCounter = 0
lazyPoints = {}
lazyCriticalPoints = []
test = []
radiusDetails = {}
zz = 0
Trig = 0
obje = 0
prev_sol = {}
sel = {}
fractionalSol_counter = 0
u_upper_bound = 0
radiusMax = 100
radiusMin = 0
bigM = 100
initialRadius = []
criRad = []
maxRadius = 0
xcoord = {}
ycoord = {}
rad = {}
zzz = {}
x_inc_v = {}
for i in range(no_circle):
    prev_sol[i + 1] = []
# min_area=7500


# Directories
data_dir = "C:/Users/shrey/OneDrive - Arizona State University/NRES2020/Forest-Conservation-Modified/Datasets/Carvajal Datasets/"
results_dir = "C:/Users/shrey/IdeaProjects/untitled3/Results/"
# data_dir = "./Datasets/"
# results_dir = "./Results/"
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
    area_p = 0.2
    bud_p = 0.5


class getData:
    A_max = {"FLG9A": 48.6, "Hardwicke": 40,
             "Shulkell": 40, "NBCL5A": 80,
             "ElDorado": 120, "Buttercreek": 120}

    XL = {"FLG9A": 100, "Hardwicke": 175,
          "ElDorado": 430}  ## 100 each for flg9a, (175,75) for hardwicke, (430,80) for el dorado

    YL = {"FLG9A": 100, "Hardwicke": 75,
          "ElDorado": 80}

    def __init__(self, type, small, ax=0, bx=0, ay=0, by=0):
        self.type = type[0]
        self.small = small
        if self.type == "Random":
            self.rows = type[1]
            self.columns = type[2]
            self.width = type[3]
            self.nodes = list(range(1, self.rows * self.columns + 1))
            self.min_area = area_p * len(self.nodes) * self.width * self.width
            self.generate_random_data(ID)
            self.xlim = self.columns * self.width
            self.ylim = self.rows * self.width
            self.budget = bud_p * sum(list(self.cost.values()))
            self.circ = [i for i in range(1, no_circle + 1)]

        else:
            self.ax = ax
            self.ay = ay
            self.bx = bx
            self.by = by
            print(f"Getting data for {self.type}")
            self.get_actual_data(ax, bx, ay, by)
            self.circ = [i for i in range(1, no_circle + 1)]

    def generate_random_data(self, ID):
        rows = self.rows
        columns = self.columns
        width = self.width

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

        self.centers = centers.copy()
        self.vertices = vertices.copy()
        self.distance = distance.copy()
        self.neighbours = neighbours.copy()
        self.area = area.copy()
        self.circ = circ.copy()
        self.cost = cost.copy()

    def get_actual_data(self, ax, bx, ay, by):
        flgid = 5
        A_max = self.A_max[self.type]
        T = 3  # T: Time periods: The value is 3 for the given model.
        Yt = 5  # Yt: Number of years in each time period. Set this to 10 for El Dorado instance, 5 otherwise.
        if self.type == "FLG9A":
            Data = lfd.FLG9A(A_max, T, flgid)
        elif self.type == "Buttercreek":
            Data = lfd.ButterCreek(A_max, T)
        elif self.type == "NBCL5A":
            Data = lfd.NBCL5A(A_max, T)
        elif self.type == "ElDorado":
            Data = lfd.ElDorado(A_max, T)
        elif self.type == "ShulKell":
            Data = lfd.ShulKellA(A_max, T)
        elif self.type == "Hardwicke":
            Data = lfd.Hardwicke(A_max, T)

        self.DATA = Data

        min_area = Data.A_min  ### - in same format
        nodes = Data.node_set  ### - in same format
        area = Data.area_set  ### - in same format

        min_area = area_p * sum(area.values())
        print("Minimum Area ", min_area)
        print("Total ", sum(area.values()))
        print("Ratio ", min_area / sum(area.values()))

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
            print("Getting subset data")
            area, cost, nodes, neighbours, vertices, centers = self.use_subset_for_landscape(nodes, vertices,
                                                                                             neighbours, area, cost,
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
        self.cost = cost.copy()
        self.area = area.copy()
        self.nodes = nodes[:]
        self.xlim = xlim
        self.ylim = ylim
        self.min_area = min_area

    def use_subset_for_landscape(self, nodes, vertices, adj, area, cost, centers, ax, bx, ay, by):
        tbr = []
        for i in nodes:
            trueval = True
            for j in vertices[i]:
                # print("hi")
                if (j[0] < ax or j[0] > bx) or (j[1] < ay or j[1] > by):
                    trueval = False
                    # print("Excluding ", i)
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
        for v in range(1, self.V + 1):
            # print(v)
            if not visited[v - 1]:
                temp = []
                cc.append(self.DFSUtil(temp, v, visited))
            # print(v,visited)
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


def plot_random_solution(centers, width, selected, rad, centre, cost):
    print("hi")


def solve_newtons_method(data):
    centers = data.centers
    cost = data.cost
    area = data.area
    budget = data.budget
    circ = data.circ
    vertices = data.vertices
    xlim = data.xlim
    ylim = data.ylim
    min_area = data.min_area
    distance = data.distance
    nodes = data.nodes
    neighbours = data.neighbours
    filename = "holdeData.csv"
    global maxRadius
    try:
        val = pd.read_csv(filename)
        flag = True
    except:
        print("File not present")
        flag = False
        ## create the data to write
        data_to_write = {"NID": [], "R": []}
        for i in nodes:
            data_to_write["NID"].append(i)
            points = []
            for vertex in vertices[i]:
                points.append(Point(vertex[0], vertex[1]))
            circleDetails = welzl(points)
            data_to_write["R"].append(circleDetails.R)
            radiusDetails[i] = circleDetails.R
        data = pd.DataFrame(data_to_write)
        data.to_csv(filename)

    if flag:
        val = val.drop(['Unnamed: 0'], axis=1)
        for i in nodes:
            radiusDetails[i] = val.loc[val.NID == i].R.values[0]

    ## create variables names
    x_bin = {}
    xVars = []
    Radius = {}
    X = {}
    Y = {}
    radi = []
    for i in nodes:
        radi.append(radiusDetails[i])
    print("radius details", min(radi))
    print("Minimum area is", min_area)
    u_upper_bound = 1 / min_area
    for i in data.nodes:
        x_bin[i] = {}
    for j in circ:
        for i in data.nodes:
            x_bin[i][j] = "x" + str(i) + str(j)
            xVars.append(x_bin[i][j])
        Radius[j] = "r" + str(j)
        X[j] = "x" + str(j)
        Y[j] = "y" + str(j)
    circlePoints = []
    for node in nodes:
        for vertex in vertices[node]:
            circlePoints.append(Point(vertex[0], vertex[1]))
    MDetails = welzl(circlePoints)
    radiusMax = MDetails.R + 1
    bigM = MDetails.R
    maxRadius = MDetails.R
    print("Maximum radius is,", maxRadius)
    print("min radius is", (maxRadius + min(radi)) / 2)
    compact, sol, maxAreaa = Alternative_CCP_Knapsack_optimization((MDetails.C.X, MDetails.C.Y, MDetails.R), nodes, [],
                                                                   u_upper_bound)
    maxArea = maxAreaa
    u_lower_bound = 1 / (math.pi * maxRadius * maxRadius)
    calcRad = maxRadius

    onal = {}
    sel[1], minRadius = Radius_subProblem()
    onal[1] = sel[1]
    u_upper_bound = (1 / math.pi * minRadius * minRadius)
    maxComp = (maxArea) / (math.pi * minRadius * minRadius)
    print("Max compactness possible is", maxComp)

    class IncumbentAnalysis(ICC):
        def __call__(self):
            global maxRadius
            x = {}
            y = {}
            RR = {}
            selected = {}
            print("Incumbent Callback")
            zz = self.get_values('z')
            if zz > 0:
                maxRadius = ((maxArea) / (round(math.pi, 8) * zz)) ** 0.5
            for j in circ:
                x[j] = self.get_values("x" + str(j))
                y[j] = self.get_values("y" + str(j))
                RR[j] = self.get_values("r" + str(j))
                w = {}
                for i in data.nodes:
                    w[(i, j)] = round(self.get_values(x_bin[i][j]))
                selected[j] = []
                for i in data.nodes:
                    if w[(i, j)] == 1:
                        selected[j].append(i)
                if use_random_data:
                    plot_solution_new(centers, width, selected, RR, (x, y), cost)
                else:
                    data.DATA.plot_soln_data(selected, [], "Newtons Method", data.nodes, (x, y), RR)

                total_cost = 0
                total_area = {}
                for l in circ:
                    total_area[l] = 0
                for j in circ:
                    for i in data.nodes:
                        if w[(i, j)] == 1:
                            total_cost += cost[i]
                            total_area[j] += area[i]
                    print("Total Cost, Percent of budget : ", round(total_cost, 2), round(total_cost / budget, 2))
                    compact[j] = total_area[j] / ((math.pi * RR[j] * RR[j]) + 1e-20)
                tot_compact = 0
                for i in circ:
                    tot_compact += compact[i]
                print("Compactness : ", tot_compact)
                print("Radius of the circles are:", RR)
                counter = True
                pa_ar = {}
                for j in circ:
                    for i in range(len(selected[j])):
                        pa_ar[selected[j][i]] = area[selected[j][i]]
                print("Selected patch : Area", pa_ar)

    class DistanceCuts(LCC):
        compactness = min_comp

        def __init__(self, env):
            LCC.__init__(self, env)

        def neighbourhood(self, z, adj):
            S = set()
            for i in z:
                for j in adj[i]:
                    S.add(j)
            return list(S - set(z))

        def get_patches(self, z, adj):
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
            return patch

        def patch_areas(self, patch, area_set):
            area = {}
            for i, j in enumerate(patch):
                a = 0
                for k in j:
                    a += area_set[k]

                area[i] = a
            return area

        def __call__(self):
            global Trig
            global counter
            global x_inc_v
            global rad
            global xcoord
            global ycoord
            global zz
            global incumbentCenterRadius
            global incumbentCounter
            global incumbentSol
            global incumbents
            global incumbentCriticalPoints
            global ccpIntegerOptimizeProblemsol
            global ccpIntegerOptimizeProblemobj
            global lazyCounter
            global lazyPoints
            global lazyCriticalPoints
            global radiusMax
            global bigM
            global ccc
            global maxRadius
            # print("Lazy Callback")
            z = self.get_values("z")
            x = {}
            y = {}
            RR = {}
            for k in circ:
                x[k] = self.get_values("x" + str(k))
                y[k] = self.get_values("y" + str(k))
                RR[k] = self.get_values("r" + str(k))
                ask = {}
                for i in data.nodes:
                    for ind, j in enumerate(vertices[i][:-1]):
                        ask[i, ind] = ((j[0] - x[k]) ** 2 + (j[1] - y[k]) ** 2) ** 0.5

                w = {}
                for i in data.nodes:
                    w[i, k] = round(self.get_values(x_bin[i][k]))

                selected = [i for i in nodes if w[i, k] > 0]

                obj_cuts = True
                if z < round(math.pi, 4) * RR[k] * RR[k] and use_obj_cuts:
                    # print("Adding objective cuts")
                    obj_cuts = False
                    var = ["z", "r" + str(k)]
                    coeff = [1, -2 * round(math.pi, 4) * RR[k]]
                    rhs = -round(math.pi, 4) * RR[k] * RR[k]
                    self.add(constraint=cplex.SparsePair(var, coeff), sense="G", rhs=rhs)

                if RR[k] > maxRadius:
                    self.add(constraint=cplex.SparsePair(["r" + str(k)], [1]), sense="L",
                             rhs=maxRadius)
                for node in nodes:
                    for node1 in nodes:
                        if distance[node, node1] > (2 * maxRadius):
                            self.add(constraint=cplex.SparsePair([x_bin[node][k], x_bin[node1][k]], [1, 1]),
                                     sense="L", rhs=1)

                lets = True
                M = max(xlim, ylim)
                for i in selected:
                    for ind, j in enumerate(vertices[i][:-1]):
                        if ask[i, ind] > RR[k] + dtol:
                            lets = False
                            xt = x[k] - j[0]
                            yt = y[k] - j[1]
                            u1t = xt / ask[i, ind]
                            u2t = yt / ask[i, ind]
                            var = []
                            coeff1 = []
                            coeff2 = []
                            var.append("x" + str(k))
                            var.append("y" + str(k))
                            coeff1.append(u1t)
                            coeff1.append(u2t)
                            coeff2.append(-u1t)
                            coeff2.append(-u2t)
                            var.append(x_bin[i][k])
                            coeff1.append(M)
                            coeff2.append(M)
                            var.append("r" + str(k))
                            coeff1.append(-1)
                            coeff2.append(-1)
                            rhs1 = M + j[0] * u1t + j[1] * u2t
                            rhs2 = M - j[0] * u1t - j[1] * u2t
                            self.add(constraint=cplex.SparsePair(var, coeff1), sense="L", rhs=rhs1)
                            self.add(constraint=cplex.SparsePair(var, coeff2), sense="L", rhs=rhs2)

                ring = True
                patches = self.get_patches(selected, neighbours)
                if lets:
                    if len(patches) > 1:
                        # print("Adding ring inequalities")
                        self.add(constraint=cplex.SparsePair([x_bin[i][1] for i in selected], [1 for i in selected]),
                                 sense='L', rhs=len(selected) - 1)
                        area_set = self.patch_areas(patches, area)
                        for i in range(len(patches)):
                            if area_set[i] < min_area:
                                neigh = self.neighbourhood(patches[i], neighbours)
                                if neigh != []:
                                    ring = False
                                    var = []
                                    coeff = []
                                    for j in patches[i]:
                                        var.append(x_bin[j][k])
                                        coeff.append(1)
                                    for j in neigh:
                                        var.append(x_bin[j][k])
                                        coeff.append(-1)
                                    rhs = len(patches[i]) - 1
                                    self.add(constraint=cplex.SparsePair(var, coeff), sense='L', rhs=rhs)
    class GenericCB:
        def neighbourhood(self, z, adj):
            S = set()
            for i in z:
                for j in adj[i]:
                    S.add(j)
            return list(S - set(z))

        def get_patches(self, z, adj):
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
            return patch

        def patch_areas(self, patch, area_set):
            area = {}
            for i, j in enumerate(patch):
                a = 0
                for k in j:
                    a += area_set[k]

                area[i] = a
            return area

        def invoke(self, context):
            # print("Lazy Callback")
            if context.is_candidate_point():
                z = context.get_candidate_point("z")
                x = {}
                y = {}
                RR = {}
                for k in circ:
                    x[k] = context.get_candidate_point("x" + str(k))
                    y[k] = context.get_candidate_point("y" + str(k))
                    RR[k] = context.get_candidate_point("r" + str(k))
                    ask = {}
                    for i in data.nodes:
                        for ind, j in enumerate(vertices[i][:-1]):
                            ask[i, ind] = ((j[0] - x[k]) ** 2 + (j[1] - y[k]) ** 2) ** 0.5

                    w = {}
                    for i in data.nodes:
                        w[i, k] = round(context.get_candidate_point(x_bin[i][k]))

                    selected = [i for i in nodes if w[i, k] > 0]

                    obj_cuts = True
                    if z < round(math.pi, 4) * RR[k] * RR[k] and use_obj_cuts:
                        # print("Adding objective cuts")
                        obj_cuts = False
                        var = ["z", "r" + str(k)]
                        coeff = [1, -2 * round(math.pi, 4) * RR[k]]
                        rhs = -round(math.pi, 4) * RR[k] * RR[k]
                        context.reject_candidate(constraints=[cplex.SparsePair(var, coeff)], senses=["G"], rhs=[rhs])

                    lets = True
                    M = max(xlim, ylim)
                    for i in selected:
                        for ind, j in enumerate(vertices[i][:-1]):
                            if ask[i, ind] > RR[k] + dtol:
                                lets = False
                                xt = x[k] - j[0]
                                yt = y[k] - j[1]
                                u1t = xt / ask[i, ind]
                                u2t = yt / ask[i, ind]
                                var = []
                                coeff1 = []
                                coeff2 = []
                                var.append("x" + str(k))
                                var.append("y" + str(k))
                                coeff1.append(u1t)
                                coeff1.append(u2t)
                                coeff2.append(-u1t)
                                coeff2.append(-u2t)
                                var.append(x_bin[i][k])
                                coeff1.append(M)
                                coeff2.append(M)
                                var.append("r" + str(k))
                                coeff1.append(-1)
                                coeff2.append(-1)
                                rhs1 = M + j[0] * u1t + j[1] * u2t
                                rhs2 = M - j[0] * u1t - j[1] * u2t
                                context.reject_candidate(constraints=[cplex.SparsePair(var, coeff1)], senses=["L"],
                                                         rhs=[rhs1])
                                context.reject_candidate(constraints=[cplex.SparsePair(var, coeff2)], senses=["L"],
                                                         rhs=[rhs2])

                    ring = True
                    patches = self.get_patches(selected, neighbours)
                    if lets:
                        if len(patches) > 1:
                            context.reject_candidate()
                            # print("Adding ring inequalities")
                            area_set = self.patch_areas(patches, area)
                            for i in range(len(patches)):
                                if area_set[i] < min_area:
                                    neigh = self.neighbourhood(patches[i], neighbours)
                                    if neigh != []:
                                        ring = False
                                        var = []
                                        coeff = []
                                        for j in patches[i]:
                                            var.append(x_bin[j][k])
                                            coeff.append(1)
                                        for j in neigh:
                                            var.append(x_bin[j][k])
                                            coeff.append(-1)
                                        rhs = len(patches[i]) - 1
                                        context.reject_candidate(constraints=[cplex.SparsePair(var, coeff)],
                                                                 senses=['L'], rhs=[rhs])



    try:
        prob = cplex.Cplex()
        prob.parameters.mip.display.set(2)
        prob.objective.set_sense(prob.objective.sense.maximize)
        prob.parameters.clocktype.set(2)
        prob.parameters.timelimit.set(7200)
        prob.parameters.mip.tolerances.mipgap.set(0)
        prob.parameters.workmem.set(1024)
        prob.parameters.emphasis.memory.set(1)
        prob.parameters.workdir.set(
            "C:/Users/shrey/OneDrive - Arizona State University/NRES2020/Forest-Conservation-Modified/Datasets/Carvajal Datasets/")
        prob.parameters.mip.strategy.file.set(2)
        prob.parameters.mip.strategy.nodeselect.set(option)
        prob.parameters.mip.strategy.dive.set(2)  ### 0:auto, 1:never probe, 2: always probe, 3: focus more on probing

        ## using this function - distance constraints
        def add_compactness_by_rad_const(inst):
            M = max(xlim, ylim)
            inst.variables.add(obj=[0], lb=[min_area], ub=[cplex.infinity], types="C", names=["z"])
            if data.type != "Random":
                ax = data.ax
                ay = data.ay
                bx = data.bx
                by = data.by
            for k in circ:
                if not data.small or data.type == "Random":
                    inst.variables.add(obj=[0], lb=[0.0], ub=[int(((xlim) ** 2 + (ylim) ** 2) ** 0.5)], types="C",
                                       names=[Radius[k]])
                    inst.variables.add(obj=[0], lb=[0.0], ub=[xlim], types="C", names=[X[k]])
                    inst.variables.add(obj=[0], lb=[0.0], ub=[ylim], types="C", names=[Y[k]])
                else:
                    inst.variables.add(obj=[0], lb=[0.0], ub=[int(((ax - bx) ** 2 + (ay - by) ** 2) ** 0.5)], types="C",
                                       names=[Radius[k]])
                    inst.variables.add(obj=[0], lb=[ax], ub=[bx], types="C", names=[X[k]])
                    inst.variables.add(obj=[0], lb=[ay], ub=[by], types="C", names=[Y[k]])

                for i in data.nodes:
                    inst.variables.add(obj=[area[i]], lb=[0], ub=[1], types="B", names=[x_bin[i][k]])
                # print("im here 1")
                for i in data.nodes:
                    for j in vertices[i][:len(vertices[i]) - 1]:
                        for u in U:
                            var = []
                            coeff = []
                            var.append(X[k])
                            var.append(Y[k])
                            coeff.append(u1[u])
                            coeff.append(u2[u])
                            var.append(x_bin[i][k])
                            coeff.append(M)
                            var.append(Radius[k])
                            coeff.append(-1)
                            rhs = M + j[0] * u1[u] + j[1] * u2[u]
                            inst.linear_constraints.add(lin_expr=[cplex.SparsePair(var, coeff)], senses="L", rhs=[rhs],
                                                        names=["radius constraint " + str(i) + "," + str(j) + "," + str(
                                                            u) + "," + str(k)])

        ## using this function - obj func approximation
        def add_obj_cutting_planes(inst, nop, min_R):
            r0 = {}
            for k in circ:
                r0[k] = np.linspace(min_R, max(xlim, ylim), nop)
                PI = round(math.pi, 8)
                for i in r0:
                    var = ["z", "r" + str(k)]
                    coeff = [1, -2 * PI * i]
                    rhs = -PI * i * i
                    inst.linear_constraints.add(lin_expr=[cplex.SparsePair(var, coeff)], senses="G", rhs=[rhs],
                                                names=["objective_cut_" + str(i) + str(k)])

        ##
        def add_lower_upper_parametric(inst, q1, q2):
            var = ["z"]
            coe1 = [-q1]
            coe2 = [-q2]
            for k in circ:
                for i in data.nodes:
                    var.append(x_bin[i][k])
                for i in data.nodes:
                    coe1.append(area[i])
                    coe2.append(area[i])
            inst.linear_constraints.add(lin_expr=[cplex.SparsePair(var, coe1)], senses="G", rhs=[0], names=["bin 1"])
            inst.linear_constraints.add(lin_expr=[cplex.SparsePair(var, coe2)], senses="L", rhs=[0], names=["bin 2"])

        ## adding VI to strengthen model
        def add_patch_distance_const(inst):
            for k in circ:
                for i in data.nodes:
                    for j in data.nodes:
                        if j > i:
                            var = [x_bin[i][k], x_bin[j][k], "r" + str(k)]
                            coeff = [np.floor(distance[i, j]), np.floor(distance[i, j]), -4]
                            rhs = 0
                            inst.linear_constraints.add(lin_expr=[cplex.SparsePair(var, coeff)], senses="L", rhs=[rhs],
                                                        names=["VI_" + str(i) + "," + str(j) + str(k)])

        print("Adding distance constraints")
        add_compactness_by_rad_const(prob)

        print("Adding budget constraint")
        var = []
        coeff = []
        for j in circ:
            for i in data.nodes:
                var.append(x_bin[i][j])
                coeff.append(cost[i])
        prob.linear_constraints.add(lin_expr=[cplex.SparsePair(var, coeff)], senses="L", rhs=[budget], names=["budget"])

        print("Adding minimum area constraint")
        var = []
        coeff = []
        for j in circ:
            for i in nodes:
                var.append(x_bin[i][j])
                coeff.append(area[i])
        prob.linear_constraints.add(lin_expr=[cplex.SparsePair(var, coeff)], senses="G", rhs=[min_area],
                                    names=["min area"])

        print("Upper bound to total area")
        var.append('z')
        coeff.append(-1)
        prob.linear_constraints.add(lin_expr=[cplex.SparsePair(var, coeff)], senses="L", rhs=[0],
                                    names=["total area UB"])
        # print("im here 1")

        #### minimum radius constraint
        print("Adding minimum radius constraint")
        min_R = round((min_area / math.pi) ** 0.5, 8)
        for j in circ:
            prob.linear_constraints.add(lin_expr=[cplex.SparsePair([f"r{j}"], [1])], senses="G", rhs=[min_R],
                                        names=["minimum radius possible"])

        print("Adding patch dist VI")
        add_patch_distance_const(prob)
        change = 0.1 * bud_p

        ## testing new constraint based on convex hull
        var = []
        coeff11 = []
        coeff12 = []
        if no_circle <= 1 and use_ch:
            print("Testing convex hull constraints")
            for i in nodes:
                for k, j in enumerate(vertices[i][:-1]):
                    prob.variables.add(obj=[0], lb=[0.0], ub=[1], types="C", names=[f"lamb_{i},{k}"])
                    var.append(f"lamb_{i},{k}")
                    coeff11.append(j[0])
                    coeff12.append(j[1])

                    prob.linear_constraints.add(lin_expr=[cplex.SparsePair([f"lamb_{i},{k}", x_bin[i][1]], [1, -1])],
                                                senses="L", rhs=[0],
                                                names=[f"chb_{i},v_{k}"])
            var1 = var[:]
            var2 = var[:]
            var1.append(X[1])
            var2.append(Y[1])
            coeff11.append(-1)
            coeff12.append(-1)
            prob.linear_constraints.add(lin_expr=[cplex.SparsePair(var1, coeff11)], senses="E", rhs=[0], names=["ch_x"])
            prob.linear_constraints.add(lin_expr=[cplex.SparsePair(var2, coeff12)], senses="E", rhs=[0], names=["ch_y"])
            coeff = [1 for i in var]
            prob.linear_constraints.add(lin_expr=[cplex.SparsePair(var, coeff)], senses="E", rhs=[1], names=["ch_sum"])

        def wasta01(inst, nodes):
            global incumbentCounter
            global incumbentSol
            global incumbents
            global x_inc_v
            global xcoord
            global ycoord
            global rad
            global zzz
            xxx = []
            nonal = []
            nvonal = []
            onalArea = 0
            k = 0
            pointsForCircle = []
            # print("Onalllllllllllllllllllllllllllllllllllllllllllllllllllllll", onal)
            for i in onal:
                for j in onal[i]:
                    nonal.append(x_bin[j][i])
                    nvonal.append(1)
                    onalArea += area[j]
            for circle in circ:
                for node in nodes:
                    if x_bin[node][circle] not in nonal:
                        nonal.append(x_bin[node][circle])
                        nvonal.append(0)
            print("Onal Solution Area is", onalArea)
            print("In verification model")
            varr, vall, sol = CCP_onal_veri(nonal, nvonal, radiusMax, minRadius, True)
            var = []
            val = []
            axi = {}
            j = 0
            for i in varr:
                if i in xVars:
                    var.append(i)
                    xxx.append(i)
                    val.append(vall[j])
                    j += 1
                elif i == "k11":
                    axi[i] = vall[j]
                    j += 1
                elif i == "k21":
                    axi[i] = vall[j]
                    j += 1
                elif i == "k31":
                    axi[i] = vall[j]
                    j += 1
                elif i == "u1":
                    axi[i] = vall[j]
                    j += 1
                elif i == "z":
                    axi['z'] = vall[j]
                    j += 1
            var.append("x1")
            val.append((axi['k11'] / axi['u1']))
            var.append("y1")
            val.append((axi['k21'] / axi['u1']))
            var.append("r1")
            val.append((axi['k31'] / axi['u1']))
            var.append("z")
            val.append(axi['z'])
            if var != False:
                incumbentCounter += 1
                incumbentSol[incumbentCounter] = sol
                incumbents[incumbentCounter] = (1, onal[1])
                j = 0
                for i in var:
                    if i == "x1":
                        xcoord[i] = val[j]
                        j += 1
                    elif i == "y1":
                        ycoord[i] = val[j]
                        j += 1
                    elif i == "r1":
                        rad[i] = val[j]
                        j += 1
                    elif i in xVars:
                        x_inc_v[i] = val[j]
                        j += 1
                    elif i == "z":
                        zzz[i] = val[j]
                        j += 1
                inst.MIP_starts.add([var, val], prob.MIP_starts.effort_level.no_check, "Warmstart1")
            else:
                print(
                    "Violating some constraint!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")

        if use_lazy:
            prob.register_callback(DistanceCuts)
        if use_inc:
            prob.register_callback(IncumbentAnalysis)
        if use_obj_cuts:
            print("Adding obj func approximation")
            add_obj_cutting_planes(prob, 50, min_R)
        if use_context:
            print("Using context callback")
            cm = CB.id.candidate
            prob.set_callback(GenericCB(), cm)

        ### iterative method starts here
        Q = []
        F = []
        q = 0.0
        prob.objective.set_linear('z', -q)
        t1 = time.time()
        gap = 0.2
        tgap = 0.01
        prob.parameters.mip.tolerances.mipgap.set(gap)
        wasta01(prob, nodes)
        prob.solve()
        Feval = prob.solution.get_objective_value()
        # prob.parameters.mip.limits.solutions.set(1)
        print("Iteration 1 : Feval,Q = ", Feval, q)
        print(f"Adjusted gap {gap}")
        ctr = 0
        itol = 0.0001
        expr = True

        ## while loop iterates until a problem is infeasible.
        while expr:  # Feval>=itol or Feval<=-itol:# and ctr<=10:
            prevq = q
            ctr += 1
            oldSol = prob.solution.get_values()
            oldProb = prob
            OldF = prob.solution.get_objective_value()
            Q.append(q)
            F.append(OldF)
            selected = {}
            for k in circ:
                selected[k] = []
                for i in data.nodes:
                    # print("here1111111")
                    count = prob.solution.get_values(x_bin[i][k])
                    if round(count) > 0:
                        selected[k].append(i)
            total_area = 0
            radius = {}
            x = {}
            y = {}
            for k in circ:
                for i in selected[k]:
                    total_area += area[i]
                radius[k] = prob.solution.get_values("r" + str(k))
                x[k] = prob.solution.get_values("x" + str(k))
                y[k] = prob.solution.get_values("y" + str(k))
            q = total_area / prob.solution.get_values('z')
            if use_random_data:
                plot_solution_new(centers, width, selected[1], radius[1], (x[1], y[1]), cost)
            else:
                data.DATA.plot_soln_data(selected, [], f"NEW_{q}", data.nodes, (x, y), radius)

            print("Starting next iteration")
            prob.objective.set_linear('z', -q)

            if not try_traditional_newtons:
                print(
                    "Adding cut that should force lower bound to be strictly positive and the objective to be atleast "
                    "tgap % better than last objective")
                varo = []
                coeffo = []
                for j in circ:
                    for i in nodes:
                        varo.append(x_bin[i][j])
                        coeffo.append(area[i])
                var = varo[:]
                coeff = coeffo[:]
                var.append('z')
                coeff.append(-q - tgap * q)
                prob.linear_constraints.add(lin_expr=[cplex.SparsePair(var, coeff)], senses="G", rhs=[0],
                                            names=["Objective cut"])
                prob.parameters.mip.limits.solutions.set(1)  ## stopping at first solution itself
            else:
                prob.parameters.mip.tolerances.mipgap.set(0.05)

            prob.solve()
            if try_traditional_newtons:
                Feval = prob.solution.get_objective_value()
                expr = (Feval >= itol) or (Feval <= -itol)
            else:
                expr = "infeas" not in prob.solution.get_status_string()

        t2 = time.time()
        print("Time Elapsed : ", round(t2 - t1, 4))
        print(f"Q (compactness) = {Q}")
        print(f"F (diff in den num) = {F}")

        #################--------------------- final solution ------------------##########################
        prob = oldProb
        selected = {}
        x = {}
        y = {}
        radius = {}
        flag = 0
        for k in circ:
            selected[k] = []
            for i in data.nodes:
                count = prob.solution.get_values(x_bin[i][k])
                if round(count) > 0:
                    selected[k].append(i)

                radius[k] = prob.solution.get_values(Radius[k])
                x[k] = prob.solution.get_values(X[k])
                y[k] = prob.solution.get_values(Y[k])

        if use_random_data:
            plot_solution_new(centers, width, selected[1], radius, (x[1], y[1]), cost)
        else:
            data.DATA.plot_soln_data(selected, [], "NEW", data.nodes, (x, y), radius)

        total_area = 0
        total_cost = 0
        for k in circ:
            for i in selected[k]:
                total_area += area[i]
                total_cost += cost[i]
            radius[k] = prob.solution.get_values("r" + str(k))
            Z = prob.solution.get_values("z")
            compact = total_area / (math.pi * radius[k] * radius[k])
            print("This is my solution :")
            print("Area of circle : ", (math.pi * radius[k] * radius[k]))
            print("value of z : ", Z)
            print(f"Total area is {total_area}")
            print("Selected nodes : ", selected[k])
            print("Compactness : ", compact)
            print("Budget utilized : {0}% of allowed".format(100 * total_cost / budget))
            print("Radius is : ", radius[k])
            print("Time Elapsed : ", round(t2 - t1, 4))


    except CplexError as e:
        print(e)


if __name__ == "__main__":
    print("Starting code")

    use_lazy = True  ## switch to implement lazy callbacks
    use_context = not use_lazy  ## switch to implement context callback (implements all constraints in lazy)
    use_random_data = False  ## switch to use random data
    use_inc = False and use_lazy  ## switch to use incumbent callback
    use_obj_cuts = True  ## switch to use objective cutting planes
    use_small_area = True  ## switch to use small area
    use_ch = True  ## switch to use new convex hull based constraints
    option = 1  ### 1: Best Bound, 0: DFS (default is 1)
    no_circle = 1
    try_traditional_newtons = True  ## switch to use tradional newtons method. If false, the problem uses
    ## an alternative iterative method to get a feasible solution in each iteration
    ## instead of finding the optimal solution

    rows = 6
    columns = 7
    width = 10
    task = ["Random", rows, columns, width]
    task = ["Hardwicke"]
    data = getData(type=task, small=use_small_area, ax=00, bx=80, ay=0, by=80)
    solve_newtons_method(data)