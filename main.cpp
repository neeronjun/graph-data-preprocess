//=======================================================================
// Copyright 2007 Aaron Windsor
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//=======================================================================
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <math.h>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/boyer_myrvold_planar_test.hpp>
#include <bits.h>

#define PI 22/7

using namespace std;

double calculateDistance(double lat1, double lon1, double lat2, double lon2, char unit){
    double radlat1 = PI * lat1 / 180;
    double radlat2 = PI * lat2 / 180;
    double radlon1 = PI * lon1 / 180;
    double radlon2 = PI * lon2 / 180;
    double theta = lon1 - lon2;
    double radtheta = PI * theta / 180;
    double dist = sin(radlat1) * sin(radlat2) + cos(radlat1) * cos(radlat2) * cos(radtheta);
    dist = acos(dist);
    dist = dist * 180 / PI;
    dist = dist * 60 * 1.1515;
    if (unit == 'K') { dist = dist * 1.609344; }
    if (unit == 'N') { dist = dist * 0.8684; }
    return dist;
}

struct Location
{
    double x, y;
};

// Compares two intervals according to staring times. 
int compareInterval (Location *l1, Location *l2)
{
    return (l1->x < l2->x);
}

int compare(const void* a, const void* b)
{
    return (*(double*)a - *(double*)b);
}

void sort_array() {

}

void nearest_neighbor(double full_array, int v_id, int* rslt, int nbr_cnt) {
    for (size_t i = 0; i < nbr_cnt; i++)
    {
       //calculateDistance()
    }

}

/* C++ implementation of QuickSort */
using namespace std;

// A utility function to swap two elements  
void swap(double* a, double* b)
{
    double t = *a;
    *a = *b;
    *b = t;
}

/* This function takes last element as pivot, places
the pivot element at its correct position in sorted
array, and places all smaller (smaller than pivot)
to left of pivot and all greater elements to right
of pivot */
double partition(double arr[], int low, int high)
{
    double pivot = arr[2*high]; // pivot  
    int i = (low - 1); // Index of smaller element  

    for (int j = low; j <= high - 1; j++)
    {
        // If current element is smaller than the pivot  
        if (arr[2*j] < pivot)
        {
            i++; // increment index of smaller element  
            swap(&arr[2 * i], &arr[2 * j]);
            swap(&arr[2 * i + 1], &arr[2 * j + 1]);
        }
    }
    swap(&arr[2*(i + 1)], &arr[2*high]);
    swap(&arr[2 * (i + 1)+1], &arr[2 * high+1]);
    return (i + 1);
}

/* The main function that implements QuickSort
arr[] --> Array to be sorted,
low --> Starting index,
high --> Ending index */
void quickSort(double arr[], int low, int high)
{
    if (low < high)
    {
        /* pi is partitioning index, arr[p] is now
        at right place */
        double pi = partition(arr, low, high);

        // Separately sort elements before  
        // partition and after partition  
        quickSort(arr, low, pi - 1);
        quickSort(arr, pi + 1, high);
    }
}

/* Function to print an array */
void printArray(double arr[], int size)
{
    int i;
    for (i = 0; i < size; i++)
        cout << arr[i] << " ";
    cout << endl;
}

//cantor pairing function that generates a unique id
int get_id(int prev, int curr) {
    return (prev + curr) * (prev + curr + 1) / 2 + prev;
}

bool e_info[9000000];
bool visited(int prev, int next) {
    bool vis = false;
    if (e_info[get_id(prev, next)] == 1)		//|| e_info[next * 10 + prev] == 1)
        vis = true;
    return vis;
}

int get_current(int prev, vector<int> adj[]) {
    int curr = -1;
    auto n = adj[prev];
    for (int i = 0; i < n.size(); i++)
    {
        //cout << "n[i]: " << n[i] << endl;
        if (!(visited(prev, n[i]) || visited(n[i], prev))) {
            curr = n[i];
            break;
        }
    }
    return curr;
}

int possible_path(int v, vector<int> adj[]) {
    int path = -1;
    auto n = adj[v];
    for (int i = 0; i < n.size(); i++)
    {
        //cout << "n[i]: " << n[i] << endl;
        if (!(visited(v, n[i]))) {
            path = n[i];
            break;
        }
    }
    return path;
}

int depair(int z) {
    long t = (int)(floor((sqrt(8 * z + 1) - 1) / 2));
    long x = t * (t + 3) / 2 - z;
    return (int)x;
}

bool edge_exist(vector<int> adj[], int u, int v) {
    bool res = false;
    for (auto x : adj[u])
    {
        printf("x: %d\t", x);
        if (x == v)
            res = true;
    }
    printf("\n");
    return res;
}

// add the edges to the graph 
void addEdge(vector<int> adj[], int u, int v)
{
    if (u == v) {
        return;
    }
    else if (!edge_exist(adj, u, v)) {
        adj[u].push_back(v);
        adj[v].push_back(u);
    }
}

// A utility function to print the adjacency list 
// representation of graph 
void printGraph(vector<int> adj[], int V)
{
    for (int v = 0; v < V; ++v)
    {
        cout << "\n Adjacency list of vertex "
            << v << "\n head ";
        for (auto x : adj[v])
            cout << "-> " << x;
        printf("\n");
    }
}

//returns true if the point v3 lies left of line (v1,v2)
bool isLeft(int v1, int v2, int v3, double v_xy[]) {
    double angle = 0, x[3], y[3], a[2], b[2], a1 = 0, b1 = 0;
    x[0] = v_xy[v1 * 2], y[0] = v_xy[v1 * 2 + 1];
    x[1] = v_xy[v2 * 2], y[1] = v_xy[v2 * 2 + 1];
    x[2] = v_xy[v3 * 2], y[2] = v_xy[v3 * 2 + 1];
    return ((x[1] - x[0]) * (y[2] - y[0]) - (y[1] - y[0]) * (x[2] - x[0])) > 0;
}

double get_angle(int v1, int v2, int v3, double v_xy[10]) {
    double angle = 0, x[3], y[3], a[2], b[2], a1 = 0, b1 = 0;
    x[0] = v_xy[v1 * 2], y[0] = v_xy[v1 * 2 + 1];
    x[1] = v_xy[v2 * 2], y[1] = v_xy[v2 * 2 + 1];
    x[2] = v_xy[v3 * 2], y[2] = v_xy[v3 * 2 + 1];
    //cout << "x[0], y[0], x[1], y[1], x[2], y[2]" << x[0]<< y[0]<< x[1]<< y[1]<< x[2] << y[2] << endl;
    a[0] = (x[1] - x[0]), a[1] = (y[1] - y[0]);
    b[0] = (x[1] - x[2]), b[1] = (y[1] - y[2]);
    a1 = sqrt(a[0] * a[0] + a[1] * a[1]);
    b1 = sqrt(b[0] * b[0] + b[1] * b[1]);
    angle = acosf((a[0] * b[0] + a[1] * b[1]) / (a1 * b1)) * 180.0 / PI;
    if (isLeft(v1, v2, v3, v_xy)) {
        angle = 360 - angle;
    }

    return angle;
}

void find_graphUnit(vector<int> adj[], double *v_xy)
{
    int start = 0, hop = 0, min[2], max[2], prev, next, curr, e_id;
    int ucomm_arr[5000], comm_arr[5000], ucomm_id = -1, comm_id = -1, ucount = 0;

    for (int j = 0; j < 8; j++) {
        start = j;
        //cout << "start: " << start << endl;
        while (possible_path(start, adj) != -1) {
            max[0] = max[1] = 0;
            min[0] = min[1] = 100;
            int temp_e[10];
            hop = 0;
            prev = start;
            next = -1;
            if (possible_path(prev, adj) == -1) {
                //cout << "prev= " << prev << endl;
                //cout << "no next vertex" << endl;
                continue;
            }
            else {
                curr = possible_path(prev, adj);
                //cout << "heres current: " << curr << endl;
            }
            cout << "visited: " << prev << "," << curr << endl;
            e_id = get_id(prev, curr);
            cout << "added: " << e_id << " at " << hop << endl;
            e_info[e_id] = 1;
            temp_e[hop] = e_id;
            //cout << "eid = " << get_id(prev, curr) << endl;
            auto n1 = adj[curr];

            while (next != start) { //&& hop <= 6) {
                double min_ang = 360;
                for (int i = 0; i < n1.size(); i++) {
                    double ang = 360;
                    //cout << "-> " << n1[i] << endl;
                    //cout << "prev= " << prev << ", curr= " << curr << ", next = " << n1[i] << endl;
                    if (visited(curr, n1[i])) {
                        ang = 361;
                        //cout << "already visited.." << endl;
                    }
                    else
                        ang = get_angle(prev, curr, n1[i], v_xy);
                    //cout << ang << endl;
                    if (ang < min_ang && ang != 0) {
                        min_ang = ang;
                        next = n1[i];
                    }

                }

                if (min_ang > 180) {
                    cout << "no right vertex to go" << endl;
                    hop = 0;
                    //cout << "prev= " << prev << ", curr= " << curr << ", next = " << next << endl;
                    //e_info[get_id(prev, curr)] = 0;
                    cout << min_ang << endl;
                    break;
                }
                prev = curr;
                curr = next;
                //cout << "prev= " << prev << ", curr= " << curr << ", next = " << next << endl;
                e_id = get_id(prev, curr);
                cout << "visited: " << prev << "," << curr << endl;
                e_info[e_id] = 1;
                temp_e[++hop] = e_id;
                cout << "added: " << e_id << " at " << hop << endl;
                //cout << "eid = " << e_id << endl;
                //cout << next << endl;
                n1 = adj[curr];
                //hop++;
            }
            cout << "hops = " << hop << endl;
            if (hop > 1) {
                ucomm_id++;
                ucomm_arr[ucount++] = ucomm_id;
                ucomm_arr[ucount++] = hop + 1;
                cout << "vertices: " << hop + 1 << endl;
                for (int k = 0; k < hop + 1; k++)
                {
                    int v_id = depair(temp_e[k]);
                    int tx = v_xy[v_id * 2], ty = v_xy[v_id * 2 + 1];
                    cout << tx << ", " << ty << endl;
                    if (tx > max[0])
                        max[0] = tx;
                    if (ty > max[1])
                        max[1] = ty;
                    if (tx < min[0])
                        min[0] = tx;
                    if (ty < min[1])
                        min[1] = ty;
                    cout << depair(temp_e[k]) << endl;
                    ucomm_arr[ucount++] = temp_e[k];

                    cout << depair(temp_e[k]) << "x: " << v_xy[depair(temp_e[k]) * 2] << ", y: " << v_xy[2 * depair(temp_e[k]) + 1] << " sss ";
                }
                cout << "minx: " << min[0] << ", miny: " << min[1] << "maxx: " << max[0] << ", maxy: " << max[1];
                cout << endl;
            }


        }
    }
}

// Driver Code 
/*int main()
{
    double arr[] = { 22, 4, 12, 33, 5, 2, 2, 1, 9.1, 5, 3, 5, 5, 63, 23, 55, 3, 0};
    int n = sizeof(arr) / sizeof(arr[0]);
    int count, max, min;
    //random_generate(50);
    quickSort(arr, 0, n / 2 -1);
    for (size_t j = 0; j < n; j++)
    {
        count = 0;
        int deg = 1 + rand() % 4;
        printf("edge_count: %d\n",deg);
        for (size_t i = 0; i < deg; i++)
        {
            max = deg + j;
            min = j - deg;
            if (min < 0)
                min = 0;
            int randNum = rand() % (max - min + 1) + min;
            if(randNum < n)
                printf("start: %d, end: %d\n", j, randNum);
        }
    }
    
    
    //quickSort(arr, 0, n - 1);
    cout << "Sorted array: \n";
    printArray(arr, n);
    return 0;
}*/

/*int main(int argc, char** argv) {
    int vertex_count = 100, v_id;
    double v_x, v_y, *full_array = new double[vertex_count*2];
    Location *loc = new Location[vertex_count];

    for (size_t i = 0; i < vertex_count; i++)
    {
        v_id = i;
        float x = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        float y = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        printf("x:%lf, y:%lf\n", x * 1000, y * 1000);
        full_array[i * 2] = x * 1000;
        full_array[i * 2 + 1] = y * 1000;
        loc[i].x = x * 1000;
        loc[i].y = y * 1000;
        //qsort(full_array, vertex_count, sizeof(double), compare);
    }
    qsort(full_array, vertex_count, sizeof(double), compare);
    
    for (size_t i = 0; i < vertex_count; i++)
    {
        printf("a: %lf\n", full_array[i]);
    }
    delete[] full_array, loc;
    return 0;
}*/

int main(int argc, char** argv)
{
    const int V = 15;
    fstream fin, fin1, fout;
    fin.open("Data topk/cali_node_norm.txt", ios::in);
    fin1.open("Data topk/Cali_Edge_info.txt", ios::in);

    double *spat_arr;
    // This program illustrates a simple use of boyer_myrvold_planar_embedding
    // as a simple yes/no test for planarity.

    using namespace boost;

    typedef adjacency_list<vecS,
        vecS,
        undirectedS,
        property<vertex_index_t, int>
    > graph;
    int n = 15, min = 0 , max = 0;
    spat_arr = new double[n * 2];
    
    graph GEN(n);
    vector<int> adj[V];

    for (size_t i = 0; i < n; i++)
    {
        float x = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        float y = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        printf("x:%lf, y:%lf\n", x * 1000, y * 1000);
        spat_arr[i * 2] = x * 1000;
        spat_arr[i * 2 + 1] = y * 1000;
    }

    //qsort(spat_arr, n, sizeof(double), compare);
    quickSort(spat_arr, 0, n - 1);

    for (size_t i = 0; i < n; i++)
    {
        printf("v_id: %d, x: %lf, y: %lf \n",i, spat_arr[i * 2], spat_arr[i * 2 + 1]);
    }
    for (size_t j = 0; j < n; j++)
    {
       
        int deg = 1 + rand() % 3;
        printf("edge_count: %d\n", deg);
        for (size_t i = 0; i < deg; i++)
        {
            max = deg + j;
            min = j - deg;
            if (min < 0)
                min = 0;
            int randNum = rand() % (max - min + 1) + min;
            if (randNum < n) {
                printf("start: %d, end: %d\n", j, randNum);
                add_edge(j, randNum, GEN); //boost add edge function
                addEdge(adj, j, randNum);  //vector add edge function
                printGraph(adj, n);
            }
        }
    }
    printGraph(adj, n);
    edge_exist(adj, 0, 2);
    
    if (boyer_myrvold_planarity_test(GEN))
        std::cout << "GEN is planar." << std::endl;
    else
        std::cout << "ERROR! K_4 should have been recognized as planar!"
        << std::endl;

    return 0;
}