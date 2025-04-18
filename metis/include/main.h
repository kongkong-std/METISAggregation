#ifndef MAIN_H_
#define MAIN_H_

#define BUF_MAX_SIZE 2048

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <metis.h>
#include <GKlib.h>
#include <cjson/cJSON.h>
#include <assert.h>

// struct
/*
 * mesh struct
 */
typedef struct
{
    /* data */
    int nn, ne;     // number of elements, nodes
    int *idx_ele;   // element index
    int *nn_ele;    // number of nodes in element
    int **idx_node; // node index in element
} DataMesh;

/*
 * gmsh struct
 */
typedef struct
{
    /* data */
    int nn, ne;               // number of nodes, elements
    int ne_bd, ne_in;         // number of elements of boundary, inner, ne = ne_bd + ne_in
    int nne_bd, nne_in;       // number of nodes in each element of boudary, inner
    idx_t nparts;             // number of partitions
    double *coordinates;      // coordinates of nodes [x1, y1, z1, x2, y2, z2, ...]
    idx_t *eptr_bd, *eind_bd; // csr mesh connectivity of boundary element
    idx_t *eptr_in, *eind_in; // csr mesh connectivity of inner element
    idx_t *npart_in;          // nodes partition of inner element
    idx_t *epart_in;          // elements partition of inner element
} DataGmsh;

typedef enum
{
    NONE,           // 0
    MESH_FORMAT,    // 1
    PHYSICAL_NAMES, // 2
    NODES,          // 3
    ELEMENTS        // 4
} Flag_Data_Block;

// function prototype
void TestMetis(void);    // metis demo test function
int TestMetisFunctionGraph(DataGmsh data /*gmsh data*/);
void GmshCoarseLevelGenerator(DataGmsh *coarse_data /*gmsh coarse level data pointer*/,
                              DataGmsh *fine_data /*gmsh fine level data pointer*/);
int TestMetisFunctionGmsh(DataGmsh data /*gmsh data*/);
int NumNodeEleTypeMap(int type /*element type*/);
void FileProcessGmsh(const char *path /*path to gmsh file*/,
                     DataGmsh *data /*gmsh data pointer*/);
void FileProcessMesh(const char *path /*path to mesh file*/,
                     DataMesh *data /*mesh data pointer*/);
int TestFunctionMetis(DataMesh data /*mesh data*/);

#endif // main.h
