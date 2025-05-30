#ifndef MAIN_H_
#define MAIN_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <assert.h>
#include <mpi.h>
#include <metis.h>
#include <parmetis.h>

#define MAX_BUF_SIZE 2048

// data struct
typedef struct
{
    /* data */
    int ele_idx;
    int ele_node[4]; // prescribed 4 nodes in each element
} DataMeshEle;

typedef struct
{
    /* data */
    int nn, ne;       // number of nodes, elements
    DataMeshEle *ele; // element data
} DataMesh;

// function prototype
int FileProcessMesh(const char *path /*path to mesh file*/,
                    DataMesh *data /*mesh data*/);

#endif // main.h