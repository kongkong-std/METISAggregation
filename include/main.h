#ifndef MAIN_H_
#define MAIN_H_

#define BUF_MAX_SIZE 2048

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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

// function prototype
void FileProcessMesh(const char *path /*path to mesh file*/,
                     DataMesh *data /*mesh data pointer*/);
int TestFunctionMetis(DataMesh data /*mesh data*/);

#endif // main.h
