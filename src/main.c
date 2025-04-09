#include "../include/main.h"

int main(int argc, char **argv)
{
    char *path_mesh = NULL;
    for (int index = 0; index < argc; ++index)
    {
        if (strstr("-mesh", argv[index]))
        {
            path_mesh = argv[index + 1];
        }
    }

    DataMesh data_mesh;
    FileProcessMesh(path_mesh, &data_mesh);

    puts("\n==== mesh data information ====");
    printf("number of nodes: %d\n", data_mesh.nn);
    printf("number of elements: %d\n", data_mesh.ne);
    for (int index = 0; index < data_mesh.ne; ++index)
    {
        printf("element %d: ", data_mesh.idx_ele[index]);
        for (int index_i = 0; index_i < data_mesh.nn_ele[index]; ++index_i)
        {
            printf("%d\t", data_mesh.idx_node[index][index_i]);
        }
        putchar('\n');
    }

    puts("\n==== test metis function ====");
    TestFunctionMetis(data_mesh);

    // free memory
    free(data_mesh.idx_ele);
    free(data_mesh.nn_ele);
    for (int index = 0; index < data_mesh.ne; ++index)
    {
        free(data_mesh.idx_node[index]);
    }
    free(data_mesh.idx_node);

    return 0;
}

// command line
/*
 * ./app_metis_exe -mesh </path/to/mesh/file>
 */
