#include "../include/main.h"

int main(int argc, char **argv)
{
    char *path_mesh = NULL;
    char *type = NULL;
    for (int index = 0; index < argc; ++index)
    {
        if (strstr("-mesh", argv[index]))
        {
            path_mesh = argv[index + 1];
        }
        if (strstr("-type", argv[index]))
        {
            type = argv[index + 1];
        }
    }

    if (strcmp(type, "gmsh") == 0)
    {
        DataGmsh data_gmsh;
        FileProcessGmsh(path_mesh, &data_gmsh);

        puts("\n==== gmsh data information ====");
        printf("number of nodes: %d\n", data_gmsh.nn);
        printf("number of elements: %d\n", data_gmsh.ne);
        for (int index = 0; index < data_gmsh.nn; ++index)
        {
            printf("node %d: %021.16le\t%021.16le\t%021.16le\n", index,
                   data_gmsh.coordinates[3 * index],
                   data_gmsh.coordinates[3 * index + 1],
                   data_gmsh.coordinates[3 * index + 2]);
        }
        printf("number of elements in boundary: %d\n", data_gmsh.ne_bd);
        for(int index = 0; index < data_gmsh.ne_bd + 1; ++index)
        {
            printf("data_gmsh.eptr_bd[%d] = %ld\n", index, data_gmsh.eptr_bd[index]);
        }
        printf("number of elements in inner: %d\n", data_gmsh.ne_in);
        for(int index = 0; index < data_gmsh.ne_in + 1; ++index)
        {
            printf("data_gmsh.eptr_in[%d] = %ld\n", index, data_gmsh.eptr_in[index]);
        }
        printf("nodes in inner elements:\n");
        for(int index = 0; index < data_gmsh.ne_in; ++index)
        {
            printf("element %d: ", index);
            for(int index_i = 0; index_i < data_gmsh.nne_in; ++index_i)
            {
                printf("%ld\t", data_gmsh.eind_in[index * data_gmsh.nne_in + index_i]);
            }
            putchar('\n');
        }

        // free memory
        free(data_gmsh.coordinates);
    }
    else
    {
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
    }

    return 0;
}

// command line
/*
 * ./app_metis_exe -mesh </path/to/mesh/file>
 *                 -type <mesh file type: gmsh/>
 */
