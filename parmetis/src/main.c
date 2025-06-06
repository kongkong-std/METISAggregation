#include "../include/main.h"

int main(int argc, char **argv)
{
    int my_rank, nprocs;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    MPI_Comm comm = MPI_COMM_WORLD;

    char *path_mesh = NULL;
    for (int index = 0; index < argc; ++index)
    {
        if (strstr("-mesh", argv[index]))
        {
            path_mesh = argv[index + 1];
        }
    }

    DataMesh mesh_data;

    if (my_rank == 0)
    {
        printf("mesh file path: %s\n", path_mesh);
        FileProcessMesh(path_mesh, &mesh_data);

#if 0
        for (int index = 0; index < mesh_data.nn; ++index)
        {
            printf("node %d:\t%021.16le\t%021.16le\t%021.16le\n", index,
                   mesh_data.coordinates[mesh_data.dim * index],
                   mesh_data.coordinates[mesh_data.dim * index + 1],
                   mesh_data.coordinates[mesh_data.dim * index + 2]);
        }

        for (int index = 0; index < mesh_data.ne_shell; ++index)
        {
            printf("element %d:\t", index);
            for (int index_i = 0; index_i < mesh_data.ele_shell[index].num_ele_node; ++index_i)
            {
                printf("%d\t", mesh_data.ele_shell[index].ele_node[index_i]);
            }
            putchar('\n');
        }
#endif // mesh data, node and element
    }

    MPI_Bcast(&mesh_data.nn, 1, MPI_INT, 0, comm);
    MPI_Bcast(&mesh_data.dim, 1, MPI_INT, 0, comm);
    MPI_Bcast(&mesh_data.ne_shell, 1, MPI_INT, 0, comm);

    if (my_rank != 0)
    {
        mesh_data.coordinates = (double *)malloc(mesh_data.dim * mesh_data.nn * sizeof(double));
        mesh_data.ele_shell = (DataMeshEle *)malloc(mesh_data.ne_shell * sizeof(DataMeshEle));
        assert(mesh_data.coordinates && mesh_data.ele_shell);
    }
    MPI_Bcast(mesh_data.coordinates, mesh_data.dim * mesh_data.nn, MPI_DOUBLE, 0, comm);
    MPI_Bcast(mesh_data.ele_shell, mesh_data.ne_shell * sizeof(DataMeshEle), MPI_BYTE, 0, comm);

    if (my_rank != 0)
    {
        for (int index = 0; index < mesh_data.ne_shell; ++index)
        {
            mesh_data.ele_shell[index].ele_node = (int *)malloc(mesh_data.ele_shell[index].num_ele_node * sizeof(int));
            assert(mesh_data.ele_shell[index].ele_node);
        }
    }

    for (int index = 0; index < mesh_data.ne_shell; ++index)
    {
        MPI_Bcast(mesh_data.ele_shell[index].ele_node, mesh_data.ele_shell[index].num_ele_node, MPI_INT, 0, comm);
    }

#if 1
    for (int index_p = 0; index_p < nprocs; ++index_p)
    {
        MPI_Barrier(comm);
        if (my_rank == index_p)
        {
            printf("in rank %d/%d:\n", my_rank, nprocs);
            puts("nodes information");
            for (int index = 0; index < mesh_data.nn; ++index)
            {
                printf("node %d:\t%021.16le\t%021.16le\t%021.16le\n", index,
                       mesh_data.coordinates[mesh_data.dim * index],
                       mesh_data.coordinates[mesh_data.dim * index + 1],
                       mesh_data.coordinates[mesh_data.dim * index + 2]);
            }

            puts("\nelements information");
            for (int index = 0; index < mesh_data.ne_shell; ++index)
            {
                printf("element %d:\t", index);
                for (int index_i = 0; index_i < mesh_data.ele_shell[index_i].num_ele_node; ++index_i)
                {
                    printf("%d\t", mesh_data.ele_shell[index].ele_node[index_i]);
                }
                putchar('\n');
            }

            puts("\n\n========\n");
        }
    }
#endif // print as my_rank ascending order

    // free memory
    for (int index = 0; index < mesh_data.ne_shell; ++index)
    {
        free(mesh_data.ele_shell[index].ele_node);
    }
    free(mesh_data.ele_shell);
    free(mesh_data.coordinates);

    MPI_Finalize();
    return 0;
}

#if 0
int main(int argc, char **argv)
{
    TestMetis();

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

#if 1
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
        for (int index = 0; index < data_gmsh.ne_bd + 1; ++index)
        {
            printf("data_gmsh.eptr_bd[%d] = %" PRIDX "\n", index, data_gmsh.eptr_bd[index]);
        }
        printf("number of elements in inner: %d\n", data_gmsh.ne_in);
        for (int index = 0; index < data_gmsh.ne_in + 1; ++index)
        {
            printf("data_gmsh.eptr_in[%d] = %" PRIDX "\n", index, data_gmsh.eptr_in[index]);
        }
        printf("nodes in inner elements:\n");
        for (int index = 0; index < data_gmsh.ne_in; ++index)
        {
            printf("element %d: ", index);
            for (int index_i = 0; index_i < data_gmsh.nne_in; ++index_i)
            {
                printf("%" PRIDX "\t", data_gmsh.eind_in[index * data_gmsh.nne_in + index_i]);
            }
            putchar('\n');
        }
#endif // gmsh data information

        TestMetisFunctionGmsh(data_gmsh);

#if 1
        puts("==== partition information ====");
        puts("element partitions:");
        for (int index = 0; index < data_gmsh.ne_in; ++index)
        {
            printf("data_gmsh.epart_in[%d] = %" PRIDX "\n", index, data_gmsh.epart_in[index]);
        }
        puts("\nnode partitions:");
        for (int index = 0; index < data_gmsh.nn; ++index)
        {
            printf("data_gmsh.npart_in[%d] = %" PRIDX "\n", index, data_gmsh.npart_in[index]);
        }
#endif // partition information

        // building coarse level
        DataGmsh coarse_data_gmsh;
        GmshCoarseLevelGenerator(&coarse_data_gmsh, &data_gmsh);

#if 1
        puts("\n==== coarse level mesh information ====");
        //printf("number of nodes: %d\n", coarse_data_gmsh.nn);
        puts("$Nodes");
        for (int index = 0; index < coarse_data_gmsh.nn; ++index)
        {
            printf("%d\t%021.16le\t%021.16le\t%021.16le\n", index,
                   coarse_data_gmsh.coordinates[3 * index],
                   coarse_data_gmsh.coordinates[3 * index + 1],
                   coarse_data_gmsh.coordinates[3 * index + 2]);
        }
        /*
        printf("number of elements: %d\n", coarse_data_gmsh.ne_in);
        for (int index = 0; index < coarse_data_gmsh.ne_in + 1; ++index)
        {
            printf("coarse_eptr[%d] = %ld\n", index, coarse_data_gmsh.eptr_in[index]);
        }
        */
       puts("$Adjacency");
        for (int index = 0; index < coarse_data_gmsh.nn; ++index)
        {
            idx_t index_start = coarse_data_gmsh.eptr_in[index];
            idx_t index_end = coarse_data_gmsh.eptr_in[index + 1];
            //printf("node %d: ", index);
            printf("%d\t", index);
            for (idx_t index_i = index_start; index_i < index_end; ++index_i)
            {
                printf("%" PRIDX "\t", coarse_data_gmsh.eind_in[index_i]);
            }
            putchar('\n');
        }
#endif

        TestMetisFunctionGraph(coarse_data_gmsh);
        puts("\n==== graph partition of coarse level ====");
        for (int index = 0; index < coarse_data_gmsh.nn; ++index)
        {
            printf("coarse_data_gmsh.npart_in[%d] = %" PRIDX "\n", index,
                   coarse_data_gmsh.npart_in[index]);
        }

        // free memory
        free(data_gmsh.coordinates);
        free(data_gmsh.eptr_bd);
        free(data_gmsh.eind_bd);
        free(data_gmsh.eptr_in);
        free(data_gmsh.eind_in);
        free(data_gmsh.epart_in);
        free(data_gmsh.npart_in);

        free(coarse_data_gmsh.coordinates);
        free(coarse_data_gmsh.eptr_in);
        free(coarse_data_gmsh.eind_in);
        free(coarse_data_gmsh.npart_in);
    }
    else if (strcmp(type, "default") == 0)
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
#endif // metis api, serial implementation
