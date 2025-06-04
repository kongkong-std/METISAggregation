#include "main.h"

#if 0
int main(int argc, char **argv)
{
    int my_rank, nprocs;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    char *path_mesh = NULL;
    for (int index = 0; index < argc; ++index)
    {
        if (strstr("-mesh", argv[index]))
        {
            path_mesh = argv[index + 1];
        }
    }

#if 0
    printf("rank %d/%d: "
           "mesh file: %s\n",
           my_rank, nprocs, path_mesh);
#endif // test rank/nprocs

    DataMesh mesh_data;
    int local_ne = 0;
    DataMeshEle *local_ele = NULL;

    idx_t *elmdist = NULL;
    elmdist = (idx_t *)calloc(nprocs + 1, sizeof(idx_t));
    assert(elmdist);

    MPI_Comm comm = MPI_COMM_WORLD;

    if (my_rank == 0)
    {
        // mesh file IO
        FileProcessMesh(path_mesh, &mesh_data);
#if 0
        printf("number of nodes: %d, number of elements: %d\n", mesh_data.nn, mesh_data.ne);
        for (int index = 0; index < mesh_data.ne; ++index)
        {
            printf("ele %d: \t%d\t%d\t%d\t%d\n", mesh_data.ele[index].ele_idx,
                   mesh_data.ele[index].ele_node[0],
                   mesh_data.ele[index].ele_node[1],
                   mesh_data.ele[index].ele_node[2],
                   mesh_data.ele[index].ele_node[3]);
        }
        puts("\n========\n");
#endif // mesh file information

        int base_num = mesh_data.ne / nprocs;
        int remainder_num = mesh_data.ne % nprocs;
        int offset = 0;
        for (int index_p = 0; index_p < nprocs; ++index_p)
        {
            int count = base_num + (index_p < remainder_num ? 1 : 0);

            if (index_p == 0)
            {
                local_ne = count;
                local_ele = (DataMeshEle *)malloc(local_ne * sizeof(DataMeshEle));
                assert(local_ele);

                memcpy(local_ele, mesh_data.ele + offset, local_ne * sizeof(DataMeshEle));
            }
            else
            {
                // send mesh_data
                MPI_Send(&count, 1, MPI_INT, index_p, 0, comm);
                MPI_Send(mesh_data.ele + offset, count * sizeof(DataMeshEle), MPI_BYTE, index_p, 0, comm);
            }

            offset += count;
            elmdist[index_p + 1] = elmdist[index_p] + count;
        }
    }
    else
    {
        // allocate data to other processors
        MPI_Recv(&local_ne, 1, MPI_INT, 0, 0, comm, MPI_STATUS_IGNORE);

        local_ele = (DataMeshEle *)malloc(local_ne * sizeof(DataMeshEle));
        assert(local_ele);

        MPI_Recv(local_ele, local_ne * sizeof(DataMeshEle), MPI_BYTE, 0, 0, comm, MPI_STATUS_IGNORE);
    }

    MPI_Bcast(elmdist, (nprocs + 1) * sizeof(idx_t), MPI_BYTE, 0, comm);

    // calling api
    /*
     * int __cdecl ParMETIS_V3_Mesh2Dual(
     *         idx_t *elmdist, idx_t *eptr, idx_t *eind, idx_t *numflag,
     *  idx_t *ncommonnodes, idx_t **xadj, idx_t **adjncy, MPI_Comm *comm);
     */
    idx_t *eptr = NULL, *eind = NULL;
    idx_t numflag = 0;
    idx_t ncommonnodes = 2; // quadrilateral
    // idx_t *xadj = NULL, *adjncy = NULL;

    eptr = (idx_t *)calloc(local_ne + 1, sizeof(idx_t));
    eind = (idx_t *)calloc(4 * local_ne, sizeof(idx_t));
    assert(eptr && eind);

    for (int index = 0; index < local_ne; ++index)
    {
        eptr[index + 1] = eptr[index] + 4;
    }

    for (int index = 0; index < local_ne; ++index)
    {
        for (int index_i = 0; index_i < 4; ++index_i)
        {
            eind[4 * index + index_i] = local_ele[index].ele_node[index_i];
        }
    }

    // int metis_status = ParMETIS_V3_Mesh2Dual(elmdist, eptr, eind, &numflag, &ncommonnodes, &xadj, &adjncy, &comm);
    // printf("metis_status = %d\n", metis_status);

    idx_t options[METIS_NOPTIONS];
    METIS_SetDefaultOptions(options);

    // parmetis api
    /*
     * int __cdecl ParMETIS_V3_PartMeshKway(
     *         idx_t *elmdist, idx_t *eptr, idx_t *eind, idx_t *elmwgt,
     *     idx_t *wgtflag, idx_t *numflag, idx_t *ncon, idx_t *ncommonnodes, idx_t *nparts,
     *     real_t *tpwgts, real_t *ubvec, idx_t *options, idx_t *edgecut, idx_t *part,
     *     MPI_Comm *comm);
     */
    idx_t ncon = 1;
    idx_t nparts = 4;
    idx_t edgecut;
    real_t ubvec = 1.05;
    idx_t wgtflag = 0;

    idx_t *part = NULL;
    real_t *tpwgts = NULL;
    part = (idx_t *)malloc(eptr[local_ne] * sizeof(idx_t));
    tpwgts = (real_t *)malloc(ncon * nparts * sizeof(real_t));
    assert(part && tpwgts);

    for (int index = 0; index < ncon * nparts; ++index)
    {
        tpwgts[index] = 1. / nparts;
    }

    int metis_status = ParMETIS_V3_PartMeshKway(elmdist, eptr, eind,
                                                NULL,
                                                &wgtflag, &numflag, &ncon, &ncommonnodes, &nparts,
                                                tpwgts, &ubvec, options,
                                                &edgecut, part, &comm);

#if 1
    for (int index_p = 0; index_p < nprocs; ++index_p)
    {
        MPI_Barrier(comm);
        if (my_rank == index_p)
        {
            printf("rank %d: local number of elements %d\n", my_rank, local_ne);
            for (int index = 0; index < local_ne; ++index)
            {
                printf("rank %d: ele[%d] \t%d\t%d\t%d\t%d\n", my_rank,
                       local_ele[index].ele_idx,
                       local_ele[index].ele_node[0],
                       local_ele[index].ele_node[1],
                       local_ele[index].ele_node[2],
                       local_ele[index].ele_node[3]);
            }
            putchar('\n');

            for (int index = 0; index < nprocs + 1; ++index)
            {
                printf("elmdist[%d] = %" PRIDX "\t", index, elmdist[index]);
            }
            puts("\n");

            for (int index = 0; index < local_ne + 1; ++index)
            {
                printf("eptr[%d] = %" PRIDX "\t", index, eptr[index]);
            }
            puts("\n");

            for (int index = 0; index < local_ne; ++index)
            {
                for (int index_i = 0; index_i < 4; ++index_i)
                {
                    printf("eind[%d] = %" PRIDX "\t", 4 * index + index_i, eind[4 * index + index_i]);
                }
                putchar('\n');
            }
            printf("\nmetis_status = %d\n", metis_status);

            puts("\npartition result:");
            for(int index = 0; index < eptr[local_ne]; ++index)
            {
                printf("node %" PRIDX " belongs to partition %" PRIDX "\n", eind[index], part[index]);
            }

            puts("\n========\n");
        }
    }
#endif // print rank information of mesh

    // free memory
    free(part);
    free(tpwgts);
    free(eptr);
    free(eind);
    free(elmdist);
    free(local_ele);
    if (my_rank == 0)
    {
        free(mesh_data.ele);
    }

    MPI_Finalize();
    return 0;
}
#endif // test api

#if 0
int main(int argc, char **argv)
{
    int my_rank, nprocs;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    printf("rank %d/%d\n", my_rank, nprocs);

    // parmetis api
    /*
    * int __cdecl ParMETIS_V3_PartKway(
             idx_t *vtxdist, idx_t *xadj, idx_t *adjncy, idx_t *vwgt,
         idx_t *adjwgt, idx_t *wgtflag, idx_t *numflag, idx_t *ncon, idx_t *nparts,
         real_t *tpwgts, real_t *ubvec, idx_t *options, idx_t *edgecut, idx_t *part,
         MPI_Comm *comm);
    */
    idx_t vtxdist[4] = {0, 5, 10, 15};
    idx_t *xadj = NULL, *adjncy = NULL;
    idx_t wgtflag = 0;
    idx_t numflag = 0, ncon = 1, nparts = 4;

    real_t *tpwgts = NULL;
    real_t ubvec = 1.05;

    idx_t options[METIS_NOPTIONS];
    METIS_SetDefaultOptions(options);

    idx_t edgecut = 0;
    idx_t *part = NULL;

    MPI_Comm comm = MPI_COMM_WORLD;

    if (my_rank == 0)
    {
        idx_t xadj_array[6] = {0, 2, 5, 8, 11, 13};
        idx_t adjncy_array[13] = {1, 5, 0, 2, 6, 1, 3, 7, 2, 4, 8, 3, 9};

        xadj = xadj_array;
        adjncy = adjncy_array;

        tpwgts = (real_t *)malloc(ncon * nparts * sizeof(real_t));
        assert(tpwgts);
        for (int index = 0; index < ncon * nparts; ++index)
        {
            tpwgts[index] = 1. / ncon / nparts;
        }

        part = (idx_t *)malloc((vtxdist[my_rank + 1] - vtxdist[my_rank]) * sizeof(idx_t));
    }

    if (my_rank == 1)
    {
        idx_t xadj_array[6] = {0, 3, 7, 11, 15, 18};
        idx_t adjncy_array[18] = {0, 6, 10, 1, 5, 7, 11, 2, 6, 8, 12, 3, 7, 9, 13, 4, 8, 14};

        xadj = xadj_array;
        adjncy = adjncy_array;

        tpwgts = (real_t *)malloc(ncon * nparts * sizeof(real_t));
        assert(tpwgts);
        for (int index = 0; index < ncon * nparts; ++index)
        {
            tpwgts[index] = 1. / ncon / nparts;
        }

        part = (idx_t *)malloc((vtxdist[my_rank + 1] - vtxdist[my_rank]) * sizeof(idx_t));
    }

    if (my_rank == 2)
    {
        idx_t xadj_array[6] = {0, 2, 5, 8, 11, 13};
        idx_t adjncy_array[13] = {5, 11, 6, 10, 12, 7, 11, 13, 8, 12, 14, 9, 13};

        xadj = xadj_array;
        adjncy = adjncy_array;

        tpwgts = (real_t *)malloc(ncon * nparts * sizeof(real_t));
        assert(tpwgts);
        for (int index = 0; index < ncon * nparts; ++index)
        {
            tpwgts[index] = 1. / ncon / nparts;
        }

        part = (idx_t *)malloc((vtxdist[my_rank + 1] - vtxdist[my_rank]) * sizeof(idx_t));
    }

    int metis_status = ParMETIS_V3_PartKway(vtxdist, xadj, adjncy,
                                            NULL, NULL, &wgtflag,
                                            &numflag, &ncon, &nparts,
                                            tpwgts, &ubvec, options,
                                            &edgecut, part, &comm);

#if 1
    for (int index_p = 0; index_p < nprocs; ++index_p)
    {
        MPI_Barrier(comm);
        if (my_rank == index_p)
        {
            printf("====in rank %d:\n", my_rank);
            printf("metis_status = %d\n", metis_status);

            for (int index = 0; index < vtxdist[my_rank + 1] - vtxdist[my_rank]; ++index)
            {
                printf("global node %" PRIDX ", part[%d] = %" PRIDX "\n", vtxdist[my_rank] + index, index, part[index]);
            }
            puts("\n--------\n");
        }
    }
#endif // print partition result

    // free
    free(part);
    free(tpwgts);

    MPI_Finalize();
    return 0;
}
#endif // test-mesh api ParMETIS_V3_PartKway

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
    idx_t num_local_node = 0;
    idx_t *vtxdist = NULL;

    vtxdist = (idx_t *)calloc(nprocs + 1, sizeof(idx_t));

    if (my_rank == 0)
    {
        printf("path to mesh file: %s\n", path_mesh);

        FileProcessMesh(path_mesh, &mesh_data);
#if 0
        for (int index = 0; index < mesh_data.ne; ++index)
        {
            printf("ele %d:\t%d\t%d\t%d\n", mesh_data.ele[index].ele_idx,
                   mesh_data.ele[index].ele_node[0],
                   mesh_data.ele[index].ele_node[1],
                   mesh_data.ele[index].ele_node[2]);
        }
        puts("\n========\n");
#endif // print mesh data
        int base_num = mesh_data.nn / nprocs;
        int remainder_num = mesh_data.nn % nprocs;
        for (int index_p = 0; index_p < nprocs; ++index_p)
        {
            int count = base_num + (index_p < remainder_num ? 1 : 0);
            vtxdist[index_p + 1] = vtxdist[index_p] + count;
        }
    }

    MPI_Bcast(vtxdist, (nprocs + 1) * sizeof(idx_t), MPI_BYTE, 0, comm);
    num_local_node = vtxdist[my_rank + 1] - vtxdist[my_rank];

    MPI_Bcast(&mesh_data.nn, sizeof(int), MPI_BYTE, 0, comm);
    MPI_Bcast(&mesh_data.ne, sizeof(int), MPI_BYTE, 0, comm);

    if (my_rank != 0)
    {
        mesh_data.ele = (DataMeshEle *)malloc(mesh_data.ne * sizeof(DataMeshEle));
        assert(mesh_data.ele);
    }

    MPI_Bcast(mesh_data.ele, mesh_data.ne * sizeof(DataMeshEle), MPI_BYTE, 0, comm);

    idx_t *xadj = NULL, *adjncy = NULL;

    CSRAdjGenerator(mesh_data.ele, mesh_data.ne, mesh_data.nn, vtxdist, my_rank, &xadj, &adjncy);

    // parmetis api
    /*
    * int __cdecl ParMETIS_V3_PartKway(
             idx_t *vtxdist, idx_t *xadj, idx_t *adjncy, idx_t *vwgt,
         idx_t *adjwgt, idx_t *wgtflag, idx_t *numflag, idx_t *ncon, idx_t *nparts,
         real_t *tpwgts, real_t *ubvec, idx_t *options, idx_t *edgecut, idx_t *part,
         MPI_Comm *comm);
    */
    idx_t wgtflag = 0;
    idx_t numflag = 0, ncon = 1, nparts = 4;
    real_t *tpwgts = NULL, ubvec = 1.05;

    idx_t options[METIS_NOPTIONS];
    METIS_SetDefaultOptions(options);

    idx_t edgecut = 0;
    idx_t *part = (idx_t *)malloc(num_local_node * sizeof(idx_t));
    tpwgts = (real_t *)malloc(ncon * nparts * sizeof(real_t));
    assert(part && tpwgts);

    for (int index = 0; index < ncon * nparts; ++index)
    {
        tpwgts[index] = 1. / ncon / nparts;
    }

    int metis_status = ParMETIS_V3_PartKway(vtxdist, xadj, adjncy,
                                            NULL, NULL, &wgtflag,
                                            &numflag, &ncon, &nparts,
                                            tpwgts, &ubvec, options,
                                            &edgecut, part, &comm);

    for (int index_p = 0; index_p < nprocs; ++index_p)
    {
        MPI_Barrier(comm);
        if (my_rank == index_p)
        {
            printf("global nn = %d, global ne = %d\n", mesh_data.nn, mesh_data.ne);
            printf("in rank %d, mesh information\n", my_rank);
#if 0
            for (int index = 0; index < mesh_data.ne; ++index)
            {
                printf("ele %d:\t%d\t%d\t%d\n", mesh_data.ele[index].ele_idx,
                       mesh_data.ele[index].ele_node[0],
                       mesh_data.ele[index].ele_node[1],
                       mesh_data.ele[index].ele_node[2]);
            }
            // puts("\n========\n");
#endif // print mesh data
            printf("number of local nodes is %" PRIDX "\n", num_local_node);
            printf("vtxdist value:\t");
            for (int index = 0; index < nprocs + 1; ++index)
            {
                printf("%" PRIDX "\t", vtxdist[index]);
            }
            puts("\n");

            printf("value of xadj:\t");
            for (int index = 0; index < num_local_node + 1; ++index)
            {
                printf("%" PRIDX "\t", xadj[index]);
            }
            puts("\n");

            printf("value of adjncy:\t");
            for (int index = 0; index < xadj[num_local_node]; ++index)
            {
                printf("%" PRIDX "\t", adjncy[index]);
            }
            puts("\n");

            printf("metis_status = %d\n", metis_status);
            for (int index = 0; index < num_local_node; ++index)
            {
                printf("global node %" PRIDX ", part[%d] = %" PRIDX "\n", vtxdist[my_rank] + index, index, part[index]);
            }

            puts("\n----\n");
        }
    }

    // free memory
    free(part);
    free(tpwgts);
    free(adjncy);
    free(xadj);
    free(vtxdist);
    free(mesh_data.ele);

    MPI_Finalize();
    return 0;
}
