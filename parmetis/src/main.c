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
    AdjDataMesh fine_graph_data;

    fine_graph_data.vtxdist = (idx_t *)calloc(nprocs + 1, sizeof(idx_t));
    assert(fine_graph_data.vtxdist);

    // parmetis api
    /*
    int __cdecl ParMETIS_V3_PartKway(
             idx_t *vtxdist, idx_t *xadj, idx_t *adjncy, idx_t *vwgt,
         idx_t *adjwgt, idx_t *wgtflag, idx_t *numflag, idx_t *ncon, idx_t *nparts,
         real_t *tpwgts, real_t *ubvec, idx_t *options, idx_t *edgecut, idx_t *part,
         MPI_Comm *comm);
     */
    idx_t *local_xadj = NULL, *local_adjncy = NULL;
    idx_t wgtflag = 0;
    idx_t numflag = 0, ncon = 1;

    real_t *tpwgts = NULL, ubvec = 1.05;

    idx_t options[METIS_NOPTIONS];
    METIS_SetDefaultOptions(options);

    idx_t local_index_start, local_index_end;

    if (my_rank == 0)
    {
        printf("mesh file path: %s\n", path_mesh);
        FileProcessMesh(path_mesh, &mesh_data);

        int base_num = mesh_data.nn / nprocs;
        int remainder_num = mesh_data.nn % nprocs;

        for (int index_p = 0; index_p < nprocs; ++index_p)
        {
            int count = base_num + (index_p < remainder_num ? 1 : 0);
            fine_graph_data.vtxdist[index_p + 1] = fine_graph_data.vtxdist[index_p] + count;
        }

        fine_graph_data.nn = mesh_data.nn;
        fine_graph_data.dim = mesh_data.dim;
        fine_graph_data.nparts = fine_graph_data.nn / 4;

        fine_graph_data.coordinates = (double *)malloc(fine_graph_data.dim * fine_graph_data.nn * sizeof(double));
        fine_graph_data.part = (idx_t *)malloc(fine_graph_data.nn * sizeof(idx_t));
        assert(fine_graph_data.coordinates && fine_graph_data.part);

        memcpy(fine_graph_data.coordinates, mesh_data.coordinates, mesh_data.dim * mesh_data.nn * sizeof(double));

        GlobalGmshCSRAdjGenerator(mesh_data.ele_shell, mesh_data.ne_shell, mesh_data.nn, &fine_graph_data);

#if 1
        for (int index = 0; index < fine_graph_data.nn; ++index)
        {
            printf("node %d:\t", index);
            for (int index_i = 0; index_i < fine_graph_data.dim; ++index_i)
            {
                printf("%021.16le\t", fine_graph_data.coordinates[fine_graph_data.dim * index + index_i]);
            }
            putchar('\n');
        }

        puts("\n\nfine_graph_data xadj value:");
        for (int index = 0; index < fine_graph_data.nn + 1; ++index)
        {
            printf("%" PRIDX "\t", fine_graph_data.xadj[index]);
        }
        putchar('\n');

        puts("\n\nfine_graph_data adjncy value:");
        for (int index = 0; index < fine_graph_data.xadj[fine_graph_data.nn]; ++index)
        {
            printf("%" PRIDX "\t", fine_graph_data.adjncy[index]);
        }
        putchar('\n');
#endif // fine_graph_data coordinates

        local_index_start = fine_graph_data.xadj[fine_graph_data.vtxdist[my_rank]];
        local_index_end = fine_graph_data.xadj[fine_graph_data.vtxdist[my_rank + 1]];

        for (int index_p = 1; index_p < nprocs; ++index_p)
        {
            MPI_Send(fine_graph_data.xadj + fine_graph_data.vtxdist[index_p], sizeof(idx_t), MPI_BYTE, index_p, 0, comm);     // index_start
            MPI_Send(fine_graph_data.xadj + fine_graph_data.vtxdist[index_p + 1], sizeof(idx_t), MPI_BYTE, index_p, 1, comm); // index_end
        }
    }
    else
    {
        MPI_Recv(&local_index_start, sizeof(idx_t), MPI_BYTE, 0, 0, comm, MPI_STATUS_IGNORE);
        MPI_Recv(&local_index_end, sizeof(idx_t), MPI_BYTE, 0, 1, comm, MPI_STATUS_IGNORE);
    }

    MPI_Bcast(&fine_graph_data.nn, 1, MPI_INT, 0, comm);
    MPI_Bcast(fine_graph_data.vtxdist, (nprocs + 1) * sizeof(idx_t), MPI_BYTE, 0, comm);
    // MPI_Bcast(fine_graph_data.xadj + fine_graph_data.vtxdist[my_rank], sizeof(idx_t), MPI_BYTE, 0, comm);     // index_start
    // MPI_Bcast(fine_graph_data.xadj + fine_graph_data.vtxdist[my_rank + 1], sizeof(idx_t), MPI_BYTE, 0, comm); // index_end

    idx_t nparts = fine_graph_data.nn / 4;
    idx_t edgecut = 0;
    idx_t num_local_node = fine_graph_data.vtxdist[my_rank + 1] - fine_graph_data.vtxdist[my_rank];
    // idx_t local_index_start = fine_graph_data.xadj[fine_graph_data.vtxdist[my_rank]];
    // idx_t local_index_end = fine_graph_data.xadj[fine_graph_data.vtxdist[my_rank + 1]];
    idx_t nnz_local_node = local_index_end - local_index_start;
    idx_t *local_part = (idx_t *)malloc(num_local_node * sizeof(idx_t));
    tpwgts = (real_t *)malloc(ncon * nparts * sizeof(real_t));
    assert(local_part && tpwgts);

    for (int index = 0; index < ncon * nparts; ++index)
    {
        tpwgts[index] = 1. / ncon / nparts;
    }

    local_xadj = (idx_t *)malloc((num_local_node + 1) * sizeof(idx_t));
    local_adjncy = (idx_t *)malloc(nnz_local_node * sizeof(idx_t));
    assert(local_xadj && local_adjncy);

    if (my_rank == 0)
    {
        memcpy(local_xadj, fine_graph_data.xadj, (num_local_node + 1) * sizeof(idx_t));
        memcpy(local_adjncy, fine_graph_data.adjncy, nnz_local_node * sizeof(idx_t));

        for (int index_p = 1; index_p < nprocs; ++index_p)
        {
            idx_t tmp_num_local_node = fine_graph_data.vtxdist[index_p + 1] - fine_graph_data.vtxdist[index_p];
            idx_t tmp_local_index_start = fine_graph_data.xadj[fine_graph_data.vtxdist[index_p]];
            idx_t tmp_local_index_end = fine_graph_data.xadj[fine_graph_data.vtxdist[index_p + 1]];
            idx_t tmp_nnz_local_node = tmp_local_index_end - tmp_local_index_start;

            MPI_Send(fine_graph_data.xadj + fine_graph_data.vtxdist[index_p], (tmp_num_local_node + 1) * sizeof(idx_t), MPI_BYTE,
                     index_p, 1, comm); // xadj data
            MPI_Send(fine_graph_data.adjncy + tmp_local_index_start, tmp_nnz_local_node * sizeof(idx_t), MPI_BYTE,
                     index_p, 2, comm); // adjncy data
        }
    }
    else
    {
        MPI_Recv(local_xadj, (num_local_node + 1) * sizeof(idx_t), MPI_BYTE, 0, 1, comm, MPI_STATUS_IGNORE);
        MPI_Recv(local_adjncy, nnz_local_node * sizeof(idx_t), MPI_BYTE, 0, 2, comm, MPI_STATUS_IGNORE);
    }

    // local xadj start from 0
    idx_t tmp_shift = local_xadj[0];
    for (int i = 0; i <= num_local_node; ++i)
    {
        local_xadj[i] -= tmp_shift;
    }

    // partition
    int metis_status = ParMETIS_V3_PartKway(fine_graph_data.vtxdist,
                                            local_xadj, local_adjncy,
                                            NULL, NULL, &wgtflag, &numflag, &ncon, &nparts,
                                            tpwgts, &ubvec, options,
                                            &edgecut, local_part, &comm);

    // gather to root rank
    int *recvcounts = NULL, *displs = NULL;
    if (my_rank == 0)
    {
        recvcounts = (int *)malloc(nprocs * sizeof(int));
        displs = (int *)malloc(nprocs * sizeof(int));
        assert(recvcounts && displs);

        for (int index = 0; index < nprocs; ++index)
        {
#if 1
            recvcounts[index] = fine_graph_data.vtxdist[index + 1] - fine_graph_data.vtxdist[index];
            displs[index] = fine_graph_data.vtxdist[index];
#endif // int
#if 0
            recvcounts[index] = (fine_graph_data.vtxdist[index + 1] - fine_graph_data.vtxdist[index]) * sizeof(idx_t);
            displs[index] = (fine_graph_data.vtxdist[index]) * sizeof(idx_t);
#endif // byte
        }
    }

#if 1
    MPI_Gatherv(local_part, num_local_node, MPI_INT,
                fine_graph_data.part, recvcounts, displs, MPI_INT, 0, comm);
#endif // int
#if 0
    MPI_Gatherv(local_part, num_local_node, MPI_BYTE,
                fine_graph_data.part, recvcounts, displs, MPI_BYTE, 0, comm);
#endif // byte

#if 1
    if (my_rank == 0)
    {
        puts("\nglobal partition value:");
        for (int index = 0; index < fine_graph_data.nn; ++index)
        {
            printf("part[%d] = %" PRIDX "\n", index, fine_graph_data.part[index]);
        }

        puts("\n========\n");
    }
#endif // print global partition value

    // coarse level
    AdjDataMesh coarse_graph_data;
    coarse_graph_data.vtxdist = (idx_t *)calloc(nprocs + 1, sizeof(idx_t));
    assert(coarse_graph_data.vtxdist);

    if (my_rank == 0)
    {
        CoarseLevelGenerator(&fine_graph_data,
                             &coarse_graph_data);

#if 1
        printf("\n>>>>coarse level has %d nodes:\n", coarse_graph_data.nn);
        for (int index = 0; index < coarse_graph_data.nn; ++index)
        {
            printf("coarse node %d:\t", index);
            for (int index_i = 0; index_i < coarse_graph_data.dim; ++index_i)
            {
                printf("%021.16le\t", coarse_graph_data.coordinates[coarse_graph_data.dim * index + index_i]);
            }
            putchar('\n');
        }

        puts("\n>>>>coarse level xadj:");
        for (int index = 0; index < coarse_graph_data.nn + 1; ++index)
        {
            printf("%" PRIDX "\t", coarse_graph_data.xadj[index]);
        }
        putchar('\n');

        puts("\n>>>>coarse level adjncy:");
        for (int index = 0; index < coarse_graph_data.xadj[coarse_graph_data.nn]; ++index)
        {
            printf("%" PRIDX "\t", coarse_graph_data.adjncy[index]);
        }
        putchar('\n');
#endif // print coarse_graph_data information

        int base_num = coarse_graph_data.nn / nprocs;
        int remainder_num = coarse_graph_data.nn % nprocs;

        for (int index_p = 0; index_p < nprocs; ++index_p)
        {
            int count = base_num + (index_p < remainder_num ? 1 : 0);
            coarse_graph_data.vtxdist[index_p + 1] = coarse_graph_data.vtxdist[index_p] + count;
        }

        local_index_start = coarse_graph_data.xadj[coarse_graph_data.vtxdist[my_rank]];
        local_index_end = coarse_graph_data.xadj[coarse_graph_data.vtxdist[my_rank + 1]];
        for (int index_p = 1; index_p < nprocs; ++index_p)
        {
            MPI_Send(coarse_graph_data.xadj + coarse_graph_data.vtxdist[index_p],
                     sizeof(idx_t), MPI_BYTE, index_p, 10, comm); // coarse level local index_start
            MPI_Send(coarse_graph_data.xadj + coarse_graph_data.vtxdist[index_p + 1],
                     sizeof(idx_t), MPI_BYTE, index_p, 11, comm); // coarse level local index_end
        }
    }
    else
    {
        MPI_Recv(&local_index_start, sizeof(idx_t), MPI_BYTE,
                 0, 10, comm, MPI_STATUS_IGNORE);
        MPI_Recv(&local_index_end, sizeof(idx_t), MPI_BYTE,
                 0, 11, comm, MPI_STATUS_IGNORE);
    }

    MPI_Bcast(&coarse_graph_data.nn, 1, MPI_INT, 0, comm);
    MPI_Bcast(coarse_graph_data.vtxdist, (nprocs + 1) * sizeof(idx_t), MPI_BYTE, 0, comm);

    nparts = coarse_graph_data.nn / 4;
    edgecut = 0;
    num_local_node = coarse_graph_data.vtxdist[my_rank + 1] - coarse_graph_data.vtxdist[my_rank];
    nnz_local_node = local_index_end - local_index_start;
    ncon = 1;
    numflag = 0;
    wgtflag = 0;
    idx_t *coarse_local_part = (idx_t *)malloc(num_local_node * sizeof(idx_t));
    real_t *coarse_tpwgts = (real_t *)malloc(ncon * nparts * sizeof(real_t));
    ubvec = 1.05;
    METIS_SetDefaultOptions(options);
    assert(coarse_local_part && coarse_tpwgts);

    for (int index = 0; index < ncon * nparts; ++index)
    {
        coarse_tpwgts[index] = 1. / ncon / nparts;
    }

    idx_t *coarse_local_xadj = (idx_t *)malloc((num_local_node + 1) * sizeof(idx_t));
    idx_t *coarse_local_adjncy = (idx_t *)malloc(nnz_local_node * sizeof(idx_t));
    assert(coarse_local_xadj && coarse_local_adjncy);

    if (my_rank == 0)
    {
        memcpy(coarse_local_xadj, coarse_graph_data.xadj, (num_local_node + 1) * sizeof(idx_t));
        memcpy(coarse_local_adjncy, coarse_graph_data.adjncy, nnz_local_node * sizeof(idx_t));

        for (int index_p = 1; index_p < nprocs; ++index_p)
        {
            idx_t tmp_num_local_node = coarse_graph_data.vtxdist[index_p + 1] - coarse_graph_data.vtxdist[index_p];
            idx_t tmp_local_index_start = coarse_graph_data.xadj[coarse_graph_data.vtxdist[index_p]];
            idx_t tmp_local_index_end = coarse_graph_data.xadj[coarse_graph_data.vtxdist[index_p + 1]];
            idx_t tmp_nnz_local_node = tmp_local_index_end - tmp_local_index_start;

            MPI_Send(coarse_graph_data.xadj + coarse_graph_data.vtxdist[index_p], (tmp_num_local_node + 1) * sizeof(idx_t), MPI_BYTE,
                     index_p, 1, comm); // xadj data
            MPI_Send(coarse_graph_data.adjncy + tmp_local_index_start, tmp_nnz_local_node * sizeof(idx_t), MPI_BYTE,
                     index_p, 2, comm); // adjncy data
        }
    }
    else
    {
        MPI_Recv(coarse_local_xadj, (num_local_node + 1) * sizeof(idx_t), MPI_BYTE, 0, 1, comm, MPI_STATUS_IGNORE);
        MPI_Recv(coarse_local_adjncy, nnz_local_node * sizeof(idx_t), MPI_BYTE, 0, 2, comm, MPI_STATUS_IGNORE);
    }

    tmp_shift = coarse_local_xadj[0];
    for (int index = 0; index <= num_local_node; ++index)
    {
        coarse_local_xadj[index] -= tmp_shift;
    }

    int coarse_metis_status = ParMETIS_V3_PartKway(coarse_graph_data.vtxdist,
                                                   coarse_local_xadj, coarse_local_adjncy,
                                                   NULL, NULL, &wgtflag, &numflag, &ncon, &nparts,
                                                   coarse_tpwgts, &ubvec, options,
                                                   &edgecut, coarse_local_part, &comm);

    int *coarse_recvcounts = NULL, *coarse_displs = NULL;
    if (my_rank == 0)
    {
        coarse_recvcounts = (int *)malloc(nprocs * sizeof(int));
        coarse_displs = (int *)malloc(nprocs * sizeof(int));
        assert(coarse_displs && coarse_displs);

        for (int index = 0; index < nprocs; ++index)
        {
#if 1
            coarse_recvcounts[index] = coarse_graph_data.vtxdist[index + 1] - coarse_graph_data.vtxdist[index];
            coarse_displs[index] = coarse_graph_data.vtxdist[index];
#endif // int
#if 0
            recvcounts[index] = (fine_graph_data.vtxdist[index + 1] - fine_graph_data.vtxdist[index]) * sizeof(idx_t);
            displs[index] = (fine_graph_data.vtxdist[index]) * sizeof(idx_t);
#endif // byte
        }
    }

#if 1
    MPI_Gatherv(coarse_local_part, num_local_node, MPI_INT,
                coarse_graph_data.part, coarse_recvcounts, coarse_displs, MPI_INT, 0, comm);
#endif // int

#if 1
    if (my_rank == 0)
    {
        puts("\ncoarse global partition value:");
        for (int index = 0; index < coarse_graph_data.nn; ++index)
        {
            printf("part[%d] = %" PRIDX "\n", index, coarse_graph_data.part[index]);
        }

        puts("\n========\n");
    }
#endif // print global partition value, coarse level

#if 0
    for (int index_p = 0; index_p < nprocs; ++index_p)
    {
        MPI_Barrier(comm);
        if (my_rank == index_p)
        {
            printf("in rank %d/%d\n\n", my_rank, nprocs);

            printf("\nfine_graph_data vtxdist value:\t");
            for (int index = 0; index < nprocs + 1; ++index)
            {
                printf("%" PRIDX "\t", fine_graph_data.vtxdist[index]);
            }

            puts("\n\nfine_graph_data local_xadj value:");
            for (int index = 0; index < num_local_node + 1; ++index)
            {
                printf("%" PRIDX "\t", local_xadj[index]);
            }
            putchar('\n');

            puts("\nfine_graph_data local_adjncy value:");
            for (int index = 0; index < local_xadj[num_local_node]; ++index)
            {
                printf("%" PRIDX "\t", local_adjncy[index]);
            }
            putchar('\n');

            printf("\nmetis_status = %d\n", metis_status);
            for (int index = 0; index < num_local_node; ++index)
            {
                printf("global node %" PRIDX ", local_part[%d] = %" PRIDX "\n",
                       fine_graph_data.vtxdist[index_p] + index,
                       index,
                       local_part[index]);
            }

            puts("\n\n========\n\n");
        }
    }
#endif // print information in rank ascending order

    // free memory
    free(coarse_graph_data.vtxdist);

    free(tpwgts);
    free(local_part);
    free(local_xadj);
    free(local_adjncy);
    free(fine_graph_data.vtxdist);
    if (my_rank == 0)
    {
        // coarse graph
        free(coarse_graph_data.part);
        free(coarse_graph_data.adjncy);
        free(coarse_graph_data.xadj);
        free(coarse_graph_data.coordinates);

        // fine graph
        free(recvcounts);
        free(displs);
        free(fine_graph_data.part);
        free(fine_graph_data.adjncy);
        free(fine_graph_data.xadj);
        free(fine_graph_data.coordinates);

        // mesh
        free(mesh_data.coordinates);

        for (int index = 0; index < mesh_data.ne_shell; ++index)
        {
            free(mesh_data.ele_shell[index].ele_node);
        }
        free(mesh_data.ele_shell);
    }

    MPI_Finalize();
    return 0;
}
