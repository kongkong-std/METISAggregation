#include "../include/main.h"

int CoarseLevelGenerator(const AdjDataMesh *fine_graph_data /*fine level graph data*/,
                         AdjDataMesh *coarse_graph_data /*coarse level graph data*/)
{
    coarse_graph_data->nn = fine_graph_data->nparts;
    coarse_graph_data->dim = fine_graph_data->dim;
    coarse_graph_data->nparts = coarse_graph_data->nn / 4;

    coarse_graph_data->coordinates = (double *)calloc(coarse_graph_data->dim * coarse_graph_data->nn, sizeof(double));
    assert(coarse_graph_data->coordinates);

    int *cnt_node_partition = (int *)calloc(coarse_graph_data->nn, sizeof(int));
    assert(cnt_node_partition);

    for (int index = 0; index < fine_graph_data->nn; ++index)
    {
        idx_t id_part = fine_graph_data->part[index];
        for (int index_i = 0; index_i < coarse_graph_data->dim; ++index_i)
        {
            coarse_graph_data->coordinates[coarse_graph_data->dim * id_part + index_i] += fine_graph_data->coordinates[fine_graph_data->dim * index + index_i];
        }
        ++(cnt_node_partition[id_part]);
    }

    for (int index = 0; index < coarse_graph_data->nn; ++index)
    {
        for (int index_i = 0; index_i < coarse_graph_data->dim; ++index_i)
        {
            coarse_graph_data->coordinates[coarse_graph_data->dim * index + index_i] /= cnt_node_partition[index];
        }
    }

    // adjacency list generator
    coarse_graph_data->xadj = (idx_t *)calloc(coarse_graph_data->nn + 1, sizeof(idx_t));

    bool **mat_adj = NULL; // adjacency matrix, coarse.nn x coarse.nn
    mat_adj = (bool **)malloc(coarse_graph_data->nn * sizeof(bool *));
    assert(mat_adj);

    for (int index = 0; index < coarse_graph_data->nn; ++index)
    {
        mat_adj[index] = (bool *)calloc(coarse_graph_data->nn, sizeof(bool));
        assert(mat_adj[index]);
    }

    for (int index = 0; index < fine_graph_data->nn; ++index)
    {
        idx_t index_start = fine_graph_data->xadj[index];
        idx_t index_end = fine_graph_data->xadj[index + 1];
        for (idx_t index_i = index_start; index_i < index_end; ++index_i)
        {
            idx_t fine_node_i = fine_graph_data->adjncy[index_i];
            idx_t coarse_node_i = fine_graph_data->part[fine_node_i];
            for (idx_t index_j = index_i + 1; index_j < index_end; ++index_j)
            {
                idx_t fine_node_j = fine_graph_data->adjncy[index_j];
                idx_t coarse_node_j = fine_graph_data->part[fine_node_j];

                mat_adj[coarse_node_i][coarse_node_j] = true;
                mat_adj[coarse_node_j][coarse_node_i] = true;
            }
        }
    }

    // set mat_adj diagonal to 0
    for (int index = 0; index < coarse_graph_data->nn; ++index)
    {
        mat_adj[index][index] = false;
    }

    // assign coarse xadj
    for (int index_i = 0; index_i < coarse_graph_data->nn; ++index_i)
    {
        int cnt_tmp = 0;
        for (int index_j = 0; index_j < coarse_graph_data->nn; ++index_j)
        {
            if (mat_adj[index_i][index_j])
            {
                ++cnt_tmp;
            }
        }
        coarse_graph_data->xadj[index_i + 1] = coarse_graph_data->xadj[index_i] + cnt_tmp;
    }

    // assign coarse adjncy
    coarse_graph_data->adjncy = (idx_t *)malloc(coarse_graph_data->xadj[coarse_graph_data->nn] * sizeof(idx_t));
    assert(coarse_graph_data->adjncy);
    int pos_tmp = 0;
    for (int index_i = 0; index_i < coarse_graph_data->nn; ++index_i)
    {
        for (int index_j = 0; index_j < coarse_graph_data->nn; ++index_j)
        {
            if (mat_adj[index_i][index_j])
            {
                coarse_graph_data->adjncy[pos_tmp] = index_j;
                ++pos_tmp;
            }
        }
    }

    coarse_graph_data->part = (idx_t *)malloc(coarse_graph_data->nn * sizeof(idx_t));
    assert(coarse_graph_data->part);

    // free memeory
    for (int index = 0; index < coarse_graph_data->nn; ++index)
    {
        free(mat_adj[index]);
    }
    free(mat_adj);
    free(cnt_node_partition);

    return 0;
}

static int CountTrueMatAdj(const bool *a, int n)
{
    int value = 0;

    for (int index = 0; index < n; ++index)
    {
        if (a[index])
        {
            ++value;
        }
    }

    return value;
}

int GlobalGmshCSRAdjGenerator(const DataMeshEle *ele_data /*mesh information*/,
                              int ne /*number of elements*/,
                              int nn /*number of nodes*/,
                              AdjDataMesh *graph_data /*csr graph data*/)
{
    graph_data->xadj = (idx_t *)calloc(nn + 1, sizeof(idx_t));
    bool *mat_adj_tmp = (bool *)calloc(nn, sizeof(bool));
    int **adjncy_tmp = (int **)malloc(nn * sizeof(int *));
    int *adjncy_tmp_len = (int *)malloc(nn * sizeof(int)); // record length of adjncy_tmp[i]

    assert(graph_data->xadj && mat_adj_tmp && adjncy_tmp && adjncy_tmp_len);

    for (int index = 0; index < nn; ++index)
    {
        for (int index_e = 0; index_e < ne; ++index_e)
        {
            // triangle element
            if (ele_data[index_e].ele_node[0] == index ||
                ele_data[index_e].ele_node[1] == index ||
                ele_data[index_e].ele_node[2] == index)
            {
                mat_adj_tmp[ele_data[index_e].ele_node[0]] = true;
                mat_adj_tmp[ele_data[index_e].ele_node[1]] = true;
                mat_adj_tmp[ele_data[index_e].ele_node[2]] = true;
            }
        }
        mat_adj_tmp[index] = false;

        int cnt_tmp = CountTrueMatAdj(mat_adj_tmp, nn);

        graph_data->xadj[index + 1] = graph_data->xadj[index] + cnt_tmp;

        adjncy_tmp[index] = (int *)malloc(cnt_tmp * sizeof(int));
        adjncy_tmp_len[index] = cnt_tmp;
        int cnt_adjncy_tmp = 0;
        for (int index_mat_adj_tmp = 0; index_mat_adj_tmp < nn; ++index_mat_adj_tmp)
        {
            if (mat_adj_tmp[index_mat_adj_tmp])
            {
                adjncy_tmp[index][cnt_adjncy_tmp] = index_mat_adj_tmp;
                ++cnt_adjncy_tmp;
            }
        }

        memset(mat_adj_tmp, 0, nn * sizeof(bool));
    }

    graph_data->adjncy = (idx_t *)malloc(graph_data->xadj[nn] * sizeof(idx_t));
    assert(graph_data->adjncy);

    int index_adjncy = 0;
    for (int index_i = 0; index_i < nn; ++index_i)
    {
        for (int index_j = 0; index_j < adjncy_tmp_len[index_i]; ++index_j)
        {
            graph_data->adjncy[index_adjncy] = adjncy_tmp[index_i][index_j];
            ++index_adjncy;
        }
    }

    // free memory
    free(adjncy_tmp_len);
    for (int index = 0; index < nn; ++index)
    {
        free(adjncy_tmp[index]);
    }
    free(adjncy_tmp);
    free(mat_adj_tmp);

    return 0;
}

int GmshCSRAdjGenerator(const DataMeshEle *ele_data /*mesh information*/,
                        int ne /*number of elements*/,
                        int nn /*number of nodes*/,
                        const idx_t *vtxdist /*global nodes array*/,
                        int my_rank /*current rank*/,
                        idx_t **xadj /*csr row pointer*/,
                        idx_t **adjncy /*adjacency nodes list*/)
{
    idx_t local_nn = vtxdist[my_rank + 1] - vtxdist[my_rank];

    *xadj = (idx_t *)calloc(local_nn + 1, sizeof(idx_t));
    bool *mat_adj_tmp = (bool *)calloc(nn, sizeof(bool));
    int **adjncy_tmp = (int **)malloc(local_nn * sizeof(int *));
    int *adjncy_tmp_len = (int *)malloc(local_nn * sizeof(int)); // record length of adjncy_tmp[i]

    assert(*xadj && mat_adj_tmp && adjncy_tmp && adjncy_tmp_len);

    for (int index = vtxdist[my_rank]; index < vtxdist[my_rank + 1]; ++index)
    {
        for (int index_e = 0; index_e < ne; ++index_e)
        {
            // triangle element
            if (ele_data[index_e].ele_node[0] == index ||
                ele_data[index_e].ele_node[1] == index ||
                ele_data[index_e].ele_node[2] == index)
            {
                mat_adj_tmp[ele_data[index_e].ele_node[0]] = true;
                mat_adj_tmp[ele_data[index_e].ele_node[1]] = true;
                mat_adj_tmp[ele_data[index_e].ele_node[2]] = true;
            }
        }
        mat_adj_tmp[index] = false;

        int cnt_tmp = CountTrueMatAdj(mat_adj_tmp, nn);
        (*xadj)[index - vtxdist[my_rank] + 1] = (*xadj)[index - vtxdist[my_rank]] + cnt_tmp;

        adjncy_tmp[index - vtxdist[my_rank]] = (int *)malloc(cnt_tmp * sizeof(int));
        adjncy_tmp_len[index - vtxdist[my_rank]] = cnt_tmp;
        int cnt_adjncy_tmp = 0;
        for (int index_mat_adj_tmp = 0; index_mat_adj_tmp < nn; ++index_mat_adj_tmp)
        {
            if (mat_adj_tmp[index_mat_adj_tmp])
            {
                adjncy_tmp[index - vtxdist[my_rank]][cnt_adjncy_tmp] = index_mat_adj_tmp;
                ++cnt_adjncy_tmp;
            }
        }

        memset(mat_adj_tmp, 0, nn * sizeof(bool));
    }

    *adjncy = (idx_t *)malloc((*xadj)[local_nn] * sizeof(idx_t));
    assert(*adjncy);

    int index_adjncy = 0;
    for (int index_i = 0; index_i < local_nn; ++index_i)
    {
        for (int index_j = 0; index_j < adjncy_tmp_len[index_i]; ++index_j)
        {
            (*adjncy)[index_adjncy] = adjncy_tmp[index_i][index_j];
            ++index_adjncy;
        }
    }

    // free memory
    free(adjncy_tmp_len);
    for (int index = 0; index < local_nn; ++index)
    {
        free(adjncy_tmp[index]);
    }
    free(adjncy_tmp);
    free(mat_adj_tmp);

    return 0;
}

int FileProcessMesh(const char *path /*path to mesh file*/, DataMesh *mesh_data /*mesh data*/)
{
    FILE *fp = fopen(path, "rb");
    assert(fp);

    char buffer[BUF_MAX_SIZE];
    FlagDataBlockGmsh current_status = NONE;

    int line_node = 0, line_ele = 0;
    DataMeshEle *ele_tmp = NULL;

    while (fgets(buffer, BUF_MAX_SIZE, fp))
    {
        char tmp_buffer[BUF_MAX_SIZE];
        if (buffer[0] == '$')
        {
            if (strstr(buffer, "$MeshFormat"))
            {
                current_status = MESH_FORMAT;
            }
            else if (strstr(buffer, "$PhysicalNames"))
            {
                current_status = PHYSICAL_NAMES;
            }
            else if (strstr(buffer, "$Nodes"))
            {
                current_status = NODES;
            }
            else if (strstr(buffer, "$Elements"))
            {
                current_status = ELEMENTS;
            }
            else
            {
                current_status = NONE;
            }
            continue;
        }

        switch (current_status)
        {
        case NODES:
        {
            /* code */
            // printf("%s", buffer);
            memcpy(tmp_buffer, buffer, BUF_MAX_SIZE);
            tmp_buffer[strcspn(tmp_buffer, "\n")] = '\0';
            int tmp_cnt = 0;
            char *tmp_token = strtok(tmp_buffer, " \t");
            while (tmp_token != NULL)
            {
                ++tmp_cnt;
                tmp_token = strtok(NULL, " \t");
            }

            if (tmp_cnt == 1)
            {
                sscanf(buffer, "%d", &mesh_data->nn);
                mesh_data->dim = 3;
                mesh_data->coordinates = (double *)malloc(mesh_data->dim * mesh_data->nn * sizeof(double));
                assert(mesh_data->coordinates);
            }
            else
            {
                sscanf(buffer, " %*d %lf %lf %lf ", mesh_data->coordinates + mesh_data->dim * line_node,
                       mesh_data->coordinates + mesh_data->dim * line_node + 1,
                       mesh_data->coordinates + mesh_data->dim * line_node + 2);
                ++line_node;
            }

            break;
        }

        case ELEMENTS:
        {
            // printf("%s", buffer);
            memcpy(tmp_buffer, buffer, BUF_MAX_SIZE);
            tmp_buffer[strcspn(tmp_buffer, "\n")] = '\0';
            int tmp_cnt = 0;
            int tmp_array[256];
            char *tmp_token = strtok(tmp_buffer, " \t");
            while (tmp_token != NULL)
            {
                tmp_array[tmp_cnt] = atoi(tmp_token);
                ++tmp_cnt;
                tmp_token = strtok(NULL, " \t");
            }

            if (tmp_cnt == 1)
            {
                sscanf(buffer, "%d", &mesh_data->ne);
                ele_tmp = (DataMeshEle *)malloc(mesh_data->ne * sizeof(DataMeshEle));
                assert(ele_tmp);
            }
            else
            {
                sscanf(buffer, " %*d %d ", &ele_tmp[line_ele].ele_type);
                ele_tmp[line_ele].num_ele_node = NumNodeEleTypeMap(ele_tmp[line_ele].ele_type);
                ele_tmp[line_ele].ele_node = (int *)malloc(ele_tmp[line_ele].num_ele_node * sizeof(int));
                assert(ele_tmp[line_ele].ele_node);

                for (int index = 0; index < ele_tmp[line_ele].num_ele_node; ++index)
                {
                    ele_tmp[line_ele].ele_node[index] = tmp_array[tmp_cnt - ele_tmp[line_ele].num_ele_node + index] - 1; // 0-base
                }

                ++line_ele;
            }

            break;
        }

        default:
            break;
        }
    }

    fclose(fp);

#if 0
    for (int index = 0; index < mesh_data->ne; ++index)
    {
        printf("element %d:\t%d\t", index, ele_tmp[index].ele_type);
        for (int index_i = 0; index_i < ele_tmp[index].num_ele_node; ++index_i)
        {
            printf("%d\t", ele_tmp[index].ele_node[index_i]);
        }
        putchar('\n');
    }
#endif // element data

    int cnt_shell = 0;
    for (int index = 0; index < mesh_data->ne; ++index)
    {
        if (ele_tmp[index].ele_type == 2)
        {
            // triangle mesh
            ++cnt_shell;
        }
    }

    mesh_data->ne_shell = cnt_shell;
    mesh_data->ele_shell = (DataMeshEle *)malloc(mesh_data->ne_shell * sizeof(DataMeshEle));
    assert(mesh_data->ele_shell);

    cnt_shell = 0;
    for (int index = 0; index < mesh_data->ne; ++index)
    {
        if (ele_tmp[index].ele_type == 2)
        {
            mesh_data->ele_shell[cnt_shell].ele_type = ele_tmp[index].ele_type;
            mesh_data->ele_shell[cnt_shell].num_ele_node = ele_tmp[index].num_ele_node;

            mesh_data->ele_shell[cnt_shell].ele_node = (int *)malloc(mesh_data->ele_shell[cnt_shell].num_ele_node * sizeof(int));
            assert(mesh_data->ele_shell[cnt_shell].ele_node);

            memcpy(mesh_data->ele_shell[cnt_shell].ele_node, ele_tmp[index].ele_node, mesh_data->ele_shell[cnt_shell].num_ele_node * sizeof(int));

            ++cnt_shell;
        }
    }

    // free memory
    for (int index = 0; index < mesh_data->ne; ++index)
    {
        free(ele_tmp[index].ele_node);
    }
    free(ele_tmp);

    return 0;
}

#if 0
void TestMetis(void)
{
    int status = 0;

    idx_t nvtxs = 15, ncon = 1, nparts = 4;
    idx_t xadj[] = {0, 2, 5, 8, 11, 13, 16, 20, 24, 28, 31, 33, 36, 39, 42, 44};
    idx_t adjncy[] = {1, 5, 0, 2, 6, 1, 3, 7, 2, 4, 8, 3, 9, 0, 6, 10, 1, 5, 7, 11, 2, 6, 8, 12, 3, 7, 9, 13, 4, 8, 14, 5, 11, 6, 10, 12, 7, 11, 13, 8, 12, 14, 9, 13};
    idx_t ne = 8, nn = 15;
    idx_t eptr[] = {0, 4, 8, 12, 16, 20, 24, 28, 32};
    idx_t eind[] = {0, 5, 6, 1, 1, 6, 7, 2, 2, 7, 8, 3, 3, 8, 9, 4, 5, 10, 11, 6, 6, 11, 12, 7, 1, 12, 13, 8, 8, 13, 14, 9};
    idx_t options[METIS_NOPTIONS];
    idx_t edgecut, objval;

    puts("\n==== test METIS_PartGraphKway ====");
    METIS_SetDefaultOptions(options);
    idx_t part[15];
    status = METIS_PartGraphKway(&nvtxs, &ncon,
                                 xadj, adjncy,
                                 NULL, NULL, NULL,
                                 &nparts,
                                 NULL, NULL,
                                 options,
                                 &edgecut, part);
    printf("status = %d, METIS_OK = %d\n", status, METIS_OK);
    for(int index = 0; index < 15; ++index)
    {
        printf("part[%d] = %" PRIDX "\n", index, part[index]);
    }

    puts("\n==== test METIS_PartMeshNodal ====");
    status = METIS_SetDefaultOptions(options);
    idx_t epart[8], npart[15];
    METIS_PartMeshNodal(&ne, &nn,
                        eptr, eind,
                        NULL, NULL,
                        &nparts,
                        NULL,
                        options,
                        &objval, epart, npart);
    printf("status = %d, METIS_OK = %d\n", status, METIS_OK);
    for(int index = 0; index < 8; ++index)
    {
        printf("epart[%d] = %" PRIDX "\n", index, epart[index]);
    }
    for(int index = 0; index < 15; ++index)
    {
        printf("npart[%d] = %" PRIDX "\n", index, npart[index]);
    }
}

int TestFunctionMetis(DataMesh data)
{
    int status = 0;

    idx_t ne = data.ne, nn = data.nn;
    idx_t *eptr, *eind;
    idx_t nparts;
    idx_t options[METIS_NOPTIONS];
    idx_t objval, *epart, *npart;

    //ncommon = 2;
    nparts = 4;

    eptr = (idx_t *)malloc((ne + 1) * sizeof(idx_t));
    assert(eptr);

    eptr[0] = 0;
    for (int index = 0; index < ne; ++index)
    {
        eptr[index + 1] = eptr[index] + data.nn_ele[index];
    }
#if 0
    for (int index = 0; index < ne + 1; ++index)
    {
        printf("eptr[%d] = %d\n", index, eptr[index]);
    }
#endif

    eind = (idx_t *)malloc(eptr[ne] * sizeof(idx_t));
    assert(eind);
    int cnt_eind = 0;
    for (int index = 0; index < ne; ++index)
    {
        for (int index_i = 0; index_i < data.nn_ele[index]; ++index_i)
        {
            eind[cnt_eind] = data.idx_node[index][index_i];
            ++cnt_eind;
        }
    }
#if 0
    for (int index = 0; index < eptr[ne]; ++index)
    {
        printf("eind[%d] = %d\n", index, eind[index]);
    }
#endif

    epart = (idx_t *)malloc(ne * sizeof(idx_t));
    npart = (idx_t *)malloc(nn * sizeof(idx_t));
    assert(epart && npart);

    METIS_SetDefaultOptions(options);

#if 0
    status = METIS_PartMeshDual(&ne, &nn,
                                eptr, eind,
                                NULL, NULL,
                                &ncommon, &nparts,
                                NULL, options,
                                objval, epart, npart);
#endif

    status = METIS_PartMeshNodal(&ne, &nn,
                                 eptr, eind,
                                 NULL, NULL,
                                 &nparts,
                                 NULL, options,
                                 &objval, epart, npart);

    printf("status = %d, METIS_OK = %d\n", status, METIS_OK);
    printf("partition objective value = %" PRIDX "\n", objval);
    printf("element partition:\n");
    for (int index = 0; index < ne; ++index)
    {
        printf("epart[%d] = %" PRIDX "\n", index, epart[index]);
    }
    printf("node partition:\n");
    for (int index = 0; index < nn; ++index)
    {
        printf("npart[%d] = %" PRIDX "\n", index, npart[index]);
    }

    // free memory
    free(eptr);
    free(eind);
    free(epart);
    free(npart);

    return status;
}

void FileProcessMesh(const char *path, DataMesh *data)
{
    FILE *fp = fopen(path, "rb");
    assert(fp);

    int nn = 0, ne = 0;
    char buffer[BUF_MAX_SIZE];

    fgets(buffer, BUF_MAX_SIZE, fp);
    sscanf(buffer, " %d %d ", &nn, &ne);
    data->nn = nn;
    data->ne = ne;

    data->idx_ele = (int *)malloc((ne + 1) * sizeof(int));
    data->nn_ele = (int *)malloc(ne * sizeof(int));
    data->idx_node = (int **)malloc(ne * sizeof(int *));
    assert(data->idx_ele && data->nn_ele);

    int cnt = 0;
    while (fgets(buffer, BUF_MAX_SIZE, fp))
    {
        char tmp_buffer[BUF_MAX_SIZE];
        memcpy(tmp_buffer, buffer, BUF_MAX_SIZE);
        tmp_buffer[strcspn(tmp_buffer, "\n")] = '\0';
        int count = 0;
        char *token = strtok(tmp_buffer, " \t");
        while (token != NULL)
        {
            ++count;
            token = strtok(NULL, " \t");
        }

        data->nn_ele[cnt] = count - 1;

        data->idx_node[cnt] = (int *)malloc((count - 1) * sizeof(int));
        assert(data->idx_node[cnt]);

        // First read the element index
        char *p = buffer;
        sscanf(p, "%d", &data->idx_ele[cnt]);

        // Skip the element index
        while (*p != ' ' && *p != '\t' && *p != '\0')
            p++;
        while (*p == ' ' || *p == '\t')
            p++;

        // Read the node indices
        for (int i = 0; i < data->nn_ele[cnt]; i++)
        {
            sscanf(p, "%d", &data->idx_node[cnt][i]);

            // Move to next number
            while (*p != ' ' && *p != '\t' && *p != '\0')
                p++;
            while (*p == ' ' || *p == '\t')
                p++;
        }

        ++cnt;
    }

    fclose(fp);
}
#endif // metis test mesh, serial implementation
