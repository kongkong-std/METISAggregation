#include "../include/main.h"

int TestMetisFunctionGraph(DataGmsh data /*gmsh data*/)
{
    int status = 0;

    idx_t nvtxs = data.nn, ncon = 1;
    idx_t nparts = data.nparts;
    idx_t edgecut;
    idx_t options[METIS_NOPTIONS];
    METIS_SetDefaultOptions(options);

    status = METIS_PartGraphKway(&nvtxs, &ncon,
                                 data.eptr_in, data.eind_in,
                                 NULL, NULL, NULL, &nparts,
                                 NULL, NULL, options,
                                 &edgecut, data.npart_in);

    return status;
}

void GmshCoarseLevelGenerator(DataGmsh *coarse_data /*gmsh coarse level data pointer*/,
                              DataGmsh *fine_data /*gmsh fine level data pointer*/)
{
    coarse_data->nn = fine_data->nparts;
    coarse_data->ne_in = coarse_data->nn; // adjacency list size, number of elements equal to number of nodes
    coarse_data->nparts = coarse_data->nn / 4;
    coarse_data->nne_bd = NumNodeEleTypeMap(1); // segment line element
    coarse_data->nne_in = NumNodeEleTypeMap(1); // segment line element

    // coarse level node coordinates
    coarse_data->coordinates = (double *)malloc(coarse_data->nn * 3 * sizeof(double));
    assert(coarse_data->coordinates);
    memset(coarse_data->coordinates, 0, 3 * coarse_data->nn * sizeof(double));

    int *cnt_node_partition = NULL; // nodes in current partition
    cnt_node_partition = (int *)malloc(coarse_data->nn * sizeof(int));
    assert(cnt_node_partition);
    memset(cnt_node_partition, 0, coarse_data->nn * sizeof(int));

    for (int index = 0; index < fine_data->nn; ++index)
    {
        idx_t id_part = fine_data->npart_in[index]; // 0-base
        coarse_data->coordinates[id_part * 3] += fine_data->coordinates[index * 3];
        coarse_data->coordinates[id_part * 3 + 1] += fine_data->coordinates[index * 3 + 1];
        coarse_data->coordinates[id_part * 3 + 2] += fine_data->coordinates[index * 3 + 2];
        ++(cnt_node_partition[id_part]);
    }

    for (int index = 0; index < coarse_data->nn; ++index)
    {
        coarse_data->coordinates[3 * index] /= cnt_node_partition[index];
        coarse_data->coordinates[3 * index + 1] /= cnt_node_partition[index];
        coarse_data->coordinates[3 * index + 2] /= cnt_node_partition[index];
    }

    // coarse level adjacency list
    coarse_data->eptr_in = (idx_t *)malloc((coarse_data->ne_in + 1) * sizeof(idx_t));
    assert(coarse_data->eptr_in);
    memset(coarse_data->eptr_in, 0, (coarse_data->ne_in + 1) * sizeof(idx_t));

    // adjacency matrix
    bool **mat_adj = NULL;
    mat_adj = (bool **)malloc(coarse_data->nn * sizeof(bool *));
    assert(mat_adj);
    for (int index = 0; index < coarse_data->nn; ++index)
    {
        mat_adj[index] = (bool *)calloc(coarse_data->nn, sizeof(bool));
        assert(mat_adj[index]);
        // mat_adj[index][index] = true;
    }

    // elements in fine level traversal
    for (int index = 0; index < fine_data->ne_in; ++index)
    {
        idx_t index_start = fine_data->eptr_in[index];
        idx_t index_end = fine_data->eptr_in[index + 1];

        for (idx_t index_i = index_start; index_i < index_end; ++index_i)
        {
            idx_t fine_node_i = fine_data->eind_in[index_i];
            idx_t coarse_node_i = fine_data->npart_in[fine_node_i];

            for (idx_t index_j = index_i + 1; index_j < index_end; ++index_j)
            {
                idx_t fine_node_j = fine_data->eind_in[index_j];
                idx_t coarse_node_j = fine_data->npart_in[fine_node_j];

                mat_adj[coarse_node_i][coarse_node_j] = true;
                mat_adj[coarse_node_j][coarse_node_i] = true;
            }
        }
    }

    // diagonal = 0
    for (int index = 0; index < coarse_data->nn; ++index)
    {
        mat_adj[index][index] = 0;
    }

    // assign value to coarse_data->eptr_in
    for (int index_i = 0; index_i < coarse_data->nn; ++index_i)
    {
        int cnt_tmp = 0;
        for (int index_j = 0; index_j < coarse_data->nn; ++index_j)
        {
            if (mat_adj[index_i][index_j])
            {
                ++cnt_tmp;
            }
        }
        coarse_data->eptr_in[index_i + 1] = coarse_data->eptr_in[index_i] + cnt_tmp;
    }

    // assign value to coarse_data->eind_in
    coarse_data->eind_in = (idx_t *)malloc(coarse_data->eptr_in[coarse_data->nn] * sizeof(idx_t));
    assert(coarse_data);
    int pos_tmp = 0;
    for (int index_i = 0; index_i < coarse_data->nn; ++index_i)
    {
        for (int index_j = 0; index_j < coarse_data->nn; ++index_j)
        {
            if (mat_adj[index_i][index_j])
            {
                coarse_data->eind_in[pos_tmp] = index_j;
                ++pos_tmp;
            }
        }
    }

    coarse_data->npart_in = (idx_t *)malloc(coarse_data->nn * sizeof(idx_t));
    assert(coarse_data->npart_in);

    // free memory
    free(cnt_node_partition);
    for (int index = 0; index < coarse_data->nn; ++index)
    {
        free(mat_adj[index]);
    }
    free(mat_adj);
}

int TestMetisFunctionGmsh(DataGmsh data)
{
    int status = 0;

    idx_t ne = data.ne_in, nn = data.nn;
    idx_t nparts = data.nparts;
    idx_t options[METIS_NOPTIONS];
    idx_t objval;
    // idx_t *epart, *npart;

#if 0
    epart = (idx_t *)malloc(ne * sizeof(idx_t));
    npart = (idx_t *)malloc(nn * sizeof(idx_t));
    assert(epart && npart);
#endif

    METIS_SetDefaultOptions(options);

#if 0
    for (int index = 0; index < ne + 1; ++index)
    {
        printf("data.eptr_in[%d] = %ld\n", index, data.eptr_in[index]);
    }
    for (int index = 0; index < data.eptr_in[ne]; ++index)
    {
        printf("data.eind_in[%d] = %ld\n", index, data.eind_in[index]);
    }
#endif // mesh csr format information

    status = METIS_PartMeshNodal(&ne, &nn,
                                 data.eptr_in, data.eind_in,
                                 NULL, NULL,
                                 &nparts,
                                 NULL, options,
                                 &objval, data.epart_in, data.npart_in);

#if 0
    puts("\n==== metis function result ====");
    printf("status = %d, METIS_OK = %d\n", status, METIS_OK);
    printf("partition objective value = %ld\n", objval);
    printf("element partition:\n");
    for (int index = 0; index < ne; ++index)
    {
        printf("epart[%d] = %ld\n", index, data.epart_in[index]);
    }
    printf("node partition:\n");
    for (int index = 0; index < nn; ++index)
    {
        // printf("npart[%d] = %ld\n", index, npart[index]);
        printf("%ld\n", data.npart_in[index]);
    }
#endif // metis function result

#if 0
    // free memory
    free(epart);
    free(npart);
#endif

    return status;
}

int NumNodeEleTypeMap(int ele_type /*element type*/)
{
    int value = 0;

    switch (ele_type)
    {
    case 1:
        // 2-node line
        value = 2;
        break;

    case 2:
        // 3-node triangle
        value = 3;
        break;

    case 3:
        // 4-node quadrangle
        value = 4;
        break;

    case 4:
        // 4-node tetrahedron
        value = 4;
        break;

    case 5:
        // 8-node hexahedron
        value = 8;
        break;

    case 6:
        // 6-node prism
        value = 6;
        break;

    case 7:
        // 5-node pyramid
        value = 5;
        break;

    default:
        break;
    }

    return value;
}

void FileProcessGmsh(const char *path /*path to gmsh file*/,
                     DataGmsh *data /*gmsh data pointer*/)
{
    FILE *fp = fopen(path, "rb");
    assert(fp);

    char buffer[BUF_MAX_SIZE];
    Flag_Data_Block current_status = NONE;
    int line_node = 0, line_ele_bd = 0, line_ele_in = 0;
    int nn_ele_bd = 0, nn_ele_in = 0;

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
                sscanf(buffer, "%d", &(data->nn));
                data->coordinates = (double *)malloc(data->nn * 3 * sizeof(double));
                assert(data->coordinates);
            }
            else
            {
                sscanf(buffer, "%*d%lf%lf%lf", data->coordinates + 3 * line_node,
                       data->coordinates + 3 * line_node + 1,
                       data->coordinates + 3 * line_node + 2);
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
            char *tmp_token = strtok(tmp_buffer, " \t");
            while (tmp_token != NULL)
            {
                ++tmp_cnt;
                tmp_token = strtok(NULL, " \t");
            }

            if (tmp_cnt == 1)
            {
                sscanf(buffer, "%d", &(data->ne));
                data->eptr_bd = (idx_t *)malloc(data->ne * sizeof(idx_t));
                data->eptr_in = (idx_t *)malloc(data->ne * sizeof(idx_t));
                nn_ele_bd = NumNodeEleTypeMap(1); // segment element
                nn_ele_in = NumNodeEleTypeMap(2); // triangle element
                data->eind_bd = (idx_t *)malloc(data->ne * nn_ele_bd * sizeof(idx_t));
                data->eind_in = (idx_t *)malloc(data->ne * nn_ele_in * sizeof(idx_t));
                assert(data->eptr_bd &&
                       data->eptr_in &&
                       data->eind_bd &&
                       data->eind_in);
                data->eptr_bd[0] = 0;
                data->eptr_in[0] = 0;
            }
            else
            {
                int ele_type = 0;
                sscanf(buffer, "%*d%d", &ele_type);
                if (ele_type == 1)
                {
                    data->eptr_bd[line_ele_bd + 1] = data->eptr_bd[line_ele_bd] + nn_ele_bd;

                    // last nn_ele_bd number in buffer
                    char *token;
                    char *last2[2] = {NULL, NULL}; // 保存最后两个token

                    token = strtok(buffer, " \t");
                    while (token != NULL)
                    {
                        // 滚动保存最后两个
                        last2[0] = last2[1];
                        last2[1] = token;
                        token = strtok(NULL, " \t");
                    }

                    int second_last = atoi(last2[0]);
                    int last = atoi(last2[1]);

                    // change to 0-base
                    data->eind_bd[2 * line_ele_bd] = second_last - 1;
                    data->eind_bd[2 * line_ele_bd + 1] = last - 1;

                    ++line_ele_bd;
                }
                else if (ele_type == 2)
                {
                    data->eptr_in[line_ele_in + 1] = data->eptr_in[line_ele_in] + nn_ele_in;

                    // last nn_ele_bd number in buffer
                    char *token;
                    char *last3[3] = {NULL, NULL, NULL}; // 保存最后两个token

                    token = strtok(buffer, " \t");
                    while (token != NULL)
                    {
                        // 滚动保存最后两个
                        last3[0] = last3[1];
                        last3[1] = last3[2];
                        last3[2] = token;
                        token = strtok(NULL, " \t");
                    }

                    int third_last = atoi(last3[0]);
                    int second_last = atoi(last3[1]);
                    int last = atoi(last3[2]);

                    // change to 0-base
                    data->eind_in[3 * line_ele_in] = third_last - 1;
                    data->eind_in[3 * line_ele_in + 1] = second_last - 1;
                    data->eind_in[3 * line_ele_in + 2] = last - 1;

                    ++line_ele_in;
                }
            }
            break;
        }
        default:
        {
            break;
        }
        }
    }

    data->ne_bd = line_ele_bd;
    data->ne_in = line_ele_in;
    data->nne_bd = nn_ele_bd;
    data->nne_in = nn_ele_in;
    data->nparts = data->nn / 4;

    data->epart_in = (idx_t *)malloc(data->ne_in * sizeof(idx_t));
    data->npart_in = (idx_t *)malloc(data->nn * sizeof(idx_t));
    assert(data->epart_in && data->npart_in);

    fclose(fp);
}
