#include "main.h"

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

int CSRAdjGenerator(const DataMeshEle *ele_data /*mesh topology data*/, int ne /*number of elements*/, int nn /*number of nodes*/,
                    const idx_t *vtxdist /*node list*/,
                    int my_rank /*rank id*/,
                    idx_t **xadj /*csr row pointer*/, idx_t **adjncy /*adjacency nodes*/)
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

int FileProcessMesh(const char *path /*path to mesh file*/,
                    DataMesh *data /*mesh data*/)
{
    FILE *fp = fopen(path, "rb");
    assert(fp);

    int nn = 0, ne = 0;
    fscanf(fp, " %d %d ", &nn, &ne);
    data->nn = nn;
    data->ne = ne;

    data->ele = (DataMeshEle *)malloc(ne * sizeof(DataMeshEle));
    assert(data->ele);

    for (int index = 0; index < ne; ++index)
    {
        fscanf(fp, " %d %d %d %d ", &data->ele[index].ele_idx,
               data->ele[index].ele_node,
               data->ele[index].ele_node + 1,
               data->ele[index].ele_node + 2);
    }

    fclose(fp);

    return 0;
}

#if 0
int FileProcessMesh(const char *path /*path to mesh file*/,
                    DataMesh *data /*mesh data*/)
{
    FILE *fp = fopen(path, "rb");
    assert(fp);

    char buffer[MAX_BUF_SIZE];

    int nn = 0, ne = 0;
    fgets(buffer, MAX_BUF_SIZE, fp);
    sscanf(buffer, " %d %d ", &nn, &ne);
    data->nn = nn;
    data->ne = ne;

    data->ele = (DataMeshEle *)malloc(ne * sizeof(DataMeshEle));
    assert(data->ele);

    for (int index = 0; index < ne; ++index)
    {
        fgets(buffer, MAX_BUF_SIZE, fp);
        sscanf(buffer, " %d %d %d %d %d ", &((data->ele + index)->ele_idx),
               (data->ele + index)->ele_node,
               (data->ele + index)->ele_node + 1,
               (data->ele + index)->ele_node + 2,
               (data->ele + index)->ele_node + 3);
    }

    fclose(fp);

    return 0;
}
#endif // 4 nodes in per element
