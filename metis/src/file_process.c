#include "../include/main.h"

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
