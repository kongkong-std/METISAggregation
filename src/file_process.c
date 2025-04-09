#include "../include/main.h"

int TestFunctionMetis(DataMesh data)
{
    int status = 0;

    idx_t ne = data.ne, nn = data.nn;
    idx_t *eptr, *eind;
    idx_t ncommon, nparts;
    idx_t options[METIS_NOPTIONS];
    idx_t *objval, *epart, *npart;

    ncommon = 2;
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

    status = METIS_PartMeshDual(&ne, &nn,
    eptr, eind,
    NULL, NULL,
    &ncommon, &nparts,
    NULL, options,
    objval, epart, npart);

    // free memory
    free(eptr);
    free(eind);

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
