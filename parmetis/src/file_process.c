#include "main.h"

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
