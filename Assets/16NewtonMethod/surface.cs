using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class surface : MonoBehaviour
{
    public int surfacenum = 0;
    private const int Dimension = 20;
    Vector3[] verts = new Vector3[Dimension * Dimension];//上一帧的高度
    protected MeshFilter meshFilter;
    protected Mesh mesh;
    int NowFrame = 0;
    public Material Mat;
    // Start is called before the first frame update
    void Start()
    {
        Application.targetFrameRate = 10;
        mesh = new Mesh();
        MeshRenderer meshrender = GetComponent<MeshRenderer>();
        meshrender.material = Mat;
        mesh.name = gameObject.name;

        mesh.vertices = GenerateVerts();
        mesh.triangles = GenerateTries();
        mesh.RecalculateBounds();

        meshFilter = gameObject.AddComponent<MeshFilter>();
        meshFilter.mesh = mesh;
    }
    private Vector3[] GenerateVerts()
    {

        for (int i = 0; i < Dimension; i++)
        {
            for (int j = 0; j < Dimension; j++)
            {
                float x = i;
                float y = j;
                float h = (x + 2 * y) * (5 * x + 6 * y) + x - 432;
                if (surfacenum == 1) h = -2 * (x + 2 * y) * (5 * x + 6 * y) + y + 854;
                verts[j * Dimension + i] = new Vector3(x, h * 0.001f, y);
            }
        }
        return verts;
    }

    private int[] GenerateTries()
    {
        var tries = new int[mesh.vertices.Length * 6];
        for (int i = 0; i < Dimension - 1; i++)
        {
            for (int j = 0; j < Dimension - 1; j++)
            {
                int TriIdx = (j * (Dimension - 1) + i) * 6;
                int VerIdx = j * Dimension + i;
                tries[TriIdx] = VerIdx;
                tries[TriIdx + 1] = VerIdx + Dimension;
                tries[TriIdx + 2] = VerIdx + Dimension + 1;

                /*
                 1 ------  2/4
                 |             |
                 |             |
               0/3 ------- 5
                 
                 */

                tries[TriIdx + 3] = VerIdx;
                tries[TriIdx + 4] = VerIdx + Dimension + 1;
                tries[TriIdx + 5] = VerIdx + 1;
            }
        }
        return tries;
    }
}
