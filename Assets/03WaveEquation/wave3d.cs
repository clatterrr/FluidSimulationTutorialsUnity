using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class wave3d : MonoBehaviour
{
    // Start is called before the first frame update
    //第一步 : Create Empty
    //第二步 ： 加上这个脚本
    //第三步 : 加一个 MeshRender 并添加材质
    //第四步 ： 按下播放键观看动画
    private const int Dimension = 100;
    Vector3[] verts = new Vector3[Dimension * Dimension];//上一帧的高度
    Vector3[] TimeMius2Verts = new Vector3[Dimension * Dimension];//上两帧的高度
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
                TimeMius2Verts[j * Dimension + i] = verts[j * Dimension + i] = new Vector3(i / 10.0f, 0, j / 10.0f);
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
    // Update is called once per frame
    void Update()
    {
        float LeftHeight = 0, RightHeight = 0, TopHeight = 0, BottomHeight = 0, CenterHeight = 0, TimeMius2CenterHeight = 0;
        NowFrame += 1;
        Vector3[] NewVerts = new Vector3[Dimension * Dimension];
        for (int i = 0; i < Dimension; i++)
        {
            for (int j = 0; j < Dimension; j++)
            {
                int VerIdx = j * Dimension + i;

                if (i > 0) LeftHeight = verts[VerIdx - 1].y;
                if (i < Dimension - 1) RightHeight = verts[VerIdx + 1].y;
                if (j > 0) BottomHeight = verts[VerIdx - Dimension].y;
                if (j < Dimension - 1) TopHeight = verts[VerIdx + Dimension].y;
                CenterHeight = verts[VerIdx].y;

                float Fsource = 0.0f;
<<<<<<< HEAD
                if (i == Dimension / 2 && j == Dimension / 2)   Fsource = Mathf.Sin(Time.time * 2.0f)/2;
                TimeMius2CenterHeight = TimeMius2Verts[VerIdx].y;
                float _Speed = 0.4f;
=======
                if (i == Dimension / 2 && j == Dimension / 2)   Fsource = Mathf.Sin(Time.time * 2.0f) * 2;
                TimeMius2CenterHeight = TimeMius2Verts[VerIdx].y;
                float _Speed = 0.6f;
>>>>>>> 33c11bdfd6ce5b4669b5c37773f6f57b8d7bceb2
                float _Height = 2 * CenterHeight - TimeMius2CenterHeight + (LeftHeight + RightHeight + BottomHeight + TopHeight - 4 * CenterHeight) * _Speed * _Speed + Fsource;
                NewVerts[j * Dimension + i] = new Vector3(i / 10.0f, _Height, j / 10.0f);
            }
        }

        mesh.vertices = NewVerts;
        TimeMius2Verts = verts;
        verts = NewVerts;
        mesh.RecalculateNormals();
    }
}
