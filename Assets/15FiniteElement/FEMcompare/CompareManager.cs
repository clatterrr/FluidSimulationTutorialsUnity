using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class CompareManager : MonoBehaviour
{

    public GameObject nodeprefab;
    GameObject[] nodes;
    Vector3[] Forces;
    Vector3[] Position;
    Vector3[] Velocity;
    Vector3[] Acceleration;
    int[,] Straint1;
    int[,] Straint2;
    int[,] Straint3;


    public float Stiffness1 = 50.0f;
    public float Stiffness2 = 0.0f;
    public float Stiffness3 = 50.0f;

    public float damp = 0.1f;
    public enum TestCases
    {
        TopCornerFixed,
        LeftFixed,
        BottomFixed
    }
    public TestCases testcases;
    public int ClothWidth = 8;
    public int ClothHeight = 3;
    public float mass = 1.0f;

    bool[] FixedNodes;
    float restLength = 0.8f;
    float dt = 0.05f;
    Vector3 gap = new Vector3(0.8f, 0.8f, 0.8f);
    float[,] Kmat; // 有限元刚度矩阵
    Vector3[] res;// 雅可比迭代结果
    int total;
    void Start()
    {
        total = ClothHeight * ClothWidth;
        Kmat = new float[total, total];
        res = new Vector3[total];
        nodes = new GameObject[total];
        Forces = new Vector3[total];
        Position = new Vector3[total];
        Velocity = new Vector3[total];
        Acceleration = new Vector3[total];
        FixedNodes = new bool[total];
        Straint1 = new int[total, 4];
        Straint2 = new int[total, 4];
        Straint3 = new int[total, 4];
        for (int i = 0; i < total; i++)
        {
            Velocity[i] = new Vector3(0, 0, 0);
            res[i] = new Vector3(0, 0, 0);
            FixedNodes[i] = false;
            Straint1[i, 0] = Straint1[i, 1] = Straint1[i, 2] = Straint1[i, 3] = -1;
            Straint2[i, 0] = Straint2[i, 1] = Straint2[i, 2] = Straint2[i, 3] = -1;
            Straint3[i, 0] = Straint3[i, 1] = Straint3[i, 2] = Straint3[i, 3] = -1;
        }
        for (int i = 0; i < total; i++)
        {
            if ((i % ClothWidth == 0) && (testcases == TestCases.LeftFixed))FixedNodes[i] = true;
            if ((i < ClothWidth) && (testcases == TestCases.BottomFixed)) FixedNodes[i] = true;
        }
        if(testcases == TestCases.TopCornerFixed)
        {
            FixedNodes[total - 1] = true;
            FixedNodes[total - ClothWidth] = true;
        }
        if(testcases == TestCases.BottomFixed)
        {
            gap.y = 1.2f;
        }
        for (int i = 0; i < ClothWidth; i++)
        {
            for (int j = 0; j < ClothHeight; j++)
            {
                Init(i, j);
            }
        }
    }

    void Init(int i, int j)
    {
        int idx = i + j * ClothWidth;
        nodes[idx] = GameObject.Instantiate(nodeprefab, new Vector3(i * gap.x + 1.5f * ClothWidth - 10 + GetComponent<Transform>().position.x, j * gap.y + 5, 0), Quaternion.identity);
        Position[idx] = nodes[idx].transform.position;
        if (i > 0) Straint1[idx, 0] = idx - 1;
        if (i < ClothWidth - 1) Straint1[idx, 1] = idx + 1;
        if (j > 0) Straint1[idx, 2] = idx - ClothWidth;
        if (j < ClothHeight - 1) Straint1[idx, 3] = idx + ClothWidth;

        if (i > 1) Straint2[idx, 0] = idx - 2;
        if (i < ClothWidth - 2) Straint2[idx, 1] = idx + 2;
        if (j > 1) Straint2[idx, 2] = idx - ClothWidth * 2;
        if (j < ClothHeight - 2) Straint2[idx, 3] = idx + ClothWidth * 2;

        if (i > 0 && j > 0) Straint3[idx, 0] = idx - 1 - ClothWidth;
        if (i < ClothWidth - 1 && j > 0) Straint3[idx, 1] = idx + 1 - ClothWidth;
        if (i > 0 && j < ClothHeight - 1) Straint3[idx, 2] = idx - 1 + ClothWidth;
        if (i < ClothWidth - 1 && j < ClothHeight - 1) Straint3[idx, 3] = idx + 1 + ClothWidth;
    }
    Vector3 FirstStrain(int me, int him)
    {
        Vector3 vec = Position[him] - Position[me];
        float dis = Vector3.Distance(Position[him], Position[me]) - restLength;
        Vector3 nor = Vector3.Normalize(vec);
        return nor * dis * Stiffness1;
    }
    Vector3 SecondStrain(int me, int him)
    {
        Vector3 vec = Position[him] - Position[me];
        float dis = Vector3.Distance(Position[him], Position[me]) - restLength * 2;
        Vector3 nor = Vector3.Normalize(vec);
        return nor * dis * Stiffness2;
    }

    Vector3 ThirdStrain(int me, int him)
    {
        Vector3 vec = Position[him] - Position[me];
        float dis = Vector3.Distance(Position[him], Position[me]) - restLength * Mathf.Sqrt(2);
        Vector3 nor = Vector3.Normalize(vec);
        return nor * dis * Stiffness3;
    }
    void Gauss()
    {
        int tmax = 100;
        while (tmax > 0)
        {
            tmax -= 1;
            for (int i = 0; i < total; i++)
            {
                float sumx = 0;
                float sumy = 0;
                float sumz = 0;
                for (int j = 0; j < total; j++)
                {
                    if (i == j) continue;
                    sumx += Kmat[i, j] * res[j].x;
                    sumy += Kmat[i, j] * res[j].y;
                    sumz += Kmat[i, j] * res[j].z;
                }
                res[i].x = (Forces[i].x / mass - sumx) / Kmat[i, i];
                res[i].y = (Forces[i].y / mass - sumy) / Kmat[i, i];
                res[i].z = (Forces[i].z / mass - sumz) / Kmat[i, i];
            }
        }
    }

    void Jacobi()
    {
        int tmax = 100;
        Vector3[] res2 = new Vector3[total];
        while (tmax > 0)
        {
            for (int i = 0; i < total; i++)
            {
                res2[i] = res[i];
            }
            tmax -= 1;
            for (int i = 0; i < total; i++)
            {
                float sumx = 0;
                float sumy = 0;
                float sumz = 0;
                for (int j = 0; j < total; j++)
                {
                    if (i == j) continue;
                    sumx += Kmat[i, j] * res2[j].x;
                    sumy += Kmat[i, j] * res2[j].y;
                    sumz += Kmat[i, j] * res2[j].z;
                }
                res[i].x = (Forces[i].x / mass - sumx) / Kmat[i, i];
                res[i].y = (Forces[i].y / mass - sumy) / Kmat[i, i];
                res[i].z = (Forces[i].z / mass - sumz) / Kmat[i, i];
            }
        }
    }
    void FiniteElement()
    {
        for (int i = 0; i < total; i++)
        {
            for (int j = 0; j < total; j++)
                Kmat[i, j] = 0;
            Forces[i] = new Vector3(0.0f, -0.98f * mass, 0.0f);
        }
        for (int i = 0; i < total; i++)
        {
            Vector3 oriVel = Velocity[i];
            if (FixedNodes[i] == true)
            {
                Kmat[i, i] = 1;
                Forces[i] = new Vector3(0.0f, 0.0f, 0.0f);
                continue;
            }
            for (int k = 0; k < 4; k++)
            {
                if (Straint1[i, k] != -1)
                {
                    Kmat[i, i] += 1;
                    Kmat[i, Straint1[i, k]] -= 1;
                    Forces[i] += FirstStrain(i, Straint1[i, k]);
                }

                if (Straint2[i, k] != -1)
                {
                    Kmat[i, i] += 1;
                    Kmat[i, Straint2[i, k]] -= 1;
                    Forces[i] += SecondStrain(i, Straint2[i, k]);
                }

                if (Straint3[i, k] != -1)
                {
                    Kmat[i, i] += 1;
                    Kmat[i, Straint3[i, k]] -= 1;
                    Forces[i] += ThirdStrain(i, Straint3[i, k]);
                }

            }
        }
        Jacobi();
        // Gauss();
        for (int i = 0; i < total; i++)
        {
            Velocity[i] = (dt * res[i] + Velocity[i]) * (1 - damp);
            Position[i] = Position[i] + dt * Velocity[i];
            nodes[i].transform.position = Position[i];
        }

    }
    void Update()
    {
        FiniteElement();
    }
}
