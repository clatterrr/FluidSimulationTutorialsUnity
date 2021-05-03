using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class jellyManager : MonoBehaviour
{
    // 用于画三角形的变量
    protected MeshFilter meshFilter;
    protected Mesh mesh;
    const int nmax = 10000;
    int NowFrame = 0;
    public GameObject sphereprefab;
    public Material Mat;
    List<Vector3> point = new List<Vector3>();
    int[] PointFront = new int[nmax];//这个点与几条前线相连
    Vector3[] verties = new Vector3[3];//三角形顶点
    int[] tries = new int[3];//三角形索引
    // 用于链式前向星的变量
    int[] EdgeFrom = new int[nmax]; // 边的起点
    int[] EdgeTo = new int[nmax]; // 边的终点
    int[] EdgeNext = new int[nmax]; // 相同起点的上一条边
    int[] Head = new int[nmax]; // 方括号里填起点编号，可得到起点最新的边
    int Edgenum = 0;// 总边数

    // 波前推进法的一些变量
    bool[] EdgeFront = new bool[nmax]; // 当前边是否为前边(front)
    int EdgeIndex = 0;// 目前索引数
    List<int> Triangles = new List<int>();
    float TriangleRadius = 1.0f;//等边三角形的边长，(1-0.5^2)^(1/2) = 0.866
    float height = 0.86f;//从边的中点生出新点的距离，尽量保持为等边三角形

    //有限元法的一些变量
    bool[] FixedNodes = new bool[nmax];
    int FixedNodesNum = 32;
    int[,] Constrain;
    float[,] RestLength;
    float[,] Kmat;
    Vector3[] Force;
    int pointnum = -1;
    float dt = 0.1f;
    float mass = 0.1f;
    float damp = 0.1f;
    GameObject[] nodes;
    Vector3[] Position;
    Vector3[] Velocity;
    Vector3[] res;

    float GroundStiffness = 1.0f;
    float JellyStiffness = 1.0f;
    public enum TestCases
    {
        square1,
        square2,
        square3,
        square4,
        circle1,
        circle2,
        circle3,
        circle4
    }
    public TestCases testcase;

    void Start()
    {

        for (int i = 0; i < nmax; i++)
        {
            Head[i] = -1;
            EdgeFront[i] = false;
            PointFront[i] = 0;
            FixedNodes[i] = false;
        }
        switch(testcase)
        {
            case TestCases.square1: Init1();GroundStiffness = 1.0f;JellyStiffness = 8.0f; break;
            case TestCases.square2: Init1(); GroundStiffness = 2.0f; JellyStiffness = 20.0f; break;
            case TestCases.square3: Init1(); GroundStiffness = 0.1f; JellyStiffness = 1.0f; break;
            case TestCases.square4: Init1(); GroundStiffness = 0.01f; JellyStiffness = 0.01f; break;
            case TestCases.circle1: Init2(); GroundStiffness = 10.0f; JellyStiffness = 10.0f; break;
            case TestCases.circle2: Init2(); GroundStiffness = 2.0f; JellyStiffness = 20.0f; break;
            case TestCases.circle3: Init2(); GroundStiffness = 0.1f; JellyStiffness = 1.0f; break;
            case TestCases.circle4: Init2(); GroundStiffness = 0.01f; JellyStiffness = 0.2f; break;
        }
        for (int i = 0; i < point.Count; i++)
        {
            PointFront[i] = 2;
            // Instantiate(sphereprefab, point[i], Quaternion.identity);
        }
        mesh = new Mesh();
        MeshRenderer meshrender = GetComponent<MeshRenderer>();
        meshrender.material = Mat;
        mesh.name = gameObject.name;
        meshFilter = gameObject.AddComponent<MeshFilter>();
        /*
        verties[0] = new Vector3(0, 1, 0);
        verties[1] = new Vector3(-1, 0, 0);
        verties[2] = new Vector3(1, 0, 0);
        tries[0] = 0;
        tries[1] = 1;
        tries[2] = 2;
        mesh.vertices = verties;
        mesh.triangles = tries;*/
        meshFilter.mesh = mesh;


    }
    void Init1()
    {
        TriangleRadius = 0.6f;
        height = TriangleRadius * 0.9f;
        int squarelen = 10; // 平板边缘的顶点数
        float distan = TriangleRadius; // 顶点之间的距离
        for (int i = 0; i < squarelen; i++)
        {
            point.Add(new Vector3(i * distan, 0.0f, 0));
            AddEdge(i, i + 1);
        }
        for (int i = 0; i < squarelen; i++)
        {
            point.Add(new Vector3(squarelen * distan, 0.0f, i * distan));
            AddEdge(i + squarelen, i + 1 + squarelen);
        }
        for (int i = 0; i < squarelen; i++)
        {
            point.Add(new Vector3((squarelen - i) * distan, 0.0f, squarelen * distan));
            FixedNodes[point.Count - 1] = true;
            AddEdge(i + squarelen * 2, i + 1 + squarelen * 2);
        }
        for (int i = 0; i < squarelen; i++)
        {
            point.Add(new Vector3(0, 0.0f, (squarelen - i) * distan));
            int ed = i + 1 + squarelen * 3;
            if (i == squarelen - 1) ed = 0;
            AddEdge(i + squarelen * 3, ed);
        }
        // 绘制中间那个洞
        Vector2 CircleCenter = new Vector2(3.0f, 3.0f);
        float CircleRadius = 1.0f;
        int CircleParticle = 12;
        for (int i = 0; i < CircleParticle; i++)
        {
            point.Add(new Vector3(CircleCenter.x + CircleRadius * Mathf.Cos(2 * Mathf.PI / CircleParticle * i), 0.0f, CircleCenter.y - CircleRadius * Mathf.Sin(2 * Mathf.PI / CircleParticle * i)));
            int ed = i + 1 + squarelen * 4;
            if (i == CircleParticle - 1) ed = 4 * squarelen;
            AddEdge(i + squarelen * 4, ed);
        }

    }
    void Init2()
    {
        TriangleRadius = 0.4f;
        height = TriangleRadius * 0.9f;
        Vector2 CircleCenter = new Vector2(4.0f, 4.0f);
        float CircleRadius = 3.0f;
        int CircleParticle = 32;
        for (int i = 0; i < CircleParticle; i++)
        {
            point.Add(new Vector3(CircleCenter.x + CircleRadius * Mathf.Cos(2 * Mathf.PI / CircleParticle * i), 0.0f, CircleCenter.y + CircleRadius * Mathf.Sin(2 * Mathf.PI / CircleParticle * i)));
            int ed = i + 1;
            if (i == CircleParticle - 1) ed = 0;
            AddEdge(i, ed);
        }

        Vector2 CircleCenter2 = new Vector2(3.0f, 3.0f);
        float CircleRadius2 = 1.0f;
        int CircleParticle2 = 12;
        for (int i = 0; i < CircleParticle2; i++)
        {
            point.Add(new Vector3(CircleCenter2.x + CircleRadius2 * Mathf.Cos(2 * Mathf.PI / CircleParticle2 * i), 0.0f, CircleCenter2.y - CircleRadius2 * Mathf.Sin(2 * Mathf.PI / CircleParticle2 * i)));
            int ed = i + 1 + CircleParticle;
            if (i == CircleParticle2 - 1) ed = CircleParticle;
            AddEdge(i + CircleParticle, ed);
        }
    }
    void AddEdge(int st, int ed)
    {
        EdgeFrom[Edgenum] = st;
        EdgeTo[Edgenum] = ed;
        EdgeNext[Edgenum] = Head[st];
        Head[st] = Edgenum;
        EdgeFront[Edgenum] = true;
        Edgenum++;
        // Debug.Log("new Edge " + st + " to " +  ed);
    }
    // Update is called once per frame
    void GenNewTriangle()
    {

        int st = EdgeFrom[EdgeIndex];
        int ed = EdgeTo[EdgeIndex];
        for (int i = Head[ed]; i != -1; i = EdgeNext[i])
        {
            if (EdgeFront[i] == false) continue;
            int to = EdgeTo[i];
            for (int j = Head[to]; j != -1; j = EdgeNext[j])
            {
                if (EdgeFront[j] == false) continue;
                if (EdgeTo[j] == st)
                {
                    Triangles.Add(st);
                    Triangles.Add(to);
                    Triangles.Add(ed);
                    PointFront[st] -= 2;
                    PointFront[to] -= 2;
                    PointFront[ed] -= 2;
                    EdgeFront[i] = false;
                    EdgeFront[j] = false;
                    EdgeFront[EdgeIndex] = false;
                    EdgeIndex++;
                    return;
                }
            }
        }

        Vector3 StartVec = point[st];
        Vector3 EndVec = point[ed];
        Vector3 MidVec = (StartVec + EndVec) / 2;
        Vector3 Normal = new Vector3(-(EndVec.z - StartVec.z), 0.0f, EndVec.x - StartVec.x);
        float div = Mathf.Sqrt(Normal.x * Normal.x + Normal.z * Normal.z);
        Normal.x /= div;
        Normal.z /= div;
        Vector3 NewVec = MidVec + Normal * height;//新点的位置
        int newindex = -1;
        float mindis = 9999;
        bool leftedge = true, rightedge = true;
        for (int i = 0; i < point.Count; i++)
        {
            if (i == st || i == ed) continue;
            // 那么之后我们只要禁止有新的三角形顶点为这些PointFront值为零的顶点就行了
            if (PointFront[i] == 0) continue;
            float dis = Vector3.Distance(NewVec, point[i]);
            // 选择离新形成顶点最近的已有顶点，防止形成过于细长的三角形
            if (dis < TriangleRadius && dis < mindis)
            {
                mindis = dis;
                newindex = i;
            }
        }
        if (newindex == -1)
        {
            newindex = point.Count;
            point.Add(NewVec);
            //  Instantiate(sphereprefab, NewVec, Quaternion.identity);
            PointFront[newindex] += 2;
        }
        else
        {
            for (int i = Head[newindex]; i != -1; i = EdgeNext[i])
            {
                if (EdgeTo[i] == st)
                {
                    EdgeFront[i] = false;
                    leftedge = false;//无需再新建三角形的左边
                    PointFront[st] -= 2;
                    break;
                }
            }
            if (leftedge == true)
                for (int i = Head[ed]; i != -1; i = EdgeNext[i])
                {
                    if (EdgeTo[i] == newindex)
                    {
                        EdgeFront[i] = false;
                        rightedge = false;//无需再新建三角形的右边
                        PointFront[ed] -= 2;
                        break;
                    }
                }
            if (leftedge == true && rightedge == true) PointFront[newindex] += 2;
        }
        EdgeFront[EdgeIndex] = false;
        Triangles.Add(st);
        Triangles.Add(newindex);
        Triangles.Add(ed);
        if (leftedge) AddEdge(st, newindex);
        if (rightedge) AddEdge(newindex, ed);
        EdgeIndex++;

    }
    void Jacobi()
    {
        int tmax = 100;
        Vector3[] res2 = new Vector3[pointnum];
        while (tmax > 0)
        {
            for (int i = 0; i < pointnum; i++)
            {
                res2[i] = res[i];
            }
            tmax -= 1;
            for (int i = 0; i < pointnum; i++)
            {
                float sumx = 0;
                float sumy = 0;
                float sumz = 0;
                for (int j = 0; j < pointnum; j++)
                {
                    if (i == j) continue;
                    sumx += Kmat[i, j] * res2[j].x;
                    sumy += Kmat[i, j] * res2[j].y;
                    sumz += Kmat[i, j] * res2[j].z;
                }
                res[i].x = (Force[i].x / mass - sumx) / Kmat[i, i];
                res[i].y = (Force[i].y / mass - sumy) / Kmat[i, i];
                res[i].z = (Force[i].z / mass - sumz) / Kmat[i, i];
            }
        }
    }
    void FemInit()
    {

        pointnum = point.Count + FixedNodesNum;
        Kmat = new float[pointnum, pointnum];
        Force = new Vector3[pointnum];
        Constrain = new int[20, pointnum];
        RestLength = new float[20, pointnum];
        Position = new Vector3[pointnum];
        Velocity = new Vector3[pointnum];
        nodes = new GameObject[pointnum];
        res = new Vector3[pointnum];
        for (int i = 0; i < pointnum; i++)
        {
            if (i < point.Count) Position[i] = point[i];
            else Position[i] = new Vector3((i - point.Count) * 0.3f - 0.4f, 0.0f, -2.0f);
        }
        for (int i = 0; i < pointnum; i++)
        {
            Constrain[0, i] = 0;
            for (int j = 1; j < 20; j++)
            {
                Constrain[j, i] = -1;
            }
            nodes[i] = Instantiate(sphereprefab, Position[i], Quaternion.identity);
        }
        for (int i = 0; i < Edgenum; i++)
        {
            int p1 = EdgeFrom[i];
            int p2 = EdgeTo[i];
            Constrain[Constrain[0, p1] + 1, p1] = p2;
            Constrain[Constrain[0, p2] + 1, p2] = p1;
            float dis = Vector3.Distance(Position[p1], Position[p2]);
            RestLength[Constrain[0, p1] + 1, p1] = dis;
            RestLength[Constrain[0, p2] + 1, p2] = dis;
            Constrain[0, p1] += 1;
            Constrain[0, p2] += 1;
        }
    }
    void Fem()
    {
        for (int i = 0; i < pointnum; i++)
        {
            for (int j = 0; j < pointnum; j++)
            {
                Kmat[i, j] = 0.0f;
            }
            Force[i] = new Vector3(0.0f, 0.0f, 0.0f);
        }
        for (int i = 0; i < pointnum; i++)
        {
            if (i < pointnum - FixedNodesNum)
            {
                for (int j = pointnum - FixedNodesNum; j < pointnum; j++)
                {
                    float dis = Vector3.Distance(Position[i], Position[j]);
                    if (dis < 0.8f)
                    {
                        Vector3 vec = Position[i] - Position[j];
                        Vector3 nor = Vector3.Normalize(vec);
                        Force[i] += nor * (1 - dis) * GroundStiffness; // 果冻粒子与地面粒子间的刚度，越大被反弹得越厉害
                        //Velocity[i] *= 0.1f;// 果冻粒子撞击地面后，有部分动能被吸收，表现为速度减少
                        Kmat[i, i] += 0.5f;
                        Kmat[i, j] += -1f;
                    }
                    if(dis < 1.5f)
                    {
                        Velocity[i] *= 0.1f;
                    }
                }
                for (int j = 0; j < Constrain[0, i]; j++)
                {
                    int k = Constrain[j + 1, i];
                    Kmat[i, i] += 1;
                    Kmat[i, k] += -1;
                    Vector3 vec = Position[i] - Position[k];
                    float dis = Vector3.Distance(Position[i], Position[k]);
                    Vector3 nor = Vector3.Normalize(vec);
                    Force[i] += nor * (RestLength[j + 1, i] - dis) * JellyStiffness; // 果冻内部的刚度，数字越大，果冻结构保持得越好，难以变形
                }
                Force[i].z -= 0.001f;
            }
            else
            {
                Kmat[i, i] = 1;
            }
        }

        Jacobi();
        for (int i = 0; i < pointnum; i++)
        {
            if (i < pointnum - FixedNodesNum)
            {
                Velocity[i] = (res[i] * 0.05f + Velocity[i]) * 0.5f;
                Position[i] = 0.1f * Velocity[i] + Position[i];
                verties[i] = Position[i];
                nodes[i].transform.position = Position[i];
            }
        }
        for (int i = 0; i < Triangles.Count; i++)
        {
            verties[i] = Position[Triangles[i]];
        }
        mesh.vertices = verties;
        mesh.RecalculateNormals();
        mesh.RecalculateBounds();
    }
    void Update()
    {
        NowFrame++;
        if (NowFrame < 200) return;
        if (EdgeIndex < Edgenum)
        {
            if (EdgeFront[EdgeIndex] == true)
            {
                GenNewTriangle();
                verties = new Vector3[Triangles.Count];
                tries = new int[Triangles.Count];
                for (int i = 0; i < Triangles.Count; i++)
                {
                    verties[i] = point[Triangles[i]];
                    tries[i] = i;
                }
                mesh.vertices = verties;
                mesh.triangles = tries;
                mesh.RecalculateNormals();
                mesh.RecalculateBounds();
            }
            else
            {
                EdgeIndex++;
            }

        }
        else
        {
            if (pointnum == -1) FemInit();
            Fem();
        }
    }
}
