using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class advancingMananger : MonoBehaviour
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
    public enum TestCases
    {
        SqureAndHole,
        CirlceAndHole,
        TriangleAnd3Hole,
        CircleAndCrack

    }
    public TestCases testcases;

    void Start()
    {

        for (int i = 0; i < nmax; i++)
        {
            Head[i] = -1;
            EdgeFront[i] = false;
            PointFront[i] = 0;
        }

        if(testcases == TestCases.SqureAndHole)
        {
            Init1();
        }else if (testcases == TestCases.CirlceAndHole)
        {
            Init2();
        }else if(testcases == TestCases.TriangleAnd3Hole)
        {
            Init3();
        }else if (testcases == TestCases.CircleAndCrack)
        {
            Init4();
        }

        for (int i = 0; i < point.Count; i++)
        {
            PointFront[i] = 2;
            Instantiate(sphereprefab, point[i], Quaternion.identity);
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
        TriangleRadius = 0.9f;
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
        float CircleRadius = 2.0f;
        int CircleParticle = 16;
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
        TriangleRadius = 0.5f;
        height = TriangleRadius * 0.9f;
        Vector2 CircleCenter = new Vector2(4.0f, 4.0f);
        float CircleRadius = 5.0f;
        int CircleParticle = 64;
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
    void Init3()
    {
        TriangleRadius = 0.8f;
        height = TriangleRadius;
        int squarelen = 32; // 平板边缘的顶点数
        float distan = TriangleRadius; // 顶点之间的距离
        for (int i = 0; i < squarelen; i++)
        {
            point.Add(new Vector3(i * distan, 0.0f, 0));
            AddEdge(i, i + 1);
        }
        for (int i = 0; i < squarelen; i++)
        {
            point.Add(new Vector3((squarelen - (float)(i) / 2) * distan, 0.0f, i * distan));
            AddEdge(i + squarelen, i + 1 + squarelen);
        }
        for (int i = 0; i < squarelen; i++)
        {
            point.Add(new Vector3((squarelen - i) * distan / 2, 0.0f, (squarelen - i) * distan));
            int ed = i + 1 + squarelen * 2;
            if (i == squarelen - 1) ed = 0;
            AddEdge(i + squarelen * 2, ed);
        }

        Vector2 CircleCenter = new Vector2(12.0f, 16.0f);
        float CircleRadius = 1.0f;
        int CircleParticle = 12;
        for (int i = 0; i < CircleParticle; i++)
        {
            point.Add(new Vector3(CircleCenter.x + CircleRadius * Mathf.Cos(2 * Mathf.PI / CircleParticle * i), 0.0f, CircleCenter.y - CircleRadius * Mathf.Sin(2 * Mathf.PI / CircleParticle * i)));
            int ed = i + 1 + squarelen * 3;
            if (i == CircleParticle - 1) ed = squarelen * 3;
            AddEdge(i + squarelen * 3, ed);
        }

        Vector2 CircleCenter2 = new Vector2(15.0f, 10.0f);
        float CircleRadius2 = 2.5f;
        int CircleParticle2 = 30;
        for (int i = 0; i < CircleParticle2; i++)
        {
            point.Add(new Vector3(CircleCenter2.x + CircleRadius2 * Mathf.Cos(2 * Mathf.PI / CircleParticle2 * i), 0.0f, CircleCenter2.y - CircleRadius2 * Mathf.Sin(2 * Mathf.PI / CircleParticle2 * i)));
            int ed = i + 1 + CircleParticle + squarelen * 3;
            if (i == CircleParticle2 - 1) ed = CircleParticle + squarelen * 3;
            AddEdge(i + CircleParticle + squarelen * 3, ed);
        }

        Vector2 CircleCenter3 = new Vector2(10.0f, 5.0f);
        float CircleRadius3 = 1.5f;
        int CircleParticle3 = 18;
        for (int i = 0; i < CircleParticle3; i++)
        {
            point.Add(new Vector3(CircleCenter3.x + CircleRadius3 * Mathf.Cos(2 * Mathf.PI / CircleParticle3 * i), 0.0f, CircleCenter3.y - CircleRadius3 * Mathf.Sin(2 * Mathf.PI / CircleParticle3 * i)));
            int ed = i + 1 + CircleParticle + squarelen * 3 + CircleParticle2;
            if (i == CircleParticle3 - 1) ed = CircleParticle + squarelen * 3 + CircleParticle2;
            AddEdge(i + CircleParticle + squarelen * 3 + CircleParticle2, ed);
        }
    }

    void Init4()
    {

        Vector2 CircleCenter = new Vector2(0.0f, 0.0f);
        float CircleRadius = 4.0f;
        int CircleParticle = 32;
        for (int i = 0; i < CircleParticle; i++)
        {
            point.Add(new Vector3(CircleCenter.x + CircleRadius * Mathf.Cos(2 * Mathf.PI / CircleParticle * i), 0.0f, CircleCenter.y + CircleRadius * Mathf.Sin(2 * Mathf.PI / CircleParticle * i)));
            int ed = i + 1;
            if (i == CircleParticle - 1) ed = 0;
            AddEdge(i, ed);
        }

        TriangleRadius = 0.4f;
        height = TriangleRadius;
        int squarelen = 4; // 平板边缘的顶点数
        float distan = TriangleRadius; // 顶点之间的距离
        for (int i = 0; i < squarelen; i++)
        {
            point.Add(new Vector3(i * distan, 0.0f, 0));
            int st = i - 1 + CircleParticle;
            if (st < CircleParticle) st = squarelen * 3 + CircleParticle - 1;
            AddEdge(i + CircleParticle, st);
        }
        for (int i = 0; i < squarelen; i++)
        {
            point.Add(new Vector3((squarelen - (float)(i) / 2) * distan, 0.0f, i * distan));
            AddEdge(i + squarelen + CircleParticle, i + squarelen - 1 + CircleParticle);
        }
        for (int i = 0; i < squarelen; i++)
        {
            point.Add(new Vector3((squarelen - i) * distan / 2, 0.0f, (squarelen - i) * distan));
            AddEdge(i + squarelen * 2 + CircleParticle, i + squarelen * 2 - 1 + CircleParticle);
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
            Instantiate(sphereprefab, NewVec, Quaternion.identity);
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
    void Update()
    {
        NowFrame++;
        if ((NowFrame % 10 != 0) || (NowFrame < 200)) return;
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
    }
}
