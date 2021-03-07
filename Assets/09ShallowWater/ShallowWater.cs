using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class ShallowWater : MonoBehaviour
{
    private const int nmax = 64;
    Vector3[] verts = new Vector3[nmax * nmax];//上一帧的高度
    protected MeshFilter meshFilter;
    protected Mesh mesh;
    int NowFrame = 0;
    public Material Mat;

    float[,] h, u, v, U1, U2, U3, F1, F2, F3, G1, G2, G3, f1, f2, f3, g1, g2, g3, u1pred, u2pred, u3pred, z;
    float dx, dy, dt, cfl = 0.8f;
    int scheme = 1;//Lax Frieshs == 1,MacCormack == 2
    int watercondition = 2;//1代表Dam Break，2代表下雨，3代表有坡度的DamBreak,4代表设置好的反射边界条件
    int Boundary = 1;//1代表不反射，2代表反射
    int[,] drop;//用于下雨的情况
	float grad = 1.0f;
    // Start is called before the first frame update
    void Start()
    {

        h = new float[nmax, nmax]; u = new float[nmax, nmax]; v = new float[nmax, nmax];
        U1 = new float[nmax, nmax]; U2 = new float[nmax, nmax]; U3 = new float[nmax, nmax];
        F1 = new float[nmax, nmax]; F2 = new float[nmax, nmax]; F3 = new float[nmax, nmax];
        G1 = new float[nmax, nmax]; G2 = new float[nmax, nmax]; G3 = new float[nmax, nmax];
        z = new float[nmax + 1, nmax + 1]; drop = new int[nmax, nmax]; 
        if (scheme == 1)
        {
            dx = dy = 10.0f / nmax;
            f1 = new float[nmax + 1, nmax]; f2 = new float[nmax + 1, nmax]; f3 = new float[nmax + 1, nmax];
            g1 = new float[nmax, nmax + 1]; g2 = new float[nmax, nmax + 1]; g3 = new float[nmax, nmax + 1];
        }
        else if (scheme == 2)
        {
            dx = dy = 100.0f / (nmax - 3);
            u1pred = new float[nmax, nmax]; u2pred = new float[nmax, nmax]; u3pred = new float[nmax, nmax];
        }


        for (int i = 0; i <= nmax; i++)
        {
            for (int j = 0; j <= nmax; j++)
            {
                z[i, j] = 0;
                if (watercondition == 3)
                {
                    if (i < nmax / 2 && i > nmax * 3 / 8) z[i, j] = z[i - 1, j] + grad;
                    if (i >= nmax / 2 && i < nmax * 5 / 8) z[i, j] = z[i - 1, j] - grad;
                }
            }
        }

        for (int i = 0; i < nmax; i++)
        {
            for (int j = 0; j < nmax; j++)
            {
                
                if (i < nmax / 2) h[i, j] = 2.0f;
                else h[i, j] = 1.0f;
                if (watercondition == 2) h[i, j] = 2.0f;//仅下雨，而无Dam Break
				if(watercondition == 4)
				{
					if (i < nmax / 2 && j < nmax / 2) h[i, j] = 2.0f;
                    else h[i, j] = 1.0f;
					Boundary = 2;
				}
                verts[j * nmax + i] = new Vector3((float)i / (float)nmax * 10.0f - 5.0f, h[i, j], (float)j / (float)nmax * 10.0f - 5.0f);
                u[i, j] = v[i, j] = 0.0f;
                U1[i, j] = h[i, j];
                U2[i, j] = h[i, j] * u[i, j];
                U3[i, j] = h[i, j] * v[i, j];
                F1[i, j] = h[i, j] * u[i, j];
                F2[i, j] = h[i, j] * u[i, j] * u[i, j] + 0.5f * 10.0f * h[i, j] * h[i, j];
                F3[i, j] = h[i, j] * u[i, j] * v[i, j];
                G1[i, j] = h[i, j] * v[i, j];
                G2[i, j] = h[i, j] * u[i, j] * v[i, j];
                G3[i, j] = h[i, j] * v[i, j] * v[i, j] + 0.5f * 10.0f * h[i, j] * h[i, j];
                drop[i, j] = 1;
            }
        }
        mesh = new Mesh();
        Vector2[] uvs = new Vector2[nmax * nmax];
        for (int i = 0; i < nmax; i++)
        {
            for (int j = 0; j < nmax; j++)
            {
                uvs[i * nmax + j] = new Vector2(i, j);
            }
        }
        MeshRenderer meshrender = GetComponent<MeshRenderer>();
        meshrender.material = Mat;
        mesh.name = gameObject.name;

        mesh.vertices = verts;
        mesh.uv = uvs;
        mesh.triangles = GenerateTries();
        mesh.RecalculateBounds();

        meshFilter = gameObject.AddComponent<MeshFilter>();
        meshFilter.mesh = mesh;

        Application.targetFrameRate = 60;
    }
    private int[] GenerateTries()
    {
        var tries = new int[mesh.vertices.Length * 6];
        for (int i = 0; i < nmax - 1; i++)
        {
            for (int j = 0; j < nmax - 1; j++)
            {
                int TriIdx = (j * (nmax - 1) + i) * 6;
                int VerIdx = j * nmax + i;
                tries[TriIdx] = VerIdx;
                tries[TriIdx + 1] = VerIdx + nmax;
                tries[TriIdx + 2] = VerIdx + nmax + 1;

                /*
                 1 ------  2/4
                 |             |
                 |             |
               0/3 ------- 5
                 
                 */

                tries[TriIdx + 3] = VerIdx;
                tries[TriIdx + 4] = VerIdx + nmax + 1;
                tries[TriIdx + 5] = VerIdx + 1;
            }
        }
        return tries;
    }
    void Update()
    {
        if (scheme == 1)
        {
            Lax_Friedrichs();
        }
        else if (scheme == 2)
        {
            MacCormack();
        }
        mesh.vertices = verts;
        mesh.RecalculateNormals();
    }
    void Lax_Friedrichs()
    {
        Vector3[] Verts = new Vector3[nmax * nmax];
        float tempx = 0, tempy = 0;
        for (int i = 0; i < nmax; i++)
        {
            for (int j = 0; j < nmax; j++)
            {
                float temp = 0;
                temp = Mathf.Abs(u[i, j]) + Mathf.Sqrt(10 * h[i, j]);
                if (temp > tempx) tempx = temp;
                temp = Mathf.Abs(v[i, j]) + Mathf.Sqrt(10 * h[i, j]);
                if (temp > tempy) tempy = temp;
            }
        }
        dt = 0.5f * cfl * Mathf.Min(dx / tempx, dy / tempy);
        for (int i = 0; i < nmax + 1; i++)
        {
            for (int j = 0; j < nmax; j++)
            {
                if (i == 0)
                {
                    f1[i, j] = F1[i, j];
                    f2[i, j] = F2[i, j];
                    f3[i, j] = F3[i, j];
                }
                else if (i == nmax)
                {
                    f1[i, j] = F1[i - 1, j];
                    f2[i, j] = F2[i - 1, j];
                    f3[i, j] = F3[i - 1, j];
                }
                else
                {
                    f1[i, j] = 0.5f * (F1[i, j] + F1[i - 1, j]) - dx / dt * 0.25f * cfl * (U1[i, j] - U1[i - 1, j]);
                    f2[i, j] = 0.5f * (F2[i, j] + F2[i - 1, j]) - dx / dt * 0.25f * cfl * (U2[i, j] - U2[i - 1, j]);
                    f3[i, j] = 0.5f * (F3[i, j] + F3[i - 1, j]) - dx / dt * 0.25f * cfl * (U3[i, j] - U3[i - 1, j]);
                }
            }
        }
        for (int j = 0; j < nmax + 1; j++)
        {
            for (int i = 0; i < nmax; i++)
            {
                if (j == 0)
                {
                    g1[i, j] = G1[i, j];
                    g2[i, j] = G2[i, j];
                    g3[i, j] = G3[i, j];
                }
                else if (j == nmax)
                {
                    g1[i, j] = G1[i, j - 1];
                    g2[i, j] = G2[i, j - 1];
                    g3[i, j] = G3[i, j - 1];
                }
                else
                {
                    g1[i, j] = 0.5f * (G1[i, j] + G1[i, j - 1]) - dy / dt * 0.25f * cfl * (U1[i, j] - U1[i, j - 1]);
                    g2[i, j] = 0.5f * (G2[i, j] + G2[i, j - 1]) - dy / dt * 0.25f * cfl * (U2[i, j] - U2[i, j - 1]);
                    g3[i, j] = 0.5f * (G3[i, j] + G3[i, j - 1]) - dy / dt * 0.25f * cfl * (U3[i, j] - U3[i, j - 1]);
                }
            }
        }
        for (int i = 0; i < nmax; i++)
        {
            for (int j = 0; j < nmax; j++)
            {
                    U1[i, j] = U1[i, j] - dt / dx * (f1[i + 1, j] - f1[i, j]) - dt / dy * (g1[i, j + 1] - g1[i, j]);
                    U2[i, j] = U2[i, j] - dt / dx * (f2[i + 1, j] - f2[i, j]) - dt / dy * (g2[i, j + 1] - g2[i, j]) - dt * 10.0f * h[i, j] * (z[i + 1, j] - z[i, j]);
                    U3[i, j] = U3[i, j] - dt / dx * (f3[i + 1, j] - f3[i, j]) - dt / dy * (g3[i, j + 1] - g3[i, j]) - dt * 10.0f * h[i, j] * (z[i, j + 1] - z[i, j]);


                if (watercondition == 2)//下雨的情况
                {
                    if (Random.Range(0f, 1f) >= 0.999995f)
                    {
                        int posx = i, posy = j;
                        if (posx < 3) posx = 3;
                        if (posy < 3) posy = 3;
                        for (int pi = posx - 3; pi <= posx; pi++)
                        {
                            for (int pj = posy - 3; pj <= posy; pj++)
                            {
                                drop[pi, pj] = 10;
                            }
                        }

                    }

                    if (drop[i, j] > 0)
                    {
                        U1[i, j] = Mathf.Sin(drop[i, j] / 3.0f) * 0.5f + 1;
                        drop[i, j]--;
                    }
                }
				if(watercondition == 3)//保证源头一直有水流出
				{
					if(i == 0)U1[i,j] = 2.0f;
				}
				h[i, j] = U1[i, j];
                u[i, j] = U2[i, j] / U1[i, j];
                v[i, j] = U3[i, j] / U1[i, j];
            }
        }
        BoundaryCondition();
        for (int i = 0; i < nmax; i++)
        {
            for (int j = 0; j < nmax; j++)
            {
                verts[j * nmax + i] = new Vector3((float)i / (float)nmax * 10.0f - 5.0f, h[i, j], (float)j / (float)nmax * 10.0f - 5.0f);
                F1[i, j] = h[i, j] * u[i, j];
                F2[i, j] = h[i, j] * u[i, j] * u[i, j] + 0.5f * 10.0f * h[i, j] * h[i, j];
                F3[i, j] = h[i, j] * u[i, j] * v[i, j];
                G1[i, j] = h[i, j] * v[i, j];
                G2[i, j] = h[i, j] * u[i, j] * v[i, j];
                G3[i, j] = h[i, j] * v[i, j] * v[i, j] + 0.5f * 10.0f * h[i, j] * h[i, j];
            }
        }
    }

    void MacCormack()
    {
        Vector3[] Verts = new Vector3[nmax * nmax];
        float tempx = 0, tempy = 0;
        for (int i = 0; i < nmax; i++)
        {
            for (int j = 0; j < nmax; j++)
            {
                float temp = 0;
                temp = Mathf.Abs(u[i, j]) + Mathf.Sqrt(10 * h[i, j]);
                if (temp > tempx) tempx = temp;
                temp = Mathf.Abs(v[i, j]) + Mathf.Sqrt(10 * h[i, j]);
                if (temp > tempy) tempy = temp;
            }
        }
        dt = cfl * Mathf.Min(dx / tempx, dy / tempy);
        float c = dt / dx;
        for (int i = 1; i < nmax - 1; i++)
        {
            for (int j = 1; j < nmax - 1; j++)
            {
                u1pred[i, j] = U1[i, j] - c * (F1[i + 1, j] - F1[i, j]) - c * (G1[i, j + 1] - G1[i, j]);
                u2pred[i, j] = U2[i, j] - c * (F2[i + 1, j] - F2[i, j]) - c * (G2[i, j + 1] - G2[i, j]);
                u3pred[i, j] = U3[i, j] - c * (F3[i + 1, j] - F3[i, j]) - c * (G3[i, j + 1] - G3[i, j]);
                h[i, j] = u1pred[i, j];
                u[i, j] = u2pred[i, j] / u1pred[i, j];
                v[i, j] = u3pred[i, j] / u1pred[i, j];
            }
        }
        BoundaryCondition();
        for (int i = 0; i < nmax; i++)
        {
            for (int j = 0; j < nmax; j++)
            {
                F1[i, j] = h[i, j] * u[i, j];
                F2[i, j] = h[i, j] * u[i, j] * u[i, j] + 0.5f * 10.0f * h[i, j] * h[i, j];
                F3[i, j] = h[i, j] * u[i, j] * v[i, j];
                G1[i, j] = h[i, j] * v[i, j];
                G2[i, j] = h[i, j] * u[i, j] * v[i, j];
                G3[i, j] = h[i, j] * v[i, j] * v[i, j] + 0.5f * 10.0f * h[i, j] * h[i, j];
            }
        }
        for (int i = 1; i < nmax - 1; i++)
        {
            for (int j = 1; j < nmax - 1; j++)
            {
                    U1[i, j] = 0.5f * (U1[i, j] + u1pred[i, j] - c * (F1[i, j] - F1[i - 1, j]) - c * (G1[i, j] - G1[i, j - 1]));
                    U2[i, j] = 0.5f * (U2[i, j] + u2pred[i, j] - c * (F2[i, j] - F2[i - 1, j]) - c * (G2[i, j] - G2[i, j - 1]));
                    U3[i, j] = 0.5f * (U3[i, j] + u3pred[i, j] - c * (F3[i, j] - F3[i - 1, j]) - c * (G3[i, j] - G3[i, j - 1]));


                if (watercondition == 2)//下雨的情况
                {
                    if (Random.Range(0f, 1f) >= 0.999995f)
                    {
                        int posx = i, posy = j;
                        if (posx < 3) posx = 3;
                        if (posy < 3) posy = 3;
                        for (int pi = posx - 3; pi <= posx; pi++)
                        {
                            for (int pj = posy - 3; pj <= posy; pj++)
                            {
                                drop[pi, pj] = 10;
                            }
                        }

                    }

                    if (drop[i, j] > 0)
                    {
                        U1[i, j] = Mathf.Sin(drop[i, j] / 3.0f) * 0.5f + 1;
                        drop[i, j]--;
                    }
                }
				if(watercondition == 3)//保证源头一直有水流出
				{
					if(i == 0)U1[i,j] = 2.0f;
				}
                h[i, j] = U1[i, j];
                u[i, j] = U2[i, j] / U1[i, j];
                v[i, j] = U3[i, j] / U1[i, j];
            }
        }
        BoundaryCondition();
        for (int i = 0; i < nmax; i++)
        {
            for (int j = 0; j < nmax; j++)
            {
                verts[j * nmax + i] = new Vector3((float)i / (float)nmax * 10.0f - 5.0f, h[i, j], (float)j / (float)nmax * 10.0f - 5.0f);
                F1[i, j] = h[i, j] * u[i, j];
                F2[i, j] = h[i, j] * u[i, j] * u[i, j] + 0.5f * 10.0f * h[i, j] * h[i, j];
                F3[i, j] = h[i, j] * u[i, j] * v[i, j];
                G1[i, j] = h[i, j] * v[i, j];
                G2[i, j] = h[i, j] * u[i, j] * v[i, j];
                G3[i, j] = h[i, j] * v[i, j] * v[i, j] + 0.5f * 10.0f * h[i, j] * h[i, j];
            }
        }
    }

    void BoundaryCondition()
    {
        if (Boundary == 1)
        {
            for (int i = 0; i < nmax; i++)
            {
                u[i, 0] = u[i, 1];
                u[i, nmax - 1] = u[i, nmax - 2];
                v[i, 0] = v[i, 1];
                v[i, nmax - 1] = v[i, nmax - 2];
                h[i, 0] = h[i, 1];
                h[i, nmax - 1] = h[i, nmax - 2];
            }
            for (int i = 0; i < nmax; i++)
            {
                u[0, i] = u[1, i];
                u[nmax - 1, i] = u[nmax - 2, i];
                v[0, i] = v[1, i];
                v[nmax - 1, i] = v[nmax - 2, i];
                h[0, i] = h[1, i];
                h[nmax - 1, i] = h[nmax - 2, i];
            }
        }
        if (Boundary == 2)
        {
            for (int i = 0; i < nmax; i++)
            {
                u[i, 0] = -u[i, 1];
                u[i, nmax - 1] = -u[i, nmax - 2];
                v[i, 0] = -v[i, 1];
                v[i, nmax - 1] = -v[i, nmax - 2];
                h[i, 0] = h[i, 1];
                h[i, nmax - 1] = h[i, nmax - 2];
            }
            for (int i = 0; i < nmax; i++)
            {
                u[0, i] = -u[1, i];
                u[nmax - 1, i] = -u[nmax - 2, i];
                v[0, i] = -v[1, i];
                v[nmax - 1, i] = -v[nmax - 2, i];
                h[0, i] = h[1, i];
                h[nmax - 1, i] = h[nmax - 2, i];
            }
        }
    }
}
