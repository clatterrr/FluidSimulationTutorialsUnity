using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class ComputeManager : MonoBehaviour
{
    private const int nmax = 128;
    Vector3[] verts = new Vector3[nmax * nmax];
    protected MeshFilter meshFilter;
    protected Mesh mesh;
    public Material Mat;

    int kernel = 0;
    ComputeBuffer[] ABuffer = new ComputeBuffer[12];
    public ComputeShader ComputeUShader;
    public ComputeShader ComputeFGShader;
    public ComputeShader ComputefgHalfShader;
    public ComputeShader UpdateUShader;

    float dx = 0.1f, dt = 0.01f;
    struct huvStruct
    {
        public float h, u, v;
    }
    huvStruct[] huv;
    ComputeBuffer huvBuffer;
    struct UStruct
    {
        public float U1, U2, U3;
    }
    UStruct[] U;
    ComputeBuffer UBuffer;
    struct FGStruct
    {
        public float F1, F2, F3, G1, G2, G3;
    }
    FGStruct[] FG = new FGStruct[nmax * nmax];
    ComputeBuffer FGBuffer;

    struct fgHalfStruct
    {
        public float f1, f2, f3, g1, g2, g3;
    }
    fgHalfStruct[] fgHalf;
    ComputeBuffer fgHalfBuffer;

    struct Data
    {
        public float A;
        public float B;
        public float C;
    }
    void Start()
    {
        huv = new huvStruct[nmax * nmax];
        U = new UStruct[nmax * nmax];
        FG = new FGStruct[nmax * nmax];
        fgHalf = new fgHalfStruct[(nmax + 1) * nmax];

        huvBuffer = new ComputeBuffer(nmax * nmax, 3 * sizeof(float));
        UBuffer = new ComputeBuffer(nmax * nmax, 3 * sizeof(float));
        FGBuffer = new ComputeBuffer(nmax * nmax, 6 * sizeof(float));
        fgHalfBuffer = new ComputeBuffer((nmax + 1) * nmax, 6 * sizeof(float));

        for (int j = 0; j < nmax; j++)
        {
            for (int i = 0; i < nmax; i++)
            {
                int idx = j * nmax + i;
                if (i < nmax / 2) huv[idx].h = 2.0f;
                else huv[idx].h = 1.0f;
                verts[idx] = new Vector3((float)i / (float)nmax * 10.0f - 5.0f, huv[idx].h, (float)j / (float)nmax * 10.0f - 5.0f);
                huv[idx].u = huv[idx].v = 0.0f;
                //uvs[i * nmax + j] = new Vector2(i, j);
            }
        }
        mesh = new Mesh();
        MeshRenderer meshrender = GetComponent<MeshRenderer>();
        meshrender.material = Mat;
        mesh.name = gameObject.name;

        mesh.vertices = verts;
        mesh.triangles = GenerateTries();
        mesh.RecalculateBounds();

        meshFilter = gameObject.AddComponent<MeshFilter>();
        meshFilter.mesh = mesh;

        huvBuffer.SetData(huv);
        kernel = ComputeUShader.FindKernel("CSMain");
        ComputeUShader.SetBuffer(kernel, "huv", huvBuffer);
        ComputeUShader.SetBuffer(kernel, "U", UBuffer);
        ComputeUShader.SetInt("nmax", nmax);
        ComputeUShader.Dispatch(kernel, nmax * nmax / 32 + 1, 1, 1);
        UBuffer.GetData(U);

        kernel = ComputeFGShader.FindKernel("CSMain");
        ComputeFGShader.SetBuffer(kernel, "huv", huvBuffer);
        ComputeFGShader.SetBuffer(kernel, "FG", FGBuffer);
        ComputeFGShader.SetInt("nmax", nmax);
        ComputeFGShader.Dispatch(kernel, nmax * nmax / 32 + 1, 1, 1);
        FGBuffer.GetData(FG);

        Application.targetFrameRate = 100;


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
    // Update is called once per frame
    void Update()
    {
        UBuffer.SetData(U);
        FGBuffer.SetData(FG);
        fgHalfBuffer.SetData(fgHalf);
        kernel = ComputefgHalfShader.FindKernel("CSMain");
        ComputefgHalfShader.SetBuffer(kernel, "U", UBuffer);
        ComputefgHalfShader.SetBuffer(kernel, "FG", FGBuffer);
        ComputefgHalfShader.SetBuffer(kernel, "fgHalf", fgHalfBuffer);
        ComputefgHalfShader.SetInt("nmax", nmax);
        ComputefgHalfShader.SetFloat("dx_dt", dx / dt);
        ComputefgHalfShader.Dispatch(kernel, (nmax + 1) * nmax / 32 + 1, 1, 1);
        fgHalfBuffer.GetData(fgHalf);

        fgHalfBuffer.SetData(fgHalf);
        kernel = UpdateUShader.FindKernel("CSMain");
        UpdateUShader.SetBuffer(kernel, "fgHalf", fgHalfBuffer);
        UpdateUShader.SetBuffer(kernel, "U", UBuffer);
        UpdateUShader.SetBuffer(kernel, "huv", huvBuffer);
        UpdateUShader.SetInt("nmax", nmax);
        UpdateUShader.SetFloat("dt_dx", dt / dx);
        UpdateUShader.Dispatch(kernel, (nmax + 1) * nmax / 32 + 1, 1, 1);
        UBuffer.GetData(U);
        huvBuffer.GetData(huv);


        huvBuffer.SetData(huv);
        kernel = ComputeFGShader.FindKernel("CSMain");
        ComputeFGShader.SetBuffer(kernel, "huv", huvBuffer);
        ComputeFGShader.SetBuffer(kernel, "FG", FGBuffer);
        ComputeFGShader.SetInt("nmax", nmax);
        ComputeFGShader.Dispatch(kernel, nmax * nmax / 32 + 1, 1, 1);
        FGBuffer.GetData(FG);

        for (int j = 0; j < nmax; j++)
        {
            for (int i = 0; i < nmax; i++)
            {
                int idx = j * nmax + i;
                verts[idx] = new Vector3((float)i / (float)nmax * 10.0f - 5.0f, huv[idx].h, (float)j / (float)nmax * 10.0f - 5.0f);
            }
        }

        mesh.vertices = verts;
        mesh.RecalculateNormals();
    }
}
