using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Cavity2D : MonoBehaviour
{
    /*
     Lid Driven Cavity Flow 使用unity在CPU上的实现
    这个问题的代码在github上一搜一大把，大部分都是python或matlab，比unity调试起来舒服多了
     */
    private Texture2D PressureTex;
    private Texture2D UTex;
    private Texture2D VTex;
    const int NumberGrid = 10;
    const int TexelNumber = NumberGrid;
    float[,] u = new float[NumberGrid, NumberGrid];//x轴方向的速度
    float[,] v = new float[NumberGrid, NumberGrid];//y轴方向的速度
    float[,] p = new float[NumberGrid, NumberGrid];//压力
    float[,] b = new float[NumberGrid, NumberGrid];//辅助压力的项，类似于散度

    float[,] un = new float[NumberGrid, NumberGrid];
    float[,] vn = new float[NumberGrid, NumberGrid];
    float[,] pn = new float[NumberGrid, NumberGrid];

    float dt = 0.001f, dx = 2 / (NumberGrid - 1), dy = 2 / (NumberGrid - 1), nu = 0.01f, rho = 1;
    void Start()
    {
        Application.targetFrameRate = 100;
        PressureTex = new Texture2D(TexelNumber, TexelNumber, TextureFormat.RGB24, false);
        UTex = new Texture2D(TexelNumber, TexelNumber, TextureFormat.RGB24, false);
        VTex = new Texture2D(TexelNumber, TexelNumber, TextureFormat.RGB24, false);
        PressureTex.filterMode = UTex.filterMode = VTex.filterMode = FilterMode.Point;
        dx = 2.0f / (NumberGrid - 1);
        dy = 2.0f / (NumberGrid - 1);
        for (int j = 0; j < TexelNumber; j++)
        {
            for (int i = 0; i < TexelNumber; i++)
            {
                PressureTex.SetPixel(i, j, new Color(i / (float)TexelNumber, 0, 0));
            }
        }
        PressureTex.Apply();
        for (int j = 0; j < NumberGrid; j++)
            for (int i = 0; i < NumberGrid; i++)
                u[j, i] = v[i, j] = 0;
        for (int j = 0; j < NumberGrid; j++)
            for (int i = 0; i < NumberGrid; i++)
                p[i, j] = 0;
    }

    private void Update()
    {
        for (int j = 0; j < NumberGrid; j++)
            for (int i = 0; i < NumberGrid; i++)
            {
                un[i, j] = u[i, j];
                vn[i, j] = v[i, j];
            }
        for (int j = 1; j < NumberGrid - 1; j++)
        {
            for (int i = 1; i < NumberGrid - 1; i++)
            {
                b[i, j] = (u[i + 1, j] - u[i - 1, j]) / (2 * dx) + (v[i, j + 1] - v[i, j - 1]) / (2 * dy);
                b[i, j] -=  ((((u[i + 1, j] - u[i - 1, j]) / (2 * dx)) * ((u[i + 1, j] - u[i - 1, j]) / (2 * dx)))
                           + 2 * (((u[i, j + 1] - u[i, j - 1]) / (2 * dy)) * ((v[i + 1, j] - v[i - 1, j]) / (2 * dy)))
                           + (((v[i, j + 1] - v[i, j - 1]) / (2 * dy)) * ((v[i, j + 1] - v[i, j - 1]) / (2 * dy))));
            }
        }
        for (int k = 0; k < 20; k++)
        {

            for (int j = 0; j < NumberGrid; j++)
                for (int i = 0; i < NumberGrid; i++) pn[i, j] = p[i, j];
            for (int j = 1; j < NumberGrid - 1; j++)
            {
                for (int i = 1; i < NumberGrid - 1; i++)
                {

                    p[i, j] = (((pn[i, j - 1] + pn[i, j + 1]) * (dx * dx) + (pn[i - 1, j] + pn[i + 1, j]) * (dy * dy)) / (
                                2 * (dx * dx + dy * dy)) - rho * (b[i, j] * ((dy * dx) * (dy * dx)) / (2 * (dx * dx + dy * dy))));

                }
            }
            for (int i = 0; i < NumberGrid; i++) p[0, i] = p[1, i];
            for (int i = 0; i < NumberGrid; i++) p[i, 0] = p[i, 1];
            for (int i = 0; i < NumberGrid; i++) p[NumberGrid - 1, i] = p[NumberGrid - 2, i];
            for (int i = 0; i < NumberGrid; i++) p[i, NumberGrid - 1] = 0;
        }


        for (int j = 1; j < NumberGrid - 1; j++)
        {
            for (int i = 1; i < NumberGrid - 1; i++)
            {
                float advection = (un[i, j] - (un[i, j] * (dt / dx) * (un[i, j] - un[i - 1, j])) - vn[i, j] * (dt / dy) * (un[i, j] - un[i, j - 1]));//平流项
                float diffuse = (nu * dt / (dx * dx)) * (un[i + 1, j] - 2 * un[i, j] + un[i - 1, j]) + (nu * dt / (dy * dx)) * (un[i, j - 1] - 2 * un[i, j] + un[i, j + 1]);//粘性项
                float density = (dt / (2 * rho * dx)) * (p[i + 1, j] - p[i - 1, j]);//密度项
                u[i, j] = advection + diffuse - density;
                advection = (vn[i, j] - (un[i, j] * dt / dx * (vn[i, j] - vn[i - 1, j])) - vn[i, j] * dt / dy * (vn[i, j] - vn[i, j - 1]));
                diffuse = (nu * dt / (dx * dx)) * (vn[i + 1, j] - 2 * vn[i, j] + vn[i - 1, j]) + (nu * dt / (dx * dx)) * (vn[i, j - 1] - 2 * vn[i, j] + vn[i, j + 1]);
                density = (dt / (2 * rho * dy)) * (p[i, j + 1] - p[i - 1, j - 1]);
                v[i, j] = advection + diffuse - density;
            }
        }

        for (int i = 0; i < NumberGrid; i++)
        {
            u[0, i] = 0.0f;
            u[NumberGrid - 1, i] = 0.0f;
            v[0, i] = 0;
            v[i, 0] = 0;
        }
        for (int i = 0; i < NumberGrid; i++)
        {
            u[i, 0] = 0.0f;
            u[i, NumberGrid - 1] = 1.0f;
            v[i, NumberGrid - 1] = 0;
            v[NumberGrid - 1, i] = 0;
        }
        for (int j = 0; j < NumberGrid; j++)
        {
            for (int i = 0; i < NumberGrid; i++)
            {
                if (u[i, j] >= 0)
                    UTex.SetPixel(i, j, new Color(u[i, j], 0, 0));
                else
                    UTex.SetPixel(i, j, new Color(0, 0, -u[i, j] * 100.0f));
                if (v[i, j] >= 0)
                    VTex.SetPixel(i, j, new Color(v[i, j] * 100.0f, 0, 0));
                else
                    VTex.SetPixel(i, j, new Color(0, 0, -v[i, j] * 100.0f));

            }
        }
        UTex.Apply();
        VTex.Apply();
    }
    private void OnRenderImage(RenderTexture source, RenderTexture destination)
    {
        Graphics.Blit(UTex, destination);
    }
}
