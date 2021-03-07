using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Manager : MonoBehaviour
{
    public GameObject prefab;
    private float PlaneSize = 10.0f;
    private const int QuadNumber = 16; //箭头在一个方向上的数量
    private GameObject[,] Arrow = new GameObject[QuadNumber, QuadNumber];
    private float[] ArrowRotation = new float[QuadNumber * QuadNumber];

    public Material Advection;//对流
    public Material InitDensity;//用于密度
    public Material UpdateVelocity;//使用箭头画出速度
    public Material UpdateLeftDensity;//更新左侧的密度图

    public RenderTexture RtPreFrame;//前一帧图像
    public RenderTexture RtNowFrame;//本帧图像
    public RenderTexture RTVelocity;//速度图
    public RenderTexture RTDensityLeft;//左侧的密度图，用于更新
    void Start()
    {
        Application.targetFrameRate = 10;
        float QuadSize = PlaneSize / QuadNumber;

        for (int j = 0; j < QuadNumber; j++)
        {
            for (int i = 0; i < QuadNumber; i++)
            {
                float posy = j * QuadSize - PlaneSize / 2.0f + QuadSize / 2.0f;
                float posx = i * QuadSize - PlaneSize / 2.0f + QuadSize / 2.0f;
                Arrow[i, j] = Instantiate(prefab, new Vector3(posx, 0.1f, posy), Quaternion.identity);
                Arrow[i, j].transform.localScale = new Vector3(2.0f, 2.0f, 2.0f);
                Arrow[i, j].transform.Rotate(new Vector3(90.0f, 0.0f, 0.0f));
                Arrow[i, j].transform.Rotate(new Vector3(0.0f, 0.0f, 90.0f));
                ArrowRotation[j * QuadNumber + i] = 0.0f;
                int HalfNumber = QuadNumber / 2;
                if (i <= HalfNumber + 2 && i >= HalfNumber + 1 && j <= HalfNumber + 2 && j >= HalfNumber + 1)//右上角
                {
                    ArrowRotation[j * QuadNumber + i] = 315.0f;
                }
                if (i <= HalfNumber && i >= HalfNumber - 1 && j <= HalfNumber + 2 && j >= HalfNumber + 1)//左上角
                {
                    ArrowRotation[j * QuadNumber + i] = 45.0f;
                }
                if (i <= HalfNumber && i >= HalfNumber - 1 && j <= HalfNumber && j >= HalfNumber - 1)//左下角
                {
                    ArrowRotation[j * QuadNumber + i] = 135.0f;
                }
                if (i <= HalfNumber + 2 && i >= HalfNumber + 1 && j <= HalfNumber && j >= HalfNumber - 1)//右下角
                {
                    ArrowRotation[j * QuadNumber + i] = 225.0f;
                }
                Arrow[i, j].transform.Rotate(new Vector3(0.0f, 0.0f, ArrowRotation[j * QuadNumber + i]));
            }
        }
        Graphics.Blit(null, RtPreFrame, InitDensity);//绘制初始的黑白棋盘格纹理
        Graphics.Blit(RtPreFrame, RTDensityLeft);
    }
    private void Update()
    {
        Advection.SetTexture("_MainTex", RtPreFrame);
        Advection.SetTexture("_VelocityTex", RTVelocity);
        Advection.SetTexture("_DensityLeft", RTDensityLeft);
        UpdateVelocity.SetFloatArray("_ArrowRotation", ArrowRotation); 
    }
    private void OnRenderImage(RenderTexture source, RenderTexture destination)
    {

        Graphics.Blit(null, RTVelocity, UpdateVelocity);//绘制速度纹理

        Graphics.Blit(RTDensityLeft, RtPreFrame);//移动左侧的纹理
        Graphics.Blit(RtPreFrame, RTDensityLeft, UpdateLeftDensity);

        Graphics.Blit(RtNowFrame, RtPreFrame);
        Graphics.Blit(RtPreFrame, RtNowFrame, Advection);//移动之前的纹理
        Graphics.Blit(source, destination);
    }
}
