using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class wave2d : MonoBehaviour
{
    public Material _mat;
    public RenderTexture _rtNull;//初始什么都没有
    public RenderTexture _rtm2;//上两帧图形
    public RenderTexture _rtm1;//上一帧图像
    public RenderTexture _rt0;//本帧图形
    [Range(2, 100)]
    public int _TexelNumber;//设定网格的数量
    [Range(0.5f, 2.0f)]
    public float _Speed;//设定波传播的速度

    private void Start()
    {
        Graphics.Blit(_rtNull, _rtm2);
        Graphics.Blit(_rtNull, _rtm1);
        Application.targetFrameRate = 10;
        _TexelNumber = 51;
        _Speed = 0.4f;
    }
    private void Update()
    {
        _mat.SetInt("_TexelNumber", _TexelNumber);
        _mat.SetFloat("_Speed", _Speed);
        _mat.SetTexture("_MainTex", _rtm2);
        _mat.SetTexture("_TimeTexM1", _rtm1);
    }

    private void OnRenderImage(RenderTexture source, RenderTexture destination)
    {
        Graphics.Blit(_rtm2, _rt0, _mat);//将上一帧的图像经过迭代计算后放到这一帧的图像上
        Graphics.Blit(_rt0, destination);//将这一帧图像输出的屏幕上
        Graphics.Blit(_rtm1, _rtm2);//将上帧图像变为上两帧图像
        Graphics.Blit(_rt0, _rtm1);//将一帧图像变为上一帧图像
    }
}
