using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class VortexStreetManager : MonoBehaviour
{
    public Material DivergenceMat;
    public Material PressureMat;
    public Material SubtractMat;
    public Material AdvectionDyeMat;
    public Material AdvectionVelocityMat;
    public Material InitDyeMat;
    public Material BlockMat;
    public Material DisplayMat;
    public Material DisplayRainbowMat;
    public Material ViscosityMat;
    private int TexWidth = Screen.width;
    private int TexHeight = Screen.height;

    private RenderTexture DivergenceRT;
    private RenderTexture DyeRT;
    private RenderTexture DyeRT2;
    private RenderTexture VelocityRT;
    private RenderTexture VelocityRT2;
    private RenderTexture PressureRT;
    private RenderTexture PressureRT2;
    private RenderTexture InitDyeRT;
    private RenderTexture BlockRT;
    

    private float dt = 0.1f;

    void Start()
    {
        DivergenceRT = new RenderTexture(TexWidth, TexHeight, 0, RenderTextureFormat.RHalf); DivergenceRT.Create();
        DyeRT = new RenderTexture(TexWidth, TexHeight, 0, RenderTextureFormat.ARGBHalf); DyeRT.Create();
        DyeRT2 = new RenderTexture(TexWidth, TexHeight, 0, RenderTextureFormat.ARGBHalf); DyeRT2.Create();
        InitDyeRT = new RenderTexture(TexWidth, TexHeight, 0, RenderTextureFormat.ARGBHalf); InitDyeRT.Create();
        VelocityRT = new RenderTexture(TexWidth, TexHeight, 0, RenderTextureFormat.RGHalf); VelocityRT.Create();
        VelocityRT2 = new RenderTexture(TexWidth, TexHeight, 0, RenderTextureFormat.RGHalf); VelocityRT2.Create();
        PressureRT = new RenderTexture(TexWidth, TexHeight, 0, RenderTextureFormat.RHalf); PressureRT.Create();
        PressureRT2 = new RenderTexture(TexWidth, TexHeight, 0, RenderTextureFormat.RHalf); PressureRT2.Create();
        BlockRT = new RenderTexture(TexWidth, TexHeight, 0, RenderTextureFormat.ARGBHalf); BlockRT.Create();
        Graphics.Blit(null, InitDyeRT, InitDyeMat);
        Graphics.Blit(null, BlockRT, BlockMat);
    }

    void OnRenderImage(RenderTexture source, RenderTexture destination)
    {
        //第一步：平流速度
        AdvectionVelocityMat.SetTexture("VelocityTex", VelocityRT2);
        AdvectionVelocityMat.SetTexture("BlockTex", BlockRT);
        AdvectionVelocityMat.SetFloat("dt", dt);
        Graphics.Blit(VelocityRT2, VelocityRT, AdvectionVelocityMat);
        Graphics.Blit(VelocityRT, VelocityRT2);

        //第二步，加大流体粘性可抑制边界层分离现象
        for (int i = 0; i < 0; i++)
        {
            ViscosityMat.SetTexture("_VelocityTex", VelocityRT2);
            Graphics.Blit(VelocityRT2, VelocityRT, ViscosityMat);
            Graphics.Blit(VelocityRT, VelocityRT2);
        }

        //第三步：计算散度
        DivergenceMat.SetTexture("VelocityTex", VelocityRT2);
        Graphics.Blit(VelocityRT2, DivergenceRT, DivergenceMat);

        //第四步：计算压力
        PressureMat.SetTexture("DivergenceTex", DivergenceRT);
        for (int i = 0; i < 100; i++)
        {
            PressureMat.SetTexture("PressureTex", PressureRT2);
            Graphics.Blit(PressureRT2, PressureRT, PressureMat);
            Graphics.Blit(PressureRT, PressureRT2);
        }
        //第五步：速度场减去压力梯度，得到无散度的速度场
        SubtractMat.SetTexture("PressureTex", PressureRT2);
        SubtractMat.SetTexture("VelocityTex", VelocityRT2);
        Graphics.Blit(VelocityRT2, VelocityRT, SubtractMat);
        Graphics.Blit(VelocityRT, VelocityRT2);

        //第六步：用最终速度去平流密度
        Graphics.Blit(DyeRT, DyeRT2);
        AdvectionDyeMat.SetTexture("VelocityTex", VelocityRT2);
        AdvectionDyeMat.SetTexture("DensityTex", DyeRT2);
        AdvectionDyeMat.SetTexture("BlockTex", BlockRT);
        AdvectionDyeMat.SetTexture("InitDyeTex", InitDyeRT);
        AdvectionDyeMat.SetFloat("dt", dt);
        Graphics.Blit(DyeRT2, DyeRT, AdvectionDyeMat);

        //第七步：显示
        DisplayMat.SetTexture("BlockTex", BlockRT);
        DisplayRainbowMat.SetTexture("BlockTex", BlockRT);
        Graphics.Blit(DyeRT, destination, DisplayMat);
        //Graphics.Blit(VelocityRT2, destination, DisplayRainbowMat);
        //Graphics.Blit(PressureRT2, destination, DisplayRainbowMat);
    }
}
