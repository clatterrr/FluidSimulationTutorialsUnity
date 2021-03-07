using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class SmokeManager : MonoBehaviour
{
	public Material AdvectionMat;
	public Material SplatMat;
	public Material CurlMat;
	public Material VorticityMat;
	public Material DivergenceMat;
	public Material PressureMat;
	public Material SubtractMat;
	private int TexWidth = Screen.width;
	private int TexHeight = Screen.height;
	
	private RenderTexture DivergenceRT;
	private RenderTexture CurlRT;
	private RenderTexture DensityRT;
	private RenderTexture DensityRT2;
	private RenderTexture VelocityRT;
	private RenderTexture VelocityRT2;
	private RenderTexture PressureRT;
	private RenderTexture PressureRT2;
	
	private float dt = 0.01f;
	private float MouseDX = 0.0f, MouseDY = 0.0f, MouseX, MouseY;

	[Range(0.95f, 1.0f)]
	public float DensityDiffusion = 0.995f;//密度消失速度，此值越大则粘度越小，越容易看到烟雾效果
	[Range(0.95f, 1.0f)]
	public float VelocityDiffusion = 0.995f;//速度扩散速度，此值越大则粘度越小，越容易看到烟雾效果
	[Range(1, 60)]
	public int Iterations = 50;//泊松方程迭代次数
	[Range(0, 60)]
	public float Vorticity = 50f;//控制漩涡缩放
	[Range(0.0001f, 0.005f)]
	public float SplatRadius = 0.001f;//鼠标点击产生烟雾的半径
	[Range(1, 15)]
	public float MouseForceScale = 10.0f;//鼠标拖动力度缩放
	
    void Start()
    {
        DivergenceRT = new RenderTexture(TexWidth, TexHeight, 0, RenderTextureFormat.RHalf); DivergenceRT.Create();
		CurlRT = new RenderTexture(TexWidth, TexHeight, 0, RenderTextureFormat.RGHalf); CurlRT.Create();
		DensityRT = new RenderTexture(TexWidth, TexHeight, 0, RenderTextureFormat.ARGBHalf); DensityRT.Create();
		DensityRT2 = new RenderTexture(TexWidth, TexHeight, 0, RenderTextureFormat.ARGBHalf); DensityRT2.Create();
		VelocityRT = new RenderTexture(TexWidth, TexHeight, 0, RenderTextureFormat.RGHalf); VelocityRT.Create();
		VelocityRT2 = new RenderTexture(TexWidth, TexHeight, 0, RenderTextureFormat.RGHalf); VelocityRT2.Create();
		PressureRT = new RenderTexture(TexWidth, TexHeight, 0, RenderTextureFormat.RHalf); PressureRT.Create();
		PressureRT2 = new RenderTexture(TexWidth, TexHeight, 0, RenderTextureFormat.RHalf); PressureRT2.Create();
    }

	void OnRenderImage(RenderTexture source, RenderTexture destination)
	{
		//鼠标输入相关
		MouseDX = (Input.mousePosition.x - MouseX) * MouseForceScale;
		MouseDY = (Input.mousePosition.y - MouseY) * MouseForceScale;
		MouseX = Input.mousePosition.x;
		MouseY = Input.mousePosition.y;

		//第一步：平流速度
		AdvectionMat.SetTexture("VelocityTex", VelocityRT2);
		AdvectionMat.SetTexture("VelocityDensityTex", VelocityRT2);
		AdvectionMat.SetFloat("dt", dt);
		AdvectionMat.SetFloat("dissipation", VelocityDiffusion);
		Graphics.Blit(VelocityRT2, VelocityRT, AdvectionMat);
		Graphics.Blit(VelocityRT, VelocityRT2);

		if (Input.GetMouseButton(0))
		{
			//第二步：鼠标拖动会影响速度
			SplatMat.SetTexture("VelocityDensityTex", VelocityRT2);
			SplatMat.SetVector("color", new Vector3(MouseDX, MouseDY, 1f));
			SplatMat.SetVector("pointerpos", new Vector2(MouseX, MouseY));
			SplatMat.SetFloat("radius", SplatRadius);
			Graphics.Blit(VelocityRT2, VelocityRT, SplatMat);
			Graphics.Blit(VelocityRT, VelocityRT2);

			//第三步：鼠标拖动也会影响密度
			SplatMat.SetTexture("VelocityDensityTex", DensityRT2);
			SplatMat.SetVector("color", new Vector3(1.0f, 0.0f, 0.0f));
			Graphics.Blit(DensityRT2, DensityRT, SplatMat);
			Graphics.Blit(DensityRT, DensityRT2);
		}

		//第四步：计算Curl
		CurlMat.SetTexture("VelocityTex", VelocityRT2);
		Graphics.Blit(VelocityRT2, CurlRT, CurlMat);

		//第五步：计算旋度，更新速度，得到有散度的速度场
		VorticityMat.SetTexture("VelocityTex", VelocityRT2);
		VorticityMat.SetTexture("CurlTex", CurlRT);
		VorticityMat.SetFloat("curl", Vorticity);
		VorticityMat.SetFloat("dt", dt);
		Graphics.Blit(VelocityRT2, VelocityRT, VorticityMat);
		Graphics.Blit(VelocityRT, VelocityRT2);

		//第六步：计算散度
		DivergenceMat.SetTexture("VelocityTex", VelocityRT2);
		Graphics.Blit(VelocityRT2, DivergenceRT, DivergenceMat);

		//第七步：计算压力
		PressureMat.SetTexture("DivergenceTex", DivergenceRT);
		for (int i = 0; i < Iterations; i++)
		{
			PressureMat.SetTexture("PressureTex", PressureRT2);
			Graphics.Blit(PressureRT2, PressureRT, PressureMat);
			Graphics.Blit(PressureRT, PressureRT2);
		}

		//第八步：速度场减去压力梯度，得到无散度的速度场
		SubtractMat.SetTexture("PressureTex", PressureRT2);
		SubtractMat.SetTexture("VelocityTex", VelocityRT2);
		Graphics.Blit(VelocityRT2, VelocityRT, SubtractMat);
		Graphics.Blit(VelocityRT, VelocityRT2);

		//第九步：用最终速度去平流密度
		AdvectionMat.SetTexture("VelocityTex", VelocityRT2);
		AdvectionMat.SetTexture("VelocityDensityTex", DensityRT2);
		AdvectionMat.SetFloat("dissipation", DensityDiffusion);
		Graphics.Blit(DensityRT2, DensityRT, AdvectionMat);
		Graphics.Blit(DensityRT, DensityRT2);

		//第十步：显示
		Graphics.Blit(DensityRT2, destination);
	}
}
