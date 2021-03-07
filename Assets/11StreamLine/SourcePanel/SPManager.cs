using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class SPManager : MonoBehaviour
{
	public Material StreamLineMat;
	public enum DisplayType
	{
		StreamLineAndObject,
		StreamLine
	}
	public enum AirfoilType
	{
		NACA1412,
		Cylinder
	}
	[Range(3,100)]
	public int PanelNumber = 8;
	[Range(0.0f,1.0f)]
	public float Scale = 0.3f;//控制圆柱和机翼缩放
	public DisplayType Displaytype;
	public AirfoilType airfoiltype;
	[Range(-10.0f,10.0f)]
	public float Vinf = 1.0f;//气流速速度
	[Range(-89.0f, 89.0f)]
	public float AngleOfAttack;//攻角

	float[] XB = new float[256];//平板端点X轴坐标
    float[] YB = new float[256];//平板端点Y轴坐标
	float[] XC = new float[256];//平板端点X轴坐标
    float[] YC = new float[256];//平板端点Y轴坐标
	float[] SourceStrength = new float[256];

	//平板长度  平板的角度  平板的法向量  受到攻角影响的法向量 切向加速度 压力
	float[] PanelLength, phi, delta, beta, Vt, Cp;

	float[,] I, J, Amatrix, AInverse;

    void Update()
    {
		for(int i = 0;i < 256;i++)
			XB[i] = YB[i] = XC[i] = YC[i] = SourceStrength[i] = 0.0f;
		if(airfoiltype == AirfoilType.Cylinder)
		{
		  for (int i = 0;i < PanelNumber;i++)
          {
            
            float theta = -Mathf.PI * 2.0f/ PanelNumber * (i + 1);//需要是顺时针
            XB[i] = -Scale * Mathf.Cos(theta) + 0.5f;
            YB[i] = -Scale * Mathf.Sin(theta) + 0.5f;
          }
		  XB[PanelNumber] = XB[0];
		  YB[PanelNumber] = YB[0];
		}else if(airfoiltype == AirfoilType.NACA1412)
		{
			NACA1412();
			for (int i = 0; i < PanelNumber + 1; i++)
			{
				float dx = XB[i] - 0.5f;
				float dy = YB[i] - 0.0f;
				float angle = Mathf.Atan2(dy, dx) - AngleOfAttack * Mathf.Deg2Rad;
				float dis = Mathf.Sqrt(dx * dx + dy * dy);
				XB[i] = dis * Mathf.Cos(angle);
				YB[i] = dis * Mathf.Sin(angle);
			}
			for (int i = 0;i < PanelNumber + 1;i++)
			{
				XB[i] = XB[i] * Scale + 0.5f;
				YB[i] = YB[i] * Scale + 0.5f;
			}
		}
		
        I = new float[PanelNumber,PanelNumber];//用于计算法向速度的辅助矩阵
        J = new float[PanelNumber,PanelNumber];//用于计算切向速度的辅助矩阵
	    PanelLength = new float[PanelNumber];//平板长度
        phi = new float[PanelNumber];//平板的角度
        delta = new float[PanelNumber];//平板的法向量
        beta = new float[PanelNumber];//受到攻角影响的法向量


        Amatrix = new float[PanelNumber, PanelNumber];//
        AInverse = new float[PanelNumber, PanelNumber];//

        Vt = new float[PanelNumber];//平板上流体的切向速度
        Cp = new float[PanelNumber];//Pressure Coffient
		
		//计算中心控制点位置与角度
		 for (int i = 0; i < PanelNumber; i++)
        {
            XC[i] = 0.5f * (XB[i] + XB[(i + 1) % PanelNumber]);
            YC[i] = 0.5f * (YB[i] + YB[(i + 1) % PanelNumber]);
            float dx = XB[(i + 1) % PanelNumber] - XB[i];
            float dy = YB[(i + 1) % PanelNumber] - YB[i];
            PanelLength[i] = Mathf.Sqrt(dx * dx + dy * dy);
            phi[i] = Mathf.Atan2(dy, dx);
            if (phi[i] < 0) phi[i] += 2.0f * Mathf.PI;
            delta[i] = phi[i] + Mathf.PI / 2.0f;//平板法向量
            beta[i] = delta[i] - AngleOfAttack * Mathf.Deg2Rad;//受到攻角影响
        }
		
		//计算参数
        for(int i = 0;i < PanelNumber;i++)
        {
            for(int j = 0;j < PanelNumber;j++)
            {
                if (i == j) continue;
                float A = -(XC[i] - XB[j]) * Mathf.Cos(phi[j]) - (YC[i] - YB[j]) * Mathf.Sin(phi[j]);
                float B = (XC[i] - XB[j]) * (XC[i] - XB[j]) + (YC[i] - YB[j]) * (YC[i] - YB[j]);
                float Cn = Mathf.Sin(phi[i] - phi[j]);
                float Dn = -(XC[i] - XB[j]) * Mathf.Sin(phi[i]) + (YC[i] - YB[j]) * Mathf.Cos(phi[i]);
                float Ct = -Mathf.Cos(phi[i] - phi[j]);
                float Dt = (XC[i] - XB[j]) * Mathf.Cos(phi[i]) + (YC[i] - YB[j]) * Mathf.Sin(phi[i]);
                float E = Mathf.Sqrt(B - A * A);
                float term1 = 0.5f * Cn * Mathf.Log((PanelLength[j] * PanelLength[j] + 2 * A * PanelLength[j] + B) / B);
                float term2 = ((Dn - A * Cn) / E) * (Mathf.Atan2((PanelLength[j] + A), E) - Mathf.Atan2(A, E));
                I[i, j] = term1 + term2;
                term1 = 0.5f * Ct * Mathf.Log((PanelLength[j] * PanelLength[j] + 2 * A * PanelLength[j] + B) / B);
                term2 = ((Dt - A * Ct) / E) * (Mathf.Atan2((PanelLength[j] + A), E) - Mathf.Atan2(A, E));
                J[i, j] = term1 + term2;
            }
        }
		
		//写方程式Ax = b
		for (int i = 0; i < PanelNumber; i++)
        {
            for (int j = 0; j < PanelNumber; j++)
            {
                if (i == j)
				{
					Amatrix[i, j] = Mathf.PI;
					AInverse[i,j] = 1;
				}
                else
				{
                    Amatrix[i, j] = I[i, j];
					AInverse[i,j] = 0;
				}
            }
        }
		
		//约旦消元法求A的逆矩阵
		for(int i = 0;i < PanelNumber;i++)
		 {
			 float temp = Amatrix[i,i];//获取A对角线上的元素
			 for(int j = 0; j < PanelNumber;j++)
			 {
				 Amatrix[i,j] = Amatrix[i,j]/temp;
				 AInverse[i,j] = AInverse[i,j]/temp;
            }
			 for(int k = 0; k < PanelNumber;k++)
			 {
                if (i == k) continue;
				 float mul = Amatrix[k,i];//其它行对应的元素，现在要将它减成零
				 for(int j = 0;j < PanelNumber;j++)
				 {
					 Amatrix[k,j] -= mul * Amatrix[i,j];
					 AInverse[k,j] -= mul * AInverse[i,j];
                }
			 }
         }
		 
		//解Ax = b，求出的x就是lambda就是源强度
        float[] BMatrix = new float[PanelNumber];
        float lambdaerror = 0.0f;
        for (int i = 0; i < PanelNumber; i++)
        {
            BMatrix[i] = -Vinf * 2.0f * Mathf.PI * Mathf.Cos(beta[i]);
        }
        for (int i = 0;i < PanelNumber;i++)
		{
			 for(int j = 0; j < PanelNumber;j++)
            {
                SourceStrength[i] += AInverse[i, j] * BMatrix[j];
            }

            lambdaerror += SourceStrength[i];
		}
		

		
		//计算其它参数
		//Debug.Log("Lambda总误差为 " + lambdaerror);//越小越好，可以通过增加PanelNumber减小误差，为零则无误差
        for (int i = 0; i < PanelNumber; i++)
        {
            float addval = 0.0f;
            for (int j = 0; j < PanelNumber; j++)
            {
                addval += (SourceStrength[j] / (2.0f * Mathf.PI)) * J[i, j];
            }
            Vt[i] = Vinf * Mathf.Sin(beta[i]) + addval;
            Cp[i] = 1 - (Vt[i] / Vinf) * (Vt[i] / Vinf);
        }
		
		StreamLineMat.SetFloatArray("XB",XB);
		StreamLineMat.SetFloatArray("YB",YB);
		StreamLineMat.SetFloatArray("SourceStrength",SourceStrength);
		StreamLineMat.SetFloatArray("XC",XC);
		StreamLineMat.SetFloatArray("YC",YC);
		StreamLineMat.SetFloat("Vinf", Vinf);
		if (Displaytype == DisplayType.StreamLine)
		StreamLineMat.SetInt("DisplayNumber",0);
		else if (Displaytype == DisplayType.StreamLineAndObject)
			StreamLineMat.SetInt("DisplayNumber", 1);
	}
	
	
    private void OnRenderImage(RenderTexture source, RenderTexture destination)
    {
        Graphics.Blit(source, destination,StreamLineMat);
    }
	void NACA1412()
	{
		
		XB[0] = 1.00000f;YB[0] = 0.00126f;
		XB[1] = 0.95025f;YB[1] = 0.00966f;
		XB[2] = 0.90040f;YB[2] = 0.01753f;
		XB[3] = 0.80058f;YB[3] = 0.03178f;
		XB[4] = 0.70061f;YB[4] = 0.04413f;
		XB[5] = 0.60051f;YB[5] = 0.05453f;
		XB[6] = 0.50029f;YB[6] = 0.06267f;
		XB[7] = 0.40000f;YB[7] = 0.06803f;
		XB[8] = 0.29925f;YB[8] = 0.06940f;
		XB[9] = 0.24889f;YB[9] = 0.06799f;
		XB[10] = 0.19857f;YB[10] = 0.06486f;
		XB[11] = 0.14833f;YB[11] = 0.05951f;
		XB[12] = 0.09824f;YB[12] = 0.05118f;
		XB[13] = 0.07330f;YB[13] = 0.04537f;
		XB[14] = 0.04845f;YB[14] = 0.03786f;
		XB[15] = 0.02378f;YB[15] = 0.02733f;
		XB[16] = 0.01158f;YB[16] = 0.01954f;
		XB[17] = 0.00000f;YB[17] = 0.00000f;
		XB[18] = 0.01342f;YB[18] = -0.01830f;
		XB[19] = 0.02622f;YB[19] = -0.02491f;
		XB[20] = 0.05155f;YB[20] = -0.03318f;
		XB[21] = 0.07670f;YB[21] = -0.03857f;
		XB[22] = 0.10176f;YB[22] = -0.04242f;
		XB[23] = 0.15167f;YB[23] = -0.04733f;
		XB[24] = 0.20143f;YB[24] = -0.04986f;
		XB[25] = 0.25111f;YB[25] = -0.05081f;
		XB[26] = 0.30075f;YB[26] = -0.05064f;
		XB[27] = 0.40000f;YB[27] = -0.04803f;
		XB[28] = 0.49971f;YB[28] = -0.04321f;
		XB[29] = 0.59949f;YB[29] = -0.03675f;
		XB[30] = 0.69939f;YB[30] = -0.02913f;
		XB[31] = 0.79942f;YB[31] = -0.02066f;
		XB[32] = 0.89960f;YB[32] = -0.01141f;
		XB[33] = 0.94975f;YB[33] = -0.00646f;
		XB[34] = 1.00000f;YB[34] = -0.00126f;
		
		XB[35] = XB[0];
		YB[35] = YB[0];
		PanelNumber = 35;
	}
}
