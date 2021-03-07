using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class SourceVortexManager : MonoBehaviour
{
    public Material Mat;

    public enum Type
    { 
	    势线和流线,
        势线,
        流线,
		x速度图,
		y速度图,
		势函数图,
		流函数图
    }

    public Type type;
	private int TypeNumber = 0;
    [Range(-2.0f, 2.0f)]
	public float UniformFlowSpeed = 1.0f;//均匀流的速度
	[Range(0.0f, 3.14f)]
	public float UniformFlowAngle = 0.0f;//均匀流的角度
	[Range(-2.0f, 2.0f)]
	public float SouceStrength1 = 0.0f;//源的强度
	[Range(0.0f, 1.0f)]
	public float SourcePosX1 = 0.3f;//源的位置
	[Range(0.0f, 1.0f)]
	public float SourcePosY1 = 0.5f;//源的位置
	[Range(-2.0f, 2.0f)]
	public float SouceStrength2 = 0.0f;//第二个源
	[Range(0.0f, 1.0f)]
	public float SourcePosX2 = 0.7f;
	[Range(0.0f, 1.0f)]
	public float SourcePosY2 = 0.5f;
	[Range(-2f, 2f)]
	public float VortexStrength = 0f;//漩涡强度
	[Range(0.0f, 1.0f)]
	public float VortexPosX = 0.5f;
	[Range(0.0f, 1.0f)]
	public float VortexPosY = 0.5f;
	[Range(-2.0f, 2.0f)]
	public float DoubletStrength = 0.0f;//Doublet强度
	[Range(0.0f, 1.0f)]
	public float DoubletPosX = 0.5f;
	[Range(.0f, 1.0f)]
	public float DoubletPosY = 0.5f;
	
	void Update()
	{
		if(type == Type.势线和流线)TypeNumber = 0;
		else if(type == Type.势线)TypeNumber = 1;
		else if(type == Type.流线)TypeNumber = 2;
		else if(type == Type.x速度图)TypeNumber = 3;
		else if(type == Type.y速度图)TypeNumber = 4;
		else if(type == Type.势函数图)TypeNumber = 5;
		else if(type == Type.流函数图)TypeNumber = 6;
	}
	
	void OnRenderImage(RenderTexture source, RenderTexture destination)
	{
		
		Mat.SetInt("type", TypeNumber);
		Mat.SetFloat("UniformFlowSpeed", UniformFlowSpeed);
		Mat.SetFloat("UniformFlowAngle", UniformFlowAngle);
		Mat.SetFloat("SouceStrength1", SouceStrength1);
		Mat.SetFloat("SourcePosX1", SourcePosX1);
		Mat.SetFloat("SourcePosY1", SourcePosY1);
		Mat.SetFloat("SouceStrength2", SouceStrength2);
		Mat.SetFloat("SourcePosX2", SourcePosX2);
		Mat.SetFloat("SourcePosY2", SourcePosY2);
		Mat.SetFloat("VortexStrength", VortexStrength);
		Mat.SetFloat("VortexPosX", VortexPosX);
		Mat.SetFloat("VortexPosY", VortexPosY);
		Mat.SetFloat("DoubletStrength", DoubletStrength);
		Mat.SetFloat("DoubletPosX", DoubletPosX);
		Mat.SetFloat("DoubletPosY", DoubletPosY);
		Graphics.Blit(source, destination,Mat);
	}
}
