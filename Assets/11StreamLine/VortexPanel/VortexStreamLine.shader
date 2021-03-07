Shader "VortexPanel/VortexStreamLine"
{
    Properties
    {
        _MainTex ("Texture", 2D) = "white" {}
    }
    SubShader
    {
        Tags { "RenderType"="Opaque" }
        LOD 100

        Pass
        {
            CGPROGRAM
            #pragma vertex vert
            #pragma fragment frag
            // make fog work
            #pragma multi_compile_fog

            #include "UnityCG.cginc"

            struct appdata
            {
                float4 vertex : POSITION;
                float2 uv : TEXCOORD0;
            };

            struct v2f
            {
                float2 uv : TEXCOORD0;
                UNITY_FOG_COORDS(1)
                float4 vertex : SV_POSITION;
            };

            sampler2D _MainTex;
            float4 _MainTex_ST;
            uniform float VortexStrength[256];
            uniform float XC[256];
            uniform float YC[256];
            uniform float XB[256];
            uniform float YB[256];
            int DisplayNumber;
            float Vinf;
            v2f vert(appdata v)
            {
                v2f o;
                o.vertex = UnityObjectToClipPos(v.vertex);
                o.uv = TRANSFORM_TEX(v.uv, _MainTex);
                UNITY_TRANSFER_FOG(o, o.vertex);
                return o;
            }
            fixed4 frag(v2f i) : SV_Target
            {
                float r,psi = 0.0f;
                for (int k = 0; k < 256; k++)
                {
                  if (VortexStrength[k] == 0.0f)break;
                  r = sqrt((i.uv.y - YC[k])*(i.uv.y - YC[k]) + (i.uv.x - XC[k])* (i.uv.x - XC[k]));
                  //因为theta在180度和-180度时不连续
                  //所以将VortexStrength的精度控制在0.1以绘制连续的流线
                  VortexStrength[k] /= 100.0f;
                  VortexStrength[k] = floor(VortexStrength[k] * 10.0f) / 10.0f;
                  psi += VortexStrength[k] / (2 * 3.14159) * log(r);
                }
                psi += Vinf * i.uv.y;//均匀流造成的势能影响
                float4 col = float4(1.0f, 1.0f, 1.0f, 1.0f);
                //VortexStrength的精度为10/100 = 0.1，如果高于这个精度则流线不连续
                //这是着色器计算theta时的问题，与数值方法无关
                if (floor(psi * 20) % 10 == 0)
                {
                    col -= float4(1.0f, 1.0f, 0.0f, 0.5f);//蓝色的线是流线
                }
                //沿着平板端点绘制圆柱/机翼本身
                if (DisplayNumber == 1)
                  for (int k = 0; k < 256; k++)
                  {
                    if (XB[k + 1] == 0.0 && YB[k + 1] == 0.0)break;
                    float p1x = XB[k] + 0.0f;
                    float p1y = YB[k] + 0.0f;
                    float p2x = XB[k + 1] + 0.0f;
                    float p2y = YB[k + 1] + 0.0f;
                    if (abs((i.uv.x - p1x) * (p2y - p1y) - (p2x - p1x) * (i.uv.y - p1y)) < 0.001f
                       && min(p1x,p2x)-0.005f <= i.uv.x && i.uv.x <= max(p1x,p2x)+0.005f
                       && min(p1y,p2y)-0.005f <= i.uv.y && i.uv.y <= max(p1y,p2y) + 0.005f)
                    {
                         col = float4(0.0f,0.0f,0.0f,1.0f);
                         break;
                    }
                  }

                return col;
                return col;
            }
            ENDCG
        }
    }
}
