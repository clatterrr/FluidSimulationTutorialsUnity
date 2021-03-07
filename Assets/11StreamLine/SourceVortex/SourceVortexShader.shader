Shader "Unlit/SourceVortexShader"
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
		
		    float UniformFlowSpeed;
			float UniformFlowAngle;
			float SouceStrength1;
			float SouceStrength2;
			float SourcePosX1;
			float SourcePosX2;
			float SourcePosY1;
			float SourcePosY2;
			float VortexStrength;
			float VortexPosX;
			float VortexPosY;
			float DoubletStrength;
			float DoubletPosX;
			float DoubletPosY;
			int type;

            v2f vert (appdata v)
            {
                v2f o;
                o.vertex = UnityObjectToClipPos(v.vertex);
                o.uv = TRANSFORM_TEX(v.uv, _MainTex);
                UNITY_TRANSFER_FOG(o,o.vertex);
                return o;
            }
			
			float4 RainbowColor(float value)//value的值在-2到2之间
            {
                float3 col;
                float x = (value + 2)/4;
                if (x < 0.25)  col = fixed3(0.0, 4.0 * x, 1.0);
                else if (x < 0.5)  col = fixed3(0.0, 1.0, 1.0 + 4.0 * (0.25 - x));
                else if (x < 0.75)    col = fixed3(4.0 * (x - 0.5), 1.0, 0.0);
                else  col = fixed3(1.0, 1.0 + 4.0 * (0.75 - x), 0.0);
                return float4(col,1.0f);
            }

            fixed4 frag (v2f i) : SV_Target
            {
                float phi = 0.0f;//势函数
                float psi = 0.0f;//流函数
                float Vx = 0.0f, Vy = 0.0f;//x轴y轴方向上的速度
				float dx,dy,r,theta;
				
				Vx += UniformFlowSpeed * cos(UniformFlowAngle);
				Vy += UniformFlowSpeed * sin(UniformFlowAngle);
				phi += Vx * i.uv.x + Vy * i.uv.y;
                psi += Vx * i.uv.y - Vy * i.uv.x;
				
				//第一个源汇
				dx = i.uv.x - SourcePosX1;
				dy = i.uv.y - SourcePosY1;
				r = sqrt(dx * dx + dy * dy);
                theta = atan2(dy,dx);
                Vx += (SouceStrength1 * dx) / (2 * 3.14159 * r * r);
                Vy += (SouceStrength1 * dy) / (2 * 3.14159 * r * r);
                phi += (SouceStrength1 / (2 * 3.14159) * log(r));//源的势能没错
                psi += (SouceStrength1 / (2 * 3.14159) * theta);  
				
				//第二个源汇
				dx = i.uv.x - SourcePosX2;
				dy = i.uv.y - SourcePosY2;
				r = sqrt(dx * dx + dy * dy);
                theta = atan2(dy,dx);
                Vx += (SouceStrength2 * dx) / (2 * 3.14159 * r * r);
                Vy += (SouceStrength2 * dy) / (2 * 3.14159 * r * r);
                phi += (SouceStrength2 / (2 * 3.14159) * log(r));
                psi += (SouceStrength2 / (2 * 3.14159) * theta);  
				
				//漩涡
				dx = i.uv.x - VortexPosX;
                dy = i.uv.y - VortexPosY;
                r = sqrt(dx * dx + dy * dy);
                theta = atan2(dy,dx);
                Vx += (VortexStrength * dy) / (2 * 3.14159 * r * r);
                Vy += (- VortexStrength * dx) / (2 * 3.14159 * r * r);
                phi -= VortexStrength / (2 * 3.14159) * (theta);
                psi += VortexStrength / (2 * 3.14159) * log(r);
				
				//Doublet
				dx = i.uv.x - DoubletPosX;
                dy = i.uv.y - DoubletPosY;
                r = sqrt(dx * dx + dy * dy);
                theta = atan2(dy,dx);
                Vx += (-DoubletStrength * (dx * dx - dy * dy)) / (3.14159 * r * r * r * r);
                Vy += (-2.0f * DoubletStrength * dx * dy) / (3.14159 * r * r * r * r);
                phi += (DoubletStrength * dx) / (2 * 3.14159 * r * r);
                psi += (-DoubletStrength * dy) / (2 * 3.14159 * r * r);
				
				float4 col = float4(1.0f, 1.0f, 1.0f, 1.0f);
				if (floor(phi * 100) % 10 == 0)
                {
				    if(type == 0 || type == 1)
                    col -= float4(0.0f, 0.2f, 1.0f, 0.5f);//橙黄色的线是等势线
                }
				if (floor(psi * 100) % 10 == 0)
                {
				    if(type == 0 || type == 2)
                    col -= float4(1.0f, 1.0f, 0.0f, 0.5f);//蓝色的线是流线
                }
				if(type == 3)col = RainbowColor(Vx);
				if(type == 4)col = RainbowColor(Vy);
				if(type == 5)col = RainbowColor(phi);
				if(type == 6)col = RainbowColor(psi);
				
				return col;
            }
            ENDCG
        }
    }
}
