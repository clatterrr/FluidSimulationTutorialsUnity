Shader "Unlit/Substract"
{
    Properties
    {
        _MainTex ("Texture", 2D) = "white" {}
		VelocityTex("VelocityTex", 2D) = "white" {}
		PressureTex ("Tex_Pressure", 2D) = "white" {}
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

			sampler2D PressureTex;
			sampler2D VelocityTex;
			float4 PressureTex_TexelSize;
			sampler2D _MainTex;
			float4 _MainTex_ST;

            v2f vert (appdata v)
            {
                v2f o;
                o.vertex = UnityObjectToClipPos(v.vertex);
                o.uv = TRANSFORM_TEX(v.uv, _MainTex);
                UNITY_TRANSFER_FOG(o,o.vertex);
                return o;
            }

			float4 frag(v2f i) :SV_Target{

				float L = tex2D(PressureTex,saturate(i.uv - float2(PressureTex_TexelSize.x,0.0f))).x;
				float R = tex2D(PressureTex,saturate(i.uv + float2(PressureTex_TexelSize.x, 0.0f))).x;
				float T = tex2D(PressureTex,saturate(i.uv + float2(0.0f, PressureTex_TexelSize.y))).x;
				float B = tex2D(PressureTex,saturate(i.uv - float2(0.0f, PressureTex_TexelSize.y))).x;

				float2 velocity = tex2D(VelocityTex, i.uv).xy;
				velocity.xy -= float2(R - L,T - B);
				return float4(velocity,0.,1.);
			}
            ENDCG
        }
    }
}
