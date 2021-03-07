Shader "Unlit/Vorticity_S"
{
    Properties
    {
        _MainTex ("Texture", 2D) = "white" {}
		VelocityTex("VelocityTex", 2D) = "white" {}
		CurlTex("CurlTex", 2D) = "white" {}
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

			sampler2D VelocityTex;
			sampler2D CurlTex;
			float4 CurlTex_TexelSize;
			uniform half curl;
			uniform half dt;
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
				float L = tex2D(CurlTex,i.uv - float2(CurlTex_TexelSize.x,0.0f)).x;
				float R = tex2D(CurlTex,i.uv + float2(CurlTex_TexelSize.x, 0.0f)).x;
				float T = tex2D(CurlTex, i.uv + float2(0.0f, CurlTex_TexelSize.y)).x;
				float B = tex2D(CurlTex, i.uv - float2(0.0f, CurlTex_TexelSize.y)).x;

				float C = tex2D(CurlTex,i.uv).x;

				float2 force = float2(abs(T)-abs(B),abs(R)-abs(L));
				force *= 1./length(force + 0.00001) * curl * C;
				float2 vel = tex2D(VelocityTex, i.uv).xy;
				return float4(vel + force * dt ,0.0f,1.0f);
			}
            ENDCG
        }
    }
}
