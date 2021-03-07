Shader "Unlit/Curl_S"
{
    Properties
    {
        _MainTex ("Texture", 2D) = "white" {}
		VelocityTex("VelocityTex", 2D) = "white" {}
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
			float4 VelocityTex_TexelSize;
			float2 TexelSize;
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
				float L = tex2D(VelocityTex,i.uv - float2(VelocityTex_TexelSize.x,0.0f)).y;
				float R = tex2D(VelocityTex,i.uv + float2(VelocityTex_TexelSize.x, 0.0f)).y;
				float T = tex2D(VelocityTex, i.uv + float2(0.0f, VelocityTex_TexelSize.y)).x;
				float B = tex2D(VelocityTex, i.uv - float2(0.0f, VelocityTex_TexelSize.y)).x;
				float vorticity = 0.5f*(R - L - T + B);
				return float4(vorticity,0.,0.,1.);
			}
            ENDCG
        }
    }
}
