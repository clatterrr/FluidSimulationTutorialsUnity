Shader "Unlit/Divergence_S"
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
			sampler2D _MainTex;
			float4 VelocityTex_TexelSize;
			float4 _MainTex_ST;
			float2 TexelSize;

            v2f vert (appdata v)
            {
                v2f o;
                o.vertex = UnityObjectToClipPos(v.vertex);
                o.uv = TRANSFORM_TEX(v.uv, _MainTex);
                UNITY_TRANSFER_FOG(o,o.vertex);
                return o;
            }

            fixed4 frag (v2f i) : SV_Target
            {
                float T = tex2D(VelocityTex, i.uv + float2(0.0f, VelocityTex_TexelSize.y)).y;
                float B = tex2D(VelocityTex, i.uv + float2(0.0f, -VelocityTex_TexelSize.y)).y;
                float R = tex2D(VelocityTex, i.uv + float2(VelocityTex_TexelSize.x, 0.0f)).x;
                float L = tex2D(VelocityTex, i.uv + float2(-VelocityTex_TexelSize.x, 0.0f)).x;
                float divergence = 0.5f * (R - L + T - B);
                return float4(divergence, 0.0f, 0.0f, 0.0f);
            }
            ENDCG
        }
    }
}
