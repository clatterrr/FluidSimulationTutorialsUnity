Shader "Unlit/Advection_S"
{
    Properties
    {
        _MainTex ("Texture", 2D) = "white" {}
		VelocityTex("VelocityTex", 2D) = "white" {}
		VelocityDensityTex("VelocityDensityTex", 2D) = "white" {}
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
		    sampler2D VelocityDensityTex;
			float4 VelocityTex_TexelSize;
		    float dt;
			float dissipation;
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

			float4 bilerp(sampler2D_half sam,float2 p) {
				float4 st;
				st.xy = floor(p - 0.5) + 0.5;
				st.zw = st.xy + 1.0;

				float4 uv = st * VelocityTex_TexelSize.xyxy;
				float4 a = tex2D(sam, uv.xy); 
				float4 b = tex2D(sam, uv.zy); 
				float4 c = tex2D(sam, uv.xw); 
				float4 d = tex2D(sam, uv.zw); 
				float2 f = p - st.xy;
				return lerp(lerp(a, b, f.x), lerp(c, d, f.x), f.y);
			}
			float4 frag(v2f i) :SV_Target{
				float2 coord = i.uv* VelocityTex_TexelSize.zw - dt * tex2D(VelocityTex, i.uv);
				float4 col = dissipation * bilerp(VelocityDensityTex, coord);
				col.a = 1.0;
				return col;
			}
            ENDCG
        }
    }
}
