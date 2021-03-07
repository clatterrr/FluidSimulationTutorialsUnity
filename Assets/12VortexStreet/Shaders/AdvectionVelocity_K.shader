Shader "VortexStreet/AdvectionVelocity_K"
{
    Properties
    {
		_MainTex("Texture", 2D) = "white" {}
		VelocityTex("VelocityTex", 2D) = "white" {}
		BlockTex("BlockTex", 2D) = "white" {}
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
			sampler2D BlockTex;
			float4 VelocityTex_TexelSize;
			float dt;
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
				float4 col = tex2D(VelocityTex, i.uv - 0.01f*tex2D(VelocityTex, i.uv));
				if(tex2D(BlockTex, i.uv).x > 0.99f)col.xy = float2(0.0f, 0.0f);
				
				//float2 dir = i.uv - float2(0.3f, 0.5f);
				//float dis = dir.x*dir.x + dir.y*dir.y;
				//if(dis < 0.001f)col.xy = float2(0.0f, 0.0f);
				if (i.uv.x < 0.01f)col.xy = float2(1.0f, 0.0f);
				return col;
			}
            ENDCG
        }
    }
}
