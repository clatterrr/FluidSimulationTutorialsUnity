Shader "VortexStreet/Display_K"
{
    Properties
    {
        _MainTex ("Texture", 2D) = "white" {}
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

            sampler2D _MainTex;
            float4 _MainTex_ST;
			sampler2D BlockTex;

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
                fixed4 col = tex2D(_MainTex, i.uv);

				if (tex2D(BlockTex, i.uv).x > 0.9f)col.xyz = float3(0.9f, 0.5f,0.0f);
				//float2 dir = i.uv - float2(0.3f, 0.5f);
				//float dis = dir.x*dir.x + dir.y*dir.y;
				//if (dis < 0.001f)col.xyz = float3(0.9f, 0.5f, 0.0f);
                return col;
            }
            ENDCG
        }
    }
}
