Shader "VortexStreet/DisplayRainow"
{
    Properties
    {
		_MainTex("Texture", 2D) = "white" {}
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

			float4 frag (v2f i) : SV_Target
            {
				float4 col = tex2D(_MainTex, i.uv);
			    float x = (col.x + 1.0) / 2;//value的值在-1到1之间，x的值在0到1之间
				if (x < 0.25)  col = float4(0.0, 4.0 * x, 1.0, 1.0f);
			    else if (x < 0.5)  col = float4(0.0, 1.0, 1.0 + 4.0 * (0.25 - x), 1.0f);
			    else if (x < 0.75)    col = float4(4.0 * (x - 0.5), 1.0, 0.0, 1.0f);
			    else  col = float4(1.0, 1.0 + 4.0 * (0.75 - x), 0.0,1.0f);
				if (tex2D(BlockTex, i.uv).x > 0.9f)col = float4(0.0f, 0.0f, 0.0f, 1.0f);
                return col;
            }
            ENDCG
        }
    }
}
