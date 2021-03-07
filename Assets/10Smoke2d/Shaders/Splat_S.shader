Shader "Unlit/Splat_S"
{
    Properties
    {
        _MainTex ("Texture", 2D) = "white" {}
		VelocityDensityTex("VelocityDensityTex", 2D) = "white" {}
		pointerpos ("pointerpos", vector) = (0,0,0,0)
		radius ("radius", float) = 0
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

			sampler2D VelocityDensityTex;
			float3 color;
			float2 pointerpos;
			float radius;
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
				float2 p = i.uv - pointerpos.xy/_ScreenParams.xy;
				p.x *= (_ScreenParams.x/_ScreenParams.y);
				half3 splat = exp(-dot(p,p)/radius)*color;
				half3 base = tex2D(VelocityDensityTex,i.uv).rgb;
				return half4(base+splat,1.);
			}
            ENDCG
        }
    }
}
