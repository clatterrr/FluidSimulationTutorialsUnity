Shader "Unlit/UpdateVelocity"
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
            uniform float _ArrowRotation[256];

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

               int _TexelNumber = 16;
               int intuvx = floor(i.uv.x * _TexelNumber);
               int intuvy = floor(i.uv.y * _TexelNumber);
               float rad = _ArrowRotation[intuvy * _TexelNumber + intuvx] * 3.1415927f / 180.0f;
               float cosRes = cos(rad);
               float sinRes = sin(rad);
               return float4(cosRes,sinRes, 0.0f, 1.0f);
            }
            ENDCG
        }
    }
}