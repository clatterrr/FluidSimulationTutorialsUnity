Shader "Unlit/Display"
{
    Properties
    {
        _MainTex ("Texture", 2D) = "white" {}
    _VelocityTex("Texture", 2D) = "white" {}
    _DyeTex("Texture", 2D) = "white" {}

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
            sampler2D _VelocityTex;
            sampler2D _DyeTex;
            float2 _VelocityTex_TexelSize;
            float4 _MainTex_ST;

            v2f vert (appdata v)
            {
                v2f o;
                o.vertex = UnityObjectToClipPos(v.vertex);
                o.uv = TRANSFORM_TEX(v.uv, _MainTex);
                UNITY_TRANSFER_FOG(o,o.vertex);
                return o;
            }

            fixed4 frag(v2f i) : SV_Target
            {
                float2 coord = i.uv -  _VelocityTex_TexelSize.x * tex2D(_VelocityTex,i.uv).xy;
                float4 col = tex2D(_DyeTex, coord);
                return col;
            }
            ENDCG
        }
    }
}
