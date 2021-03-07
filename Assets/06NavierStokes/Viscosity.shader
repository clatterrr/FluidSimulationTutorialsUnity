Shader "Unlit/Viscosity"
{
    Properties
    {
        _MainTex ("Texture", 2D) = "white" {}
       _VelocityTex("Texture", 2D) = "white" {}
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
            sampler2D _VelocityTex;
            float2 _VelocityTex_TexelSize;

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
                float2 Top = tex2D(_VelocityTex, i.uv + float2(0.0f, _VelocityTex_TexelSize.y)).xy;
                float2 Bottom = tex2D(_VelocityTex, i.uv + float2(0.0f, -_VelocityTex_TexelSize.y)).xy;
                float2 Right = tex2D(_VelocityTex, i.uv + float2(_VelocityTex_TexelSize.x, 0.0f)).xy;
                float2 Left = tex2D(_VelocityTex, i.uv + float2(-_VelocityTex_TexelSize.x, 0.0f)).xy;
                float2 Center = tex2D(_VelocityTex, i.uv ).xy;
                float nu = 0.1;//运动粘性系数
                float VelocityX = Center.x + nu * ((Right.x - 2 * Center.x + Left.x) + ( Top.x - 2 * Center.x + Bottom.x));
                float VelocityY = Center.y + nu * ((Right.y - 2 * Center.y + Left.y) + (Top.y - 2 * Center.y + Bottom.y));
                return float4(VelocityX, VelocityY, 0.0f, 0.0f);
            }
            ENDCG
        }
    }
}
