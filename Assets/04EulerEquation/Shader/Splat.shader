Shader "Unlit/Splat"
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
            sampler2D _VelocityTex;
            float PointerX;
            float PointerY;
            float PointerDX;
            float PointerDY;
            float4 _MainTex_ST;
            int MouseDown;

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
                float2 pointeruv = float2(PointerX,PointerY);
                float2 p = i.uv - pointeruv;
                float radius = 0.001;//圆的半径
                float3 color = float3(PointerDX, PointerDY, 0.0f) * 50.0f;//根据鼠标的方向产生不同方向的速度
                float3 splat = pow(2.1,-dot(p, p) / radius) * color;
                float3 base = tex2D(_VelocityTex, i.uv).xyz;
                if (MouseDown == 1)
                    base += splat;
                float4 col = float4(base, 1.0f);
                return col;
            }
            ENDCG
        }
    }
}
