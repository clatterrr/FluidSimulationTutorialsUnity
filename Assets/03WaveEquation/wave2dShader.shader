Shader "Unlit/wave2dShader"
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
            sampler2D _TimeTexM1;
            int _TexelNumber;
            float _Speed;

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
                float _Height = 0.0f;

                float FloatNumber = _TexelNumber;

                //把一个区域内的uv值都改为该区域左下角像素的uv值
                int flag = 0;//检测uv是否位于中心方块

                float error = 0.01f;//取样误差
                int intuvx = floor(i.uv.x * _TexelNumber);
                if (intuvx == (_TexelNumber - 1) / 2)flag++;
                float uvx = intuvx / FloatNumber + error;


                int intuvy = floor(i.uv.y * _TexelNumber);
                if (intuvy == (_TexelNumber - 1) / 2)flag++;
                float uvy = intuvy / FloatNumber + error;

                if (intuvx == 0 || intuvx == _TexelNumber - 1 || intuvy == 0 || intuvy == _TexelNumber - 1)
                    _Height = 0.0f;
                else
                {


                        float uvxPlus = (intuvx + 1.0f) / FloatNumber + error;
                        float uvxMinus = (intuvx - 1.0f) / FloatNumber + error;
                        float uvyPlus = (intuvy + 1.0f) / FloatNumber + error;
                        float uvyMinus = (intuvy - 1.0f) / FloatNumber + error;


                        float2 uv = float2(uvx, uvy);
                        float TimeMinus2Center = tex2D(_MainTex, uv).g;
                        float centerT = tex2D(_TimeTexM1, uv).g;
                        float FSource = 0.0f;
                        if (flag == 2)
                        {
                            FSource = sin(_Time.y * 10.0f) * 50.0f;
                        }

                        uv = float2(uvxPlus , uvy);
                        float rightT = tex2D(_TimeTexM1, uv).g;

                         uv = float2(uvxMinus, uvy);
                         float leftT = tex2D(_TimeTexM1, uv).g;

                         uv = float2(uvx, uvyPlus);
                         float topT = tex2D(_TimeTexM1, uv).g;

                         uv = float2(uvx, uvyMinus);
                         float BottomT = tex2D(_TimeTexM1, uv).g;

                         _Speed = 0.6f;
                         _Height = 2 * centerT - TimeMinus2Center + (leftT + rightT + BottomT + topT - 4 * centerT) * _Speed * _Speed + FSource;
                }
                return float4(1.0f, _Height, 0.0f, 1.0f);
            }
            ENDCG
        }
    }
}
