using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class newtonSphere : MonoBehaviour
{
    public GameObject spherefab;
    private GameObject ball;
    private int NowFrame = 0;
    float valx = 15, valy = 15;//最初猜的值
    void Start()
    {
        ball = Instantiate(spherefab, new Vector3(valx, Fun(0, valx, valy) * 0.001f, valy), Quaternion.identity);
    }
    float Fun(float n, float x, float y)
    {
        if (n == 0) return (x + 2 * y) * (5 * x + 6 * y) + x - 432;
        else return -2 * (x + 2 * y) * (5 * x + 6 * y) + y + 854;
    }
    void Update()
    {
        NowFrame++;
        if (NowFrame < 60) return;
        float a = (Fun(0, valx + 1, valy) - Fun(0, valx - 1, valy)) / 2;
        float b = (Fun(0, valx, valy + 1) - Fun(0, valx, valy - 1)) / 2;
        float c = (Fun(1, valx + 1, valy) - Fun(1, valx - 1, valy)) / 2;
        float d = (Fun(1, valx, valy + 1) - Fun(1, valx, valy - 1)) / 2;
        float e = Fun(0, valx, valy);
        float f = Fun(1, valx, valy);
        float div = 1 / (a * d - b * c);
        float speed = 0.1f;
        valx = valx - speed * (d * e - b * f) * div;
        valy = valy - speed * (-c * e + a * f) * div;
        ball.transform.position = new Vector3(valx, Fun(0, valx, valy) * 0.001f, valy);
        Debug.Log("x = " + valx + " y = " + valy);

    }
}
